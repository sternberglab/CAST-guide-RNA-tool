import json
import time
import os
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from src.genbank import retrieve_annotation, get_regions
from src.finder import get_target_region_for_gene, get_candidates_for_region, remove_offtarget_matches, order_candidates_for_region
from src.outputs import make_spacer_gen_output, make_eval_outputs
from src.bowtie import find_offtargets

def spacer_gen(args):
	# unpack the arguments
	email = args['email']
	output_path = args['output_path']
	genbank_id = args['genbank_id']
	start_pct = args['start_pct']
	end_pct = args['end_pct']
	spacers_per_region = args['spacers_per_region']
	GC_requirement = args['GC_requirement']
	overlapping_spacers = args['overlapping_spacers']
	target_gene_names = args['target_gene_names']
	coding_spacer_direction = args['coding_spacer_direction']
	region_type = args['region_type']
	noncoding_boundary = args['noncoding_boundary']
	custom_regions = args['custom_regions']

	# Check for required parameters
	if not email or '@' not in email:
		print("Please enter an email for NCBI API calls")
		return
	if region_type not in ['coding', 'noncoding', 'custom']:
		print("region_type must be either 'coding', 'noncoding', or 'custom")
		return

	# Load genome data
	print("Starting spacer search...")
	genbank_info = retrieve_annotation(genbank_id, email)
	genome = Seq(genbank_info['GBSeq_sequence'])

	# Generate coding regions for each region_type
	if region_type == 'coding':
		all_genes = get_regions(genbank_info, coding=True)
		# match gene names case insensitive, lowercase both sides
		target_gene_names = [t.lower() for t in target_gene_names]
		target_genes = [gene for gene in all_genes if gene['name'].lower() in target_gene_names]
		if not len(target_genes):
			print("You must enter at least one valid gene name")
			return
		regions = target_genes

	if region_type == 'noncoding':
		all_noncoding = get_regions(genbank_info, coding=False)
		regions = [r for r in all_noncoding if r['end'] > noncoding_boundary[0] and r['start'] < noncoding_boundary[1]]

	if region_type == 'custom':
		regions = [{'name': f'custom-{index}', 'start': c[0], 'end': c[1], 'direction': 'fw'} for (index, c) in enumerate(custom_regions)]
		start_pct = 0
		end_pct = 100

	for region in regions:
		start = time.perf_counter()
		print(region['name'])
		print(f"Finding gRNA for \"{region['name']}\"")
		[start_mark, end_mark] = get_target_region_for_gene(region, start_pct, end_pct)
		
		candidates = get_candidates_for_region(genome, start_mark, end_mark, region['name'], GC_requirement)
		if region_type == 'coding':
			if coding_spacer_direction not in ['N_to_C', 'C_to_N']:
				print("Invalid 'coding_spacer_direction' parameter, see the valid inputs")
				return
			candidates = order_candidates_for_region(candidates, region, coding_spacer_direction)
		
		candidates = remove_offtarget_matches(genbank_id, region['name'], candidates, spacers_per_region, overlapping_spacers)
		region['candidates'] = candidates
		if len(candidates) > spacers_per_region:
			region['candidates'] = candidates[:spacers_per_region]

		elapsed_time = round(time.perf_counter() - start, 2)
		print(f"Identified {len(candidates)} spacers for {region['name']} in {elapsed_time} seconds")
	make_spacer_gen_output(regions, os.path.join(output_path, 'spacer_gen_output'))

def spacer_eval(args):
	genbank_id = args['genbank_id']
	output_path = args['output_path']
	email = args['email']
	user_spacers = args['spacers']

	if not email or '@' not in email:
		print("Please enter an email for NCBI API calls")
		return

	genbank_info = retrieve_annotation(genbank_id, email)
	genome = Seq(genbank_info['GBSeq_sequence'])
	
	spacers = [s for s in user_spacers if len(s) == 32]
	print(f"Starting evaluation for {len(spacers)} spacers...")
	spacer_batch = []
	for (index, spacer) in enumerate(spacers):
		flexible_seq = spacer[:]
		for i in range(5, 32, 6):
			flexible_seq = flexible_seq[:i] + 'N' + flexible_seq[i+1:]
		spacer_record = SeqRecord(Seq(flexible_seq), id=f'spacer_{index+1}', description=genbank_id)  # use 'description' to store ref_genome info
		spacer_batch.append(spacer_record)

	# Write the batch of candidate sequences to a fasta for bowtie2 to use
	fasta_name = 'offtarget-check.fasta'
	root_dir = Path(__file__).parent.parent
	fasta_name = os.path.join(root_dir, 'assets', 'bowtie', genbank_id, fasta_name)
	with open(fasta_name, 'w') as targets_file:
		SeqIO.write(spacer_batch, targets_file, 'fasta')

	output_location = find_offtargets(genbank_id, fasta_name)
	make_eval_outputs(spacer_batch, output_location, genome, output_path)

	os.remove(fasta_name)
	os.remove(output_location)