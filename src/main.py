import json
import time
import os
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from src.genbank import retrieve_annotation, get_regions
from src.finder import get_target_region_for_gene, get_candidates_for_region, remove_offtarget_matches
from src.outputs import make_spacer_gen_output, make_eval_outputs
from src.bowtie import find_offtargets
from src.advanced_parameters import SPACER_LENGTH, flex_base, flex_spacing


def spacer_gen(args):
	email = args['email']
	outputPath = args['outputPath']
	genbankId = args['genbankId']
	startPct = args['startPct']
	endPct = args['endPct']
	spacersPerRegion = args['spacersPerRegion']
	target_gene_names = args['target_gene_names']
	regionType = args['regionType']
	noncoding_boundary = args['noncoding_boundary']
	custom_regions = args['custom_regions']

	if not email or '@' not in email:
		print("Please enter an email for NCBI API calls")
		return
	if regionType not in ['coding', 'noncoding', 'custom']:
		print("regionType must be either 'coding', 'noncoding', or 'custom")
		return

	print("Starting spacer search...")
	genbank_info = retrieve_annotation(genbankId, email)
	genome = Seq(genbank_info['GBSeq_sequence'])

	if regionType == 'coding':
		all_genes = get_regions(genbank_info, coding=True)
		# match gene names case insensitive, lowercase both sides
		target_gene_names = [t.lower() for t in target_gene_names]
		target_genes = [gene for gene in all_genes if gene['name'].lower() in target_gene_names]
		if not len(target_genes):
			print("You must enter at least one valid gene name")
			return
		regions = target_genes

	if regionType == 'noncoding':
		all_noncoding = get_regions(genbank_info, coding=False)
		regions = [r for r in all_noncoding if r['end'] > noncoding_boundary[0] and r['start'] < noncoding_boundary[1]]

	if regionType == 'custom':
		regions = [{'name': f'custom-{index}', 'start': c[0], 'end': c[1], 'direction': 'fw'} for (index, c) in enumerate(custom_regions)]
		startPct = 0
		endPct = 100

	for region in regions:
		start = time.perf_counter()
		print(region['name'])
		print(f"Finding crRNA for \"{region['name']}\"")
		[start_mark, end_mark] = get_target_region_for_gene(region, startPct, endPct)
		
		candidates = get_candidates_for_region(genome, start_mark, end_mark, region['name'])
		# print("{} candidates with a valid PAM targeting the region".format(len(candidates)))
		# candidates = filter_non_unique_fingerprints(candidates)
		# print("{} candidates with a unique fingerprint".format(len(candidates)))
		
		candidates = remove_offtarget_matches(genbankId, region['name'], candidates, spacersPerRegion)
		region['candidates'] = candidates
		if len(candidates) > spacersPerRegion:
			region['candidates'] = candidates[:spacersPerRegion]

		elapsed_time = round(time.perf_counter() - start, 2)
		print(f"Identified {len(candidates)} spacers for {region['name']} in {elapsed_time}")
	make_spacer_gen_output(regions, os.path.join(outputPath, 'spacer_gen_output'))

def spacer_eval(args):
	genbankId = args['genbankId']
	outputPath = args['outputPath']
	email = args['email']
	user_spacers = args['spacers']

	if not email or '@' not in email:
		print("Please enter an email for NCBI API calls")
		return

	genbank_info = retrieve_annotation(genbankId, email)
	genome = Seq(genbank_info['GBSeq_sequence'])
	
	spacers = [s for s in user_spacers if len(s) == SPACER_LENGTH]
	print(f"Starting evaluation for {len(spacers)} spacers...")
	spacer_batch = []
	spacer_batch_unmod = []  # make unmodified copy of spacer recs for output
	if flex_base:
		for (index, spacer) in enumerate(spacers):
			flexible_seq = spacer[:]
			for i in range(flex_spacing-1, SPACER_LENGTH, flex_spacing):
				flexible_seq = flexible_seq[:i] + 'N' + flexible_seq[i+1:]
			spacer_record = SeqRecord(Seq(flexible_seq), id=f'spacer_{index+1}', description=genbankId)  # use 'description' to store ref_genome info
			spacer_batch.append(spacer_record)
			spacer_record_unmod = SeqRecord(Seq(spacer.upper()), id=f'spacer_{index + 1}',
									  description=genbankId)  # unmodified spacer for output
			spacer_batch_unmod.append(spacer_record_unmod)
	else:
		for (index, spacer) in enumerate(spacers):
			spacer_record = SeqRecord(Seq(spacer), id=f'spacer_{index + 1}',
									  description=genbankId)  # use 'description' to store ref_genome info
			spacer_batch.append(spacer_record)
			spacer_batch_unmod.append(spacer_record)

	# Write the batch of candidate sequences to a fasta for bowtie2 to use
	fasta_name = 'offtarget-check.fasta'
	root_dir = Path(__file__).parent.parent
	fasta_name = os.path.join(root_dir, 'assets', 'bowtie', genbankId, fasta_name)
	with open(fasta_name, 'w') as targets_file:
		SeqIO.write(spacer_batch, targets_file, 'fasta')

	output_location = find_offtargets(genbankId, fasta_name)
	make_eval_outputs(spacer_batch_unmod, output_location, genome, outputPath)

	os.remove(fasta_name)
	# os.remove(output_location)