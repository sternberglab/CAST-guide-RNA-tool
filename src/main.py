import json
import time
import os
import csv
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from src.genbank import retrieve_annotation, get_regions
from src.finder import get_target_region_for_gene, get_candidates_for_region, remove_offtarget_matches, order_candidates_for_region
from src.outputs import make_spacer_gen_output, make_eval_outputs
from src.parse import extract_column_from_csv
from src.bowtie import find_offtargets
from src.advanced_parameters import SPACER_LENGTH, flex_base, flex_spacing

def spacer_gen(args):
	# unpack the arguments
	email = args['email']
	nonessential_only = args['nonessential_only']
	output_path = args['output_path']
	genbank_id = args['genbank_id']
	genbank_file = args['genbank_file']
	genome_fasta_file = args['genome_fasta_file']
	start_pct = args['start_pct']
	end_pct = args['end_pct']
	spacers_per_region = args['spacers_per_region']
	GC_requirement = args['GC_requirement']
	overlapping_spacers = args['overlapping_spacers']
	target_locus_tags_csv = args['target_locus_tags_csv']
	coding_spacer_direction = args['coding_spacer_direction']
	region_type = args['region_type']
	noncoding_boundary = args['noncoding_boundary']
	custom_regions_csv = args['custom_regions_csv']

	# Check for required parameters
	if region_type not in ['coding', 'noncoding', 'custom']:
		print("region_type must be either 'coding', 'noncoding', or 'custom")
		return

	# Load genome data
	if genbank_id:
		if not email or '@' not in email:
			print("Please enter an email for NCBI API calls")
			return
		record = retrieve_annotation(genbank_id, email)
	elif genbank_file:
		if not Path(genbank_file).exists():
			print("Invalid genbank_file provided. Either leave it empty or provide a valid file path. ")
			return
		with open(genbank_file, 'r') as f:
			record = SeqIO.read(Path(genbank_file), 'genbank')
			genbank_id = record.id
	elif genome_fasta_file:
		if not Path(genome_fasta_file).exists():
			print("Invalid genome_fasta_file provided. Either leave it empty or provide a valid file path. ")
			return
		if not region_type == 'custom':
			print("Cannot do coding or noncoding spacer generation against a fasta file. Use a genbank input or use custom regions against this fasta. ")
			return
		record =SeqIO.read(Path(genome_fasta_file), "fasta")
		genbank_id = record.id
	genome = record.seq.upper()

	print("Starting spacer search...")
	# Generate coding regions for each region_type
	if region_type == 'coding':
		all_genes = get_regions(record, region='coding')

		# match gene names case insensitive, lowercase both sides
		target_locus_tag_ids = extract_column_from_csv(target_locus_tags_csv, 'locus_tags')
		target_genes = [gene for gene in all_genes if gene['name'] in target_locus_tag_ids]
		if not len(target_genes):
			print(f"You must enter at least one valid locus tag identifier, ex. {all_genes[0]['name']}")
			return
		regions = target_genes

	if region_type == 'noncoding':
		if len(noncoding_boundary) != 2:
			print("You must specify the boundaries around which to search for noncoding regions")
			return
		region_name = 'nonessential' if nonessential_only else 'noncoding'
		all_noncoding = get_regions(record, region=region_name)
		regions = [r for r in all_noncoding if r['end'] > noncoding_boundary[0] and r['start'] < noncoding_boundary[1]]
		start_pct = 0
		end_pct = 100

	if region_type == 'custom':
		if not Path(custom_regions_csv).exists():
			print("Must input a valid filepath for the custom regions csv")
			return
		with open(Path(custom_regions_csv), 'r', encoding='utf-8-sig') as csv_file:
			reader = csv.DictReader(csv_file)
			custom_regions = [[int(r['start_ref']), int(r['end_ref'])] for r in reader]
		regions = [{'name': f'custom-{index}', 'start': c[0], 'end': c[1], 'direction': 'fw'} for (index, c) in enumerate(custom_regions)]
		start_pct = 0
		end_pct = 100

	for region in regions:
		start = time.perf_counter()
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

	record = retrieve_annotation(genbank_id, email)
	genome = record.seq.upper()
	
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
	fasta_name = os.path.join(root_dir, 'assets', 'bowtie', genbank_id, fasta_name)
	with open(fasta_name, 'w') as targets_file:
		SeqIO.write(spacer_batch, targets_file, 'fasta')

	output_location = find_offtargets(genbank_id, fasta_name)
	make_eval_outputs(spacer_batch_unmod, output_location, genome, output_path)

	os.remove(fasta_name)
	# os.remove(output_location)