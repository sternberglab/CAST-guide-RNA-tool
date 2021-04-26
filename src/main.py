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
from src.filters import filter_non_unique_fingerprints, filter_re_sites, filter_homopolymers
from src.advanced_parameters import SPACER_LENGTH, flex_base, flex_spacing

def spacer_gen(args):
	# unpack the arguments
	email = args['email']
	nonessential_only = args['nonessential_only']
	output_path = args['output_path']
	genbank_ids = args['genbank_ids']
	genbank_files = args['genbank_files']
	genome_fasta_files = args['genome_fasta_files']
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
	custom_sequences = args['custom_sequences']

	# Check for required parameters
	if region_type not in ['coding', 'noncoding', 'custom']:
		print("region_type must be either 'coding', 'noncoding', or 'custom")
		return

	# Load genome data
	genome_input_type = None
	if len(genbank_ids):
		if not email or '@' not in email:
			print("Please enter an email for NCBI API calls")
			return
		genome_input_type = 'genbank_ids'
		default_genome_id = genbank_ids[0]
		record = retrieve_annotation(default_genome_id, email)
		for genbank_id in genbank_ids:
			retrieve_annotation(genbank_id, email, return_record=False)
		genome_ids = genbank_ids

	elif len(genbank_files):
		for genbank_file in genbank_files:
			if not Path(genbank_file).exists():
				print("Invalid genbank_file provided. Either leave it empty or provide a valid file path. ")
				return
		genome_input_type = 'genbank_files'
		default_genome_id = genbank_files[0]
		record = SeqIO.read(Path(default_genome_id), 'genbank')
		genome_ids = []
		for genbank_file in genbank_files:
			record = SeqIO.read(Path(default_genome_id), 'genbank')
			genome_ids += [record.id]
	elif genome_fasta_files:
		for genome_fasta_file in genome_fasta_files:
			if not Path(genome_fasta_file).exists():
				print("Invalid genome_fasta_file provided. Either leave it empty or provide a valid file path. ")
				return
		if not region_type == 'custom':
			print("Cannot do coding or noncoding spacer generation against a fasta file. Use a genbank input or use custom regions against this fasta. ")
			return
		genome_input_type = 'fasta_files'
		genome_ids = []
		for genome_fasta_file in genome_fasta_files:
			record = SeqIO.read(Path(genome_fasta_file), "fasta")
			genome_ids += [record.id]

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
		for r in regions:
			r['genome_input_type'] = genome_input_type
			r['genome_id'] = default_genome_id

	if region_type == 'noncoding':
		if len(noncoding_boundary) != 2:
			print("You must specify the boundaries around which to search for noncoding regions")
			return
		region_name = 'nonessential' if nonessential_only else 'noncoding'
		all_noncoding = get_regions(record, region=region_name)
		regions = [r for r in all_noncoding if r['end'] > noncoding_boundary[0] and r['start'] < noncoding_boundary[1]]
		for r in regions:
			r['genome_input_type'] = genome_input_type
			r['genome_id'] = default_genome_id
		start_pct = 0
		end_pct = 100

	if region_type == 'custom':
		if len(custom_regions_csv):
			if not Path(custom_regions_csv).exists():
				print("Must input a valid filepath for the custom regions csv")
				return
			with open(Path(custom_regions_csv), 'r', encoding='utf-8-sig') as csv_file:
				reader = csv.DictReader(csv_file)
				custom_regions = [r for r in reader]
			regions = [{'name': f'custom-{index}', 'start': int(c['start_ref']), 'end': int(c['end_ref']), 'direction': 'fw', 'genome_id': c['genome_id']} for (index, c) in enumerate(custom_regions)]
		elif len(custom_sequences) > 0:
			regions = [{'name': f'custom-{index}', 'start': 0, 'end': len(c), 'direction': 'fw', 'genome_id': None, 'sequence': c.upper()} for (index, c) in enumerate(custom_sequences)]
		else:
			print("No custom regions csv or custom sequences provided. Please fill in one of those fields")
			return
		start_pct = 0
		end_pct = 100

	for idx, region in enumerate(regions):
		start = time.perf_counter()
		print(f"Finding gRNA for \"{region['name']}\"")

		if 'sequence' in region and len(region['sequence']) >= 20:
			genome = Seq(region['sequence'])
			start_mark = 0
			end_mark = len(genome)
		else:
			[start_mark, end_mark] = get_target_region_for_gene(region, start_pct, end_pct)
			if genome_input_type == 'genbank_ids':
				record = retrieve_annotation(region['genome_id'], email)
			elif genome_input_type == 'genbank_files':
				record = SeqIO.read(Path(region['genome_id']), 'genbank')
			elif genome_input_type == 'fasta_files':
				record = SeqIO.read(Path(region['genome_id']), "fasta")
			genome = record.seq.upper()

		candidates = get_candidates_for_region(genome, start_mark, end_mark, region['name'], GC_requirement)
		print(f"Identified {len(candidates)} candidates by PAM and GC content, filtering them now...")
		if region_type == 'coding':
			if coding_spacer_direction not in ['N_to_C', 'C_to_N']:
				print("Invalid 'coding_spacer_direction' parameter, see the valid inputs")
				return
			candidates = order_candidates_for_region(candidates, region, coding_spacer_direction)
		
		candidates = filter_re_sites(candidates)
		candidates = filter_homopolymers(candidates)
		for genome_id in genome_ids:
			candidates = filter_non_unique_fingerprints(candidates, genome_id)

		print(f"{len(candidates)} remain after filtering for re_sites, homopolymers, and non unique fingerprints")
		for genome_id in genome_ids:
			candidates = remove_offtarget_matches(genome_id, region['name'], candidates, spacers_per_region, overlapping_spacers, check_all=len(genome_ids) > 1)
		region['candidates'] = candidates
		if len(candidates) > spacers_per_region:
			region['candidates'] = candidates[:spacers_per_region]

		elapsed_time = round(time.perf_counter() - start, 2)
		print(f"Identified {len(candidates)} spacers for {region['name']} in {elapsed_time} seconds")
		make_spacer_gen_output(region, os.path.join(output_path, 'spacer_gen_output.csv'))

def spacer_eval(args):
	genbank_ids = args['genbank_ids']
	output_path = args['output_path']
	email = args['email']
	user_spacers = args['spacers']

	if not email or '@' not in email:
		print("Please enter an email for NCBI API calls")
		return
	
	if type(user_spacers) is dict:
		spacers = [v for k,v in user_spacers.items() if len(v) == SPACER_LENGTH]
	else:
		spacers = [s for s in user_spacers if len(s) == SPACER_LENGTH]
	print(f"Starting evaluation for {len(spacers)} spacers...")
	spacer_batch = []
	spacer_batch_unmod = []  # make unmodified copy of spacer recs for output
	if flex_base:
		for (index, spacer) in enumerate(spacers):
			flexible_seq = spacer.upper()[:]
			for i in range(flex_spacing-1, SPACER_LENGTH, flex_spacing):
				flexible_seq = flexible_seq[:i] + 'N' + flexible_seq[i+1:]
			spacer_record = SeqRecord(Seq(flexible_seq), id=f'spacer_{index+1}', description=f'spacer_{index+1}')  # use 'description' to store ref_genome info
			spacer_batch.append(spacer_record)
			spacer_record_unmod = SeqRecord(Seq(spacer.upper()), id=f'spacer_{index + 1}',
									  description=f'spacer_{index+1}')  # unmodified spacer for output
			spacer_batch_unmod.append(spacer_record_unmod)
	else:
		for (index, spacer) in enumerate(spacers):
			spacer_record = SeqRecord(Seq(spacer), id=f'spacer_{index + 1}',
									  description=f'spacer_{index+1}')  # use 'description' to store ref_genome info
			spacer_batch.append(spacer_record)
			spacer_batch_unmod.append(spacer_record)

	# Write the batch of candidate sequences to a fasta for bowtie2 to use
	fasta_name = 'eval-spacers.fasta'
	root_dir = Path(__file__).parent.parent
	fasta_name = os.path.join(root_dir, 'assets', 'bowtie', fasta_name)
	with open(fasta_name, 'w') as targets_file:
		SeqIO.write(spacer_batch, targets_file, 'fasta')

	output_locations = []
	for genbank_id in genbank_ids:
		record = retrieve_annotation(genbank_id, email, return_record=False)
		output_locations.append(find_offtargets(genbank_id, fasta_name))

	make_eval_outputs(spacer_batch_unmod, output_locations, email, output_path, user_spacers)

	os.remove(fasta_name)
	for loc in output_locations:
		os.remove(loc)