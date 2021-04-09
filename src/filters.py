from Bio import Restriction, Seq, SeqIO
import multiprocessing
import subprocess
import os
from pathlib import Path
from simplesam import Reader as samReader
from src.advanced_parameters import restriction_enzymes, homopolymer_length

def filter_re_sites(candidates):
	rb = Restriction.RestrictionBatch(restriction_enzymes)
	filtered_candidates = []
	for c in candidates:
		rbsearch = rb.search(c['seqrec'].seq)
		matched = any([match for re in rbsearch.keys() for match in rbsearch[re]])
		if not matched:
			filtered_candidates.append(c)
	return filtered_candidates

def filter_non_unique_fingerprints(candidates, genbank_id):
	root_dir = Path(__file__).parent.parent
	temp_fasta_name = 'fp_check.fasta'
	with open(temp_fasta_name, 'w') as temp_fasta:
		SeqIO.write([c['fp_seq'] for c in candidates], temp_fasta, 'fasta')
	cores = multiprocessing.cpu_count()
	output_name = 'fp_check_out.sam'
	index_location = os.path.join(root_dir, 'assets', 'bowtie', genbank_id, 'index')
	align_command = 'bowtie2 -x {} -k 2 -f {} -p {} -S {}'.format(index_location, temp_fasta_name, cores-1, output_name)
	subprocess.run(align_command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

	filtered_candidates = []
	sam_reads = {}
	with open(output_name, 'r') as sam_file:
		reader = samReader(sam_file)
		for r in reader:
			if r.safename not in sam_reads:
				sam_reads[r.safename] = 1
			else:
				sam_reads[r.safename] += 1


	for c in candidates:
		if sam_reads[c['fp_seq'].id] == 1:
			filtered_candidates.append(c)
	os.remove(temp_fasta_name)
	os.remove(output_name)
	return filtered_candidates

def filter_homopolymers(candidates, G_only=True):
	filtered_candidates = []
	bases_to_avoid = ['G']
	if not G_only:
		bases_to_avoid = ['G', 'A', 'T', 'C']
	for c in candidates:
		should_filter = False
		seq = c['seqrec'].seq.upper()
		for base in bases_to_avoid:
			homopolymer = base*homopolymer_length
			if seq.find(homopolymer) != -1:
				should_filter = True
		if not should_filter:
			filtered_candidates.append(c)

	return filtered_candidates
