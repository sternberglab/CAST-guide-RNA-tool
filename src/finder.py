import subprocess
import sys
import tempfile
import os
import csv
import time
import multiprocessing
from collections import Counter
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from simplesam import Reader as samReader

from src.bowtie import find_offtargets
from src.advanced_parameters import PAM_SEQ, INTEGRATION_SITE_DISTANCE, SPACER_LENGTH

def candidates_for_seq(seq, descriptor):
	candidates = []
	i=0
	while i < len(seq):
		nextPAM = seq[i:].find(PAM_SEQ)
		if nextPAM == -1 or (i + nextPAM + 34) > len(seq):
			i += 100000
			break
			
		targetSeq = seq[i+nextPAM+2:i+nextPAM+34]
		name = descriptor + str(i+nextPAM+2)
		target = SeqRecord(targetSeq, id=name, name=name, description=name)
		candidate = {'name': target.id, 'seqrec': target, 'location': i+nextPAM+2}
		candidates.append(candidate)
		i += nextPAM + 1
	return candidates
	
def get_target_region_for_gene(gene, startPct, endPct):
	# Calculate the target region of a gene by pct of the coding region, from N to C
	gene_length = gene['end'] - gene['start']
	
	fifty_pct_mark = gene['start'] + int(gene_length/2)
	if gene['direction'] == 'fw':
		start_mark = gene['start'] + int(gene_length*startPct/100)
		end_mark = gene['start'] + int(gene_length*endPct/100)
	else:
		start_mark = gene['start'] + int(gene_length*(100-endPct)/100)
		end_mark = gene['start'] + int(gene_length*(100-startPct)/100)
	return [start_mark, end_mark]
	
def get_candidates_for_region(genome, start_mark, end_mark, name):
	genome_seq = genome.upper()

	search_offset = len(PAM_SEQ) + SPACER_LENGTH + INTEGRATION_SITE_DISTANCE
	
	fw_search_seq = genome_seq[start_mark-search_offset:end_mark-search_offset]
	rv_search_seq = genome_seq[start_mark+search_offset:end_mark+search_offset].reverse_complement()
	
	candidates = candidates_for_seq(fw_search_seq, name+'--fw')
	for c in candidates:
		# the initial "location" here is the location in the search sequence, needs to be 
		# placed in the genome location
		c['location'] = start_mark - search_offset + c['location']
	rv_candidates = candidates_for_seq(rv_search_seq, name+'--rv')
	for c in rv_candidates:
		# move back additional 32 (spacer length) for the reverse oriented spacers
		c['location'] = end_mark + search_offset - c['location'] - 32
	candidates.extend(rv_candidates)
	
	for candidate in candidates:
		if 'fw' in candidate['name']:
			# for fw strand inserts, the fingerprint is downstream
			fp_start = candidate['location'] + SPACER_LENGTH + INTEGRATION_SITE_DISTANCE - 20
			fp_end = candidate['location'] + SPACER_LENGTH + INTEGRATION_SITE_DISTANCE
		else:
			# for rv strand inserts, the fingerprint is upstream, but location is the END of the spacer
			fp_start = candidate['location'] - INTEGRATION_SITE_DISTANCE - 20
			fp_end = candidate['location'] - INTEGRATION_SITE_DISTANCE
		name = candidate['name']
		candidate['fp_seq'] = SeqRecord(genome_seq[fp_start:fp_end], id=name, name=name, description=name)
	return candidates

def filter_non_unique_fingerprints(candidates):
	temp_fasta_name = 'fp_check.fasta'
	with open(temp_fasta_name, 'w') as temp_fasta:
		SeqIO.write([c['fp_seq'] for c in candidates], temp_fasta, 'fasta')
	cores = multiprocessing.cpu_count()
	output_name = 'fp_check_out.sam'
	align_command = 'bowtie2 -x {} -f {} -p {} -S {}'.format(GENOME_NAME, temp_fasta_name, cores-1, output_name)
	subprocess.run(align_command, shell=True)
	
	filtered_candidates = []
	sam_reads = []
	with open(output_name, 'r') as sam_file:
		reader = samReader(sam_file)
		sam_reads = [r for r in reader]
	
	for c in candidates:
		reads = [r for r in sam_reads if r.safename == c['name']]
		if len(reads) == 1:
			filtered_candidates.append(c)
	os.remove(temp_fasta_name)
	os.remove(output_name)
	return filtered_candidates


def remove_offtarget_matches(genbankId, name, candidates, minMatches):
	print("Finding potential off-target matches...", flush=True)
	
	no_offtargets = []
	# go through 10 at a time since this is slow, stop once we have the specified # of spacers
	for i in range(0, len(candidates), 10):
		test_candidates = candidates[i:i+10]
		
		# Get the candidate sequences to use, and
		# make every 6th bp an N to allow for ambiguous matches
		match_seqs = [c['seqrec'].upper() for c in test_candidates]
		for seq in match_seqs:
			flexible_seq = seq.seq[:]
			for i in range(5, 32, 6):
				flexible_seq = flexible_seq[:i] + 'N' + flexible_seq[i+1:]
			seq.seq = flexible_seq

		# Write the batch of candidate sequences to a fasta for bowtie2 to use
		fasta_name = name+'-candidates.fasta'
		root_dir = Path(__file__).parent.parent
		fasta_name = os.path.join(root_dir, 'assets', 'bowtie', genbankId, fasta_name)
		with open(fasta_name, 'w') as targets_file:
			SeqIO.write(match_seqs, targets_file, 'fasta')

		output_location = find_offtargets(genbankId, fasta_name)

		filtered_candidates = []
		sam_reads = []
		with open(output_location) as sam_file:
			reader = samReader(sam_file)
			sam_reads = [r for r in reader]

		for c in candidates:
			reads = [r for r in sam_reads if r.safename == c['name']]
			if len(reads) == 1:
				no_offtargets.append(c)
		os.remove(fasta_name)
		os.remove(output_location)
		# return once at least minMatches are found without off-targets
		if len(no_offtargets) >= minMatches:
			return no_offtargets
	return no_offtargets
