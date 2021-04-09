import random
import time
import csv
import os
from pathlib import Path
from simplesam import Reader as samReader
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from src.filters import filter_homopolymers, filter_re_sites
from src.bowtie import find_offtargets
from src.genbank import retrieve_annotation
from src.advanced_parameters import SPACER_LENGTH, flex_spacing

barcode_length = 12
min_hamming_distance = 4

def hamming_dist(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def remove_offtargets(spacers, genbank_id):
	flex_seqs = []
	for (index, spacer) in enumerate(spacers):
		flexible_seq = spacer["seqrec"].seq[:]
		for i in range(flex_spacing-1, SPACER_LENGTH, flex_spacing):
			flexible_seq = flexible_seq[:i] + 'N' + flexible_seq[i+1:]
		spacer_record = SeqRecord(flexible_seq, id=f'spacer_{index+1}', description=genbank_id)  # use 'description' to store ref_genome info
		flex_seqs.append(spacer_record)

	fasta_name = 'offtarget-check.fasta'
	root_dir = Path(__file__).parent.parent
	fasta_name = os.path.join(root_dir, 'assets', fasta_name)
	with open(fasta_name, 'w') as targets_file:
		SeqIO.write(flex_seqs, targets_file, 'fasta')

	output_location = find_offtargets(genbank_id, fasta_name)
	
	filtered_spacers = []
	sam_reads = []
	with open(output_location, 'r') as sam_file:
		reader = samReader(sam_file)
		sam_reads = [r for r in reader]

	for i, s in enumerate(spacers):
		reads = [r for r in sam_reads if r.safename == f'spacer_{i+1}']
		reads = [r for r in reads if r.mapped]
		if len(reads) == 0:
			filtered_spacers.append(s)
	os.remove(fasta_name)
	os.remove(output_location)
	return filtered_spacers

def make_barcodes():
	barcodes = []
	start = time.perf_counter()
	i = 0
	while len(barcodes) < 35000:
		if i % 5000 == 0:
			print(f"{len(barcodes)} barcodes created out of {i} potential sequences. {round(time.perf_counter() - start, 2)} seconds")
		i+=1
		potential = ''.join(random.choices('ATGC', k=barcode_length))
		can_obj = {' seq': SeqRecord(Seq(potential))}
		if len(filter_homopolymers([can_obj], G_only=False)) != 1:
			continue
		if len(filter_re_sites([can_obj])) != 1:
			continue
		if potential in barcodes:
			continue
		for crosscheck in barcodes:
			if hamming_dist(potential, crosscheck) < min_hamming_distance:
				continue

		barcodes.append(potential)
	print("done generating, writing to csv")
	with open('barcodes_output.csv', 'w', newline='') as csvf:
		writer = csv.writer(csvf)
		writer.writerow(['id', 'barcode_sequence'])
		for (index, bc) in enumerate(barcodes):
			writer.writerow([index+1, bc])
	print("DONE")


def make_nontargeting(genbank_id, number=100):
	genome = retrieve_annotation(genbank_id, 'acreechristopher@gmail.com')

	nt_spacers = []
	while len(nt_spacers) < number:
		needed = number - len(nt_spacers)
		potentials = []
		while len(potentials) < needed:
			seq = ''.join(random.choices('ATGC', k=SPACER_LENGTH))
			potential = {'seqrec': SeqRecord(Seq(seq), id=seq)}
			if len(filter_homopolymers([potential], G_only=False)) != 1:
				continue
			if len(filter_re_sites([potential])) != 1:
				continue
			potentials.append(potential)
		print(f"Screening {len(potentials)} potentials")
		potentials = remove_offtargets(potentials, genbank_id)
		nt_spacers += potentials
		print(f"{len(potentials)} potentials passed a screening round, now have {len(nt_spacers)} total")

	with open('NT_spacers.csv', 'w', newline='') as csvf:
		writer = csv.writer(csvf)
		writer.writerow('spacer_seq')
		for nt in nt_spacers:
			writer.writerow([nt['seqrec'].seq.upper()])