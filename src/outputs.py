from Bio import SeqIO, SeqRecord, Seq
import subprocess
import sys
from simplesam import Reader as samReader
import tempfile
import os
import csv
import time
import multiprocessing
from collections import Counter
from pathlib import Path

def make_spacer_gen_output(regions, output_filename):
	fieldnames = ['spacer_id', 'region', 'sequence', 'genomic_coordinate', 'GC_content', 'PAM', 'strand']
	gc_counter = Counter(['G', 'C'])

	spacer_id = 0
	with open(f'{output_filename}.csv', 'w', newline='') as out_file:
		writer = csv.DictWriter(out_file, fieldnames=fieldnames)
		writer.writeheader()
		for region in regions:
			if 'candidates' in region:
				for c in region['candidates']:
					seq = c['seqrec'].seq
					counter = Counter(seq)
					GC_content = (counter['G'] + counter['C'])/len(seq)

					row = {
						'spacer_id': spacer_id,
						'region': region['name'],
						'sequence': c['seqrec'].seq, 
						'genomic_coordinate': c['location'],
						'GC_content': GC_content,
						'PAM': 'CC',
						'strand': c['name'].split('--')[1]
					}
					writer.writerow(row)
					spacer_id += 1
	return


def make_eval_outputs(spacers, output_sam, genome, output_path):
	genome_seq = genome.upper()

	with open(output_sam) as sam_file:
		reader = samReader(sam_file)
		sam_reads = [r for r in reader]

	spacer_output = []
	for spacer in spacers:
		reads = [r for r in sam_reads if r.safename == spacer.id]
		mismatch_list = []  # to get min/max mismatches
		if len(reads) == 1:
			print(f"Spacer '{spacer.id}' has no potential off-targets")
			spacer_loc = reads[0].coords[0] if reads[0].reverse else reads[0].coords[-1]

		elif len(reads) > 1:
			print(f"Spacer '{spacer.id}' has {len(reads)} potential off-targets - see output files for details")
			with open(os.path.join(output_path, f'crrna_{spacer.id}_off_target.txt'), 'w') as log:
				log.write(f"Potential off-target sites for {spacer.id}")
				log.write("\n(Top row is spacer sequence, following rows are potential off-targets; closest matches listed first)")
				log.write("\n-------------")
				proto_list = []
				for i in reads:
					if i.reverse:
						pam = genome_seq[i.coords[-1]:i.coords[-1] + 2].upper()
						pam = pam.reverse_complement()
						protospacer = genome_seq[i.coords[0] - 1:i.coords[-1]].upper()
						protospacer = protospacer.reverse_complement()
					else:
						pam = genome_seq[i.coords[0] - 3:i.coords[0] - 1]
						protospacer = genome_seq[i.coords[0] - 1:i.coords[-1]]
					proto_list.append(protospacer)
					log.write(f"\n{protospacer.upper()} "
							  f"(Coordinates: {i.coords[0]}; PAM: {pam.upper()}; RevCom = {i.reverse})")
				for i in reads[1:]:
					mismatch_count = abs(i.tags['AS'])//6  # use AS instead of NM since NM counts wobble bases
					mismatch_list.append(mismatch_count)
				log.write("\n-------------")
				log.write("\nAlignment: (ambiguous base mismatches shown as lower case")
				log.write((f"\n{proto_list[0]}"))
				for target in range(1, len(proto_list)):
					log.write("\n")
					pos = 0
					for a, b in zip(proto_list[0], proto_list[target]):
						pos += 1
						if a == b:
							log.write("-")
						else:
							if pos%6 == 0:  # check if ambiguous base
								log.write(str(b).lower())
							else:
								log.write(str(b).upper())
				log.write(f"\nTotal potential off-targets = {len(reads) - 1}")
			spacer_loc = reads[0].coords[0] if reads[0].reverse else reads[0].coords[-1]

		else:
			print(f"Warning - no protospacer found for spacer '{spacer.id}'")
			spacer_loc = 'No protospacer found'

		spacer_dict = {'name': spacer.id, 'sequence': spacer.seq, 'refseq': spacer.description, 'location': spacer_loc,
					   'offtar_count': len(reads)-1 if len(reads)>0 else 'N/A',
					   'mismatch_min': min(mismatch_list) if len(mismatch_list)>0 else 'N/A',
					   'mismatch_max': max(mismatch_list) if len(mismatch_list)>0 else 'N/A'}
		spacer_output.append(spacer_dict)

	print('------------')

	fieldnames = ['Spacer Name', 'Spacer Sequence', 'Reference Genome',
				  'Protospacer Location', 'Number of Off-targets', 'Fewest Mismatches', 'Most Mismatches']
	with open(os.path.join(Path(output_path), 'spacer_eval_output.csv'), 'w', newline='') as out_file:
		writer = csv.DictWriter(out_file, fieldnames=fieldnames)
		writer.writeheader()
		for spacer in spacer_output:
			row = {'Spacer Name': spacer['name'], 'Spacer Sequence': spacer['sequence'],
				   'Reference Genome': spacer['refseq'], 'Protospacer Location': spacer['location'],
				   'Number of Off-targets': spacer['offtar_count'],
				   'Fewest Mismatches': spacer['mismatch_min'], 'Most Mismatches': spacer['mismatch_max']}
			writer.writerow(row)
	return

def spacer_eval():
	make_eval_output_files(spacers, output_name)

	os.remove(fasta_name)
	os.remove(output_name)
