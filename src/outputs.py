import os
import csv
from collections import Counter
from pathlib import Path
from Bio.pairwise2 import format_alignment as f_align

from Bio import SeqIO, SeqRecord, Seq, pairwise2
from simplesam import Reader as samReader

from src.advanced_parameters import PAM_SEQ, SPACER_LENGTH, flex_base, flex_spacing

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
						'region': f"{region['name']} - ({region['start']}-{region['end']})",
						'sequence': c['seqrec'].seq, 
						'genomic_coordinate': c['location'],
						'GC_content': GC_content,
						'PAM': PAM_SEQ,
						'strand': c['name'].split('--')[1][:2]
					}
					writer.writerow(row)
					spacer_id += 1
	return

def align_offtargets(spacer_seq, offtar_info):  # uses Bio.pairwise2, offtar_info is a list of 2 things
	if not offtar_info[1]:  # ungapped alignment, penalize gaps more as it usually looks better
		alignments = pairwise2.align.globalms(spacer_seq, offtar_info[0], 2, -1, -5, -.1) # opening gaps = at least 5 mismatches
	else:  # probably has gaps that need to be shown
		alignments = pairwise2.align.globalxx(spacer_seq, offtar_info[0])
	align_out = f_align(*alignments[0])
	items = align_out.split('\n')
	del items[-1]  # manually remove the 'score=XX' line
	del items[-1]
	return '\n'.join(items)+'\n'

def make_eval_outputs(spacers, output_sam, genome, output_path):
	genome_seq = genome.upper()

	with open(output_sam) as sam_file:
		reader = samReader(sam_file)
		sam_reads = [r for r in reader]

	spacer_output = []
	for spacer in spacers:
		reads = [r for r in sam_reads if r.safename == spacer.id]
		mismatch_list = []  # to get min/max mismatches
		# if len(reads) == 1:
		#     print(f"Spacer '{spacer.id}' has no potential off-targets")
		#     spacer_match = reads[0].coords[0] if reads[0].reverse else reads[0].coords[-1]
		spacer_match = 'No perfect match found'  # update if perfect match (true protospacer) is found
		if len(reads) >= 1:
			print(f"Spacer '{spacer.id}' has {len(reads)} potential match(es) - see output files for details")
			with open(os.path.join(output_path, f'{spacer.id}_off_target.txt'), 'w') as text_out:
				text_out.write(f"Potential off-target sites for {spacer.id}")
				text_out.write("\n(Closer matches to protospacer listed first)")
				text_out.write("\n-------------")
				proto_list = []
				perfect_match = False
				for i in reads:
					if i.reverse:
						pam = genome_seq[i.coords[-1]:i.coords[-1] + 2].upper()  # todo allow different PAM lengths
						pam = pam.reverse_complement()
						protospacer = genome_seq[i.coords[0] - 1:i.coords[-1]].upper()
						protospacer = protospacer.reverse_complement()
					else:
						pam = genome_seq[i.coords[0] - 3:i.coords[0] - 1]
						protospacer = genome_seq[i.coords[0] - 1:i.coords[-1]]
					gapped_align = i.tags['XO'] > 0
					proto_list.append([protospacer, gapped_align])
					if protospacer == spacer.seq:
						perfect_match = True
					flex_count = len(spacer.seq)//flex_spacing if flex_base else 0
					mismatch_count = abs(i.tags['XM']) - flex_count if len(
						protospacer) >= len(spacer.seq) else 'N/A'  # use XM, taking into account flexible bases
					if perfect_match:
						spacer_match = 'Perfect match(es) found'
					text_out.write(f"\n{protospacer.upper()} "
								   f"  Coordinates: {i.coords[0]}; PAM: {pam.upper()}; RevCom = {i.reverse}; "
								   f"Perfect Match = {str(perfect_match)}; Mismatches = {mismatch_count}; Gapped Alignment = {gapped_align}")
					if mismatch_count != 'N/A':
						mismatch_list.append(mismatch_count)
				text_out.write("\n-------------")
				text_out.write(
					"\nPairwise Alignments: (ruler for ambiguous base positions shown above each alignment pair for ungapped alignments)\n")
				# text_out.write('\n')
				for target in proto_list:
					if not target[1]:  # ungapped; print out a ruler
						text_out.write('\n')
						for _ in range(0, len(spacer.seq)):
							if (_ + 1) % 6 == 0:
								text_out.write('X')
							else:
								text_out.write('-')
					text_out.write('\n' + align_offtargets(spacer.seq, target))

				text_out.write(f"\nTotal potential matches = {len(reads)}")
				if not perfect_match:
					text_out.write('\nNo identical protospacer found')
				else:
					text_out.write('\nIdentical protospacer(s) found')

		else:
			print(f"Warning - no matches, including protospacer, found for spacer '{spacer.id}'")

		spacer_dict = {'name': spacer.id, 'sequence': spacer.seq, 'refseq': spacer.description,
					   'match_found': spacer_match,
					   'offtar_count': len(reads),
					   'mismatch_min': min(mismatch_list) if len(mismatch_list) > 0 else 'N/A',
					   'mismatch_max': max(mismatch_list) if len(mismatch_list) > 0 else 'N/A'}
		spacer_output.append(spacer_dict)

	print('------------')

	fieldnames = ['Spacer Name', 'Spacer Sequence', 'Reference Genome',
				  'Perfect Match', 'Number of Matches', 'Fewest Mismatches', 'Most Mismatches']
	with open(os.path.join(Path(output_path), 'spacer_eval_output.csv'), 'w', newline='') as out_file:
		writer = csv.DictWriter(out_file, fieldnames=fieldnames)
		writer.writeheader()
		for spacer in spacer_output:
			row = {'Spacer Name': spacer['name'], 'Spacer Sequence': spacer['sequence'],
				   'Reference Genome': spacer['refseq'], 'Perfect Match': spacer['match_found'],
				   'Number of Matches': spacer['offtar_count'],
				   'Fewest Mismatches': spacer['mismatch_min'], 'Most Mismatches': spacer['mismatch_max']}
			writer.writerow(row)
	return

