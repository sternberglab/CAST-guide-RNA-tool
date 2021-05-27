import sys
import subprocess
import multiprocessing
import os
from pathlib import Path

from src.advanced_parameters import mismatch_threshold, allow_gaps

def build(genbank_id, fasta_file=None):
	root_dir = Path(__file__).parent.parent
	bowtie_genome_dir = os.path.join(root_dir, 'assets', 'bowtie', genbank_id)
	os.makedirs(bowtie_genome_dir, exist_ok=True)

	build_output_name = os.path.join(bowtie_genome_dir, 'index')
	# if index exists no need to rebuild
	if Path(f'{build_output_name}.1.bt2').exists():
		return
	
	if not fasta_file:
		fasta_build_file = os.path.join(root_dir, 'assets', 'genbank', f'{genbank_id}.fasta')
	else:
		fasta_build_file = fasta_file
	build_command = f'bowtie2-build {fasta_build_file} {build_output_name}'
	try:
		subprocess.run(build_command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
	except Exception as e:
		raise Exception(f"Error running bowtie, see Installation notes in the README: \n{e}")

def find_offtargets(genbank_id, fasta_name):
	root_dir = Path(__file__).parent.parent
	index_location = os.path.join(root_dir, 'assets', 'bowtie', genbank_id, 'index')
	# Run the bowtie2 alignment command
	# -x {} : the name of the genome index file (already built by bowtie2-build)
	# -a : return all results, not just highest match
	# -f : read is a fasta file
	# -t : include time in the command line output
	# -p {} : use multiple cores for faster processing
	# -S {} : output reads to SAM file at this path
	# --no-1mm-upfront : "no 1 mismatch upfront", normally bowtie2 will return if it finds an exact or 1 mismatch read before doing a deeper search
	# --np 0 --n-ceil 5 : No penalty for ambiguous base pairs in the sequence ("N"), and allow up to 5
	# --score-min : scoring equation, for us just flat min score of (-6 per mismatch X mismatch_threshold) or greater is a pass
	# -N 0 -L 5 -i S,6,0 -D 10: Sets the multiseed alignment rules. 0 mismatches allowed, seed length 5, and skip 6bp between seeds. Search up to 10 seeds before failing
	# (This means there must be one fully matching 5bp sequence between a set of ambiguous characters to pass the seed filtering)
	# --rdg XX,1 : read gap-open penalty of XX and gap-extension penalty of 1. Set XX to scale with mismatch_threshold
	# --rfg XX,1 : reference gap-open penalty of 50 and gap-extension penalty of 1
	gap_option = f'--rdg {mismatch_threshold*100},1 --rfg {mismatch_threshold*100},1' if not allow_gaps else ''

	cores = multiprocessing.cpu_count()
	output_name = f"{fasta_name.split('.')[0]}-{genbank_id}-offtarget.sam"
	output_location = os.path.join(root_dir, 'assets', 'bowtie', genbank_id, output_name)

	root_dir = Path(__file__).parent.parent
	index_location = os.path.join(root_dir, 'assets', 'bowtie', genbank_id, 'index')
	align_command = f'bowtie2 -x {index_location} -a -f {fasta_name} -t -p {cores - 1} {gap_option} -S {output_location} --no-1mm-upfront --np 0 --n-ceil 5 --score-min L,-{6*mismatch_threshold+1},0 -N 1 -L 11 -i S,6,0 -D 6 --no-unal'
	try:
		subprocess.run(align_command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
	except Exception as e:
		raise Exception(f'Error running bowtie, see Installation notes in the README: \n{e}')
	return output_location