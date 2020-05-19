import sys
import json
import os
from pathlib import Path

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from src.bowtie import build
from src.advanced_parameters import minimum_intergenic_region_length

def get_from_cache(genbankId):
	root_dir = Path(__file__).parent.parent
	genbank_assets_path = os.path.join(root_dir, 'assets', 'genbank')
	os.makedirs(genbank_assets_path, exist_ok=True)

	genbank_file = os.path.join(genbank_assets_path, f'{genbankId}.json')
	if Path(genbank_file).exists():
		with open(genbank_file, 'r') as f:
			return json.load(f)
	return False

def save_to_cache(genbankId, genbank_info):
	root_dir = Path(__file__).parent.parent
	genbank_assets_path = os.path.join(root_dir, 'assets', 'genbank')
	os.makedirs(genbank_assets_path, exist_ok=True)

	genbank_file = os.path.join(genbank_assets_path,  f'{genbankId}.json')
	with open(genbank_file, 'w') as f:
		json.dump(genbank_info, f)
	genbank_fasta = os.path.join(genbank_assets_path,  f'{genbankId}.fasta')
	with open(genbank_fasta, 'w') as f:
		seqrec = SeqRecord(Seq(genbank_info['GBSeq_sequence']), id=genbankId, name=genbankId, description=genbankId)
		SeqIO.write(seqrec, f, 'fasta')
	build(genbankId)

def retrieve_annotation(genbankId, email):
	# *Always* tell NCBI who you are
	Entrez.email = email

	cached = get_from_cache(genbankId)
	if cached:
		save_to_cache(genbankId, cached)
		return cached	
	"""
	Annotates Entrez Gene IDs using Bio.Entrez, in particular epost (to
	submit the data to NCBI) and esummary to retrieve the information.
	Returns a list of dictionaries with the annotations.
	"""
	print(f"Fetching genbank info for {genbankId}")
	handle = Entrez.efetch("nucleotide", id=genbankId, retmode="xml")
	try:
		result = Entrez.read(handle)[0]
	except RuntimeError as e:
		print(e)
		sys.exit(-1)
	save_to_cache(genbankId, result)
	return result

def basic_gene_info(gene):
	# Returns a structure in which direction is explicitly stated, and
	# start is always before end in the refseq
	genbank_start = int(gene['GBFeature_intervals'][0]['GBInterval_from'])
	genbank_end = int(gene['GBFeature_intervals'][0]['GBInterval_to'])
	direction = 'fw' if genbank_start < genbank_end else 'rv'
	try:
		name = next(i for i in gene['GBFeature_quals'] if i['GBQualifier_name'] == 'gene')['GBQualifier_value']
	except:
		name = f'Unknown_Name_{min(genbank_start, genbank_end)}'
	return {
		"name": name,
		"start": min(genbank_start, genbank_end),
		"end":  max(genbank_start, genbank_end),
		"direction": direction
	}

def get_noncoding_regions_from_genes(genes, genome):
	genome_end = len(genome)
	noncoding_regions = []
	
	gene_index = 0
	prev_gene = {'end': -1}
	while gene_index < len(genes):
		next_gene = genes[gene_index]
		if (int(next_gene['start']) - int(prev_gene['end'])) > minimum_intergenic_region_length:
			new_region = {
				'name': 'noncoding_'+str(len(noncoding_regions)+1),
				'start': int(prev_gene['end'])+1,
				'end': int(next_gene['start'])-1,
				'direction': 'fw'
			}
			noncoding_regions.append(new_region)
		prev_gene = next_gene
		gene_index += 1
	return noncoding_regions

def get_regions(annotation_result, coding=True):
	genes = [feat for feat in annotation_result['GBSeq_feature-table'] if feat['GBFeature_key'] == 'gene']
	genes = [basic_gene_info(g) for g in genes]
	if coding is False:
		genome = annotation_result['GBSeq_sequence']
		noncoding_regions = get_noncoding_regions_from_genes(genes, genome)
		return noncoding_regions
	return genes