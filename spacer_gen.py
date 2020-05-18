from src.main import spacer_gen
parameters = {
	# Folder to put outputs in
	# Can be relative or absolute
	# Ex. 'C:\\Users\\Me\\Downloads', '../../Downloads'
	outputPath = ''

	# Genbank Accession number
	# Ex. "CP001509.3"
	genbankId = ''

	# Email: only used for NCBI API calls, be polite!
	email = ''

	# Percent of each region to target, from N terminus to C terminus
	startPct = 10
	endPct = 50

	# The number of spacers to generate per region
	spacersPerRegion = 5

	# Target regions
	# Ex. ['ThrA', 'ThrB', 'htgA']
	target_gene_names = []
}


if __name__ == "main":
	spacer_gen(parameters)
