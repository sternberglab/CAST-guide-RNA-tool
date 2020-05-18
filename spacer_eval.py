from src.main import spacer_eval
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
	
	# List of spacers to evaluate for off-target matches
	# Each must be 32 basepairs in length
	spacers = ['']

}

if __name__ == "main":
	spacer_eval(parameters)
