from src.main import spacer_eval

# Folder to put outputs in
# Can be relative or absolute
# Ex. 'C:\\Users\\Me\\Downloads', '../../Downloads'
outputPath = './'

# Genbank Accession number
# Ex. "CP001509.3"
genbankId = ""

# Email: only used for NCBI API calls, be polite!
email = ''

# List of spacers to evaluate for off-target matches
# Each must be 32 basepairs in length, invalid lengths will be filtered out
spacers = ['AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA']







# Do not modify, this calls the function when run with 'python spacer_eval.py'
if __name__ == "__main__":
	spacer_eval({"outputPath": outputPath, "genbankId": genbankId, "email": email, "spacers": spacers})
