from src.main import spacer_gen

# Folder to put outputs in
# Can be relative or absolute
# Ex. 'C:\\Users\\Me\\Downloads', '../../Downloads'
outputPath = ''

# Genbank Accession number
# Ex. "CP001509.3"
genbankId = ""

# Email: only used for NCBI API calls, be polite!
email = ''

# The number of spacers to generate per region
spacersPerRegion = 5

# ---------- REGION TYPES - CODING, NONCODING, OR CUSTOM --------------#
# regionType must be 'coding', 'noncoding', or 'custom'
# see the specific settings for each below
regionType = ''

# --------  CODING REGION SETTINGS ------------#
# Percent of each region to target, from N terminus to C terminus
startPct = 10
endPct = 50

# Target genes: gene names as known from the genbank accession number, referenced above
# Ex. ['ThrA', 'ThrB', 'htgA']
target_gene_names = []


# -------- NONCODING REGION SETTINGS ------------#
# genome_boundary takes a single pair of coordinates, and will only target
# intergenic regions that lie between the two values
# Ex. [13555, 20000]
noncoding_boundary = []


# -------- CUSTOM REGION SETTINGS ------------#
# A list of coordinate pairs in the reference sequence to generate spacers for, 
# as a pair
# Ex. [[13555, 20000], [25000,26000]]
custom_regions = []


# Do not modify, this calls the function when run with 'python spacer_gen.py'
if __name__ == "__main__":
	spacer_gen({"outputPath": outputPath, "spacersPerRegion": spacersPerRegion, "genbankId": genbankId, "email": email, "regionType": regionType, "startPct": startPct, "endPct": endPct, "target_gene_names": target_gene_names, "noncoding_boundary": noncoding_boundary, "custom_regions": custom_regions})
