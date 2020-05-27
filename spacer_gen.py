from src.main import spacer_gen

# Folder to put outputs in
# Can be relative or absolute
# Ex. 'C:\\Users\\Me\\Downloads', '../../Downloads'
output_path = './'

# Genbank Accession number
# The downloaded genbank files are cached locally in the assets folder, and will only
# be downloaded from NCBI if not present
# Ex. "CP001509.3"
genbank_id = "CP001509.3"

# Email: only used for NCBI API calls, be polite!
email = 'acreechristopher@gmail.com'

# The number of spacers to generate per region
spacers_per_region = 5

# GC Content Filter - as a min and max percent of GC content allowed
# Ex. GC_requirement = [40,60] ensures only spacers with GC content of 40% - 60% will be generated
GC_requirement = [46,54]

# Preferences on overlapping spacers - can be "allowed", "avoid", "forbidden"
# "allowed" will not manipulate ordering apart from other settings
# "avoid" will permit overlapping spacers only if sufficient non-overlapping spacers cannot be generated
# "forbidden" will not return overlapping spacers, even if it means the specified number per region cannot be generated
overlapping_spacers = "avoid"

# -------------- REGION TYPES - CODING, NONCODING, OR CUSTOM --------------#
# region_type must be 'coding', 'noncoding', or 'custom'
# see the specific settings for each below
region_type = 'coding'

# -------------- CODING REGION SETTINGS -----------------------------------#
# Percent of each region to target, from N terminus to C terminus
start_pct = 10
end_pct = 50

# How to prioritize choosing spacers in each region:
# Must be either "N_to_C" search or "C_to_N" search
coding_spacer_direction = "N_to_C"

# Target genes: gene names as known from the genbank accession number, referenced above
# Ex. ['ThrA', 'ThrB', 'htgA']
target_gene_names = ['ThrA']


# --------------- NONCODING REGION SETTINGS -------------------------------#
# genome_boundary takes a single pair of coordinates, and will only target
# intergenic regions that lie between the two values
# Ex. [13555, 20000]
noncoding_boundary = []


# --------------- CUSTOM REGION SETTINGS ----------------------------------#
# A list of coordinate pairs in the reference sequence to generate spacers for, 
# as a pair
# Ex. [[13555, 20000], [25000,26000]]
custom_regions = []


# Do not modify, this calls the function when run with 'python spacer_gen.py'
if __name__ == "__main__":
	spacer_gen({"output_path": output_path, "spacers_per_region": spacers_per_region, "GC_requirement": GC_requirement, "genbank_id": genbank_id, "email": email, "region_type": region_type, "overlapping_spacers": overlapping_spacers, "coding_spacer_direction": coding_spacer_direction, "start_pct": start_pct, "end_pct": end_pct, "target_gene_names": target_gene_names, "noncoding_boundary": noncoding_boundary, "custom_regions": custom_regions})
