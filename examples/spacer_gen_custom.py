from src.main import spacer_gen

# Folder to put outputs in
# Can be relative or absolute
# Ex. 'C:\\Users\\Me\\Downloads', '../../Downloads'
# Note that outputs to the 'spacer_gen_output.csv' will be APPENDED to existing entries
# in that file, not overridden, in the case of repeat runs. 
output_path = './'

# The number of spacers to generate per region
# The script will stop searching early if possible once it reaches this number, 
# and only output this many at the end, even if more were identified
# Make this very high to see all matching spacers
spacers_per_region = 100

# GC Content Filter - as a min and max percent of GC content allowed
# Ex. GC_requirement = [40,60] ensures only spacers with GC content of 40% - 60% will be generated
GC_requirement = [35,65]

# Preferences on overlapping spacers - one of "allowed", "avoid", "forbidden"
# "allowed" will not manipulate ordering apart from other settings
# "avoid" will permit overlapping spacers only if sufficient non-overlapping spacers cannot be generated
# "forbidden" will not return overlapping spacers, even if it means the specified number per region cannot be generated
overlapping_spacers = 'avoid'

# --------------- GENOME INPUT --------------------------------------------#
# Enter ONE of the following:
# 1. genbank_id and email, to download and use cached genbank files from NCBI
# 2. genbank_file to use a locally saved genbank file
# 3. genome_fasta_file to use a locally saved fasta file with the target genome (only for 'custom' region type)

# Genbank Accession number
# The downloaded genbank files are cached locally in the assets folder, and will only
# be downloaded from NCBI if not present
#
# Currently the coding and noncoding modes will only use FIRST file
# of genbank id listed to identify coding and noncoding regions to generate
# spacers for, though all files and genbank ids will be checked for 
# off-target activity and bad ones filtered out
# Ex. ["CP001509.3"]
genbank_ids = ["CP001509.3"]
# Email: only used for NCBI API calls, be polite!
email = ''

# Genbank file is a path to a local genbank file
# Must be a FULL genbank file for coding and noncoding modes, containing both the 
# full genome sequence and the feature definitions
# Ex. genbank_file = 'C:\Users\Me\Documents\experiments\integrate2010\genome.gb'
genbank_files = []

# genome_fasta_file is a path to a local fasta file containing the genome for 
# spacer generation with custom boundaries
# Ex. genome_fasta_file = 'C:\Users\Me\Documents\experiments\integrate2010\genome.fasta'
genome_fasta_files = []


# -------------- REGION TYPES - CODING, NONCODING, OR CUSTOM --------------#
# region_type must be 'coding', 'noncoding', or 'custom'
# see the specific additional settings for each below
region_type = 'custom'

# -------------- CODING REGION SETTINGS -----------------------------------#
# Percent of each region to target, from N terminus to C terminus
start_pct = 10
end_pct = 50

# How to prioritize choosing spacers in each coding region:
# Must be either 'N_to_C' or 'C_to_N'
coding_spacer_direction = 'N_to_C'

# Target locus tags: locus tags as known from the genbank file, referenced above
# File *must* be a csv with a column named "locus_tags"
# Ex. target_locus_tags_csv = 'C:\Users\Me\Documents\experiments\integrate2010\target_locus_tags.csv'
target_locus_tags_csv = './needed_regions.csv'


# --------------- NONCODING REGION SETTINGS -------------------------------#
# genome_boundary takes a single pair of coordinates, and will only target
# intergenic regions that lie between the two values
# Ex. noncoding_boundary = [13555, 20000]
noncoding_boundary = [0,20000]

# Restricts regions to sites that are "nonessential": those with the 
# C terminus of a coding region on both ends. Either False or True
nonessential_only = True


# --------------- CUSTOM REGION SETTINGS ----------------------------------#
# Custom regions should be defined either in the csv file, specifying what
# regions (by refseq) in the input files should be used, 
# OR
# as custom sequences input in the 'custom_sequences' field
# Leave the one not used empty


# A filepath to a csv containing pairs of corrdinates for the reference sequence to generate spacers for
# The CSV should have 3 columns with headers "genome_id", start_ref" and "end_ref"
# The "genome_id" should exactly match one of the ids, genbank files, or fasta files
# used in the matching field above
# Ex. custom_regions_csv = './custom_regions.csv'
custom_regions_csv = './custom_regions.csv'


# A list of sequences in which to generate spacers
# NOTE: when using this setting, you might want to turn 'offset' in 'advanced_parameters.py' to False. 
# See that setting for more details
# Ex. custom_sequences = ['CTCTCTCCTACTCTCTCGTACTCTCTCTACTCTCTCTACTCTCTCTACTCTCTCTACTCTCTCTACTCTCTCTACTCTCTCCTACTCTCTCTACTCTCTCTACTCTCTCTACTCTCTCTACTCTCTCTACTCTCTCTACTCTCTCTACTCTCTCTACTCTCTCTACTCTCTCTA']
custom_sequences = []

# Do not modify, this calls the function when run with 'python spacer_gen.py'
if __name__ == "__main__":
	spacer_gen({"output_path": output_path, "genbank_files": genbank_files, "genome_fasta_files": genome_fasta_files, "nonessential_only": nonessential_only, "spacers_per_region": spacers_per_region, "GC_requirement": GC_requirement, "genbank_ids": genbank_ids, "email": email, "region_type": region_type, "overlapping_spacers": overlapping_spacers, "coding_spacer_direction": coding_spacer_direction, "start_pct": start_pct, "end_pct": end_pct, "target_locus_tags_csv": target_locus_tags_csv, "noncoding_boundary": noncoding_boundary, "custom_regions_csv": custom_regions_csv, "custom_sequences": custom_sequences})
