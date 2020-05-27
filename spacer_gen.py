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
GC_requirement = [0, 100]

# Preferences on overlapping spacers - can be "allowed", "avoid", "forbidden"
# "allowed" will not manipulate ordering apart from other settings
# "avoid" will permit overlapping spacers only if sufficient non-overlapping spacers cannot be generated
# "forbidden" will not return overlapping spacers, even if it means the specified number per region cannot be generated
overlapping_spacers = "avoid"

# -------------- REGION TYPES - CODING, NONCODING, OR CUSTOM --------------#
# region_type must be 'coding', 'noncoding', or 'custom'
# see the specific settings for each below
region_type = 'noncoding'

# -------------- CODING REGION SETTINGS -----------------------------------#
# Percent of each region to target, from N terminus to C terminus
start_pct = 10
end_pct = 50

# How to prioritize choosing spacers in each region:
# Must be either "N_to_C" search or "C_to_N" search
coding_spacer_direction = "N_to_C"

# Target locus tags: locus tags as known from the genbank file, referenced above
# File *must* be a csv with a column named "locus_tags"
# Ex. target_locus_tags_csv = 'C:\Users\Me\Documents\experiments\integrate2010\target_locus_tags.csv'
target_locus_tags_csv = './locus_tags.csv'


# --------------- NONCODING REGION SETTINGS -------------------------------#
# genome_boundary takes a single pair of coordinates, and will only target
# intergenic regions that lie between the two values
# Ex. [13555, 20000]
noncoding_boundary = [0,10000]

# Restricts regions to sites that are "nonessential": those with the 
# C terminus of a coding region on both ends
nonessential_only = True


# --------------- CUSTOM REGION SETTINGS ----------------------------------#
# A filepath to a csv containing pairs of corrdinates for the reference sequence to generate spacers for
# The CSV should have 2 columns with headers "start_ref" and "end_ref"
# Ex. custom_regions_csv = 'C:\Users\Me\Documents\experiments\integrate2010\custom_regions.csv'
custom_regions_csv = './locus_tags.csv'


# Do not modify, this calls the function when run with 'python spacer_gen.py'
if __name__ == "__main__":
	spacer_gen({"output_path": output_path, "nonessential_only": nonessential_only, "spacers_per_region": spacers_per_region, "GC_requirement": GC_requirement, "genbank_id": genbank_id, "email": email, "region_type": region_type, "overlapping_spacers": overlapping_spacers, "coding_spacer_direction": coding_spacer_direction, "start_pct": start_pct, "end_pct": end_pct, "target_locus_tags_csv": target_locus_tags_csv, "noncoding_boundary": noncoding_boundary, "custom_regions_csv": custom_regions_csv})
