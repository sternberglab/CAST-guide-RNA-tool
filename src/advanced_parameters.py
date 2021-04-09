# Search for up to this number of mismatches (ignoring every 6th bp)
# Note: This is the primary driver of search speed

# Ex. Searching EColi for 10 sequences per coding regions
# at 5 mismatch threshold takes ~3 minutes on a normal PC with allow_gaps = True
# Not recommended for values of 8 or higher at the moment
mismatch_threshold = 5

# Only relevant for noncoding spacer generation
# Minimum basepairs length to qualify as an intergenic region
# Note: There are single digit, even 1bp gaps between some coding sequences
# that would otherwise be included and attempted to get spacers for
minimum_intergenic_region_length = 50

# Length of spacers
# VchINTEGRATE always uses a 32bp spacer
SPACER_LENGTH = 32

# PAM Sequence to use to identify candidates
# CC is known to be a good PAM for VchINTEGRATE
# Does not currently tolerate 'N'
PAM_SEQ = 'CC'

# The distance from the spacer at which transposition occurs
# Is roughly 49bp (downstream) for INTEGRATE
# This is used to identify spacers that target specific genomic regions
INTEGRATION_SITE_DISTANCE = 49

# Changes whether the algorithm to take into account INTEGRATION_SITE_DISTANCE when searching for candidates
# e.g. offset = True limits spacer candidates to spacers that will direct integration into the target search window
# offset = False will return spacers whose 3' ends are within the target search window
offset = True

# Change flex_base = True to convert flexible bases to N for bowtie alignments. Change to False to turn off
flex_base = True
# Number of bases per 1 flexible base (default = 6)
flex_spacing = 6

# Change alow_gaps = True to allow gaps for bowtie alignment between spacer and off-targets (default = False)
allow_gaps = False

# Restriction enzymes: sequences will be filtered against this list of restriction enzymes
# Ex. restriction_enzymes = ['SalI', 'BsaI', 'BamHI', 'AvrII', 'HindIII']
restriction_enzymes = []

# homopolymer length: length of consecutive bp to count as a homopolymer for filtering
homopolymer_length = 5