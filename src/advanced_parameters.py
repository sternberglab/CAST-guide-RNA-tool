# Search for up to this number of mismatches (ignoring every 6th bp)
# Note: This is the primary driver of search speed
# Ex. Searching EColi for 10 sequences at 5 mismatch threshold takes ~3 minutes on a normal PC
mismatch_threshold = 4

# Minimum basepairs length to qualify as an intergenic region
# Note: There are single digit, even 1bp gaps between some coding sequences
# that would otherwise be included and attempted to get spacers for
minimum_intergenic_region_length = 50

# Length of spacers
# VchINTEGRATE always uses a 32bp spacer
SPACER_LENGTH = 32

# PAM Sequence to use to identify candidates
# CC is known to be a good PAM for INTEGRATE
PAM_SEQ = 'CC'

# The distance from the spacer at which transposition occurs
# Is roughly 49bp (downstream) for INTEGRATE
# This is used to identify spacers that target specific genomic regions
INTEGRATION_SITE_DISTANCE = 49