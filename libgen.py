import re
import csv
from src.genbank import retrieve_annotation, get_regions
from src.advanced_parameters import PAM_SEQ, SPACER_LENGTH, INTEGRATION_SITE_DISTANCE

genbank_id="CP001509.3"
email="acreechristopher@gmail.com"

def get_duplicate_target_sites():
	# scans the genome for duplicate target sites

	# create the output csv
	with open('dupe_targets.csv', 'w', newline='') as csvfile:
		writer = csv.writer(csvfile, delimiter=',')
		writer.writerow(['number', 'target_seq', 'count', 'positions'])
	record = retrieve_annotation(genbank_id, email)
	genome_length = len(record.seq)
	match_length = len(PAM_SEQ) + SPACER_LENGTH

	# add match_length bp to end of genome to allow for matches at the end
	genome = str(record.seq.upper() + record.seq.upper()[:match_length])
	rc_genome = str(record.seq.reverse_complement() + record.seq.reverse_complement()[:match_length])
	
	candidates = dict()
	for (idx, position) in enumerate(record.seq.upper()):
		if idx % 10000 == 0:
			print(idx, len(candidates))
		if idx > 0 and idx % 50000 == 0:
			break
		match_seq = str(genome[idx:idx+match_length])
		target_seq = str(genome[idx+len(PAM_SEQ):idx+match_length])
		if match_seq[:len(PAM_SEQ)] != PAM_SEQ:
			continue

		if match_seq in candidates:
			continue

		duplicates = genome.count(match_seq) + rc_genome.count(match_seq)
		if duplicates > 1 and duplicates < 10:
			fw_dupes = [m.start() for m in re.finditer(match_seq, genome)]
			fw_dupes = [f"{str(site + len(PAM_SEQ))}_fw" for site in fw_dupes]

			rc_dupes = [m.start() for m in re.finditer(match_seq, rc_genome)]
			rc_dupes = [f"{str(genome_length - site - len(PAM_SEQ))}_rv" for site in rc_dupes]
			dupes = fw_dupes + rc_dupes
			if duplicates is not len(dupes):
				raise Exception(f"Mismatch: counted {duplicates} duplicates but found {len(dupes)} for seq {match_seq}")
			candidates[target_seq] = {"seq": target_seq, "positions": "||".join(dupes), "count": duplicates}
	with open('dupe_targets.csv', 'w', newline='') as csvfile:
		writer = csv.writer(csvfile, delimiter=',')
		i=1
		for item in candidates:
			c = candidates[item]
			writer.writerow([i, c["seq"], c["count"], c["positions"]])
			i += 1

get_duplicate_target_sites()
