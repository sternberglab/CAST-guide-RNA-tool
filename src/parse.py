import csv
from pathlib import Path

def extract_column_from_csv(filepath, column_name):
	print(Path(filepath))
	with open(Path(filepath), 'r', encoding='utf-8-sig') as csv_file:
		reader = csv.DictReader(csv_file)
		data = [r[column_name] for r in reader]
	return data
