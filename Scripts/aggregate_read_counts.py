from argparse import ArgumentParser
import pandas as pd 

def get_args():
	p = ArgumentParser()
	p.add_argument("-i", "--inputs", nargs="+", help="List of files containing read counts for each sample")
	p.add_argument("output_file", help="Destination file for stats.")
	args = vars(p.parse_args())
	return args



# Extracts read counts for a sample at different stages of processing and places them in a single-row dataframe.
# @file: A single-column csv file containing read counts at each stage of processing.
def file_to_csv_line(file):
	df = pd.read_csv(file, header = 0, index_col=0)
	return df.T

def import_basic_read_stats(files):
	basic_read_stats = pd.DataFrame()
	for file in files:
		line = file_to_csv_line(file)
		if basic_read_stats.empty:
			basic_read_stats = line 
		else:
			basic_read_stats = basic_read_stats.append(line)
	return basic_read_stats

def calculate_percentage_read_stats(read_stats):
	if "Trimmed" in read_stats.columns.values:
		if "Deduplicated" in read_stats.columns.values:
			read_stats["Duplication Rate"] = (read_stats["Trimmed"]-read_stats["Deduplicated"])/read_stats["Trimmed"]*100
			read_stats["Duplication Rate"] = read_stats["Duplication Rate"].round(decimals = 2)
			if "Normalized" in read_stats.columns.values:
				read_stats["% Reads Removed by Normalization"] = (read_stats["Deduplicated"]-read_stats["Normalized"])/read_stats["Deduplicated"]
				read_stats["% Reads Removed by Normalization"] = read_stats["% Reads Removed by Normalization"].round(decimals = 2)

	return read_stats

def write_csv_to_file(df, output_file):
	df.to_csv(output_file)


args = get_args()

basic_read_stats = import_basic_read_stats(args["inputs"])
read_stats = calculate_percentage_read_stats(basic_read_stats)

write_csv_to_file(basic_read_stats, args["output_file"])
