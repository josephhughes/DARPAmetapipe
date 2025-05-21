from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from urllib.error import HTTPError
from argparse import ArgumentParser
import datetime
import re
import os
import sys


Entrez.email = "Matej.Vucak@glasgow.ac.uk"

def get_args():
	parser = ArgumentParser()
	parser.add_argument("-o", "--output_file", required = True, help="Target file where downloaded records are stored.")
	parser.add_argument("-s","--search_term", required = True, help="The complete term to search genbank with. e.g. \"Lassa Mammarenavirus[ORGN]\"")
	parser.add_argument("-f", "--format", required = True, help="Format to store downloaded records in. e.g. fasta, gb")
	args = parser.parse_args()
	return args




def download_accessions(output_file, search_term, output_format):
	# Search for the query, save history details and number of records found.
	search = Entrez.esearch(db="nucleotide", term=search_term, usehistory="y", idtype="acc")
	res = Entrez.read(search)
	search.close()
	
	num_records_to_download = int(res["Count"])
	webenv = res["WebEnv"]
	query_key = res["QueryKey"]
	
	batch_size = 400 # Arbitrary batch size.
	max_attempts = 3


	records_to_download = []
	# Download records in batches of <batch_size> and save them to file.
	for start in range(0, num_records_to_download, batch_size):
		next_max = start + batch_size
		print("Downloading records " + str(start) + " through " + str(next_max))
		attempt=0
		while attempt < max_attempts:
			attempt += 1
			try:
				fetch_handle = Entrez.efetch(db="nucleotide", rettype = output_format, retmode = "text", retstart = start, retmax = batch_size, webenv = webenv, query_key=query_key, idtype="acc")
			except HTTPError as err:
				if 500<= err.code <=599:
					print("Received error from server %s" % err)
					print("Attempt %i of 3" % attempt)
					time.sleep(15)
				else:
					raise
		# data = fetch_handle.read()
		
		recs = SeqIO.parse(fetch_handle, output_format)
		for rec in recs:
			records_to_download.append(rec)
		
		fetch_handle.close()
	SeqIO.write(records_to_download, output_file, output_format)


args = get_args()


download_accessions(args.output_file, args.search_term, args.format) 






