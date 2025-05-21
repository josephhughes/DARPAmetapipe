from argparse import ArgumentParser
import pandas as pd 
from Bio import SeqIO
import re
from datetime import datetime


def get_args():
    parser = ArgumentParser()
    parser.add_argument("-i", "--input_files", required = True, nargs = '+',
        help = ".gbk file(s) containing genbank records to amend")
    parser.add_argument("-m", "--metadata_file", required = True,
        help = ".tsv file containing metadata to add to records")
    args = vars(parser.parse_args())
    return args 

# Parses metadata file into a pandas csv, formats column names and dates
def read_metadata_file(metadata_file):
    # Read in mouse metadata and set index to mouse ID (e.g. 00106)
    meta_df = pd.read_csv(metadata_file, sep = '\t', dtype = {"Animal ID (SL#####)":str}).set_index("Animal ID (SL#####)")

    # Simplify column names
    meta_df.rename(columns = {"Collection Date:": "Collection Date", "Species ID (final)":"Species"}, inplace = True)

    # Convert date format from YYYY-MM-DD to DD-MMM-YYYY (2019-09-22 becomes 22-Sep-2019)
    # Conforms to general GenBank date format
    meta_df["Collection Date"] = meta_df["Collection Date"].apply(lambda x: datetime.strptime(x, "%Y-%m-%d").strftime("%d-%b-%Y"))

    return meta_df

# Return source feature for a SeqRecord
def get_source_feature(record):
    for feature in record.features:
        if feature.type == "source":
            return feature 
    return None


def main():
    args = get_args()

    # Read in metadata
    meta_df = read_metadata_file(args["metadata_file"])

    # For each input file, add species and collection date for every included record.
    for file in args["input_files"]:
        updated_recs = []
        for rec in SeqIO.parse(file, "genbank"):
            source_feature = get_source_feature(rec)
            # Record name includes animal number and segment, want only animal number.
            animal_id = re.split("_", rec.name)[0] # 00106_L becomes 00106

            # Add/replace collection date and host species data
            source_feature.qualifiers["collection_date"] = [meta_df.loc[animal_id, "Collection Date"]]
            source_feature.qualifiers["host"] = [meta_df.loc[animal_id, "Species"]]
            updated_recs.append(rec)
        # Rewrite input file
        SeqIO.write(updated_recs, file, "genbank")



if __name__ == "__main__":
    main()
