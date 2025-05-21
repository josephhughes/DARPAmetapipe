from ete3 import NCBITaxa
import pandas as pd 
from argparse import ArgumentParser


def get_args():
    parser = ArgumentParser()

    parser.add_argument("-i", "--input_file", required = True, 
        help = "Input .tsv with occurence counts and TaxID for each SILVA rRNA accession.")
    parser.add_argument("-r", "--ranks", nargs = '*', type = str, default = ["superkingdom"], 
        help = "Taxonomic ranks to identify for contigs.")
    parser.add_argument("-o", "--output_file", required = True, 
        help = "Output .tsv file to write contig data updated with taxonomic rank.")

    args = vars(parser.parse_args())
    return args


def read_in_taxonomy_db():
    ncbi = NCBITaxa()
    return ncbi


def read_ribopicker_tsv(input_file):
    col_names = ["SILVA_accession", "TaxID", "Count"]
    tax_df = pd.read_csv(input_file, sep = '\t', names = col_names, index_col = 0, dtype = {"SILVA_accession": str, "TaxID": str, "Count": int})
    print(tax_df)
    return tax_df


def get_specified_parent_ranks(row, ranks):
    print(row["TaxID"])
    try:
        lineage = ncbi.get_lineage(row["TaxID"])
        lineage_ranks = ncbi.get_rank(lineage)

        for rank in ranks:
            if rank in lineage_ranks.values():
                for lineage_taxid in lineage_ranks:
                    if lineage_ranks.get(lineage_taxid) == rank:
                        row[rank] = str(lineage_taxid)
            else:
                row[rank] = "N/A"
    except ValueError:
        for rank in ranks:
            row[rank] = "N/A"

    return row


ncbi = None

def main():
    # Get input parameters.
    args = get_args()
    
    # Read NCBI Taxonomy database into memory
    global ncbi
    ncbi = read_in_taxonomy_db()
    
    # Read TaxID data for contigs in as dataframe
    tax_df = read_ribopicker_tsv(args.get("input_file"))
    
    # Identify taxonomic units for each contig based on LCA.
    tax_df = tax_df.apply(get_specified_parent_ranks, axis = 1, args = (args.get("ranks"),))
    
    # Write updated contig dataframe to file.
    tax_df.to_csv(args.get("output_file"), sep = '\t')


if __name__ == "__main__":
    main()