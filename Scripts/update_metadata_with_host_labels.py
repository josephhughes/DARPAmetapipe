from argparse import ArgumentParser 
import pandas as pd 

def get_args():
    parser = ArgumentParser()
    parser.add_argument("-i", "--input_file", help = ".tsv file containing metadata")
    parser.add_argument("-o", "--output_file", help = ".tsv file to write updated metadata to.")
    args = vars(parser.parse_args())
    return args 


args = get_args() 

metadata_df = pd.read_csv(args["input_file"], sep = '\t')

print(metadata_df)

homo_sapiens_mask = metadata_df.loc[:, "Host_Simple"] == "Homo sapiens"

metadata_df.loc[homo_sapiens_mask, "Host"] = "Human"

metadata_df.loc[~homo_sapiens_mask, "Host"] = "Rodent"

unspecified_mask = metadata_df.loc[:, "Host_Simple"] == "Unspecified"

metadata_df.loc[unspecified_mask, "Host"] = "Unspecified"

metadata_df.to_csv(args["output_file"], sep = '\t', index = False)