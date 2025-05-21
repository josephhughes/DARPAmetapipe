import re
from argparse import ArgumentParser 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

def get_args():
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", help = "File containing sequence to be condensed.")
    parser.add_argument("-l", "--output_l_seg", help = "Output file to save condensed sequence to.")
    parser.add_argument("-s", "--output_s_seg", help = "Output file to save condensed sequence to.")
    args = vars(parser.parse_args())
    return args 


def condense_sequence(input_file, l_output, s_output):
    with open(input_file, 'r') as f:
        lines = [line.rstrip('\n') for line in f]
        new_file = ''
        for line in lines:
            new_file += line

        first_seg = re.search(".*(?=(>))", new_file)[0]

        first_seg_name = re.search("[0-9]{5}_L", first_seg)[0]
        first_seg_seq = Seq(re.search("(?<=([0-9]{5}_L)).*", first_seg)[0], "ambiguous_dna") 

        first_seg_record = SeqRecord(first_seg_seq, id = first_seg_name, description = "")

        SeqIO.write([first_seg_record], l_output, "fasta")


        second_seg = re.search("(?<=(.))>.*$", new_file)[0]


        second_seg_name = re.search("[0-9]{5}_S", second_seg)[0] 
        second_seg_seq = Seq(re.search("(?<=([0-9]{5}_S)).*", second_seg)[0], "ambiguous_dna") 


        second_seg_record = SeqRecord(second_seg_seq, id = second_seg_name, description = "")

        SeqIO.write([second_seg_record], s_output, "fasta")

    return




args = get_args()

condense_sequence(args["input"], args["output_l_seg"], args["output_s_seg"])