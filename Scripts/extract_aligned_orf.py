import pandas as pd 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
import re

def read_in_blast_output(blast_out_file):
    return pd.read_csv(blast_out_file, sep = '\t', header = 0).set_index("qseqid")


def get_contig_with_blastx_hit(contig_name, contigs_file):
    for rec in SeqIO.parse(contigs_file, "fasta"):
         if rec.name == contig_name:
            return rec
    return None


def get_aligned_contig_region_as_seqfeature(row):
    contig_align_start = min(row["qstart"], row["qend"])
    contig_align_end = max(row["qstart"], row["qend"])

    # Determine if alignment is sense or antisense based on start/end positions
    if row["qstart"] < row["qend"]:
        strand = +1
    else:
        strand = -1

    return SeqFeature(FeatureLocation(contig_align_start, contig_align_end, strand = strand))


# Given alignment metrics for a subject protein sequence, calculates
# how many nucleotides are missing on either end, then adds another 90.
def _get_alignment_expansion_parameters(alignment_info):
    polymerase_align_start = alignment_info["sstart"]
    polymerase_align_end = alignment_info["send"]
    polymerase_length = alignment_info["slen"]

    upstream_expansion = (polymerase_align_start*3)+45
    downstream_expansion = ((polymerase_length - polymerase_align_end)*3)+45

    return upstream_expansion, downstream_expansion


def expand_alignment_feature(aligned_contig_region, alignment_info, contig):
    upstream_expansion, downstream_expansion = _get_alignment_expansion_parameters(alignment_info)
    contig_length = len(contig.seq)

    # Subject alignment metrics are always in sense orientation.
    # Any down/upstream adjsutments must take into account query orientation.
    strand = aligned_contig_region.location.strand
    # Sense-oriented contig
    if strand == 1:
        # Start location represents 5' end, subtract upstream expansion to expand
        new_start = aligned_contig_region.location.start - upstream_expansion
        # End location represents 3' end, add downstream expansion to expand
        new_end = aligned_contig_region.location.end + downstream_expansion
    # Antisense-oriented contig
    elif strand == -1:
        # Start location represents 3' end, subtract downstream expansion to expand.
        new_start = aligned_contig_region.location.start - downstream_expansion
        # End location represents 5' end, add upstream expansion to expand.
        new_end = aligned_contig_region.location.end + upstream_expansion

    # Alignment shouldn't extend past contig bounds.
    # Correct start/end locatins if they violate this.
    if new_start<0:
        new_start = 0
    if new_end>contig_length:
        new_end = contig_length

    return SeqFeature(FeatureLocation(new_start, new_end, strand = strand))


def extract_longest_orf_from_region_in_contig(contig, aligned_contig_region):
    sequence = aligned_contig_region.extract(contig.seq)

    all_orfs = []

    longest_orf = max(re.findall(r'ATG(?:(?!TAA|TAG|TGA)...)*(?:TAA|TAG|TGA)', str(sequence)), key = len)

    return(longest_orf)


def contig_is_shorter_than_polymerase(contig, polymerase_length):
    if len(contig.seq) < polymerase_length*3:
        return True
    else:
        return False

def aligned_region_is_shorter_than_polymerase(contig, aligned_contig_region, polymerase_length):
    aligned_sequence = aligned_contig_region.extract(contig.seq)
    aligned_length = len(aligned_sequence)
    print(contig.name)
    print(str(polymerase_length*3-3))
    print(aligned_length)
    if aligned_length < polymerase_length*3-3:
        return True 
    else:
        return False


def aligned_region_to_SeqRecord(contig, contig_description, aligned_contig_region):
    aligned_sequence = str(aligned_contig_region.extract(contig.seq))
    # print(aligned_sequence)

    aligned_record = SeqRecord(Seq(aligned_sequence), id=contig.id, name=contig.name, description = contig_description)
    return aligned_record


def extract_cds_for_each_blastx_hit(blastx_df, contigs_file):
    cds_recs = []
    for contig_name, alignment_info in blastx_df.iterrows():
        blastx_hit_contig = get_contig_with_blastx_hit(contig_name, contigs_file)
        polymerase_hit = alignment_info["sseqid"]

        # If contig is too small to contain the whole ORF, no point in doing further processing
        # Just extract it as is.
        if contig_is_shorter_than_polymerase(blastx_hit_contig, alignment_info["slen"]):
            blastx_hit_contig.description = "Incomplete {pol}-like Sequence".format(pol=polymerase_hit)
            cds_recs.append(blastx_hit_contig)
            continue
        else:
            aligned_contig_region = get_aligned_contig_region_as_seqfeature(alignment_info)

            aligned_contig_region = expand_alignment_feature(aligned_contig_region, alignment_info, blastx_hit_contig)

            if aligned_region_is_shorter_than_polymerase(blastx_hit_contig, aligned_contig_region, alignment_info["slen"]):
                description = "Incomplete {pol}-like Sequence".format(pol=polymerase_hit)
                cds_recs.append(aligned_region_to_SeqRecord(blastx_hit_contig, description, aligned_contig_region))
                continue
            else:
                description = "{pol}-like Sequence".format(pol=polymerase_hit)
                longest_orf_in_aligned_region = extract_longest_orf_from_region_in_contig(blastx_hit_contig, aligned_contig_region)
                longest_orf_rec = SeqRecord(Seq(longest_orf_in_aligned_region), id=contig_name, name=contig_name, description = description)
                cds_recs.append(longest_orf_rec)
    return cds_recs

        
blastx_df = read_in_blast_output(snakemake.input.blastx_out)
polymerase_like_cds = extract_cds_for_each_blastx_hit(blastx_df, snakemake.input.herpesvirus_contigs)
SeqIO.write(polymerase_like_cds, snakemake.output.output_fasta, "fasta-2line")

# print(blastx_df)