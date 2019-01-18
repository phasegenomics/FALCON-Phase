# take a combined p and h contig fasta, split it into p and h contigs,
# and create name_mapping.txt file for falcon-phase.

from __future__ import print_function
import os
import sys
import re
import argparse
from Bio import SeqIO

FC_PAT = re.compile(r'([0-9]{6}F[_0-9]*)')

def parse_args(desc):
    '''parse command-line args

    Args:
        desc(str): program description, e.g. __file__

    Returns:
        args (dict): dict of the form {arg_name: arg_value}
    '''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-c", "--combined_fasta", required=False, default=None,
                        help="Assembly fasta with both p and h contigs to separate and map to each other.")
    parser.add_argument("-p", "--p_contig_fasta", required=False, default=None,
                        help="Assembly fasta with p contigs only. If supplied, requires also h contig fasta (-h).")
    parser.add_argument("-h", "--h_contig_fasta", required=False, default=None,
                        help="Assembly fasta with h contigs only. If supplied, requires also p contig fasta (-p).")
    parser.add_argument("-o", "--outfile_prefix", required=True,
                        help="Prefix for new files to be written.")

    args = parser.parse_args()
    return vars(args)

def is_p_contig(scrubbed_id):
    '''based on tig name, test for whether it is p
    '''
    return len(scrubbed_id.split("_")) == 1
    
def is_h_contig(scrubbed_id):
    '''based on tig name, test for whether it is h
    '''
    return len(scrubbed_id.split("_")) == 2

def scrub_name(fasta_id, fc_pat=FC_PAT):
    '''remove non-falcon-unzip stuff from contig names
    '''
    match = re.search(fc_pat, fasta_id).group()
    if match[-1] == "_":
        match = match[:-1]
    return match

def get_p_of_h_tig(scrubbed_h_id):
    '''based on h tig name, get corresponding p tig
    '''
    return scrubbed_h_id.split("_")[0]

def parse_combined_fasta(fasta_path):
    '''parse the combined fasta, collect p and h tigs separately in lists of SeqRecords.
    '''
    print("parsing combined fasta of p and h tigs.")
    handle = SeqIO.parse(fasta_path, "fasta")
    p_tigs = []
    h_tigs = []
    all_records = set()
    for record in handle:
        scrubbed_id = scrub_name(record.id)
        record.id = scrubbed_id
        record.name = None
        record.description = ""
        if scrubbed_id in all_records:
            raise ValueError("tig with name {0} appears multiple times!!\nnote that " \
                            "tig names are scrubbed back to falcon-unzip output!".format(scrubbed_id)
                            )
        all_records.add(scrubbed_id)
        if is_p_contig(scrubbed_id):
            p_tigs.append(record)
        elif is_h_contig(scrubbed_id):
            h_tigs.append(record)
        else:
            print(record, scrubbed_id)
            raise ValueError("neither p nor h contig! is this falcon-unzip-derived??")
    print("found", len(p_tigs), "p tigs and", len(h_tigs), "h tigs")
    return p_tigs, h_tigs

def make_name_mappings(p_tigs, h_tigs):
    '''create name_mappings of p to h tigs for falcon-phase, write it out to a file. 
    If p tigs are missing, omit corresponding h tigs while printing out that info.
    '''
    mappings = {}
    omitted_h_tigs = []
    for record in p_tigs:
        mappings[record.id] = []
    for record in h_tigs:
        p_tig = get_p_of_h_tig(record.id)
        if p_tig not in mappings.keys():
            print("p tig {0} doesn't exist! Omitting corresponding h tigs".format(p_tig))
            omitted_h_tigs.append(record.id)
            continue
        mappings[p_tig].append(record.id)
    
    with open("name_mappings.txt", "w") as mapping_file:
        for p_name in sorted(mappings.keys()):
            for h_name in sorted(mappings[p_name]):
                mapping_file.write("{0}\t{1}\n".format(p_name, h_name))
    print("omitted {0} h tigs without corresponding p tigs.".format(len(omitted_h_tigs)))
    return omitted_h_tigs
    
def main():
    c_args = parse_args()

    if c_args["combined_fasta"] is not None:
        p_tigs, h_tigs = parse_combined_fasta(c_args["combined_fasta"])
    elif c_args["p_contig_fasta"] is not None and c_args["h_contig_fasta"] is not None:
        p_tigs = SeqIO.parse(c_args["p_contig_fasta"], "fasta")
        h_tigs = SeqIO.parse(c_args["h_contig_fasta", "fasta"])
    else:
        ValueError("Script requires either a combined fasta (-c arg) or both p and h fastas (-p and -h args)." \
                   "If both supplied, combined FASTA is used.")

    omitted_h_tigs = make_name_mappings(p_tigs, h_tigs)
    # filter omitted h tigs
    h_tigs = [h_tig for h_tig in h_tigs if h_tig.id not in omitted_h_tigs]

    # write out 2 fastas
    SeqIO.write(p_tigs, c_args["outfile_prefix"] + ".p_tigs.fasta", "fasta")
    SeqIO.write(h_tigs, c_args["outfile_prefix"] + ".h_tigs.fasta", "fasta")
    
if __name__ == "__main__":
    main()
