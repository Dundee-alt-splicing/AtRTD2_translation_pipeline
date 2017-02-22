# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 11:25:23 2017

@author: Mark Spensley, Juan C Entizne
@email: maspensley[at]gmail.com, e.entizne[at]dundee.ac.uk
"""

import os
import argparse
import lib.translate_fix_aug
import lib.translate_longest_orf


description = "Description:\n\n" + \
              "Transfix is a translation program for nucleotide sequences. \n" \
              "It can translate the sequences either by fixing the AUG for group of transcripts " \
              "that translates into the longest ORF;  \n " \
              "or by translating each sequence in all possible ORFs and keeping the longest peptide." \

parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
subparsers = parser.add_subparsers()

parser.add_argument('-fix', '--fix-aug',
                    dest="fix",
                    action="store",
                    required=True,
                    choices=["True", "False"],
                    help="Boolean. If True, the program fix the AUG common to a group of transcripts that translates "
                         "to the longest ORF. If False, the program translated the longest ORF of each single transcript.")

parser.add_argument('-fa', '--fasta',
                    dest="fasta",
                    action="store",
                    nargs=1,
                    required=True,
                    default=None,
                    help="Path of the Fasta file.")

parser.add_argument('-o', '--output',
                    dest="outfile",
                    action="store",
                    required=True,
                    default=None,
                    help="Name of the output files.")

parser.add_argument('-g3', '--gff3',
                    dest="gff3",
                    action="store",
                    nargs=1,
                    default=None,
                    help="Path of GFF3 file.")

parser.add_argument('-g', '--gtf',
                    dest="gtf",
                    action="store",
                    nargs=1,
                    default=None,
                    help="Path of the GTF file.")

parser.add_argument('-f', '--filter',
                    dest="filter",
                    action="store",
                    required=False,
                    default=None,
                    choices=['gene', 'transcript'],
                    help="Filter fasta file by specific transcripts or by gene.")

parser.add_argument('-i', '--ids',
                    dest="ids",
                    action="store",
                    nargs=1,
                    default=[None],
                    help="Path to the file containing the list of gene/transcript ids to keep.")

parser.add_argument('-l', '--len',
                    dest="len_th",
                    action="store",
                    nargs=1,
                    type=int,
                    default=[0],
                    help="Minimum length of the translated peptide to be kept in the output file. (Default: 0)")


def create_path(lst):
    '''Check if path is absolute, if not the program use the current working path.'''

    temp_lst = []
    for fl in lst:
        if not os.path.isabs(fl):
            fl_path = os.getcwd()+"/"+fl
            temp_lst.append(fl_path)
        else:
            temp_lst.append(fl)
    return temp_lst


def main():

    args = parser.parse_args()

    fasta_file = create_path(args.fasta)

    if args.fix == "True":
        if args.gff3 and args.gtf:
            gff3_file = create_path(args.gff3)
            gtf_file = create_path(args.gtf)

            lib.translate_fix_aug.execute_translation(gff3_file[0], gtf_file[0], fasta_file[0], args.outfile)
        else:
            print("An input file (gff3 or gtf) is missing. Aborting.")
            exit()
    else:
        if args.filter:
            if not args.ids[0]:
                print("A list of the genes or transcripts to keep must be specified for filtering. Aborting.")
                exit()
            args.ids = create_path(args.ids)

        lib.translate_longest_orf.main(fasta_file[0], args.ids[0], args.filter,
                                          args.len_th[0], args.outfile)

if __name__ == "__main__":
    main()
