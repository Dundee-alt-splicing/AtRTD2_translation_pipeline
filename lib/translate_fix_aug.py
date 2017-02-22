#!/usr/bin/python

from lib.load_gff3 import *
from lib.load_gtf import *
from lib.reannotate_to_max_orf import *


def translate_seq(sequence):
    'Returns the amino acid translation of a given sequence'
    codon_table = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 
                   'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
                   'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*', 
                   'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',
                   'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 
                   'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
                   'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 
                   'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 
                   'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 
                   'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
                   'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 
                   'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
                   'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V', 
                   'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
                   'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 
                   'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}
    codon_pos = 0
    trans = ''
    while codon_pos < (len(sequence) - 2):     
        codon = sequence[codon_pos:codon_pos + 3]
        if codon in codon_table:
            aa = codon_table[codon]
        else:
            aa = 'x'
        trans = trans + aa
        if aa == '*':
            return trans
        codon_pos += 3
    return trans


def print_list(my_list, filename):
    outfile = open(filename, 'w')
    for item in my_list:
        outfile.writelines(item+'\n')
    outfile.close()


def execute_translation(gff3_file, gtf_file, fasta_file, outfile_name):

    gff3_models = get_models_gff3(gff3_file)
    gff3_max_orf_models, ara_n_list, gff3_models = characterise_max_orfs2(gff3_models)

    gtf_models = get_models_gtf(gtf_file)
    gtf_models = add_gtf_seqs(fasta_file, gtf_models)

    outfile = open(outfile_name+".fa", 'w')

    no_exon_counter, model_count = (0 for _ in range(2))
    absent_gff3, not_atg_start, rep_atg_absent, retro_transposons = ([] for _ in range(4))
    for locus_id in sorted(gtf_models.keys()):
        if locus_id not in gff3_models:
            absent_gff3.append(locus_id)
            continue
        if gff3_models[locus_id].transposon is True:
            retro_transposons.append(locus_id)
            continue

        atg_pos = gff3_models[locus_id].rep_atg
        gene_model = gtf_models[locus_id]
        transcript_dict = gene_model.transcript_dict
        for transcript_id in transcript_dict:
            model_count += 1

            transcript_model = transcript_dict[transcript_id]

            if transcript_model.sense == '+':
                exon_list = transcript_model.exon_list
            else:
                exon_list = transcript_model.exon_list[::-1]

            cds_index = 0
            z = False
            for exon in exon_list:
                # Note exon[0] < exon[1]
                if atg_pos not in range(exon[0], exon[1]+1):
                    cds_index += 1+exon[1] - exon[0]
                else:
                    if transcript_model.sense == '+':
                        cds_index += atg_pos - exon[0]
                        z = True
                        break
                    else:
                        cds_index += exon[1] - atg_pos
                        z = True
                        break
            if z is False:
                rep_atg_absent.append(transcript_id)
                continue

            seq = transcript_model.seq[cds_index:]
            if not seq.startswith('ATG'):
                not_atg_start.append(transcript_id)
                continue
            else:
                peptide = translate_seq(seq)
                outfile.writelines('>'+transcript_id+'\n')
                outfile.writelines(peptide+'\n\n')

    outfile.close()

    print("Total models considered: ", model_count)

    print("locus id not in GFF3 (%s_non-gff3_locus.txt) " % outfile_name, len(absent_gff3))
    print_list(absent_gff3, outfile_name+'_non-gff3_locus.txt')

    print("transposon loci (%s_retrotransposons.txt) " % outfile_name, len(retro_transposons))
    print_list(retro_transposons, outfile_name+'_retrotransposons.txt')

    print("rep atg not in transcript (%s_start_absent.txt) " % outfile_name, len(rep_atg_absent))
    print_list(rep_atg_absent, outfile_name+'_start_absent.txt')

    print("start position is not ATG (%s_non_atg_start.txt): " % outfile_name, len(not_atg_start))
    print_list(not_atg_start, outfile_name+'_non_atg_start.txt')
