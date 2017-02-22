from lib.gene_models import *


def get_models_gtf(gtf_filename):
    '''Build gene models from gtf file'''
    infile = open(gtf_filename, 'rU')
    gff_lines = infile.readlines()
    infile.close()
    locus_dict = {}
    for line in gff_lines:
        line = line.strip()
        if line.startswith('#'):
            continue
        line_fields = line.split('\t')
        transcript_ID = line_fields[8].split(';')[0][15:-1]
        locus_ID = line_fields[8].split(';')[1][10:19]

        if locus_ID not in locus_dict:
            gene_model = Gene(locus_ID)
            gene_model.add_sense(line_fields[6])
            locus_dict.update({locus_ID: gene_model})
        else:
            gene_model = locus_dict[locus_ID]
            
        if transcript_ID not in gene_model.transcript_dict:
            transcript_model = Transcript(transcript_ID)
            transcript_model.sense = line_fields[6]
            gene_model.transcript_dict.update({transcript_ID:transcript_model})
        transcript_model = gene_model.transcript_dict[transcript_ID]    
  
        if line_fields[2] == 'CDS':# and proteinCoding == True:
            cds_coords = [int(line_fields[3]), int(line_fields[4])]
            transcript_model.add_CDS(cds_coords)

        if line_fields[2] == 'exon':# and proteinCoding == True:
            exon_coords = [int(line_fields[3]), int(line_fields[4])]
            transcript_model.add_exon(exon_coords)
        
        gene_model.transcript_dict.update({transcript_ID:transcript_model})
        locus_dict.update({locus_ID: gene_model})
    return locus_dict


# add sequences to GTF models
def add_gtf_seqs(fasta_filename, locus_dict):
    '''Read gtf sequences from fasta file and add
    to gene models in locus_dict'''
    infile = open(fasta_filename, 'rU')
    for line in infile:
        line = line.strip()
        if line.startswith('>'):
            seq = ''
            line = line.split(' ')
            model_id = line[0][1:]
            gene_id = line[0][1:10]
            gene_model = locus_dict[gene_id]
            transcript_model = gene_model.transcript_dict[model_id]
        else:
            transcript_model.seq += line
    infile.close()
    return locus_dict
