# Functions to read GFF3 file and extract details to...
#   populate the Gene and Transcript models in gene_models.py 
from lib.gene_models import *


def add_gff3_feature(locus_model, line_fields):
    '''Called by get_models_gff3()
    instatiantiates Transcript() and adds
    features to the model
    '''
    feature = line_fields[8]
    if 'retrotransposon' in feature:
        locus_model.transposon = True
    if line_fields[2] == 'mRNA':
        transcript_ID_list = [feature[3:14]]
    else:
        if feature.startswith('ID='): # Araport11 style
            feature = feature.split(';')
            transcript_ID_list = feature[1][7:].split(',')

    for transcript_ID in transcript_ID_list:
        assert transcript_ID.startswith('AT')
        assert transcript_ID[:9] == locus_model.name
            
        if transcript_ID not in locus_model.transcript_dict:
            transcript_model = Transcript(transcript_ID)
            transcript_model.sense = line_fields[6]
            locus_model.transcript_dict.update({transcript_ID:transcript_model})
        transcript_model = locus_model.transcript_dict[transcript_ID]
       
        if line_fields[2] == 'mRNA':
            transcript_ID = feature[3:14]
            mRNA_bounds = [line_fields[3], line_fields[4]]
            transcript_model.mRNA = mRNA_bounds

        if line_fields[2] == 'CDS':
            cds_coords = [int(line_fields[3]), int(line_fields[4])]
            transcript_model.add_CDS(cds_coords)

        if line_fields[2] == 'exon':
            exon_coords = [int(line_fields[3]), int(line_fields[4])]
            transcript_model.add_exon(exon_coords)
            locus_model.transcript_dict.update({transcript_ID: transcript_model})
    return locus_model


def get_models_gff3(filename):
    '''Opens GFF3 file with <filename>, reads and
    instantiates the Gene model at the first occurrence of
    each locus ID. Additional locus features are extracted
    by calling add_feature_gff3()
    '''
    infile = open(filename, 'rU')
    gff_lines = infile.readlines()
    infile.close()
    locus_dict = {}
    for line in gff_lines:
        line = line.strip()
        line_fields = line.split('\t')
        if line.startswith('#'):
            continue
        # second index allows for 'genes' and 'transposable_element_gene' etc
        if line_fields[2][-4:] == 'gene':
            locus_ID = line_fields[8].split(';')[0][3:12]
            assert locus_ID not in locus_dict   
            gene_model = Gene(locus_ID)
            gene_coords = [int(line_fields[3]), int(line_fields[4])]
            gene_model.add_bounds(gene_coords)
            gene_model.add_sense(line_fields[6])
            locus_dict.update({locus_ID: gene_model})                       
        else:  
            locus_model = locus_dict[locus_ID]      
            locus_model = add_gff3_feature(locus_model, line_fields)
    return locus_dict



