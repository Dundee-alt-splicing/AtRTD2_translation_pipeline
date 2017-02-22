import os
from Bio.Seq import Seq
from collections import defaultdict


def filter_fasta_by_gene_id(fasta, ids_lst):

    ids = []
    with open(ids_lst, 'r') as fh:
        for ln in fh:
            ids.append(ln.strip("\n"))

    lines = []
    with open(fasta) as fh:
        for ln in fh:
            if ln.startswith(">"):
                trans_id = ln.strip(">").split(" ")[0]
                gene_id = trans_id.replace(".", " ").replace("_", " ").split(" ")[0]
                if gene_id in ids:
                    keep = True
                else:
                    keep = False
            if keep:
                lines.append(ln)

    outfile_fasta = fasta.split("/")[-1].split(".")[0]+"_filtered.fa"

    with open(outfile_fasta, "w+") as fh:
        for ln in lines:
            fh.write(ln)

    return os.getcwd()+"/"+outfile_fasta


def filter_fasta_by_trans_id(fasta, ids_lst):

    ids = []
    with open(ids_lst, 'r') as fh:
        for ln in fh:
            ids.append(ln.strip("\n"))

    lines = []
    with open(fasta) as fh:
        for ln in fh:
            if ln.startswith(">"):
                trans_id = ln.strip(">").split(" ")[0]
                if trans_id in ids:
                    keep = True
                else:
                    keep = False
            if keep:
                lines.append(ln)

    outfile_fasta = fasta.split("/")[-1].split(".")[0] + "_filtered.fa"

    with open(outfile_fasta, "w+") as fh:
        for ln in lines:
            fh.write(ln)

    return os.getcwd()+"/"+outfile_fasta


def get_fasta_seq(fasta):

    lines, nuc_seq, seq_dt = [], "", {}
    with open(fasta) as fh:
        for ln in fh:
            if ln.startswith(">"):
                key = ln.strip("\n")
                nuc_seq = ""
            else:
                nuc_seq += ln.strip("\n")
                seq_dt[key] = nuc_seq

    return seq_dt


def translate_frames(fasta_dt, outname, len_th):

    frames_dt = defaultdict(list)
    for key in fasta_dt:
        for i in range(3):
            seq = Seq(fasta_dt[key][i:])
            rec_seq = seq.reverse_complement()
            for s in [seq, rec_seq]:
                frames_dt[key].append(s.translate())

    with open("%s.fa" % outname, "w+") as fh:
        prev_len = 0
        for key in frames_dt:
            for el in frames_dt[key]:
                longest_fragment = str(sorted(el.split("*"), key=len)[-1])
                if "M" in longest_fragment:
                    longest_orf = "M"+"M".join(longest_fragment.split("M")[1:])
                else:
                    longest_orf = ""
                new_len = len(longest_orf)
                if new_len > prev_len:
                    longest_peptide = longest_orf
                    prev_len = new_len
            if longest_orf != "":
                if len(longest_peptide) >= len_th:
                    fh.write(key+"\n")
                    fh.write(longest_peptide+"*"+"\n")
                else:
                    pass
            else:
                pass
            prev_len = 0


def main(fasta, ids_lst, filter, len_th, outname):

    if not filter:
        filtered_fasta = fasta
    else:
        if filter == "gene":
            filtered_fasta = filter_fasta_by_gene_id(fasta, ids_lst)
        else:
            filtered_fasta = filter_fasta_by_trans_id(fasta, ids_lst)

    fasta_dt = get_fasta_seq(filtered_fasta)
    translate_frames(fasta_dt, outname, len_th)

    os.remove(filtered_fasta)
