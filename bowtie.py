import os
from Bio import SeqIO, Entrez
from variables import *
from run_kallisto import *


def get_fasta(accession):           # writes fasta file from accession provided
    Entrez.email = "ecrum@luc.edu"
    net_handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta")      # fetches GenBank record
    out_handle = open(str(accession) + '.fasta', "w")                               # writes record to local file
    recs = SeqIO.parse(net_handle, 'fasta')
    for r in recs:
        SeqIO.write(r, out_handle, 'fasta')
    return accession + '.fasta'


def bow_index(file, accession):     # builds bowtie index from accession fasta
    os.system('bowtie2-build ~/mini_proj/%s %s_ref' % (file, accession))
    return '%s_ref' % accession


def bow(ind, data):
    pe = []                         # groups paired end reads
    for i in data:
        for j in data:
            if i[:10] == j[:10] and i != j and i[-7] != '2':
                pe.append([i, j])
    os.system('cd ~/mini_proj/')
    for k in pe:                    # runs bowtie2 for each pair
        os.system('nohup bowtie2 -x ' + ind + ' -1 ~/mini_proj/transcr_data/' + k[0] + ' -2 ~/mini_proj/transcr_data/' +
                  k[1] + ' -S ' + k[0][:10] + '.sam' + ' --al-conc-gz ' + k[0][:10] + '_mapped.fq.gz')


# Driver
a = kb_map_accession()      # Accession provided in variables.py
f = get_fasta(a)            # Pulls fasta seq of accession

bi = bow_index(f, a)        # creates a bowtie index from fasta file

d = wanted_files(os.listdir('~/mini_proj/transcr_data/'))        # only the .fastq files
bow(bi, d)                  # runs bowtie
