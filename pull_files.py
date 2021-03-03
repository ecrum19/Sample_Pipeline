import os
from Bio import SeqIO, Entrez
from variables import *


def get_seqs(seq_addresses):        # downloads seqs via SRR addresses
    for i in seq_addresses:
        os.system("wget " + i)


def pull_map(accession):
    Entrez.email = "ecrum@luc.edu"
    net_handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb")      # fetches GenBank record
    out_handle = open(str(accession) + '.fasta', "w")                               # writes record to local file
    recs = SeqIO.parse(net_handle, 'gb')
    for r in recs:
        output_log(accession, 'The HCMV genome (EF999921) has %s CDS.' % str(len(r.features)))    # for q2
    out_handle.close()
    net_handle.close()


def output_log(name, args):                     # writes output to .log file
    o = open('miniProject.log', 'a')
    o.write('s% Fasta:\n' % name)
    o.write(args)
    o.close()


def new_dir(addresses):                         # makes a directory for transcriptome files
    os.system('mkdir /homes/ecrum/mini_proj/transcr_data')
    for i in addresses:
        c = i.split('/')
        os.system('mv /homes/ecrum/mini_proj/%s /homes/ecrum/mini_proj/transcr_data/' % c[-1])


# Driver
# internet addresses for the transcriptome data (change addresses for different transcriptome data analysis)
transcript_seqs = transcriptome_sequences()
get_seqs(transcript_seqs)
new_dir(transcript_seqs)        # makes new directory for the transcriptome files

# accession used for the Bowtie2 reference genome (change accession for different map)
map_accession = kb_map_accession()
pull_map(accession=map_accession)
