import os
import argparse
from Bio import Entrez, SeqIO

# Make sure all files in github repo are in current directory before beginning
os.system('cd ~/mini_proj')

# supplies the -t flag in terminal for a test run
run = argparse.ArgumentParser(description='Run test data or real data.')
run.add_argument('-t', action='store_true',
                 help='"-t" flag uses test data instead of full dataset (default: full dataset)')
args = run.parse_args()
test_or_full = vars(args)['t']

'''
Fill out the variables here with data if run is not with supplied data
'''
# edit sra addresses her
if test_or_full:  # if '-t' flag present --> test data used
    transcript_addresses = ['https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660030/SRR5660030.1']
else:
    transcript_addresses = ['https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660030/SRR5660030.1',
                            'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660033/SRR5660033.1',
                            'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660044/SRR5660044.1',
                            'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660045/SRR5660045.1']


# edit record used for kallisto reference
map_accession = 'EF999921'

# edit record used for kallisto map
kallisto_map_file = 'EF999921_CDS.txt'

# dictionary that pairs transcriptome data and time collected
transcriptome_conditions = {}
if test_or_full:                  # if '-t' flag present --> test data used
    transcriptome_conditions.update({'sample': 'condition'})
    transcriptome_conditions.update({'SRR5660030': '2dpi'})
else:
    transcriptome_conditions = {}
    transcriptome_conditions.update({'sample': 'condition'})
    transcriptome_conditions.update({'SRR5660030': '2dpi'})  # fill in data here
    transcriptome_conditions.update({'SRR5660033': '6dpi'})
    transcriptome_conditions.update({'SRR5660044': '2dpi'})
    transcriptome_conditions.update({'SRR5660045': '6dpi'})

# list of all time and donor conditions of transcriptome data
if test_or_full:                  # if '-t' flag present --> test data used
    conditions = ['Donor 1 (2dpi)']
else:
    conditions = ['Donor 1 (2dpi)', 'Donor 1 (6dpi)','Donor 2 (2dpi)','Donor 2 (6dpi)']

# name of the file used to produce the BLAST database
blast_db = 'Betaherpesvirinae_ref.fasta'


def for_test(file):
    f = open('~/mini_proj/transcr_data/' + file, 'r')
    o = open('~/mini_proj/transcr_data/test_' + file, 'w')
    f_s = f.read().strip().split('\n')
    for i in range(40):
        o.write(f_s[i] + '\n')
    return 'test_%s' % file


def pull_files():       # pulls transcriptome + reference CDS map data
    def get_seqs(seq_addresses):  # downloads seqs via SRR addresses
        for i in seq_addresses:
            os.system("wget " + i)

    def pull_map(accession):
        Entrez.email = "ecrum@luc.edu"
        net_handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb")  # fetches GenBank record
        out_handle = open(str(accession) + '.fasta', "w")  # writes record to local file
        recs = SeqIO.parse(net_handle, 'gb')
        for r in recs:
            output_log(accession, 'The HCMV genome (EF999921) has %s CDS.' % str(len(r.features)))  # for q2
        out_handle.close()
        net_handle.close()

    def output_log(name, argus):  # writes output to .log file
        o = open('miniProject.log', 'a')
        o.write('%s Fasta:\n' % name)
        o.write(argus)
        o.close()

    def new_dir(addresses):  # makes a directory for transcriptome files
        os.system('mkdir /homes/ecrum/mini_proj/transcr_data')
        for i in addresses:
            c = i.split('/')
            os.system('mv /homes/ecrum/mini_proj/%s /homes/ecrum/mini_proj/transcr_data/' % c[-1])

    # Driver
    # internet addresses for the transcriptome data (change addresses for different transcriptome data analysis)
    get_seqs(transcript_addresses)
    new_dir(transcript_addresses)  # makes new directory for the transcriptome files

    # accession used for the Bowtie2 reference genome (change accession for different map)
    pull_map(accession=map_accession)
#pull_files()


def split_fastq():
    def split_fq(files):  # splits transcriptome files into fwd and rev fastq reads
        for i in files:
            os.system('fastq-dump -I --split-files ~/mini_proj/transcr_data/%s' % i)

    def clean(files):  # cleans dir, moving fastq files to the transcr_data dir
        for j in files:
            if '.fastq' in j:
                os.system('mv /homes/ecrum/mini_proj/%s ~/mini_proj/transcr_data/' % j)
                if test_or_full:
                    for_test(j)

    # Driver
    split_fq(os.listdir('/mini_proj/transcr_data'))  # extracts .fastq files
    clean(os.listdir('/mini_proj'))  # moves paired end .fastq files to transr_data dir
split_fastq()


#os.system('python3 split_fastq.py')     # splits fastq files and cleans directory
#os.system('python3 run_kallisto.py')    # analyzes transcriptome data using kallisto
#os.system('python3 sleuth_tab.py')      # produces table for sleuth
#os.system('python3 sleuth_write.py')    # writes sleuth output to .log
#os.system('python3 bowtie.py')          # runs bowtie seq mapping
#os.system('python3 bowtie_charact.py')  # writes number of reads before and after mapping to .log
#os.system('python3 spades.py')          # assemples a de novo genome from the mapped bowtie transcriptome reads
#os.system('python3 contig_len.py')      # determines the number of contigs >1000 bp obtained from spades assembly
#os.system('python3 assembly_len.py')    # determines the length of assembled genome with >1000 bp scaffolds
#os.system('python3 blast.py')           # runs local BLAST with longest scaffold from assembly, writes outputs to .log


