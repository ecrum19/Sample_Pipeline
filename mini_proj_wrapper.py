import os
import argparse


# Download EF999921 CDS fasta document here: https://www.ncbi.nlm.nih.gov/nuccore/EF999921.1/
# Download Betaherpesvirinae reference genomes here: https://www.ncbi.nlm.nih.gov/nuccore
#   (by searching Betaherpesvirinae[organism] and designating only RefSeq records)


def test_or_full():
    run = argparse.ArgumentParser(description='Run test data or real data.')
    run.add_argument('-t', action='store_true',
                 help='"-t" flag uses test data instead of full dataset (default: full dataset)')
    args = run.parse_args()
    return vars(args)['t']

os.system('cd ~/mini_proj/')
os.system('python3 pull_files.py')      # pulls transcriptome + reference CDS map data
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


