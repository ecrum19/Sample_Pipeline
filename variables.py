'''
Fill out the functions here with data if run is not with supplied data
'''
from mini_proj_wrapper import *


def transcriptome_sequences():          # edit sra addresses here
    transcript_addresses = ['https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660030/SRR5660030.1',
                            'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660033/SRR5660033.1',
                            'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660044/SRR5660044.1',
                            'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660045/SRR5660045.1']
    test_t = ['https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660030/SRR5660030.1']

    if test_or_full():                  # if '-t' flag present --> test data used
        return test_t
    else:
        return transcript_addresses


def kb_map_accession():           # edit record used for kallisto reference
    map_accession = 'EF999921'
    return map_accession


def kallisto_map_file():
    file_name = 'EF999921_CDS.txt'
    return file_name


def transcr_data_col():
    # dictionary that pairs transcriptome data and time collected
    transcriptome_conditions = {}
    transcriptome_conditions.update({'sample': 'condition'})
    transcriptome_conditions.update({'SRR5660030': '2dpi'})     # fill in data here
    transcriptome_conditions.update({'SRR5660033': '6dpi'})
    transcriptome_conditions.update({'SRR5660044': '2dpi'})
    transcriptome_conditions.update({'SRR5660045': '6dpi'})

    test_con = {}
    test_con.update({'sample': 'condition'})
    test_con.update({'SRR5660030': '2dpi'})

    if test_or_full():                  # if '-t' flag present --> test data used
        return test_con
    else:
        return transcriptome_conditions


def poss_conditions():
    # list of all time and donor conditions of transcriptome data
    conditions = ['Donor 1 (2dpi)', 'Donor 1 (6dpi)','Donor 2 (2dpi)','Donor 2 (6dpi)']
    test_c = ['Donor 1 (2dpi)']

    if test_or_full():                  # if '-t' flag present --> test data used
        return test_c
    else:
        return conditions


def BLAST_subject_db():
    # name of the file used to produce the BLAST database
    file = 'Betaherpesvirinae_ref.fasta'
    return file


def for_test(file):
    f = open('~/mini_proj/transcr_data/' + file, 'r')
    o = open('~/mini_proj/transcr_data/test_' + file, 'w')
    f_s = f.read().strip().split('\n')
    for i in range(40):
        o.write(f_s[i] + '\n')
    return 'test_%s' % file
