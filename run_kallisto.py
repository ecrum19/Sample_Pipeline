import os
from variables import *
from mini_proj_wrapper import *


def wanted_files(dir):                        # limits input files to only paired end fastq reads
    w_f = os.listdir('/homes/ecrum/mini_proj/%s' % dir)
    nwf = w_f.copy()
    for n in w_f:
        if '.fastq' not in n:
            nwf.remove(n)
        if test_or_full() and 'test_' not in n:     # selects test or full data
            nwf.remove(n)
    return nwf


def prep_index(index):                          # creates kallisto index
    a = index[:-4]
    os.system('time kallisto index -i %s_index.idx %s' % (a, index))
    return '%s_index.idx' % a


def kallisto(index, data):
    def kcommand(i, d):         # concatenates the command passed to kallisto
        os.system('kallisto quant -i ~/mini_proj/'+i+' -o ~/mini_proj/kall_out/'+d[0][:10]+' -b 30 -t 2 ~/mini_proj/transcr_data/'+d[0]+' ~/mini_proj/transcr_data/'+d[1])
    pe = []
    for i in data:                      # groups paired end reads
        for j in data:
            if i[:10] == j[:10] and i != j and i[-7] != '2':
                pe.append([i, j])
    os.system('mkdir ~/mini_proj/kall_out/')
    for k in pe:                        # runs kallisto for each of the paired end read pairs
        kcommand(index, k)


# Driver
transcriptomes = wanted_files('transcr_data')       # directory containing the transcriptomes (returns only fastq files)

ref = kallisto_map_file()                 # name of the reference CDS fasta file downloaded from NCBI
ind = prep_index(ref)                     # creates kallisto index

kallisto(ind, transcriptomes)       # runs kallisto by calling index and .fastq files
