import os
from variables import *
from run_kallisto import *


def write(bnumbers, anumbers, samps):  # writes output to .log file
    o = open('miniProject.log', 'a')
    pos = 0
    o.write('\nBowtie Results:\n')
    for i in samps:
        o.write('%s had %s read pairs before Bowtie2 filtering and %s read pairs after.\n'
                % (i, bnumbers[pos], anumbers[pos]))
        pos += 1


def mapped_wf(file_lst):                # selects bowtie mapped files from directory
    wf = file_lst.copy()
    for h in file_lst:
        if '_mapped_' not in h:
            wf.remove(h)
        if test_or_full() and 'test_' not in h:
            wf.remove(h)
    return wf


def unzip_files(gzfiles):               # unzips mapped .gz files produced from bowtie
    gzfs = os.listdir('~/mini_proj/%s' % gzfiles)
    for g in gzfs:
        os.system('gzip -d %s' % g)


def calc_nums(qfiles, mfiles):          # calculates number of reads based on file length / 4 (because of fastq fmt)
    bnums = []
    anums = []
    for j in qfiles:
        binp = open('~/mini_proj/transcr_data/'+j, 'r')
        binp_s = binp.read().strip().split('\n')
        bnums.append(len(binp_s) / 4)
    for k in mfiles:
        ainp = open('~/mini_proj/'+k, 'r')
        ainp_s = ainp.read().strip().split('\n')
        anums.append(len(ainp_s) / 4)
    return bnums, anums


# Driver
conds = poss_conditions()           # conditions from variables.py
trans = wanted_files(os.listdir('~/mini_proj/transcr_data'))        # transcriptome fastq files (before bowtie)
unzip_files(mapped_wf(os.listdir('~/mini_proj')))                   # unzip mapped bowtie files
mapped_f = mapped_wf(os.listdir('~/mini_proj'))      # CHECK                # unzipped mapped bowtie files
b, c = calc_nums(trans, mapped_f)                               # count number before + after bowtie mapping
write(b, c, conds)                                              # writes nnumbers to .log file


