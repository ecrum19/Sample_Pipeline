import os
from mini_proj_wrapper import *
from variables import *


def split_fastq(files):             # splits transcriptome files into fwd and rev fastq reads
    for i in files:
        os.system('fastq-dump -I --split-files ~/mini_proj/transcr_data/%s' % i)


def clean(files):       # cleans dir, moving fastq files to the transcr_data dir
    for j in files:
        if '.fastq' in j:
            os.system('mv /homes/ecrum/mini_proj/%s ~/mini_proj/transcr_data/' % j)
            if test_or_full():
                for_test(j)


# Driver
transcriptomes = os.listdir('~/mini_proj/transcr_data')
split_fastq(transcriptomes)                     # extracts .fastq files

clean(os.listdir('~/mini_proj'))                # moves paired end .fastq files to transr_data dir
