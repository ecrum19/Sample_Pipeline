import os
from Bio import SeqIO
from bowtie_charact import *


def run_spades(files):      # runs spades and writes command used to .log
    def write(args):
        o = open('miniProject.log', 'a')
        o.write('\n\nSPAdes Command:\n' + args + '\n')
    pe = []
    for i in files:
        for j in files:
            if i[:10] == j[:10] and i != j and i[-4] != '2':
                pe.append([i, j])  # groups paired end reads
    command_string = 'nohup spades -k 127 -t 2 --only-assembler'
    count = 1
    for i in pe:
        command_string += ' --pe' + str(count) + '-1 ~/mini_proj/' + i[0] + \
                         ' --pe' + str(count) + '-2 ~/mini_proj/'+ i[1]
        count += 1
    command_string += ' -o ~/mini_proj/spades_out'
    write(command_string)
    os.system(command_string)


# Driver
f = mapped_f(os.listdir('~/mini_proj'))         # only want mapped .fastq files

run_spades(f)                           # runs spades with the mapped .fastq files
