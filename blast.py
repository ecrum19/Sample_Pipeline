import os
from variables import *


def longest_scaffold(file):     # parses scaffold.fasta file and writes new file containing the longest scaffold
    def write(args):
        o = open('~/mini_proj/spades_out/long_scaf.fasta', 'w')
        for j in args:
            o.write(j+'\n')
    f = open('~/mini_proj/spades_out/' + file, 'r')
    f_s = f.read().strip().split('\n')
    long = 0
    startp = 0
    endp = 0
    all_pos = []
    long_pos = 0
    pos = 0
    for i in f_s:
        if '_length_' in i:
            if pos == 0:
                r = i.split('_')
                if int(r[3]) > long:
                    long = int(r[3])
                    long_pos = pos
            else:
                endp = pos-1
                all_pos.append([startp, endp])
                startp = pos
                r = i.split('_')
                if int(r[3]) > long:
                    long = int(r[3])
                    long_pos = pos
        pos += 1
    for j in all_pos:
        if j[0] == long_pos:
            startp = long_pos
            endp = j[1]
    write(f_s[startp:endp+1])


def run_blast(dbf, q, o):           # runs local BLAST using longest scaffold
    os.system('makeblastdb -in ~/mini_proj/%s -out %s -title %s -dbtype nucl' % (dbf, dbf[:17], dbf[:17]))
    os.system(
        'blastn -query %s -db %s -out %s.csv -outfmt "10 sseqid pident length qstart qend sstart send bitscore evalue stitle"' % (q, dbf[:17], o))


def write_blast(outfile):           # writes top 10 BLAST hits to .log
    o = open('miniProject.log', 'a')
    c = open('/homes/ecrum/mini_proj/%s' % outfile)
    c_s = c.read().strip().split('\n')
    header = ['sacc', 'pident', 'length', 'qstart', 'qend', 'sstart', 'send', 'bitscore', 'evalue', 'stitle']
    o.write('\n\nBLAST+ Output:\n')
    pos = 0
    for i in header:
        if pos == len(header)-1:
            o.write(i + '\n')
        else:
            o.write(i + '\t')
        pos += 1
    for j in range(10):
        v = c_s[j].split(',')
        pos = 0
        for k in v:
            if pos == len(v)-1:
                o.write(k + '\n')
            else:
                o.write(k + '\t')
            pos += 1


# Driver
longest_scaffold('scaffolds.fasta')     # finds longest scaffold

db = BLAST_subject_db()
query = '~/mini_proj/spades_out/long_scaf.fasta'
out = 'blast_out'
run_blast(db, query, out)               # runs blast Query: longest scaffold, Subject: Herpes RefSeq files
write_blast('blast_out.csv')     # writes blast outputs to .log
