import os
import argparse
from Bio import Entrez, SeqIO

# Make sure all files in github repo are in current directory before beginning
os.system('cd ~/miniProject_elias_crum')

# supplies the -t flag in terminal for a test run
run = argparse.ArgumentParser(description='Run test data or real data.')
run.add_argument('-t', action='store_true',
                 help='"-t" flag uses test data instead of full dataset (default: full dataset)')
args = run.parse_args()
test_or_full = vars(args)['t']      # determines if the run is a test run or not


'''
Fill out the variables here with data if run is not with supplied data
'''
# edit sra addresses her
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
    transcriptome_conditions.update({'test_SRR56.1': '2dpi'})
    transcriptome_conditions.update({'test_SRR56.2': '6dpi'})
else:
    transcriptome_conditions.update({'sample': 'condition'})
    transcriptome_conditions.update({'SRR5660030': '2dpi'})  # fill in time aspects of data here
    transcriptome_conditions.update({'SRR5660033': '6dpi'})
    transcriptome_conditions.update({'SRR5660044': '2dpi'})
    transcriptome_conditions.update({'SRR5660045': '6dpi'})

# list of all time and donor conditions of transcriptome data
if test_or_full:                  # if '-t' flag present --> test data used
    conditions = ['Donor 1 (2dpi)', 'Donor 1 (6dpi)']
else:
    conditions = ['Donor 1 (2dpi)', 'Donor 1 (6dpi)','Donor 2 (2dpi)','Donor 2 (6dpi)']

# name of the file used to produce the BLAST database
blast_db = 'Betaherpesvirinae_ref.fasta'


cwd = os.getcwd()


def pull_files():       # pulls transcriptome + reference CDS map data
    def get_seqs(seq_addresses):  # downloads seqs via SRR addresses
        for i in seq_addresses:
            os.system("wget " + i)

    def new_dir(addresses, cd):  # moves transcriptome files to transcr_data dir
        if not test_or_full:
            for i in addresses:
                c = i.split('/')
                os.system('mv ' + cd + '/' + c[-1] + ' ' + cd + '/transcr_data/')

    # Driver
    # internet addresses for the transcriptome data (change addresses for different transcriptome data analysis)
    if not test_or_full:
        get_seqs(transcript_addresses)
        new_dir(transcript_addresses, cwd)  # makes new directory for the transcriptome files
pull_files()


def pull_map(accession):  # gets CDS info from Entrez
    def output_log(name, argus):  # writes CDS info to .log file
        o = open('miniProject.log', 'a')
        o.write('%s Fasta:\n' % name)
        o.write(argus)
        o.close()
    Entrez.email = "ecrum@luc.edu"
    net_handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb")  # fetches GenBank record
    out_handle = open(str(accession) + '.fasta', "w")  # writes record to local file
    recs = SeqIO.parse(net_handle, 'gb')
    for r in recs:
        output_log(accession, 'The HCMV genome (EF999921) has %s CDS.' % str(len(r.features)))  # for q2
    out_handle.close()
    net_handle.close()
pull_map(map_accession)


def split_fastq():        # splits fastq files and cleans directory
    def split_fq(files):  # splits transcriptome files into fwd and rev fastq reads
        for i in files:
            os.system('fastq-dump -I --split-files ' + cwd + '/transcr_data/%s' % i)

    def clean(cd, files):  # cleans dir, moving fastq files to the transcr_data dir
        for j in files:
            if '.fastq' in j:
                os.system('mv ' + cd + '/' + j + ' %s/transcr_data/' % cd)

    # Driver
    if not test_or_full:
        split_fq(os.listdir(cwd + '/transcr_data'))  # extracts .fastq files
    clean(cwd, os.listdir(cwd))  # moves paired end .fastq files to transr_data dir
split_fastq()


def wanted_files(cd):  # limits input files to only paired end fastq reads or test reads
    w_f = os.listdir(cd + '/transcr_data')
    nwf = w_f.copy()
    for n in w_f:
        if '.fastq' not in n or ('test_' not in n and test_or_full) or ('test_' in n and not test_or_full):
            nwf.remove(n)
    return nwf


def run_kallisto():        # analyzes transcriptome data using kallisto
    def prep_index(index):  # creates kallisto index
        a = index[:-4]
        os.system('time kallisto index -i %s_index.idx %s' % (a, index))
        return '%s_index.idx' % a

    def kallisto(index, data, cd):      # runs kallisto
        def kcommand(i, d, m):  # concatenates the command passed to kallisto
            if test_or_full:
                os.system('kallisto quant -i '+cwd+'/'+i+' -o ' + cwd+'/kallisto_out/'+d[0][:10]+'.'+str(m) +
                          ' -b 30 -t 2 '+cwd+'/transcr_data/'+d[0]+' '+cwd+'/transcr_data/' + d[1])
            else:
                os.system('kallisto quant -i '+cwd+'/'+i+' -o ' + cwd+'/kallisto_out/'+d[0][:10]+
                          ' -b 30 -t 2 '+cwd+'/transcr_data/'+d[0]+' '+cwd+'/transcr_data/' + d[1])
        pe = []
        for i in data:  # groups paired end reads
            for j in data:
                if not test_or_full:
                    if i[:10] == j[:10] and i != j and i[-7] != '2':
                        pe.append([i, j])
                else:
                    if i[:15] == j[:15] and i != j and i[-7] != '2':        # accounts for test_ in file name or not
                        pe.append([i, j])
        os.system('mkdir ' + cd + '/kallisto_out')
        num = 1
        for k in pe:  # runs kallisto for each of the paired end read pairs
            kcommand(index, k, num)
            num += 1

    # Driver
    transcriptomes = wanted_files(cwd)  # directory containing the transcriptomes (returns only fastq file names)
    ref = kallisto_map_file  # name of the reference CDS fasta file downloaded from NCBI
    ind = prep_index(ref)  # creates kallisto index
    kallisto(ind, transcriptomes, cwd)  # runs kallisto by calling index and .fastq files
run_kallisto()


def sleuth_tab():           # produces table for sleuth
    def table(conds_dict, pths):  # makes table of results for sleuth analysis
        o = open('sleuth_table.txt', 'w')
        count = 0
        for a in conds_dict:
            o.write('%s %s %s\n' % (str(a), str(conds_dict[a]), pths[count]))
            count += 1

    # Driver
    # dictionary denoting the transcriptome sample + time of data collection
    conds = transcriptome_conditions
    paths = ['path']
    for key in conds:  # makes list to produce table
        if key != 'sample':
            paths.append(cwd + '/kallisto_out/' + key)
    table(conds, paths)  # makes file with table
sleuth_tab()


def sleuth_write():         # writes sleuth output to .log
    def sleuth_out(file):
        os.system('Rscript ' + file)  # runs sleuth Rscript

    def sleuth_rite(f):  # writes the results of sleuth stat analysis to .log
        i = open(f, 'r')
        i_s = i.read().rstrip().split('\n')
        i_s_s = []
        for j in i_s:
            i_s_s.append(j.split(' '))
        o = open('miniProject.log', 'a')
        o.write('\n\nSLEUTH Outputs:\n')
        for k in i_s_s:
            for m in range(len(k)):
                if m == 0 and k[m] != 'target_id':
                    o.write(str(k[m][4:]) + '\t')
                elif m == len(k) - 1:
                    o.write(str(k[m]))
                    o.write('\n')
                else:
                    o.write(str(k[m]) + '\t')
        o.close()
    # Driver
    fi = cwd + '/sleuth.R'
    sleuth_out(fi)  # sleuth R analysis

    a = cwd + '/EF999921_sleuth_results.txt'
    sleuth_rite(a)  # writes results to .log file
sleuth_write()


def paird_end(inp):
    pe = []  # groups paired end reads
    for i in inp:
        for j in inp:
            if test_or_full:
                if i[:15] == j[:15] and i != j and i[-7] != '2':
                    pe.append([i, j])
            else:
                if i[:10] == j[:10] and i != j and i[-7] != '2':
                    pe.append([i, j])
    return pe


def mapped_wf(file_lst):  # selects bowtie mapped files from directory
    wf = file_lst.copy()
    for h in file_lst:
        if 'mapped' not in h or (test_or_full and 'test_' not in h) or (not test_or_full and 'test_' in h):
            wf.remove(h)
    return wf


def bowtie():           # runs bowtie seq mapping
    def get_fasta(accession):  # writes fasta file from accession provided
        Entrez.email = "ecrum@luc.edu"
        net_handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta")  # fetches GenBank record
        out_handle = open(str(accession) + '.fasta', "w")  # writes record to local file
        recs = SeqIO.parse(net_handle, 'fasta')
        for r in recs:
            SeqIO.write(r, out_handle, 'fasta')
        return accession + '.fasta'

    def bow_index(file, accession):  # builds bowtie index from accession fasta
        os.system('bowtie2-build '+cwd+'/'+file+' %s_ref' % accession)
        return '%s_ref' % accession

    def bow(ind, data):
        p = paird_end(data)
        for k in p:  # runs bowtie2 for each pair
            os.system(
                'nohup bowtie2 -x ' + ind + ' -1 '+cwd+'/transcr_data/' + k[0] + ' -2 '+cwd+'/transcr_data/' +
                k[1] + ' -S ' + k[0][:-6] + '.sam' + ' --al-conc-gz ' + k[0][:-6] + '_mapped.fq.gz')

    def rewrite(files):
        for s in files:
            if s[-17:] == '_1_mapped.fq.2.gz':
                ns = s.replace('_1_mapped.fq.2.gz', '_2_mapped.fq.gz')
            elif s[-17:] == '_1_mapped.fq.1.gz':
                ns = s.replace('_1_mapped.fq.1.gz', '_1_mapped.fq.gz')
            os.system('mv ' + cwd + '/' + s + ' ' + cwd + '/' + ns)


    # Driver
    f = get_fasta(map_accession)  # Pulls fasta seq of accession
    bi = bow_index(f, map_accession)  # creates a bowtie index from fasta file
    d = wanted_files(cwd)  # only the .fastq files
    bow(bi, d)  # runs bowtie
    rewrite(mapped_wf(os.listdir(cwd))) # rewites file so that it works with spades
bowtie()


def bowtie_charact():           # writes number of reads before and after mapping to .log
    def write(bnumbers, anumbers, samps):  # writes output to .log file
        o = open('miniProject.log', 'a')
        pos = 0
        o.write('\nBowtie Results:\n')
        for i in samps:
            o.write('%s had %s read pairs before Bowtie2 filtering and %s read pairs after.\n'
                    % (i, bnumbers[pos], anumbers[pos]))
            pos += 1

    def unzip_files(gzfiles):  # unzips mapped .gz files produced from bowtie
        for g in gzfiles:
            if '.gz' in g:
                os.system('gzip -d %s' % g)

    def calc_nums(qfiles, mfiles):  # calculates number of reads based on file length / 4 (because of fastq fmt)
        bnums = []
        fbnums = []
        anums = []
        fanums = []
        cb = 1
        ca = 1
        for j in sorted(qfiles):
            binp = open(cwd + '/transcr_data/' + j, 'r')
            binp_s = binp.read().strip().split('\n')
            bnums.append(len(binp_s) / 4)
            if cb % 2 == 0:
                if bnums[cb-2] == bnums[cb-1]:
                    fbnums.append(bnums[cb-2])
                else:
                    fbnums.append(min(bnums[cb-2],bnums[cb-1]))
            cb += 1
        for k in sorted(mfiles):
            ainp = open(cwd + '/' + k, 'r')
            ainp_s = ainp.read().strip().split('\n')
            anums.append(len(ainp_s) / 4)
            if ca % 2 == 0:
                newa = anums[ca-2] + anums[ca-1]
                fanums.append(newa)
            ca += 1
        return fbnums, anums

    # Driver
    trans = wanted_files(cwd)  # transcriptome fastq files (before bowtie)
    unzip_files(mapped_wf(os.listdir(cwd)))  # unzip mapped bowtie files
    mapped_f = mapped_wf(os.listdir(cwd))  # CHECK                # unzipped mapped bowtie files
    b, c = calc_nums(trans, mapped_f)  # count number before + after bowtie mapping
    write(b, c, conditions)  # writes numbers to .log file


def spades():           # assemples a de novo genome from the mapped bowtie transcriptome reads
    def run_spades(files):  # runs spades and writes command used to .log
        def write(argas):
            o = open('miniProject.log', 'a')
            o.write('\n\nSPAdes Command:\n' + argas + '\n')

        ge = []  # groups paired end reads
        for i in files:
            for j in files:
                if test_or_full:
                    if i[:15] == j[:15] and i != j and i[-14] != '2':
                        ge.append([i, j])
                else:
                    if i[:12] == j[:12] and i != j and i[-11] != '2':
                        ge.append([i, j])
        command_string = 'nohup spades -k 55,127 -t 2 --only-assembler'
        count = 1
        for i in ge:
            command_string += ' --pe' + str(count) + '-1 ' + cwd + '/' + i[0] + \
                              ' --pe' + str(count) + '-2 ' + cwd + '/' + i[1]
            count += 1
        command_string += ' -o ' + cwd + '/spades_out'
        write(command_string)
        os.system(command_string)

    # Driver
    f = mapped_wf(os.listdir(cwd))  # only want mapped .fastq files
    run_spades(f)  # runs spades with the mapped .fastq files
spades()
bowtie_charact()


def contig_len():           # determines the number of contigs >1000 bp obtained from spades assembly
    def num_big(file):  # writes number of contigs > 1000 bp to .log
        def write(arges):
            o = open('miniProject.log', 'a')
            o.write('\n\nSPAdes Contigs:\n' +
                    'There are %d contigs > 1000 bp in the assembly.\n' % arges)
        f = open(cwd + '/spades_out/' + file, 'r')
        f_s = f.read().strip().split('\n')
        num = 0
        for i in f_s:
            if '_length_' in i:
                r = i.split('_')
                if int(r[3]) >= 1000:
                    num += 1
        write(num)

    # Driver
    num_big('contigs.fasta')  # determines number of contigs > 1000 bp
contig_len()


def assembly_len():         # determines the length of assembled genome with >1000 bp scaffolds
    def tot_bp(file):  # prints length of assembly seq of contigs > 1000 bp
        def write(argqs):
            o = open('miniProject.log', 'a')
            o.write('\n\nSPAdes Scaffold Length:\n' +
                    'There are %d bp in the assembly.\n' % argqs)

        f = open(cwd + '/spades_out/' + file, 'r')
        f_s = f.read().strip().split('\n')
        num = 0
        for i in f_s:
            if '_length_' in i:
                r = i.split('_')
                if int(r[3]) >= 1000:
                    num += int(r[3])
        write(num)

    # Driver
    tot_bp('scaffolds.fasta')  # executes function
assembly_len()


def blast():            # runs local BLAST with longest scaffold from assembly, writes outputs to .log
    def longest_scaffold(file):  # parses scaffold.fasta file and writes new file containing the longest scaffold
        def write(args):
            o = open(cwd + '/spades_out/long_scaf.fasta', 'w')
            for j in args:
                o.write(j + '\n')

        f = open(cwd + '/spades_out/' + file, 'r')
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
                    endp = pos - 1
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
        write(f_s[startp:endp + 1])

    def run_blast(dbf, q, o):  # runs local BLAST using longest scaffold
        os.system('makeblastdb -in '+cwd+'/%s -out %s -title %s -dbtype nucl' % (dbf, dbf[:17], dbf[:17]))
        os.system(
            'blastn -query %s -db %s -out %s.csv -outfmt "10 sseqid pident length qstart qend sstart send bitscore evalue stitle"' % (q, dbf[:17], o))

    def write_blast(outfile):  # writes top 10 BLAST hits to .log
        o = open('miniProject.log', 'a')
        c = open(cwd + '/' + outfile)
        c_s = c.read().strip().split('\n')
        header = ['sacc', 'pident', 'length', 'qstart', 'qend', 'sstart', 'send', 'bitscore', 'evalue', 'stitle']
        o.write('\n\nBLAST+ Output:\n')
        pos = 0
        for i in header:
            if pos == len(header) - 1:
                o.write(i + '\n')
            else:
                o.write(i + '\t')
            pos += 1
        if test_or_full:
            w = 2
        else:
            w = 10
        for j in range(w):
            v = c_s[j].split(',')
            pos = 0
            for k in v:
                if pos == len(v) - 1:
                    o.write(k + '\n')
                else:
                    o.write(k + '\t')
                pos += 1

    # Driver
    longest_scaffold('scaffolds.fasta')  # finds longest scaffold

    db = blast_db
    query = cwd + '/spades_out/long_scaf.fasta'
    out = 'blast_out'
    run_blast(db, query, out)  # runs blast Query: longest scaffold, Subject: Herpes RefSeq files
    write_blast('blast_out.csv')  # writes blast outputs to .log
blast()


'''
Output results are summarized in the miniProject.log file
'''
