import os
from Bio import SeqIO, Entrez


def get_seqs(seq_addresses):        # downloads seqs via SRR addresses
    for i in seq_addresses:
        os.system("wget " + i)


def pull_map(accession):
    Entrez.email = "ecrum@luc.edu"
    net_handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb")      # fetches GenBank record
    out_handle = open(str(accession) + '.fasta', "w")                               # writes record to local file
    recs = SeqIO.parse(net_handle, 'gb')
    for r in recs:
        output_log('The HCMV genome (EF999921) has %s CDS.' % str(len(r.features)))    # for q2
    out_handle.close()
    net_handle.close()


def output_log(args):                   # writes output to .log file
    o = open('miniProject.log', 'a')
    o.write(args + '\n\n')
    o.close()


# Driver
# internet addresses for the transcriptome data
transcript_seqs = ['https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660030/SRR5660030.1',
                   'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660033/SRR5660033.1',
                   'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660044/SRR5660044.1',
                   'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660045/SRR5660045.1']

get_seqs(transcript_seqs)
pull_map(accession='EF999921')
