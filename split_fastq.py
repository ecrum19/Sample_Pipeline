import os


def split_fastq(files):             # moves transcriptome files to new directory
    os.system('mkdir /homes/ecrum/mini_proj/transcr_data')
    for i in files:
        os.system('mv /homes/ecrum/mini_proj/' + i + ' /homes/ecrum/mini_proj/transcr_data/' + i)
    # Run 'fastq-dump -I --split-files transcrfile' for each of the four files

    # (when used in a script the error: 'Failed to call external services.' is thrown)


# Driver
transcriptomes = ['SRR5660033.1', 'SRR5660045.1', 'SRR5660030.1', 'SRR5660044.1']
split_fastq(transcriptomes)
