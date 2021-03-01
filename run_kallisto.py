import os


def clean_dir(files):                                           # cleans directory of unnecessary files
    all_f = os.listdir('/homes/ecrum/mini_proj/transcr_data/')
    for i in all_f:
        if i in files:
            os.system('rm ~/mini_proj/transcr_data/' + i)


def run_kallisto(index, data):                      # creates index
    n = index[:-4]
    os.system('time kallisto index -i ' + n + '_index.idx ' + index)
    pe = []
    for i in data:
        for j in data:
            if i[:10] == j[:10] and i != j:
                pe.append([i, j])                   # groups paired end reads

    os.system('cd /homes/ecrum/mini_proj/')
    for j in pe:                                    # runs kallisto
        os.system('time kallisto quant -i index/'+n+'_index.idx -o kall_results/'+j[0][:10]+' -b 30 -t 2 transcr_data/'+j[0]+' transcr_data/'+j[1])


# Driver
transcriptomes = ['SRR5660033.1', 'SRR5660045.1', 'SRR5660030.1', 'SRR5660044.1']
clean_dir(transcriptomes)

ref = 'EF999921_CDS.txt'
data = ['SRR5660033.1_1.fastq', 'SRR5660044.1_2.fastq']
d = os.listdir('/homes/ecrum/mini_proj/transcr_data/')
run_kallisto(ref, data)

