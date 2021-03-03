def num_big(file):          # writes number of contigs > 1000 bp to .log
    def write(args):
        o = open('miniProject.log', 'a')
        o.write('\n\nSPAdes Contigs:\n' +
                'There are %d contigs > 1000 bp in the assembly.\n' % args)
    f = open('/homes/ecrum/mini_proj/spades_out/' + file, 'r')
    f_s = f.read().strip().split('\n')
    num = 0
    for i in f_s:
        if '_length_' in i:
            r = i.split('_')
            if int(r[3]) >= 1000:
                num += 1
    write(num)


# Driver
num_big('contigs.fasta')        # determines number of contigs > 1000 bp
