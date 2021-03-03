def tot_bp(file):           # prints length of assembly seq of contigs > 1000 bp
    def write(args):
        o = open('miniProject.log', 'a')
        o.write('\n\nSPAdes Scaffold Length:\n' +
                'There are %d bp in the assembly.\n' % args)
    f = open('/homes/ecrum/mini_proj/spades_out/' + file, 'r')
    f_s = f.read().strip().split('\n')
    num = 0
    for i in f_s:
        if '_length_' in i:
            r = i.split('_')
            if int(r[3]) >= 1000:
                num += int(r[3])
    write(num)


# Driver
tot_bp('scaffolds.fasta')       # executes function
