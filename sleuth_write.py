import os


def sleuth_out(file):
    os.system('Rscript ' + file)        # runs sleuth Rscript


def sleuth_write(file):                 # writes the results of sleuth stat analysis to .log
    i = open(file, 'r')
    i_s = i.read().rstrip().split('\n')
    i_s_s = []
    for j in i_s:
        i_s_s.append(j.split(' '))
    print(i_s_s)
    o = open('miniProject.log', 'a')
    o.write('\n\nSLEUTH Outputs:\n')
    for k in i_s_s:
        for m in range(len(k)):
            if m == 0 and k[m] != 'target_id':
                o.write(str(k[m][4:]) + '\t')
            elif m == len(k)-1:
                o.write(str(k[m]))
                o.write('\n')
            else:
                o.write(str(k[m]) + '\t')
    o.close()


# Driver
f = '~/mini_proj/sleuth.R'
sleuth_out(f)               # sleuth R analysis

a = '~/mini_proj/sleuth_table.txt'
sleuth_write(a)             # writes results to .log file
