from variables import *


def table(conds_dict, paths):               # makes table of results for sleuth analysis
    o = open('sleuth_table.txt', 'w')
    count = 0
    for a in conds_dict:
        o.write('%s %s  %s\n' % (str(a), str(conds_dict[a]), paths[count]))
        count += 1


# Driver
# dictionary denoting the transcriptome sample + time of data collection
conds = transcr_data_col()
paths = ['path']

for key in conds:               # makes list to produce table
    if key != 'sample':
        paths.append('~/mini_proj/kall_out/'+key)

table(conds, paths)             # makes file with table
