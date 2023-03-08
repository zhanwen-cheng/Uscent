import sys
import getopt
import pandas as pd
import numpy as np
from collections import Counter

def usage():
    print('''
          Usage: python3 script.py [option] [patameter]
          -p/--psl            input a psl file
          -d/--depth          input depth file
          -o/--outfile        input outfile name
          -h/--help           show possible options
          ''')


opts, args = getopt.getopt(sys.argv[1:], 'hp:d:o:', ['help', 'psl=', 'depth=', 'outfile='])
for opt, val in opts:
    if opt == '-p' or opt == '--psl':
        psl_file = val
    elif opt == '-d' or opt == '--depth':
        depth_file = val
    elif opt == '-o' or opt == '--outfile':
        out_file = val
    elif opt == '-h' or opt == '--help':
        usage()
        sys.exit(1)


psl = pd.read_table(psl_file, sep=" ",header=None, usecols=[9,13,14,20])#genome size 5199318
psl.columns = ['query','genome','genome_size','U_pos']

psl_list = []
for pos in psl.U_pos:
    psl_list.extend(pos.split(',')[:-1])#最后一个元素是,
U_pos = Counter(psl_list)
U_pos = pd.DataFrame(U_pos.items(),columns=['lable','Ucounts'])
U_pos['lable'] = U_pos[['lable']].astype(int)+1
U_pos = U_pos.sort_values(by = ['lable']).reset_index(drop=True)
all_pos = pd.DataFrame({'lable':range(1,(psl.loc[0,'genome_size']+1))})
all_pos = pd.merge(all_pos,U_pos,how='left')
all_pos.fillna(0,inplace=True)

depth = pd.read_table(depth_file,header=None, usecols=[1,2])
depth.columns = ['lable','Bcounts']
#depth['lable'] = depth['lable']
all_pos = pd.merge(all_pos,depth,how='left')#depth file is longer than genome
all_pos.fillna(0,inplace=True)
all_pos[['Ucounts','Bcounts']] = all_pos[['Ucounts','Bcounts']].astype(int)

all_pos.to_csv(out_file, index=False, sep=',')

