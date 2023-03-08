import sys
import getopt
import pandas as pd
import numpy as np
from collections import Counter

def usage():
    print('''
          Usage: python3 script.py [option] [patameter]
          -d/--depth          input depth and U counts file
          -o/--outfile        input outfile name
          -h/--help           show possible options
          ''')


opts, args = getopt.getopt(sys.argv[1:], 'hd:o:', ['help', 'depth=', 'outfile='])
for opt, val in opts:
    if opt == '-d' or opt == '--depth':
        depth_file = val
    elif opt == '-o' or opt == '--outfile':
        out_file = val
    elif opt == '-h' or opt == '--help':
        usage()
        sys.exit(1)

import numpy as np
import pandas as pd

df = pd.read_table(depth_file, sep=',')
df.loc[df.lable==1,'Ucounts']=1
df_new = df.groupby(df.index//1000).sum()
df_new.lable = df_new.index+1
df_new = df_new[['lable','Ucounts','Bcounts']]
df_new.to_csv(out_file, index=False, sep=',')
