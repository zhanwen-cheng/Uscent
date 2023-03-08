import sys
import getopt
import pandas as pd
import numpy as np
from collections import Counter
from Bio import SeqIO
from Bio import SeqUtils
import pandas as pd
import matplotlib.pyplot as plt

def usage():
    print('''
          Usage: python3 script.py [option] [patameter]
          -r/--depth          input reference fasta file
          -o/--outfile        input outfile name
          -h/--help           show possible options
          ''')


opts, args = getopt.getopt(sys.argv[1:], 'hr:o:', ['help', 'ref=', 'outfile='])
for opt, val in opts:
    if opt == '-r' or opt == '--ref':
        ref_file = val
    elif opt == '-o' or opt == '--outfile':
        out_file = val
    elif opt == '-h' or opt == '--help':
        usage()
        sys.exit(1)


genome=SeqIO.read(ref_file,"fasta")
genome_gc=[]
for i in range(0,len(genome.seq)-1000+1,1000):
    subseq=genome.seq[i:i+1000]
    genome_gc.append(SeqUtils.GC(subseq))
genome_gc=pd.DataFrame(genome_gc)
genome_gc.index=genome_gc.index+1
genome_gc.columns=['gc']
genome_gc.to_csv(out_file, index=False, sep=',')
