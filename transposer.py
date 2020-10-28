#https://stackoverflow.com/questions/54208323/how-can-i-efficiently-transpose-a-67-gb-file-dask-dataframe-without-loading-it-e

import argparse, subprocess
import numpy as np
import pandas as pd
import dask.dataframe as dd
from contextlib import closing
from os import cpu_count
from multiprocessing import Pool

parser = argparse.ArgumentParser(description='Transpose csv')
parser.add_argument('-i', '--infile', help='Path to input folder',
                    default=None)
parser.add_argument('-s', '--sep', help='input separator',
                    default='\t')

args = parser.parse_args()
infile = args.infile
sep = args.sep    
df = pd.read_csv(infile, nrows=3)    

def READ_COL(item):
    print(item)
    outfile = 'outfile{}.temp'.format(item)
    if item !=0:
        x = "awk -F ',' '{print $"+str(item)+"}' "+infile+" > "+outfile
        subprocess.check_call([x], shell=True)
        col = pd.read_csv(outfile)
        row = col.transpose()
        out = 'col_{:09d}.csv'.format(item)
        row.to_csv(out, header=False)
        subprocess.check_call(['rm '+outfile], shell=True)
    
with closing(Pool(processes=cpu_count())) as pool:
    pool.map(READ_COL, list(range(1, len(df.columns)+1)))