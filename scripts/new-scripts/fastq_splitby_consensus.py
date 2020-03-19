#!/hpc/cog_bioinf/ridder/users/lchen/miniconda3/bin/python
import pandas as pd
import sys
import numpy as np
import os

path = sys.argv[1]
csv = sys.argv[2]
ins_or_bb = path.split("_")[-1]

    os.mkdir(path+f"_{ins_or_bb}")

    mytable = pd.read_csv(csv, delimiter=',', names = ['readID',
                                                       'bb_count',
                                                       'reversed_bb_count',
                                                       'bb_len',
                                                       'ins_count',
                                                       'ins_len'])

    for n in range(1, 40):
        filename = f"con{n}.txt"
        readID = mytable[mytable[f'{ins_or_bb}_count'] == n]['readID'].values
        np.savetxt(f"{path}_{ins_or_bb}/{filename}", readID, fmt='%s' )


    filename = f"con{n+1}more.txt"
    readID = mytable[mytable[f'{ins_or_bb}_count'] > n]['readID'].values
    np.savetxt(f"{path}_{ins_or_bb}/{filename}", readID, fmt='%s' )

