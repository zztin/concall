#!/hpc/cog_bioinf/ridder/users/lchen/miniconda3/envs/medaka_env/bin/python

import os
import sys
import time

outpath = sys.argv[1]
txt_folder = sys.argv[2]
fasta_file_path = sys.argv[3]
os.mkdir(f"{outpath}")
ins_or_bb = outpath.split("_")[-1]
for n in range(1,40):
    cmd = f"fgrep -A 1 -f {txt_folder}/con{n}.txt {fasta_file_path} > {outpath}/consensus_{ins_or_bb}_{n}.fasta"
    os.system(cmd)
    time.sleep(0.2)

cmd = f"fgrep -A 1 -f {txt_folder}/con{n+1}more.txt {fasta_file_path} > {outpath}/consensus_{ins_or_bb}_{n+1}.fasta"
os.system(cmd)

