#/bin/python3
import os
import sys

# Concat gzipped fastq
path = sys.argv[1] # input dir
out_path = sys.argv[2] # output dir: must be not relative path but exact path since chdir happens in the code.
batch_size = sys.argv[3] # 10
name = sys.argv[4] # P260

os.makedirs(out_path)
os.chdir(path)

space_sep_10_filenames = []
x = 0
batch = ""
for i, file in enumerate(os.listdir(os.getcwd())):
    y = i //int(batch_size)
    if y == x :
        batch += file + " "
    else:
        space_sep_10_filenames.append(batch)
        x +=1
        batch = file + " "

        
for j in range(i//int(batch_size)):
    cmd = f"cat {space_sep_10_filenames[j]} > {out_path}/{name}_b{j}.fastq.gz"
    os.system(cmd)


