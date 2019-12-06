#/bin/env/python3
import sys
import os
in_dir = sys.argv[1]
out_dir = sys.argv[2]
for file in os.listdir(in_dir):
    print(file)
    if file.endswith("ins.fasta") or file.endswith("bb.fasta"):
        basename = file.split(".")[0]
        fa = ""
        with open(f"{in_dir}/{file}", "r") as f:
            fasta = f.readlines()
            for line in fasta:
                if line.startswith(">"):
                    readname = line.split(" ")[0]
                    idx = line.split("c_idx=")[1]
                    fa= fa + f"{readname}_{idx}"
                else:
                    fa= fa + line
        with open(f"{out_dir}/{basename}.fasta", "w") as f:
            f.write(fa)

