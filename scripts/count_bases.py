
import sys
import os
fasta_file = sys.argv[1]
#out_dir = sys.argv[2]
        
basename = fasta_file.split(".")[0]
fa = ""
length_list = []
with open(f"{fasta_file}", "r") as f:
    fasta = f.readlines()
    for line in fasta:
        if line.startswith(">"):
            readname = line.split(" ")[0]
            fa= fa + f"{readname}\t"
        else:
            fa = fa + f"{len(line)}\n"
            length_list.append(len(line))
with open(f"{basename}_length_count.txt", "w") as f:
    f.write(fa)
print(min(length_list))
