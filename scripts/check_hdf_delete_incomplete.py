#!/bin/python3
import os
import sys
import shutil
#in_path = "/Users/lchen/Projects/sam_add_tag/raw_data/05_aggregate/"
in_path = sys.argv[1]
folder_to_delete = []
for root, dirs, files in os.walk(in_path):
#    print(root, dirs, files)
    if "consensus.hdf" in files and "consensus.fasta" not in files:
        folder_to_delete.append(root)
for path in folder_to_delete:
#    pass
    shutil.rmtree(path)
print(f"{len(folder_to_delete)} incomplete reads deleted\n" , folder_to_delete)
