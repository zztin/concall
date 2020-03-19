#!/bin/python3
import sys
import pandas as pd
in_file = sys.argv[1]
out_fasta = sys.argv[2]
out_pickle = sys.argv[3]
#out_file = "/Users/lchen/Projects/sam_add_tag/raw_data/05_aggregated/new.fasta"
#in_file = "/Users/lchen/Projects/sam_add_tag/raw_data/05_aggregated/consensus_tide.tsv"

tide = pd.read_csv(in_file, sep="\t", names=['readname', 
                                          "consN", 
                                          "readLen", 
                                          "start", 
                                          "end", 
                                          "consLen", 
                                          "copyNum", 
                                          "fullLen",
                                          "subPos",
                                          "consensus"])
tp = list(zip(tide.readname, tide.consensus))
alist = [f">{i}\n{j}\n" for i,j in tp]
string = "".join(alist)
with open(out_fasta, "w") as f:
    f.write(string)
pd.to_pickle(tide, out_pickle)
print("trimming task finished.")
## for testing
# with open(out_file, "r") as f:
#     a = f.read()
#     print(a)
