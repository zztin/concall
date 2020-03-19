import pandas as pd
#import snakemake
#in_file = sys.argv[1]
#out_fasta = sys.argv[2]
#out_pickle = sys.argv[3]

tide = pd.read_csv(snakemake.input[0], sep="\t", names=['readname',
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
with open(snakemake.output[0], "w") as f:
    f.write(string)
pd.to_pickle(tide, snakemake.output[1])
