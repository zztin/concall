import pickle
import collections
from dateutil import parser
import os
import argparse

def get_timestamp(in_path, out_path, name, datype="fa"):
    time_string = collections.OrderedDict()
    time_stamp = collections.OrderedDict()

    if datype == "fq":
        for file in os.listdir(in_path):
            if file.endswith(".fastq"):
                with open(os.path.join(in_path, file), "r") as file:
                    content = file.read().splitlines()

                for i, txt in enumerate(content):

                    if i % 4 == 0:
                        readname = txt.split(" ")[0][1:]
                        time_string[readname] = txt.split(" ")[4].split("T")[1][:-1] # tstamp[readname] = time (hh:mm:ss)
                        t = parser.isoparse(txt.split(" ")[4].split("=")[1])
                        time_stamp[readname] = t.timestamp()
    if datype == "fa":
        for file in os.listdir(in_path):
            if file.endswith(".fasta"):
                with open(os.path.join(in_path, file), "r") as file:
                    content = file.read().splitlines()

                for i, txt in enumerate(content):

                    if i % 2 == 0:
                        readname = txt.split(" ")[0][1:]
                        time_string[readname] = txt.split(" ")[4].split("T")[1][
                                                :-1]  # time_string[readname] = time (hh:mm:ss)
                        t = parser.isoparse(txt.split(" ")[4].split("=")[1])
                        time_stamp[readname] = t.timestamp()

    if out_path == None:
        out_path= in_path

    with open(out_path, 'wb') as handle:
        pickle.dump(time_string, handle, protocol=pickle.HIGHEST_PROTOCOL)
    print("get_timestamp.py finished.")

if __name__ == '__main__':
    pa = argparse.ArgumentParser(description='Get timestamp of nanopore reads from raw fastq files.')
    pa.add_argument('-i', "--in_path", type=str,
                        help='provide exact (not relative) path to fq/fa folder of reads of interests.')
    pa.add_argument('-o', "--out_path", type=str,
                        help='provide path to store timestamp info as pickle.')
    pa.add_argument('-n', "--sample_name", type=str,
                        help='provide sample name.')
    pa.add_argument('-t', "--datype", type=str, default="fa",
                        help='fa or fq. default = fa')

    args = pa.parse_args()
    get_timestamp(args.in_path, args.out_path,  args.sample_name, args.datype)


