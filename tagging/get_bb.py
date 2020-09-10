import re
import collections
import pysam
import pandas as pd
import argparse
import os
import matplotlib.pyplot as plt
import ruptures as rpt
import numpy as np



def get_ref_barcode_positions(barcode_fa):
    with open("../reference_genome/bbcr.fa", 'r') as f:
        BB = f.readlines()
    BB = BB[1].strip()
    bb_pos = [match.start() for match in re.finditer(re.escape("n"), BB)]
    return bb_pos

def convert_base_to_num(value):
    if value == "T":
        out = 4
    elif value == "C":
        out = 1
    elif value == "A":
        out = -1
    elif value == "G":
        out = -4
    else:
        out = -100
    return out

def convert_num_to_base(value):
    if value == 4:
        out = "T"
    elif value == 1:
        out = "C"
    elif value == -1:
        out = "A"
    elif value == -4:
        out = "G"
    elif value == -100:
        out = "X"
    else:
        out = value
    return out


def get_barcode(bb_bam, barcode_pos=[32, 62, 79, 97], barcode_name="BB25", barcode_fa=None):
    if barcode_pos is None:
        barcode_pos = get_ref_barcode_positions(barcode_fa)
    ## Find random bases (barcodes) in backbones.
    barcode = collections.defaultdict(lambda: collections.defaultdict(list))
    # random positions, 0-based: [32, 62, 79, 97]
    start = barcode_pos[0]
    end = barcode_pos[-1]
    if os.path.isdir(bb_bam):
        for file in os.listdir(bb_bam):
            if file.endswith("sorted.bam"):
                # file path + file name
                samfile = pysam.AlignmentFile( bb_bam + file, "rb")
                for pileupcolumn in samfile.pileup(barcode_name, start, end):
                    if (pileupcolumn.pos == 32) or (pileupcolumn.pos == 62) or (pileupcolumn.pos == 79) or (pileupcolumn.pos == 97):
                        for pileupread in pileupcolumn.pileups:
                            if not pileupread.is_del and not pileupread.is_refskip:
                                read_name = (pileupread.alignment.query_name).split("_")[0]
                                base_rep = convert_base_to_num(pileupread.alignment.query_sequence[pileupread.query_position])
                                barcode[read_name][f"b{pileupcolumn.pos}"].append(base_rep)
                                # query position is None if is_del or is_refskip is set.
                            elif pileupread.is_del or pileupread.is_refskip:
                                barcode[read_name][f"b{pileupcolumn.pos}"].append(0)
                            else:
                                barcode[read_name][f"b{pileupcolumn.pos}"].append(0)
                                # barcode[read_name][f"b{pileupcolumn.pos}_range7"] = pileupread.alignment.query_sequence[
                                #                                                     pileupread.query_position - 3:pileupread.query_position + 3]

                samfile.close()
    bb_combi = pd.DataFrame.from_dict(barcode).T
    bb_combi.fillna(-100, inplace=True)
    bb_info = collections.defaultdict(dict)

    for i, row in bb_combi.iterrows():
        a = []
        most_freq = ""
        for x in [row['b32'], row['b62'], row['b79'], row['b97']]:
            if x == -100:
                x = [0, 0]
            else:
                pass
            try:
                most_common_2 = collections.Counter(x).most_common(2)
                most_common = most_common_2[0][0]
                if most_common == 0:
                    if len(most_common_2) == 1:
                        most_freq += "X"
                    else:
                        # pick the second common base rather than X
                        most_freq += convert_num_to_base(most_common_2[1][0])
                else:
                    most_freq += convert_num_to_base(most_common)
            except Exception as e:
                print(e)
                print(x)
            a.append(np.asarray(x))
        bb_info[i]['bc'] = most_freq
        max_len = max([len(x) for x in a])
        signal = [np.pad(x, (0, max_len - len(x)), 'constant') for x in a]
        new_signal = np.vstack(signal)
        break_points = get_chimeric_reads(i, new_signal.T, plot=False)
        if len(break_points) > 0:
            bb_info[i]['break_points'] = break_points
        else:
            bb_info[i]['break_points'] = False
    bb_info_df = pd.DataFrame.from_dict(bb_info).T
    return bb_info_df


def get_chimeric_reads(read_name, signal, fig_out="../data/bc_breakpoint", plot = False):
    algo = rpt.Pelt(model="rbf").fit(signal)
    result = algo.predict(pen=10)
    rpt.display(signal, [0,1], result)
    if plot:
        if not os.path.exists(fig_out):
            os.makedirs(fig_out)
        plt.savefig(f"{fig_out}/{read_name[:6]}.png")
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--bb_bam", help="bam file with the backbone sequences mapped with bowtie. directory 01_bowtie", type=str)
    parser.add_argument("--bb_type", help="Backbone type. This determines where the barcode is located",
                        type=str, default="BB25")
    args = parser.parse_args()
    # read name is the contig name. Use fasta to map to barcode again. Construct reference with only barcode
    bb_info_df = get_barcode(args.bb_bam, barcode_pos=[32, 62, 79, 97], barcode_name="BB24", barcode_fa=None)

    bb_info_df.to_pickle("../data/bb_info.pickle")

