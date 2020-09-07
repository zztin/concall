import re
import collections
import pysam
import pandas as pd
import argparse
import os

def get_ref_barcode_positions(barcode_fa):
    with open("../reference_genome/bbcr.fa", 'r') as f:
        BB = f.readlines()
    BB = BB[1].strip()
    bb_pos = [match.start() for match in re.finditer(re.escape("n"), BB)]
    return bb_pos


def get_barcode(bb_bam, barcode_pos=[32, 62, 79, 97], barcode_name="BB25", barcode_fa=None):
    if barcode_pos is None:
        barcode_pos = get_ref_barcode_positions(barcode_fa)
    ## Find random bases (barcodes) in backbones.
    barcode = collections.defaultdict(dict)
    # random positions, 0-based: [32, 62, 79, 97]
    start = barcode_pos[0]
    end = barcode_pos[-1]
    if os.path.isdir(bb_bam):
        for file in os.listdir(bb_bam):
            if file.endswith("sorted.bam"):
                samfile = pysam.AlignmentFile(file, "rb")
                print(file)
                for pileupcolumn in samfile.pileup(barcode_name, start, end):
                    if (pileupcolumn.pos == 32) or (pileupcolumn.pos == 62) or (pileupcolumn.pos == 79) or (pileupcolumn.pos == 97):
                        for pileupread in pileupcolumn.pileups:
                            if not pileupread.is_del and not pileupread.is_refskip:
                                read_name = (pileupread.alignment.query_name).split("_")[0]
                                print(read_name)
                                # query position is None if is_del or is_refskip is set.
                                barcode[read_name][f"b{pileupcolumn.pos}"].append(pileupread.alignment.query_sequence[
                                    pileupread.query_position])

                                # barcode[read_name][f"b{pileupcolumn.pos}_range7"] = pileupread.alignment.query_sequence[
                                #                                                     pileupread.query_position - 3:pileupread.query_position + 3]

                samfile.close()
    bb_combi = pd.DataFrame.from_dict(barcode).T
    bb_combi.fillna("X", inplace=True)
    bb_combi['BB'] = bb_combi.apply(lambda row: row.loc["b32"]+row.loc["b62"]+row.loc["b79"]+row.loc["b97"], axis=1)
    return bb_combi



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--bb_bam", help="bam file with the backbone sequences mapped with bowtie. directory 01_bowtie", type=str)
    parser.add_argument("--bb_type", help="Backbone type. This determines where the barcode is located",
                        type=str, default="BB25")
    args = parser.parse_args()
    bb_combi = get_barcode(args.bb_bam, barcode_pos=[32, 62, 79, 97], barcode_name="BB25", barcode_fa=None)
    print(bb_combi.head())
    pd.to_pickle(bb_combi, "../data/bb_combi.pickle")