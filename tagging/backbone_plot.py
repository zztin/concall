
def convert_base_to_num(value):
    if value == "T":
        out = 9
    elif value == "C":
        out = 7
    elif value == "A":
        out = 5
    elif value == "G":
        out = 3
    elif value == "X":
        out = 0
    else:
        out = value
    return out

def convert_num_to_base(value):
    if value == 9:
        out = "T"
    elif value == 7:
        out = "C"
    elif value == 5:
        out = "A"
    elif value == 3:
        out = "G"
    elif value == 0:
        out = "X"
    else:
        out = value
    return out

def convert_base_to_nan(value):
    if value == "T":
        out = 9
    elif value == "C":
        out = 7
    elif value == "A":
        out = 5
    elif value == "G":
        out = 3
    elif value == "X":
        out = np.nan
    else:
        out = value
    return out

if __name__ == '__main__':
    DER4622_bc_df_BB25 = get_barcode(bb_bam, barcode_pos=[32, 62, 79, 97], barcode_name="BB25", barcode_fa=None)
    df_num = DER4622_bc_df_BB25.apply(np.vectorize(convert_base_to_num), axis=0)
    sns.clustermap(bb_df,
                   col_cluster=False)

