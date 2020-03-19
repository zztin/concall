
def extract_good_reads(roy_input_filter_txt, output_good_read_list ):
    with open(roy_input_filter_txt,"r") as file:
        namelist = file.readlines()
    goodReads = []
    for x in range(len(namelist)):
        try:
            readID = namelist[x].split(" ")[1].split('_')[1]
            flag = namelist[x].split(" ")[3]
            if flag == "OK\n":
                goodReads.append(readID)
        except IndexError:
            pass
    with open(output_good_read_list, "w") as file:
        file.writelines(["%s\n" % item  for item in goodReads])

roy_input_filter_txt = "/Users/Alice/Projects/sam_add_tag/DER4126.txt"
output_good_read_list = "/Users/Alice/Projects/sam_add_tag/DER4126-filtered-good.txt"

extract_good_reads(roy_input_filter_txt, output_good_read_list )
