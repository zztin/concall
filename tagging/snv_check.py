
from pysam import VariantFile

vcf_in = VariantFile("/Users/lchen/Projects/sam_add_tag/raw_data/DER4535/HGS-1-R2_tumor_somatic_post_processed.vcf",
                    index_filename = "/Users/lchen/Projects/sam_add_tag/raw_data/DER4535/HGS-1-R2_tumor_somatic_post_processed.vcf.tbi")  # auto-detect input format
vcf_out = VariantFile('/Users/lchen/Projects/sam_add_tag/raw_data/DER4535/DER4535_chr1.vcf', 'w', header=vcf_in.header)

for rec in vcf_in.fetch('1', 10610, 14620):
    vcf_out.write(rec)


vcf_loc = []
for rec in vcf_in:
    print(rec)
    a = str(rec)
    vcf_loc += [(a.split("\t")[0], a.split("\t")[1])]
