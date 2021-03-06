#! /usr/bin/env python
import sys, os,re, subprocess
from optparse import OptionParser
from optparse import OptionGroup

if __name__ == "__main__":
        parser = OptionParser();
        group = OptionGroup(parser, "Main options")
        group.add_option("-k", dest="target_key", metavar="[PATH]", help="target_key [e.g. TP53]")
        group.add_option("-v", dest="target_value", metavar="[PATH]", help="value_key [e.g. 17:7565097-7590856]")
        group.add_option("-r", default="yes", dest="cosmic", help="calculate COSMIC position TP53 yes/no [default= yes]")
        group.add_option("-i", default=os.getcwd(), dest="input_dir", help="assign directory to find sorted.bam files, if not given, default is cwd")
        group.add_option("-o", default=None, dest="output", help="output directory")
        parser.add_option_group(group)
        (opt, args) = parser.parse_args()

dic={}
if opt.target_key and opt.target_value:
    dic[opt.target_key]=opt.target_value
else:
#    dic={"exon12":"17:7577000-7577180","BB22":"BB22:1-248","BB24":"BB24:1-248","BB25":"BB25:1-248","BBCR":"BB22:1-248","TP53":"17:7565097-7590856","PJET":"PJET:1-2974","EGFR":"7:55081714-55329313"}
    dic={"exon12":"17:7577000-7577180"}
sambamba = "sambamba"
#sambamba="/hpc/local/CentOS7/cog/software/sambamba-0.6.5/sambamba"
#cosmic= "/hpc/cog_bioinf/ridder/tools/Cyclomics_consensus_pipeline/data_files/COSMIC_mutations.bed"
if opt.output == None:
    opt.output= opt.input_dir
for f in os.listdir(opt.input_dir):
    if f.endswith((".sorted.bam")):
        for item in dic:
            sambamba_file_name = str(f).replace('.sorted.bam','_sambamba_output_'+str(item)+'.txt')

            action = f"{sambamba} depth base -L {dic[item]} --min-coverage=0 {opt.input_dir}/{f} > {opt.output}/{sambamba_file_name}"
            #action = str(sambamba)+" depth base -L "+str(dic[item])+ " --min-coverage=0 " + opt.input_dir+ str(f) + " > "+ opt.output+str(f).replace('.sorted.bam','_sambamba_output_'+str(item)+'.txt')
            os.system(action)
        if opt.cosmic=="yes":
            cosmic_file_name = str(f).replace('.sorted.bam','_sambamba_output_cosmic.txt')
            action = f"{sambamba} depth base -L {dic[item]} --min-coverage=0 {opt.input_dir}/{f} > {opt.output}/{cosmic_file_name}"            
            os.system(action)

# Create folders
move = False
if move ==True:
    os.system(f"mkdir -p {opt.output}/01_sam {opt.output}/02_bam {opt.output}/03_sorted {opt.output}/04_sambamba")
    os.system(f"mv {opt.output}/*.sorted.* {opt.output}/03_sorted")
    os.system(f"mv {opt.output}/*.sam {opt.output}/01_sam")
    os.system(f"mv {opt.output}/*.bam {opt.output}/02_bam")
    os.system(f"mv {opt.output}/*sambamba_out* {opt.output}/04_sambamba")
