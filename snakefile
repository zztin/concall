configfile: "./config_test.yaml"
# SAMPLES = ["40reads_119r10"] # <--- THIS IS WORKING
# SAMPLES = ["FAK80297_b08ac56b5a71e0628cfd2168a44680a365dc559f_301"]
#IN_PATH = "/hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/conbow2/"
#IN_PATH = "/hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/conbow2/results/"
#IN_PATH = os.path.expanduser(config["IN_PATH"])
#SUP_SAMPLES = ["CY_SS_PC_HC_0001_001_000_GU_303_HAC_20191031"]

#SEG = ["bb4", "targets"] 
#TYPE = ["bb", "ins"]
# snakemake define the first rule of the Snakefile as the target.
# Therefore, it is best practice to have a rule "all" on top of the workflow which define the target files as input files.
# end point
ruleorder: gz_fastq_get_fasta > fastq_get_fasta
#ruleorder: aggregation > aggregate_csv
# ruleorder: bowtie_wrapper_map > bowtie_map_backbone_to_read
#ruleorder: bowtie_map_backbone_to_read > bowtie_wrapper_map
SUP_SAMPLES = config['SUP_SAMPLES']
SAMPLES, = glob_wildcards(config['rawdir']+"/{sample}.fastq.gz")

print("SAMPLES", SAMPLES)

#SAMPLES = ['ABD169_9b86e52523af3f63ffea1043c200f43472e41222_19']
#print("SAMPLES:", SAMPLES)

#SAMPLES = ["FAK58127_e3b7026e6c44a11096370b0cfd31042b469e95fc_1"]
#SAMPLES = ["40reads_119r10"]


#rule all:
#    input:
#        expand("output/{SUP_SAMPLE}/05_aggregated/stats.csv", SUP_SAMPLE=SUP_SAMPLES),
#        expand("output/{SUP_SAMPLE}/05_aggregated/bb_consensus.fasta",SUP_SAMPLE=SUP_SAMPLES),
#        expand("output/{SUP_SAMPLE}/05_aggregated/ins_consensus.fasta",SUP_SAMPLE=SUP_SAMPLES),

rule all:
    input:
#        expand("output/{SUP_SAMPLE}/04_done/{sample}_bb.done", SUP_SAMPLE=SUP_SAMPLES, sample=SAMPLES),
#        expand("output/{SUP_SAMPLE}/04_done/{sample}_ins.done", SUP_SAMPLE=SUP_SAMPLES, sample=SAMPLES)
        expand("output/{SUP_SAMPLE}/02_split/stats/{sample}.csv", SUP_SAMPLE=SUP_SAMPLES, sample=SAMPLES)
#d        expand("output/03_consensus/bb/{sample}/consensus.fasta", sample=SAMPLES)
#        expand("output/011/{SUP_SAMPLE}_{sample}.done", sample=SAMPLES)
#rule print_name:
#    input:
#        "{IN_PATH}/{SUP_SAMPLE}/{sample}/{sample}.fastq"
#    output:
#        "otuput/011/{SUP_SAMPLE}_{sample}.done"
#    shell:
#        "echo {IN_PATH}/{SUP_SAMPLE}/{sample}/{sample}.fastq"


rule bedtool_getfasta:
    group: "bowtie_split"
    input:
        seg = "data/seg/{seg}.bed",
        ref = "/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta"
	#ref = config["genome"]
    output:
	"data/seg/{seg}.fa"
#        expand("data/seg/{seg}.fa",seg = SEG)
    shell:
        "bedtools getfasta -fi /hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta -bed {input.seg} -fo {output}"
rule gz_fastq_get_fasta:
    group: "bowtie_split"
    input:
        gz = config['rawdir']+"/{sample}.fastq.gz"
    output:
        touch("output/{SUP_SAMPLE}/01_bowtie/{sample}/createfolder.done"),
        fastq = temp("output/{SUP_SAMPLE}/00_fasta/{sample}.fastq"),
        fasta = "output/{SUP_SAMPLE}/00_fasta/{sample}.fasta"
    shell:
        "zcat {input.gz} > {output.fastq};\
         sed -n '1~4s/^@/>/p;2~4p' {output.fastq} > {output.fasta}"
rule fastq_get_fasta:
    group: "bowtie_split"
    input:
        fastq = "output/{SUP_SAMPLE}/00_fasta/{sample}.fastq"
    output:
        touch("output/{SUP_SAMPLE}/01_bowtie/{sample}/createfolder.done"),
        fasta = "output/{SUP_SAMPLE}/00_fasta/{sample}.fasta"
    shell:
        "sed -n '1~4s/^@/>/p;2~4p' {input.fastq} > {output.fasta}"

rule bowtie_build:
    # did not remove output while rule failed...
    # How to write output? answer here: https://www.biostars.org/p/342988/
    group: "bowtie_split"
    input:
        fasta = "output/{SUP_SAMPLE}/00_fasta/{sample}.fasta",
        done =  "output/{SUP_SAMPLE}/01_bowtie/{sample}/createfolder.done"
    params:
        re_path="output/{SUP_SAMPLE}/01_bowtie/{sample}/reference"
    log:
        "log/{SUP_SAMPLE}/bt-build_{sample}.log"
    output:
        touch('output/{SUP_SAMPLE}/01_bowtie/{sample}/bowtie_build_{sample}.done')
    shell:
        "/hpc/cog_bioinf/ridder/users/aallahyar/My_Works/Useful_Sample_Codes/Bowtie2/bowtie2-2.2.6/bowtie2-build -f {input.fasta} {params.re_path} > {log} 2>&1"

#rule bowtie_wrapper_map:
#    group:"bowtie_split"
#    input:
#        split_by = "data/seg/bb4.fa",
#        done = "output/{SUP_SAMPLE}/01_bowtie/{sample}/bowtie_build_{sample}.done"
#    output:
#        sam = "output/{SUP_SAMPLE}/01_bowtie/{sample}.sam"
#    log:
#        "log/bt-split_{SUP_SAMPLE}_{sample}.log"    
#    params:
#        index = "output/{SUP_SAMPLE}/01_bowtie/{sample}/reference",
#        extra="--local -D 20 -R 3 -N 0 -L 15 -i S,1,0.5 --rdg 2,1 --rfg 2,1 --mp 3,2 --ma 2 -a -f"
#    threads : 4
#    wrapper:
#        "0.42.0/bio/bowtie2/align"
        
    
rule bowtie_map_backbone_to_read:
    ## --very-sensitive-local -D 20 -R 3 -N 0 -L 20 -i S,1,0.5
    # we use -L 15: even more sensitive
    ## -i the interval between extracted seeds. (eg S, 1, 0.5)
    # interval function f(x) = <init1> + <init2> * sqrt(x). x = read length
    ## -rdg: read gap open, -rfg: ref gap open.
    # A gap of N gets a penalty of var1 + var2*N. default: 5,3
    ## -mp MX,MN (maximum and minimum mismatch penalties): modify alignment score considering Q(Phred score)
    ## -ma match bonus
    # http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer
    group: "bowtie_split"
    input:
        split_by = "data/seg/bb4.fa",
        done = "output/{SUP_SAMPLE}/01_bowtie/{sample}/bowtie_build_{sample}.done"
    params:
        # -x expect to find index files at current folder, then in the
        # directory specified in the BOWTIE2_INDEXES environment variable.
        # if this doesn't work, try:
        # export BOWTIE2_INDEXES=/path/to/my/bowtie2/databases/
        basename = "output/{SUP_SAMPLE}/01_bowtie/{sample}/reference"
    threads: 4
    benchmark:
        "log/benchmark/{SUP_SAMPLE}_{sample}_map_backbone_to_read_time.txt"
    log:
        "log/{SUP_SAMPLE}/bt-split_{sample}.log"
    output:
        sam = "output/{SUP_SAMPLE}/01_bowtie/{sample}.sam"
    shell: "/hpc/cog_bioinf/ridder/users/aallahyar/My_Works/Useful_Sample_Codes/Bowtie2/bowtie2-2.2.6/bowtie2 --local -D 20 -R 3 -N 0 -L 15 -i S,1,0.5 --rdg 2,1 --rfg 2,1 --mp 3,2 --ma 2 -a -p {threads} -f -x {params.basename} -U {input.split_by} -S {output.sam} > {log} 2>&1"

#rule clean_bowtie_ref:
#    input:
#        done = "output/02_split/{sample}_split.done",
#    output:
#        touch("output/02_split/{sample}_clean.done")
#    shell:
#        "rm -rf output/01_bowtie/{sample}/"

rule split_by_backbone:
    group: "bowtie_split"
    input:
        sam = "output/{SUP_SAMPLE}/01_bowtie/{sample}.sam",
        fasta = "output/{SUP_SAMPLE}/00_fasta/{sample}.fasta"
    benchmark:
        "log/benchmark/{SUP_SAMPLE}_{sample}.sam2_split_time.txt"
    output:
        bb = "output/{SUP_SAMPLE}/02_split/bb/{sample}.fasta",
        ins = "output/{SUP_SAMPLE}/02_split/ins/{sample}.fasta",
        stats = "output/{SUP_SAMPLE}/02_split/stats/{sample}.csv"
    log:
        "log/{SUP_SAMPLE}/sam2_split_{sample}.log"
    conda:
        "envs/pysam-env.yaml"
    shell:
        #"python scripts/sam2.py {input.sam} {input.fastq} {output}"
        "python3 scripts/sam2.py {input.sam} {input.fasta} {output.bb} {output.ins} {output.stats}"

rule smolecule_ins:
    group: "smolecule"
    input:
#        venv = ancient(IN_MEDAKA),
        fasta = "output/{SUP_SAMPLE}/02_split/ins/{sample}.fasta"
    params:
        path = directory("output/{SUP_SAMPLE}/03_consensus/ins/{sample}/"),
    output:
        consensus = "output/{SUP_SAMPLE}/03_consensus/ins/{sample}/consensus.fasta",
        done = touch("output/{SUP_SAMPLE}/04_done/{sample}_ins.done"),
#        consensus = "output/{SUP_SAMPLE}/03_consensus/ins/{sample}/consensus.fasta"
    log:
        "log/{SUP_SAMPLE}/smol_ins_{sample}.log"
    threads: 4
    benchmark:
        "log/benchmark/{SUP_SAMPLE}_{sample}.smolecule_ins_time.txt"
    conda:
        "envs/smolecule-env.yaml"
    shell:
         "ulimit -c 0; medaka smolecule --length 30 --depth 1 --threads {threads} {input.fasta} {params.path} > {log} 2>&1"
         #set +u; source activate snake-mdk; set -u;
         #medaka smolecule --length 50 --depth 5 --threads {threads} {input.fasta} {output.path} > {log} 2>&1

rule smolecule_bb:
    group: "smolecule"
    input:
        fasta = "output/{SUP_SAMPLE}/02_split/bb/{sample}.fasta"
    output:
        done = touch("output/{SUP_SAMPLE}/04_done/{sample}_bb.done"),
        consensus = "output/{SUP_SAMPLE}/03_consensus/bb/{sample}/consensus.fasta"
    params:
        path = directory("output/{SUP_SAMPLE}/03_consensus/bb/{sample}/")
    log:
        "log/{SUP_SAMPLE}/smol_bb_{sample}.log"
    threads: 4
    benchmark:
        "log/benchmark/{SUP_SAMPLE}_{sample}.smolecule_bb_time.txt"
    conda:
        "envs/smolecule-env.yaml"
    shell:
         "ulimit -c 0; medaka smolecule --length 30 --depth 1 --threads {threads} {input.fasta} {params.path} > {log} 2>&1"

# how to make sure all samples are done?
#ALL_BB, = glob_wildcards(output/{SUP_SAMPLE}/03_consensus/bb/{sample}/consensus.fasta")


#checkpoint clustering:
#    input:
#        expand("output/{SUP_SAMPLE}/03_consensus/bb/{sample}/consensus.fasta", SUP_SAMPLE=SUP_SAMPLES,  sample = SAMPLES)
#    output:
#        clusters=directory("output/{SUP_SAMPLE}/05_aggregated")
#    shell:
#        "mkdir {output}"

#rule intermediate:
#    # everything post processing goes in here.
#    input:
#        "output/aggregated/{SUP_SAMPLE}/{sample}.csv"
#    output:
#        "output/post_processing/{SUP_SAMPLE}/{sample}.csv"
#    shell:
#        "cp {input} {output}"

#def aggregate_input(wildcards):
#    checkpoint_output = checkpoints.clustering.get(**wildcards).output[0]
#    return expand("output/{SUP_SAMPLE}/03_consensus/bb/{sample}/consensus.fasta",
#           SUP_SAMPLE=SUP_SAMPLES,
#           sample=SAMPLES)
#           sample=glob_wildcards(os.path.join(checkpoint_output, "{sample}/consensus.fasta")).sample)

#rule aggregate_csv:
#    input:
#        aggregate_input
#    output:
#        csv = "output/{SUP_SAMPLE}/05_aggregated/stats.csv"
#    shell:
#        "cat {input} > {output.csv}"

TYPES = ["bb","ins"]

rule aggregation:
    input:
        csv = expand("output/{SUP_SAMPLE}/02_split/stats/{sample}.csv", SUP_SAMPLE=SUP_SAMPLES, sample=SAMPLES),
#        bb_fasta = aggregate_input
        bb_fasta = expand("output/{SUP_SAMPLE}/03_consensus/bb/{sample}/consensus.fasta",
           SUP_SAMPLE=SUP_SAMPLES,
           sample=SAMPLES),
        ins_fasta = expand("output/{SUP_SAMPLE}/03_consensus/ins/{sample}/consensus.fasta",
           SUP_SAMPLE=SUP_SAMPLES,
           sample=SAMPLES)
    output:
        csv = "output/{SUP_SAMPLE}/05_aggregated/stats.csv",
        bb = "output/{SUP_SAMPLE}/05_aggregated/bb_consensus.fasta",
        ins = "output/{SUP_SAMPLE}/05_aggregated/ins_consensus.fasta"
    shell:
        "cat {input.csv} > {output.csv}; cat {input.bb_fasta} > {output.bb}; cat {input.ins_fasta} > {output.ins}"

#rule count_repeat:
#    input:
#        csv = "output/{SUP_SAMPLE}/05_aggregated/bb/stats.csv",
#    output:
#        folder = directory("output/{SUP_SAMPLE}/05_aggregated/{type}/01_txt") 
#    params:
#        type = expand(type, type = TYPES)
#    script:
#        "scripts/fastq_splitby_consensus.py {output.folder} {input.csv} {params.type}"

#rule split_fasta:
    
        
#rule map_consensus:
#    input:
#
#    output:
#
#    log:
#    threads:
#    benchmark:
#    conda:
#    shell:



#rule split_by_backbone:
#    input:
#        sam = "output/01_bowtie/{sample}.sam",
#        fastq = "data/samples/{sample}.fastq"
#    output:
#        bb ="output/02_split/bb/{sample}.fastq",
#        ins = "output/02_split/ins/{sample}.fastq",
#        stats = "output/02_split/stats/{sample}.csv",
#        stats_subread = "output/02_split/stats/{sample}_subreads.csv"
#        touch("output/02_split/{sample}_split.done")
#    shell:
#        "/hpc/local/CentOS7/common/lang/python/2.7.10/bin/python scripts/sam.py --sam {input.sam} --fasta {input.fastq} --bb {output.bb} --ins {output.ins} --stats {output.stats} --st_sub {output.stats_subread}"
#rule create_read_repeats:
#    input: "data/output/02_split/{type}/{sample}.fastq" # run_name = {sample}_bb or {sample}_ins
#    output: directory("data/output/03_read_repeats/{type}/{sample}")
#    shell:
#        "sh ./scripts/split_fastq.sh {output} {input}"


#rule medaka:
#    input:
#        read = dynamic("data/output/03_read_repeats/{type}/{sample}/{readname}.fastq"),
#        draft = "data/seg/{seg}.fa",
#        ref = config["genome"]
#    params:
#        model = "r941_min_high"
#    output:
#        "data/output/04_consensus/{seg}/{sample}/{readname}"
#    shell:
#        "medaka_consensus -i {input.ref} -d {input.draft} -o {output} -m {params.model}"

#rule clean_fastq:
#    input:
#        "output/04_done/{sample}_bb.done"
#    output:
#        done = touch("output/04_done/{sample}_cleanup.done")
#    shell:
#        "rm -rf output/00_fasta/{sample}/* output/01_bowtie/{sample}*"

onsuccess:
    print("Workflow finished, no error. Success!")
#    shell("mail -s 'Workflow finished, no error!' litingchen16@gmail.com < {log}")
onerror:
    print("An notice sent to Liting by mail.")
 #   shell("mail -s 'an error occurred' litingchen16@gmail.com < {log}")
