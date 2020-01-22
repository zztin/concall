configfile: "./config_test2.yaml"
#configfile: "./config.yaml"
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
#ruleorder: gz_fastq_get_fasta > fastq_get_fasta
#ruleorder: aggregation > aggregate_csv
# ruleorder: bowtie_wrapper_map > bowtie_map_backbone_to_read
#ruleorder: bowtie_map_backbone_to_read > bowtie_wrapper_map
ruleorder: bwasw > bwa_mem> postprocessing
SUP_SAMPLES = config['SUP_SAMPLES']
TYPES = ["bb","ins"]

# gz or not gz
if config['gz'] == True:
    SAMPLES, = glob_wildcards(config['rawdir']+"/{sample}.fastq.gz")
else:
    SAMPLES, = glob_wildcards(config['rawdir']+"/{sample}.fastq")


print("SAMPLES", SAMPLES)

rule all:
    input:
#         expand("output/{SUP_SAMPLE}/07_stats_done/stats_{type}.done", SUP_SAMPLE=SUP_SAMPLES, type=TYPES)
#        expand("output/{SUP_SAMPLE}/05_aggregated/03_bwa_{type}/{count_name}.sam", SUP_SAMPLE=SUP_SAMPLES , type = TYPES, count_name = range(1,41) )
#        expand("output/{SUP_SAMPLE}/07_stats_done/postprocessing_{type}.done", SUP_SAMPLE=SUP_SAMPLES, type = TYPES)
        expand("output/{SUP_SAMPLE}/07_stats_done/postprocessing_{type}.done", SUP_SAMPLE=SUP_SAMPLES, type = TYPES)
#         expand("output/{SUP_SAMPLE}/04_done/{type}_split_fasta.done", SUP_SAMPLE=SUP_SAMPLES , type = TYPES)
#        expand("output/{SUP_SAMPLE}/05_aggregated/stats.csv", SUP_SAMPLE=SUP_SAMPLES)
#        expand("output/{SUP_SAMPLE}/04_done/{sample}_bb.done", SUP_SAMPLE=SUP_SAMPLES, sample=SAMPLES),
#        expand("output/{SUP_SAMPLE}/04_done/{sample}_ins.done", SUP_SAMPLE=SUP_SAMPLES, sample=SAMPLES)
#        expand("output/{SUP_SAMPLE}/02_split/stats/{sample}.csv", SUP_SAMPLE=SUP_SAMPLES, sample=SAMPLES)
#d        expand("output/03_consensus/bb/{sample}/consensus.fasta", sample=SAMPLES)
#        expand("output/011/{SUP_SAMPLE}_{sample}.done", sample=SAMPLES)

localrules: all, bedtool_getfasta, gz_fastq_get_fasta, fastq_get_fasta, bwasw, bwa_mem, aggregation, count_repeat, split_fasta, postprocessing

rule bedtool_getfasta:
#    group: "bowtie_split"
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
    # this rule takes > 5 sec < 60 sec to generate output files while submitted to cluster.
#    group: "bowtie_split"
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
#    group: "bowtie_split"
    input:
        fastq  = config['rawdir']+"/{sample}.fastq"
    output:
        touch("output/{SUP_SAMPLE}/01_bowtie/{sample}/createfolder.done"),
        fasta = temp("output/{SUP_SAMPLE}/00_fasta/{sample}.fasta")
    shell:
        "sed -n '1~4s/^@/>/p;2~4p' {input.fastq} > {output.fasta}"

#rule batch_fasta:
#    input:
#        fasta = "output/{SUP_SAMPLE}/00_fasta/{meta_sample}.fasta"
#    output:
#        fasta = dynamic("output/{SUP_SAMPLE}/001_split_fasta/{sample}.fasta")
#    params:
#        prefix = "output/{SUP_SAMPLE}/001_split_fasta/{meta_sample}",
#    shell:
#        "split -l 8000 --numeric-suffixes {input.fasta} {params.prefix} --additional-suffix=.fasta"

rule bowtie_build:
    # did not remove output while rule failed...
    # How to write output? answer here: https://www.biostars.org/p/342988/
#    group: "bowtie_split"
    input:
        fasta = "output/{SUP_SAMPLE}/00_fasta/{sample}.fasta",
        done =  "output/{SUP_SAMPLE}/01_bowtie/{sample}/createfolder.done"
    params:
        re_path="output/{SUP_SAMPLE}/01_bowtie/{sample}/reference"
    log:
        "log/{SUP_SAMPLE}/bt-build_{sample}.log"
    threads: 8
    conda:
        "envs/bt.yaml"
    output:
        touch('output/{SUP_SAMPLE}/01_bowtie/{sample}/bowtie_build_{sample}.done')
    shell:
        "bowtie2-build --threads {threads} -f {input.fasta} {params.re_path} > {log} 2>&1"

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
#    group: "bowtie_split"
    input:
        split_by = "data/seg/bb-C.fa",
        done = "output/{SUP_SAMPLE}/01_bowtie/{sample}/bowtie_build_{sample}.done"
    params:
        # -x expect to find index files at current folder, then in the
        # directory specified in the BOWTIE2_INDEXES environment variable.
        # if this doesn't work, try:
        # export BOWTIE2_INDEXES=/path/to/my/bowtie2/databases/
        basename = "output/{SUP_SAMPLE}/01_bowtie/{sample}/reference"
    threads: 8
    conda:
        "envs/bt.yaml"
    benchmark:
        "log/benchmark/{SUP_SAMPLE}_{sample}_map_backbone_to_read_time.txt"
    log:
        "log/{SUP_SAMPLE}/bt-split_{sample}.log"
    output:
        sam = "output/{SUP_SAMPLE}/01_bowtie/{sample}.sam"
    shell: "bowtie2 --local -D 20 -R 3 -N 0 -L 15 -i S,1,0.5 --rdg 2,1 --rfg 2,1 --mp 3,2 --ma 2 -a -p {threads} -f -x {params.basename} -U {input.split_by} -S {output.sam} > {log} 2>&1"

#rule clean_bowtie_ref:
#    input:
#        done = "output/02_split/{sample}_split.done",
#    output:
#        touch("output/02_split/{sample}_clean.done")
#    shell:
#        "rm -rf output/01_bowtie/{sample}/"

rule split_by_backbone:
#    group: "bowtie_split"
    input:
        sam = "output/{SUP_SAMPLE}/01_bowtie/{sample}.sam",
        fasta = "output/{SUP_SAMPLE}/00_fasta/{sample}.fasta"
    params:
        min_insert_length = 50, 
        max_insert_length = 500
    benchmark:
# specify wildcard!!!
        "log/benchmark/{SUP_SAMPLE}xxx{sample}.sam2_split_time.txt"
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
        "python3 scripts/sam2.py {input.sam} {input.fasta} {output.bb} {output.ins} {output.stats} {params.min_insert_length} {params.max_insert_length}"

rule smolecule_ins:
#    group: "smolecule"
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
#    group: "smolecule"
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
        bb = "output/{SUP_SAMPLE}/05_aggregated/all_consensus_bb.fasta",
        ins = "output/{SUP_SAMPLE}/05_aggregated/all_consensus_ins.fasta"
    shell:
        "cat {input.csv} > {output.csv}; cat {input.bb_fasta} > {output.bb}; cat {input.ins_fasta} > {output.ins}"


rule count_repeat:
    input:
        csv = rules.aggregation.output.csv
    output:
        folder = directory("output/{SUP_SAMPLE}/05_aggregated/01_txt_{type}"),
        done = touch("output/{SUP_SAMPLE}/04_done/{type}_bin_name.done")
    conda:
        "envs/bt.yaml"
    log:
        "log/{SUP_SAMPLE}/{type}_count_repeat.log"
    shell:
        "python scripts/fastq_splitby_consensus.py {output.folder} {input.csv} > {log};"

rule split_fasta:
    input:
         donefile = "output/{SUP_SAMPLE}/04_done/{type}_bin_name.done",
         fasta = "output/{SUP_SAMPLE}/05_aggregated/all_consensus_{type}.fasta"
    output:
        split = directory("output/{SUP_SAMPLE}/05_aggregated/02_split_{type}"),
        done = touch("output/{SUP_SAMPLE}/04_done/{type}_split_fasta.done")
    conda:
        "envs/bt-new.yaml"
    params:
        txt = "output/{SUP_SAMPLE}/05_aggregated/01_txt_{type}",
        split = "output/{SUP_SAMPLE}/05_aggregated/02_split_{type}"
    shell:
        "python scripts/split_fasta.py {output.split} {params.txt} {input.fasta}"

rule postprocessing:
    input:
        done= "output/{SUP_SAMPLE}/04_done/{type}_split_fasta.done"
    output:
        done = touch("output/{SUP_SAMPLE}/07_stats_done/postprocessing_{type}.done"),
        split_folder = directory("output/{SUP_SAMPLE}/05_aggregated/02_split_{type}/"),
        bwa_folder = directory("output/{SUP_SAMPLE}/05_aggregated/03_bwa_{type}")
#        sam = "output/{SUP_SAMPLE}/05_aggregated/03_bwa_{type}/{count_name}.sdam",
#        bam = "output/{SUP_SAMPLE}/05_aggregated/03_bwa_{type}/{count_name}.bam",
#        sorted = "output/{SUP_SAMPLE}/05_aggregated/03_bwa_{type}/{count_name}.sorted.bam"
    threads: 2
    params:
#        fasta = "output/{SUP_SAMPLE}/05_aggregated/02_split_{type, \s+[2-3]}/consensus_{type, \s+[2-3]}_{count_name, \d+}.fasta",
#        sorted = "output/{SUP_SAMPLE}/06_sorted",
        ref_genome_fasta = config["ref_genome_final"],
        type= "{type}"
        
    shell:
        "bash scripts/postprocessing.sh {output.split_folder} {output.bwa_folder} {params.type}"

rule bwasw:
    input:
        done= "output/{SUP_SAMPLE}/04_done/{type}_split_fasta.done"
    output:
        sam = "output/{SUP_SAMPLE}/05_aggregated/03_bwa_{type}/{count_name}.sam",
        bam = "output/{SUP_SAMPLE}/05_aggregated/03_bwa_{type}/{count_name}.bam",
        sorted = "output/{SUP_SAMPLE}/05_aggregated/03_bwa_{type}/{count_name}.sorted.bam"
    threads: 1
    params:
        fasta = "output/{SUP_SAMPLE}/05_aggregated/02_split_{type, \s+[2-3]}/consensus_{type, \s+[2-3]}_{count_name, \d+}.fasta",
        sorted = "output/{SUP_SAMPLE}/06_sorted",
        ref_genome_fasta = config["ref_genome_final"],
        name = "{SUP_SAMPLE}"
    shell:
        "bwasw \
        -b 5 \
        -q 2 \
        -r 1 \
        -z 10 \
        -T 15 \
        -t 4 \
        -f {output.sam} {params.ref_genome_fasta} {params.fasta}"

rule bwa_mem:
    input:
#        fasta = "output/{SUP_SAMPLE}/05_aggregated/02_split_{type, \s+[2-3]}/consensus_{type, \s+[2-3]}_{count_name, \d+}.fasta",
        done= "output/{SUP_SAMPLE}/04_done/{type}_split_fasta.done"
    output:
        sam = "output/{SUP_SAMPLE}/05_aggregated/03_bwa_{type}/{count_name}.sam",
        bam = "output/{SUP_SAMPLE}/05_aggregated/03_bwa_{type}/{count_name}.bam",
        sorted = "output/{SUP_SAMPLE}/05_aggregated/03_bwa_{type}/{count_name}.sorted.bam"
    threads: 2
    params:
        fasta = "output/{SUP_SAMPLE}/05_aggregated/02_split_{type, \s+[2-3]}/consensus_{type, \s+[2-3]}_{count_name, \d+}.fasta",
        sorted = "output/{SUP_SAMPLE}/06_sorted",
        ref_genome_fasta = config["ref_genome_final"],
        name = "{SUP_SAMPLE}"
    shell:
        "bash scripts/bwa_mem.sh {params.fasta} {output.sam} {params.ref_genome_fasta} {params.name} {threads} {output.bam} {output.sorted}"
#        "bwa mem -t {threads} -c 100 -M -R’@RG\tID:{params.name}\tSM:{params.name}\tPL:NANOPORE\tLB:{params.name}’ {params.ref_genome_fasta} {input.fasta} > {output.sam} "
#FILE_FASTA=$1
#FILE_SAM=$2
#FILE_REF=$3
#NAME=$4
#THREADS=$5
#FILE_BAM=$6
#FILE_SORTED_BAM=$7
#

rule final_stats_aggregation_bb:
    input:
        sorted = expand("output/{SUP_SAMPLE}/05_aggregated/03_bwa_bb/{count_name}.sorted.bam", SUP_SAMPLE=SUP_SAMPLES, count_name=range(1,41))
    output:
        done=touch("output/{SUP_SAMPLE}/07_stats_done/stats_bb.done"),
        out_dir = directory("output/{SUP_SAMPLE}/06_stats_bb/")
    params:
        base_folder = "output/{SUP_SAMPLE}/05_aggregated/03_bwa_bb/",
#        out_folder = "output/{SUP_SAMPLE}/06_stats_bb/",
        type = "bb"
    conda:
        "envs/bt.yaml"
    shell:
        "python scripts/calculate_depth_new_py37.py -i {params.base_folder} -o {output.out_dir}"
#        "python scripts/calculate_depth_new_py37.py -d {params.base_folder}; mkdir -p {params.base_folder}/04_sam {params.base_folder}/05_bam {params.base_folder}/06_sorted_bam {params.base_folder}/../../06_stats_{params.type};mv {params.base_folder}/*.sorted.bam {params.base_folder}/06_sorted_bam; mv {params.base_folder}/*.sam {params.base_folder}/04_sam; mv {params.base_folder}/*.bam* {params.base_folder}/05_bam; mv {params.base_folder}/*sambamba_out* {params.base_folder}/../../06_stats_{params.type}"


rule final_stats_aggregation_ins:
    input:
        sorted = expand("output/{SUP_SAMPLE}/05_aggregated/03_bwa_ins/{count_name}.sorted.bam", SUP_SAMPLE=SUP_SAMPLES, count_name=range(1,41))
    output:
        done=touch("output/{SUP_SAMPLE}/07_stats_done/stats_ins.done"),
        out_dir = directory("output/{SUP_SAMPLE}/06_stats_ins/")
    params:
        base_folder = "output/{SUP_SAMPLE}/05_aggregated/03_bwa_ins/",
#        out_folder = "output/{SUP_SAMPLE}/06_stats_ins/",
        type = "ins"
    conda:
        "envs/bt.yaml"
    shell:
        "python scripts/calculate_depth_new_py37.py -i {params.base_folder} -o {output.out_dir}"
#        "mkdir -p {params.base_folder}/../../06_stats_{params.type}; python scripts/calculate_depth_new_py37.py -i {params.base_folder} -o {params.out_folder}"
#        "python scripts/calculate_depth_new_py37.py -i {params.base_folder} -o {params.out_folder}; mkdir -p {params.base_folder}/../../06_stats_{params.type};cp {params.base_folder}/*sambamba_out* {params.base_folder}../../06_stats_{params.type}"
#        "python scripts/calculate_depth_new_py37.py -d {params.base_folder}; mkdir -p {params.base_folder}/04_sam {params.base_folder}/05_bam {params.base_folder}/06_sorted_bam {params.base_folder}/../../06_stats_{params.type};mv {params.base_folder}/*.sorted.bam {params.base_folder}/06_sorted_bam; mv {params.base_folder}/*.sam {params.base_folder}/04_sam; mv {params.base_folder}/*.bam* {params.base_folder}/05_bam; mv {params.base_folder}/*sambamba_out* {params.base_folder}/../../06_stats_{params.type}"

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
    shell("mail -s 'Workflow finished, no error!' litingchen16@gmail.com < {log}")
onerror:
    print("An notice sent to Liting by mail.")
    shell("mail -s 'an error occurred' litingchen16@gmail.com < {log}")
