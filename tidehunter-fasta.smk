#configfile:"./config-DER4535.yaml"
SUP_SAMPLES = config['SUP_SAMPLES']
TYPES = ["bb","ins"]

# gz or not gz
if config['gz'] == True:
    SAMPLES, = glob_wildcards(config['rawdir']+"/{sample}.fastq.gz")
else:
    SAMPLES, = glob_wildcards(config['rawdir']+"/{sample}.fastq")

rule all:
    input:
        expand("output/{SUP_SAMPLE}/07_stats_done/bwa_wrapper_tide.done", SUP_SAMPLE=SUP_SAMPLES)
#        expand("output/{SUP_SAMPLE}/07_stats_done/bwa-whole-{SUP_SAMPLE}-{type}.done", SUP_SAMPLE=SUP_SAMPLES, type = TYPES),
#        expand("output/{SUP_SAMPLE}/06_cut/{SUP_SAMPLE}_consensus_clean_bb_tide.fasta", SUP_SAMPLE=SUP_SAMPLES),
#        expand("output/{SUP_SAMPLE}/07_stats_done/bwa-whole-{SUP_SAMPLE}_clean_bb_tide.done", SUP_SAMPLE=SUP_SAMPLES)
       # expand("output/{SUP_SAMPLE}/07_stats_done/bwa-whole-{SUP_SAMPLE}_clean_bb.done", SUP_SAMPLE=SUP_SAMPLES),
        #expand("output/{SUP_SAMPLE}/06_cut/{SUP_SAMPLE}_cut_info.csv", SUP_SAMPLE=SUP_SAMPLES, type = TYPES),

#ruleorder: tidehunter_sing > tidehunter 
localrules: all, get_timestamp, bedtool_getfasta, gz_fastq_get_fasta, fastq_get_fasta, aggregate_python, aggregate_tide, bwasw, bwa_mem,  count_repeat, sambamba

rule get_timestamp:
    input:
        fasta = "output/{SUP_SAMPLE}/00_fasta/{sample}.fasta"
    output:
        timestamp = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_timestamp.pickle"
    conda:
        "envs/bt.yaml"
    script:
        "scripts/get_timestamp.py -i {input.fasta} -o {output.timestamp} -n {SUP_SAMPLE} --datype 'fa' " 

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
        fasta = temp("output/{SUP_SAMPLE}/00_fasta/{sample}.fasta")
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
    # How to write output? answer here: https://www.biostars.org/p/342988/
#    group: "bowtie_split"
    input:
        fasta = "output/{SUP_SAMPLE}/00_fasta/{sample}.fasta",
        done =  "output/{SUP_SAMPLE}/01_bowtie/{sample}/createfolder.done"
    params:
        re_path="output/{SUP_SAMPLE}/01_bowtie/{sample}/reference"
    log:
        "log/{SUP_SAMPLE}/bt-build_{sample}.log"
    benchmark:
        "log/benchmark/{SUP_SAMPLE}_{sample}_bowtie_build.txt"
    threads: 8
    conda:
        "envs/bt.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20000,
        runtime=lambda wildcards, attempt, input: ( attempt * 4)    
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
        split_by = config["BBTYPE"],
        done = "output/{SUP_SAMPLE}/01_bowtie/{sample}/bowtie_build_{sample}.done"
    params:
        # -x expect to find index files at current folder, then in the
        # directory specified in the BOWTIE2_INDEXES environment variable.
        # if this doesn't work, try:
        # export BOWTIE2_INDEXES=/path/to/my/bowtie2/databases/
        basename = "output/{SUP_SAMPLE}/01_bowtie/{sample}/reference",
        #runtime="6h"
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 40000,
        runtime=lambda wildcards, attempt, input: ( attempt * 1),
    conda:
        "envs/bt.yaml"
    benchmark:
        "log/benchmark/{SUP_SAMPLE}_{sample}_map_backbone_to_read_time.txt"
    log:
        "log/{SUP_SAMPLE}/bt-split_{sample}.log"
    output:
        sam = "output/{SUP_SAMPLE}/01_bowtie/{sample}.sam"
    shell: "bowtie2 --local --very-sensitive -a -p {threads} -f -x {params.basename} -U {input.split_by} -S {output.sam} > {log} 2>&1"
    #shell: "bowtie2 --local -D 20 -R 3 -N 0 -L 15 -i S,1,0.5 --rdg 2,1 --rfg 2,1 --mp 3,2 --ma 2 -a -p {threads} -f -x {params.basename} -U {input.split_by} -S {output.sam} > {log} 2>&1"

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
        min_insert_length = config['min_insert_length'],
        max_insert_length = config['max_insert_length']
    benchmark:
# specify wildcard!!!
        "log/benchmark/{SUP_SAMPLE}xxx{sample}.sam2_split_time.txt"
    output:
        bb = "output/{SUP_SAMPLE}/02_split/bb/{sample}.fasta",
        ins = "output/{SUP_SAMPLE}/02_split/ins/{sample}.fasta",
        stats = "output/{SUP_SAMPLE}/02_split/stats/{sample}.csv"
    log:
        "log/{SUP_SAMPLE}/sam2_split_{sample}.log"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20000,
        runtime=lambda wildcards, attempt, input: ( attempt * 4)
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
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 6000,
        runtime=lambda wildcards, attempt, input: ( attempt * 1)
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
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 6000,
        runtime=lambda wildcards, attempt, input: ( attempt * 1)
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

rule aggregate_python:
    input:
        mdk_bb = expand("output/{SUP_SAMPLE}/04_done/{sample}_bb.done", SUP_SAMPLE=SUP_SAMPLES, sample=SAMPLES),
        mdk_ins = expand("output/{SUP_SAMPLE}/04_done/{sample}_ins.done", SUP_SAMPLE=SUP_SAMPLES, sample=SAMPLES),
    output:
        csv = "output/{SUP_SAMPLE}/05_aggregated/stats.csv",
        bb = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_consensus_bb.fasta",
        ins = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_consensus_ins.fasta",
        done = touch("output/{SUP_SAMPLE}/04_done/aggregate.done")
    params:
        bb = "output/{SUP_SAMPLE}/03_consensus/bb",
        ins = "output/{SUP_SAMPLE}/03_consensus/ins",
        stats = "output/{SUP_SAMPLE}/02_split/stats"
    shell:
        "cat {params.bb}/*/consensus.fasta > {output.bb};"
        "cat {params.ins}/*/consensus.fasta > {output.ins};"
        "cat {params.stats}/*.csv > {output.csv};"

rule aggregate_tide:
    input:
        tide_fasta = expand("output/{SUP_SAMPLE}/05_aggregated/tide/{sample}_tide_consensus.fasta", SUP_SAMPLE=SUP_SAMPLES, sample=SAMPLES),
    output:
        tide = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_consensus_tide.fasta",
        done = touch("output/{SUP_SAMPLE}/04_done/aggregate_tide.done")
    params:
        tide = "output/{SUP_SAMPLE}/05_aggregated/tide",
    shell:
        "cat {params.tide}/*_tide_consensus.fasta > {output.tide};"

rule cutadapt_tide:
    input:
        tide="output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_consensus_tide.fasta",
        done="output/{SUP_SAMPLE}/04_done/aggregate_tide.done"
    
    output:
        cut_info="output/{SUP_SAMPLE}/06_cut/{SUP_SAMPLE}_tide_cut_info.csv",
        fasta="output/{SUP_SAMPLE}/06_cut/{SUP_SAMPLE}_consensus_clean_bb_tide.fasta",
        summary="output/{SUP_SAMPLE}/06_cut/{SUP_SAMPLE}_tide_summary.txt"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000,
        runtime=lambda wildcards, attempt, input: ( attempt * 4)
    conda:"envs/bt.yaml"
    shell:
        "cutadapt -e 0.15 -b GGGCGGTATGTCATGCACACGAATCCCGAAGAnTGTTGTCCATTCATTGAATATGAGATCTCnATGGTATGATCAATATnCGGATGCGATATTGATAnCTGATAAATCATATATGCATAATCTCACATTATATTTATTATAATAAATCATCGTAGATATACACAATGTGAATTGTATACAATGGATAGTATAACTATCCAATTTCTTTGAGCATTGGCCTTGGTGTAGATTGCATGACATACCGCCC --action=lowercase --info-file {output.cut_info} -o {output.fasta} {input.tide} > {output.summary}"


rule bwa_wrapper_tide:
    input:
        reads="output/{SUP_SAMPLE}/06_cut/{SUP_SAMPLE}_consensus_clean_bb_tide.fasta"
    output:
        bam = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}-ins-clean-tide.bam",
        done = touch("output/{SUP_SAMPLE}/07_stats_done/bwa_wrapper_tide.done")
    log:
        "log/{SUP_SAMPLE}/{SUP_SAMPLE}_wrapper_bwa.log"
    params:
        index=config["ref_genome_final"],
        extra=r"-R '@RG\tID:{SUP_SAMPLE}\tSM:{SUP_SAMPLE}'",
        sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
        sort_extra="-l 9"            # Extra args for samtools/picard.
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000,
        runtime=lambda wildcards, attempt, input: ( attempt * 4)
    wrapper:
        "0.50.0/bio/bwa/mem"



#rule tidehunter_sing:
#    input:
#       prime_3="data/seg/5_prime_bbcr.fa",
#       prime_5="data/seg/3_prime_bbcr.fa",
#       fasta="output/{SUP_SAMPLE}/00_fasta/{sample}.fasta"
#    output:
#        fasta="output/{SUP_SAMPLE}/05_aggregated/tide/{sample}_tide_consensus.fasta",
#        done=touch("output/{SUP_SAMPLE}/04_done/{sample}_tide.done")
#    threads: 4
#    singularity:
#        "tidehunter.sif"
#    resources:
#        mem_mb=lambda wildcards, attempt: attempt * 20000,
#        runtime=lambda wildcards, attempt, input: ( attempt * 1)
#    shell:
#        "/TideHunter-v1.2.2/bin/TideHunter -f 2 -t {threads} -5 {input.prime_5} -3 {input.prime_3} {input.fasta} > {output.fasta}"

rule tidehunter:
    # -f 2: output is fasta format (9 columns)
    # column 8: fullLen: if the read contains backbone (5' and 3' end found)
    input:
       prime_3="data/seg/5_prime_bbcr.fa",
       prime_5="data/seg/3_prime_bbcr.fa",
       fasta="output/{SUP_SAMPLE}/00_fasta/{sample}.fasta"
    output:
        fasta="output/{SUP_SAMPLE}/05_aggregated/tide/{sample}_tide_consensus.fasta",
        done=touch("output/{SUP_SAMPLE}/04_done/{sample}_tide.done")
    threads: 4
    conda:
        "envs/bt.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000,
        runtime=lambda wildcards, attempt, input: ( attempt * 1)
    shell:
        "TideHunter-local/bin/TideHunter -f 2 -t {threads} -5 {input.prime_5} -3 {input.prime_3} {input.fasta} > {output.fasta}"

rule trim_tide:
    # before BWA
    input:
        fasta="output/{SUP_SAMPLE}/05_aggregated/tide/{sample}_tide_consensus.fasta"
    output:
        fasta="output/{SUP_SAMPLE}/05_aggregated/tide/{sample}_tide_consensus.fasta",
        pickle="output/{SUP_SAMPLE}/05_aggregated/tide/{sample}_tide_consensus_tsv.pickle.gz"
    conda:
       "envs/bt.yaml"
    shell:
        "python3 scripts/trim_tide_fasta_long_readname_fasta.py {input.fasta} {output.fasta} {output.pickle}"
#rule trim_tide_script:
#    input:
#        "output/{SUP_SAMPLE}/05_aggregated/tide/{sample}_tide_consensus.tsv"
#    output:
#        "output/{SUP_SAMPLE}/05_aggregated/tide/{sample}_tide_consensus.fasta",
#        "output/{SUP_SAMPLE}/05_aggregated/tide/{sample}_tide_consensus_tsv.pickle.gz"
#    conda:
#       "envs/bt.yaml"
#    script:
#        "scripts/trim_tide_fasta_long_readname_smk_script.py"


onsuccess:
    print("Workflow finished, no error. Success!")
    shell("mail -s 'Workflow finished, no error!' litingchen16@gmail.com < {log}")
onerror:
    print("An notice sent to Liting by mail.")
    shell("mail -s 'an error occurred' litingchen16@gmail.com < {log}")
