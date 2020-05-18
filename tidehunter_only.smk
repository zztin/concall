#configfile:"configfiles/config-test2.yaml"
SUP_SAMPLES = config['SUP_SAMPLES']
TYPES = ["bb","ins"]

# gz or not gz
if config['gz'] == True:
    SAMPLES, = glob_wildcards(config['rawdir']+"/{sample}.fastq.gz")
else:
    SAMPLES, = glob_wildcards(config['rawdir']+"/{sample}.fastq")
print(f"There are {len(SAMPLES)} samples, {SAMPLES}.")
rule all:
    input:
        # timestamp, tide results, medaka results, cutadapt_info, bb_barcode(to be implement), tagged_bam_files (to be implement), samtools stats
        expand("output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_timestamp.pickle", SUP_SAMPLE=SUP_SAMPLES),
# bwa index for both tide and medaka
        expand("output/{SUP_SAMPLE}/07_stats_done/bwa_index_tide.done", SUP_SAMPLE=SUP_SAMPLES),
# samtool stats for tide
#        expand("output/{SUP_SAMPLE}/07_stats_done/samtools_stats.done", SUP_SAMPLE=SUP_SAMPLES),
 #       expand("output/{SUP_SAMPLE}/07_stats_done/samtools_stats_medaka.done", SUP_SAMPLE=SUP_SAMPLES),
# align bb sequence for extracting barcode (can also use tide to extract barcode)
#        expand( "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_bb.sorted.bam", SUP_SAMPLE=SUP_SAMPLES),
# per repeat bwa + sambamba result -- only for targeted
        #expand("output/{SUP_SAMPLE}/04_done/{type}_sambamba.done", SUP_SAMPLE=SUP_SAMPLES, type = TYPES)
localrules: all, bwasw, bwa_mem, get_timestamp, bedtool_getfasta, gz_fastq_get_fasta, fastq_get_fasta, aggregate_python, aggregate_tide, count_repeat, sambamba
ruleorder: tidehunter_conda > tidehunter_sing
rule get_timestamp:
    input:
        fasta = expand("output/{SUP_SAMPLE}/00_fasta/{SAMPLE}.fasta", SAMPLE=SAMPLES, SUP_SAMPLE=SUP_SAMPLES)
    output:
        timestamp = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_timestamp.pickle"
    params:
        fa = "output/{SUP_SAMPLE}/00_fasta/",
        name= "{SUP_SAMPLE}"
    conda:
        "envs/bt.yaml"
    shell:
        "python scripts/get_timestamp.py -i {params.fa} -o {output.timestamp} -n {params.name} --datype 'fa' " 

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
        gz = ancient(config['rawdir']+"/{sample}.fastq.gz")
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
        fastq  = ancient(config['rawdir']+"/{sample}.fastq")
    output:
        touch("output/{SUP_SAMPLE}/01_bowtie/{sample}/createfolder.done"),
        fasta = temp("output/{SUP_SAMPLE}/00_fasta/{sample}.fasta")
    shell:
        "sed -n '1~4s/^@/>/p;2~4p' {input.fastq} > {output.fasta}"


rule aggregate_tide:
    input:
        tide_fasta = expand("output/{SUP_SAMPLE}/05_aggregated/tide/{sample}_tide_consensus.fasta", SUP_SAMPLE=SUP_SAMPLES, sample=SAMPLES),
        tsv=expand("output/{SUP_SAMPLE}/05_aggregated/tide/{sample}_tide_consensus.tsv", SUP_SAMPLE=SUP_SAMPLES, sample=SAMPLES)
    output:
        tide = temp("output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_consensus_tide.fasta"),
        tsv = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_tide.tsv",
        done = touch("output/{SUP_SAMPLE}/04_done/aggregate_tide.done")
    params:
        tide = "output/{SUP_SAMPLE}/05_aggregated/tide",
    shell:
        "cat {params.tide}/*_tide_consensus.fasta > {output.tide};"
        "cat {params.tide}/*_tide_consensus.tsv > {output.tsv}"


rule cutadapt:
    # -e error rate
    # -b adaptor sequence (bb) can be anywhere in sequence. If in middle, all down stream of read is trimmed.
    # check cutadapt program
    input:
        ins = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_consensus_ins.fasta"
    output:
        cut_info="output/{SUP_SAMPLE}/06_cut/{SUP_SAMPLE}_cut_info.csv",
        fasta="output/{SUP_SAMPLE}/06_cut/{SUP_SAMPLE}_consensus_clean_bb.fasta",
        summary="output/{SUP_SAMPLE}/06_cut/{SUP_SAMPLE}_cutadapt_summary.txt"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000,
        runtime=lambda wildcards, attempt, input: ( attempt * 4)
    conda:"envs/bt.yaml"
    shell:
        "cutadapt -e 0.15 -b GGGCGGTATGTCATGCACACGAATCCCGAAGAnTGTTGTCCATTCATTGAATATGAGATCTCnATGGTATGATCAATATnCGGATGCGATATTGATAnCTGATAAATCATATATGCATAATCTCACATTATATTTATTATAATAAATCATCGTAGATATACACAATGTGAATTGTATACAATGGATAGTATAACTATCCAATTTCTTTGAGCATTGGCCTTGGTGTAGATTGCATGACATACCGCCC --action=lowercase --info-file {output.cut_info} -o {output.fasta} {input.ins} > {output.summary}"

rule bwa_wrapper_after_cutadapt:
    input:
        reads="output/{SUP_SAMPLE}/06_cut/{SUP_SAMPLE}_consensus_clean_bb.fasta"
    output:
        bam = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}-ins-clean.sorted.bam",
        done = touch("output/{SUP_SAMPLE}/07_stats_done/bwa_wrapper.done")
    log:
        "log/{SUP_SAMPLE}/{SUP_SAMPLE}_wrapper_bwa.log"
    params:
        index=config["ref_genome_final"],
        extra=r"-R '@RG\tID:{SUP_SAMPLE}\tSM:{SUP_SAMPLE}'",
        sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sort_extra="-l 9"            # Extra args for samtools/picard.
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 10000,
        runtime=lambda wildcards, attempt, input: ( attempt * 4)
    wrapper:
        "0.50.0/bio/bwa/mem"

rule bwa_wrapper_bb:
    input:
        reads="output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_consensus_bb.fasta"
    output:
        bam = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_bb.sorted.bam",
    log:
        "log/{SUP_SAMPLE}/{SUP_SAMPLE}_wrapper_bwa.log"
    params:
        index=config["ref_genome_final"],
        extra=r"-R '@RG\tID:{SUP_SAMPLE}\tSM:{SUP_SAMPLE}'",
        sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sort_extra="-l 9"            # Extra args for samtools/picard.
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 10000,
        runtime=lambda wildcards, attempt, input: ( attempt * 4)
    wrapper:
        "0.50.0/bio/bwa/mem"

rule bwa_index:
    input:
        tide_bam= "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}-ins-clean-tide.sorted.bam",
    output:
        done = touch("output/{SUP_SAMPLE}/07_stats_done/bwa_index_tide.done")
    conda:
        "envs/bt.yaml"
    shell:
        "samtools index {input.tide_bam};"

rule bwa_mem:
    input:
#        fasta = "output/{SUP_SAMPLE}/05_aggregated/02_split_{type, \s+[2-3]}/consensus_{type, \s+[2-3]}_{count_name, \d+}.fasta",
        done= "output/{SUP_SAMPLE}/04_done/{type}_split_fasta.done"
    output:
        sam = "output/{SUP_SAMPLE}/05_aggregated/03_bwa_{type}/{count_name}.sam",
        bam = "output/{SUP_SAMPLE}/05_aggregated/03_bwa_{type}/{count_name}.bam",
        sorted = "output/{SUP_SAMPLE}/05_aggregated/03_bwa_{type}/{count_name}.sorted.bam"
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000,
        runtime=lambda wildcards, attempt, input: ( attempt * 4)
    conda:
       "envs/bt.yaml"
    params:
        fasta = "output/{SUP_SAMPLE}/05_aggregated/02_split_{type, \s+[2-3]}/consensus_{type, \s+[2-3]}_{count_name, \d+}.fasta",
        sorted = "output/{SUP_SAMPLE}/06_sorted",
        ref_genome_fasta = config["ref_genome_final"],
        name = "{SUP_SAMPLE}"
    shell:
        "bash scripts/bwa_mem.sh {params.fasta} {output.sam} {params.ref_genome_fasta} {params.name} {threads} {output.bam} {output.sorted}"

rule tidehunter_sing:
    input:
       prime_3="data/seg/5_prime_bbcr.fa",
       prime_5="data/seg/3_prime_bbcr.fa",
       fasta="output/{SUP_SAMPLE}/00_fasta/{sample}.fasta"
    output:
        tsv="output/{SUP_SAMPLE}/05_aggregated/tide/{sample}_tide_consensus.tsv",
        done=touch("output/{SUP_SAMPLE}/04_done/{sample}_tide.done")
    threads: 4
    singularity:
        "tidehunter.sif"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20000,
        runtime=lambda wildcards, attempt, input: ( attempt * 1)
    shell:
        "/TideHunter-v1.2.2/bin/TideHunter -f 2 -t {threads} -5 {input.prime_5} -3 {input.prime_3} {input.fasta} > {output.tsv}"


rule tidehunter_conda:
    input:
       prime_3="data/seg/5_prime_bbcr.fa",
       prime_5="data/seg/3_prime_bbcr.fa",
       fasta="output/{SUP_SAMPLE}/00_fasta/{sample}.fasta"
    output:
        tsv="output/{SUP_SAMPLE}/05_aggregated/tide/{sample}_tide_consensus.tsv",
        done=touch("output/{SUP_SAMPLE}/04_done/{sample}_tide.done")
    threads: 4
    conda:
        "envs/tidehunter.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20000,
        runtime=lambda wildcards, attempt, input: ( attempt * 1)
    shell:
        "TideHunter -f 2 -t {threads} -5 {input.prime_5} -3 {input.prime_3} {input.fasta} > {output.tsv}"

rule trim_tide:
    # before BWA
    input:
        tsv="output/{SUP_SAMPLE}/05_aggregated/tide/{sample}_tide_consensus.tsv"
    output:
        fasta="output/{SUP_SAMPLE}/05_aggregated/tide/{sample}_tide_consensus.fasta",
#        pickle="output/{SUP_SAMPLE}/05_aggregated/tide/{sample}_tide_consensus_tsv.pickle.gz"
    conda:
       "envs/bt.yaml"
    shell:
        "python3 scripts/trim_tide_fasta_long_readname.py {input.tsv} {output.fasta}"

rule cutadapt_tide:
    input:
        tide="output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}_consensus_tide.fasta",
        done="output/{SUP_SAMPLE}/04_done/aggregate_tide.done"

    output:
        cut_info="output/{SUP_SAMPLE}/06_cut/{SUP_SAMPLE}_tide_cut_info.csv",
        fasta="output/{SUP_SAMPLE}/06_cut/{SUP_SAMPLE}_consensus_clean_bb_tide.fasta",
        summary="output/{SUP_SAMPLE}/06_cut/{SUP_SAMPLE}_cutadapt_tide_summary.txt"
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
        bam = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}-ins-clean-tide.sorted.bam",
        done = touch("output/{SUP_SAMPLE}/07_stats_done/bwa_wrapper_tide.done")
    log:
        "log/{SUP_SAMPLE}/{SUP_SAMPLE}_wrapper_bwa.log"
    params:
        index=config["ref_genome_final"],
        extra=r"-R '@RG\tID:{SUP_SAMPLE}\tSM:{SUP_SAMPLE}'",
        sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sort_extra="-l 9"            # Extra args for samtools/picard.
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 10000,
        runtime=lambda wildcards, attempt, input: ( attempt * 4)
    wrapper:
        "0.50.0/bio/bwa/mem"

rule plot_samtools_stats:
    input:
        bam = "output/{SUP_SAMPLE}/05_aggregated/{SUP_SAMPLE}-ins-clean-tide.sorted.bam", 
    output:
        stats = "output/{SUP_SAMPLE}/08_samtools_stats/tide/{SUP_SAMPLE}.stats",
        plot = directory("output/{SUP_SAMPLE}/08_samtools_stats/tide/{SUP_SAMPLE}_plot/"),
        done = touch("output/{SUP_SAMPLE}/07_stats_done/samtools_stats.done"),
    params:
        name = "{SUP_SAMPLE}_tide"
    conda:
        "envs/bt.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000,
    shell:
        "samtools stats {input.bam} > {output.stats};"
        "plot-bamstats -p {output.plot}{params.name} {output.stats};"
        "cat {output.stats} | grep ^RL | cut -f 2- > {output.plot}{params.name}_RL.txt;"


onsuccess:
    print("Workflow finished, no error. Success!")
    shell("mail -s 'Workflow finished, no error!' litingchen16@gmail.com < {log}")
onerror:
    print("An notice sent to Liting by mail.")
    shell("mail -s 'an error occurred' litingchen16@gmail.com < {log}")
