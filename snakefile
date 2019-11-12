configfile: "./config.yaml"
SAMPLES = ["40reads_119r10"]
SEG = ["bb4", "targets"] 
TYPE = ["bb", "ins"]
# snakemake define the first rule of the Snakefile as the target.
# Therefore, it is best practice to have a rule "all" on top of the workflow which define the target files as input files.
# end point
rule all:
    input:
        expand("output/02_split/ins/{sample}.fastq", sample=SAMPLES)
#        dynamic("data/output/04_consensus/{type}/{sample}/{readname}/consensus/consensus.fasta" )

rule bedtool_getfasta:
    input:
        seg = "data/seg/{seg}.bed",
        ref = "/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta"
	#ref = config["genome"]
    output:
	"data/seg/{seg}.fa"
#        expand("data/seg/{seg}.fa",seg = SEG)
    shell:
        "bedtools getfasta -fi /hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta -bed {input.seg} -fo {output}"

rule fastq_get_fasta:
    input:
        fastq = "data/samples/{sample}.fastq"
    output:
        fasta = "data/samples/{sample}.fasta"
    shell:
        "sed -n '1~4s/^@/>/p;2~4p' {input.fastq} > {output.fasta}"

rule bowtie_build:
    # How to write output? answer here: https://www.biostars.org/p/342988/
    input:
        fasta = expand("data/samples/{sample}.fasta", sample=SAMPLES)
    params:
        re_path="output/01_bowtie/reference"
    log:
        "log/bt-build/{sample}.log"
    output:
        touch('log/bowtie_build_{sample}.done')
    shell:
        "/hpc/cog_bioinf/ridder/users/aallahyar/My_Works/Useful_Sample_Codes/Bowtie2/bowtie2-2.2.6/bowtie2-build -f {input.fasta} {params.re_path} > {log}"

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
    input:
        split_by = "data/seg/bb4.fa",
        done = "log/bowtie_build_{sample}.done"
    params:
        # -x expect to find index files at current folder, then in the
        # directory specified in the BOWTIE2_INDEXES environment variable.
        # if this doesn't work, try:
        # export BOWTIE2_INDEXES=/path/to/my/bowtie2/databases/
        basename = "output/01_bowtie/reference"
    output:
        sam = "output/01_bowtie/{sample}.sam"
    shell: "/hpc/cog_bioinf/ridder/users/aallahyar/My_Works/Useful_Sample_Codes/Bowtie2/bowtie2-2.2.6/bowtie2 --local -D 20 -R 3 -N 0 -L 15 -i S,1,0.5 --rdg 2,1 --rfg 2,1 --mp 3,2 --ma 2 -a -p 4 -f -x {params.basename} -U {input.split_by} -S {output.sam}"

rule split_by_backbone:
    input:
        sam = "output/01_bowtie/{sample}.sam",
        fastq = "data/samples/{sample}.fastq"
    output:
        "output/02_split/bb/{sample}.fastq",
        "output/02_split/ins/{sample}.fastq"
    shell:
        "/hpc/local/CentOS7/common/lang/python/2.7.10/bin/python scripts/sam.py {input.sam} {input.fastq} {output}"
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



#onsuccess:
#    print("Workflow finished, no error")

#onerror:
#    print("An error occurred")
#    shell("mail -s 'an error occurred' litingchen16@gmail.com < {log}")
