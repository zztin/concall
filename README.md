# concall

![test-concall-on-mac](https://github.com/zztin/concall/workflows/test-concall-on-mac/badge.svg?branch=master)
![test-concall-ubuntu-linux](https://github.com/zztin/concall/workflows/test-concall-ubuntu-linux/badge.svg)

`concall` is a consensus calling pipeline based on open source deep learning based nanopore consensus sequences polishing "medaka", an experimental project developed and activately maintained by ONT Nanopore. We use the smolecule functionale to polish repeats of RCA product of short cfDNA fragments, aiming to provide high fidelity reads for clinical use. Accurate mutation calling is beneficial for early detection of cancer recurrence.

# Features
- snakemake pipeline management tool takes care of resource distribution 

# Installation

**Installation with conda**

	conda create -n snake-mdk -c conda-forge -c bioconda pip snakemake 
**Installation with pip**

	pip install snakemake

# Dependencies
You need to provide an indexed reference genome to the snakemake pipeline in the configfiles. If you do not have a reference genome yet, you need to download them. For example, download hg19 human reference genome from here:
hg19
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/ 
using 
```
rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet .
gunzip hg19.p13.plusMT.no_alt_analysis_set.fa.gz
bwa index -a bwtsw hg19.p13.plusMT.no_alt_analysis_set.fa
```
This could take up to 2 hours.

# Usage
`concall` can be run from command line after installation of snakemake.
1. Locate where your fastq/fastq.gz files are. Fastq files up to 40000 reads (from 40, 4000, to 40000) were tested. If fastq is over this range, batch them with scripts/batch_fastq_gz.sh; scripts/batch_fastq.sh. If the fastq is too small, use scripts/concat_fastq_gz.py to combine them into bigger fastq.  ## Bowtie build will increase run time exponentially when input file is bigger. Therefore smaller batch of fastq files is preferred. To retain the number of files in a managable range, limit the input files < 700 files (tested - of 20000 reads - P260).  
	
		scripts/batch_fastq_gz.sh <input_dir> <output_dir> 
		scripts/concat_fastq.py <input_dir> <output_dir: exact path> <batch-size> <prefix>

2. Secondly, prepare the relevant files including 
- reference genome
- prepare backbone sequences flanking the target sequences (optional)
`python ./backbone_processing/get_primers.py --help`
- construct a configfiles. A configfile example can be found at ./configfiles/. Locate the configfile at ./configfiles/config-XXXX.yaml

3. Lastly, decide where you want to run it (locally, or with qsub)

**run locally**

	cd concall/ # Where snakefile is located 
	snakemake --configfiles configfiles/config-test2.yaml --use-conda
	or
	snakemake --configfiles configfiles/config-test2.yaml--use-conda --use-singularity  # at the moment singularity image is stored locally at .sif file which is not included in this repo. Please contact authors for further information. This takes around 2 mins. Output will be generated at output/TEST40/
	
**run on cluster: qsub or slurm **
	
	# testing: 
	sh qsub_snakemake.sh configfiles/config-test2.yaml
	# or use this command:

	snakemake --cluster sge_wrapper.py --jobs 50 --latency-wait 240 --use-conda --use-singularity --rerun-incomplete --keep-going --restart-times 3 --configfile ./test/config/config-test2.yaml
	
	# replace sge_wrapper.py with slurm_wrapper.py to submit to slurm

**submit via bash script for submission on cluster**
	1. create config file at configfiles/config-<name>.yaml
	2. submit to cluster with:
		
		sh qsub_snakemake.sh <name>
	
	3. this will automatically generate an html report after the pipeline is finished.

# Development
1. Tidehunter program can be run locally on a Singularity container (tidehunter.sif) is used as a local file at the moment.
Recipe of the container is in the file Singularity. Currently Singularity Hub uses v2.5 but this image is built and tested on v3.5 machine. It gives errors if image is built by singularity hub and pull down to hpc. (this needs to fix. Before then, please contact the author to get the correct image of Singularity container.)
2. Tidehunter is installed via 
	
	conda install -c bioconda tidehunter
	
It runs normally on head node, compute node, but when submitted to hpc, fails occationally if running conda version. # need to check if it is resource issue? # need to check error message.

3. You can git clone the TideHunter package from https://github.com/Xinglab/TideHunter if to be run locally.



Contributions are welcome!"**


 
