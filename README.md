# concall

`concall` is a consensus calling pipeline based on open source deep learning based nanopore consensus sequences polishing "medaka", an experimental project developed and activately maintained by ONT Nanopore. We use the smolecule functionale to polish repeats of RCA product of short cfDNA fragments, aiming to provide high fidelity reads for clinical use. Accurate mutation calling is beneficial for early detection of cancer recurrence.

# Features
- snakemake pipeline management tool takes care of resource distribution 

# Installation

**Installation with conda**

	conda create -n snake-mdk -c conda-forge -c bioconda pip snakemake 
**Installation with pip**

	pip install snakemake

# Usage
`concall` can be run from command line after installation of snakemake.
1. First, change the hard coded path in config files.
2. Second, decide where you want to run it (locally, or with qsub)

**run locally**

	cd concall/ # Where snakefile is located 
	snakemake --use-conda
	or
	snakemake --use-conda --use-singularity  # at the moment singularity image is stored locally at .sif file. 
	
**run on cluster: qsub or slurm **
	
	# testing: 
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
	
	conda -c bioconda tidehunter
	
It runs normally on head node, compute node, but when submitted to hpc, fails occationally if running conda version. # need to check if it is resource issue? # need to check error message.


Comtributions are welcome!"**


 
