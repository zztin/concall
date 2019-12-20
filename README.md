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

**run locally**

	cd concall/ # Where snakefile is located 
	snakemake   # this automatically search for snakefile and execute it

**run on cluster: qsub**

(This is also recorded in qsub_snake.sh in the repo)

	snakemake --jobs 100 --latency-wait 600 --use-conda --cluster "qsub -l h_rt=02:00:00 -l h_vmem=30G  -l tmpspace=100G -o /hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/concall/log/$(date "+%Y.%m.%d-%H.%M.%S")snake_stdout.txt -e /hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/concall/log/$(date "+%Y.%m.%d-%H.%M.%S")snake_stderr.txt -pe threaded {threads}"`

# Development
**Currently the qsub functionality is still facing fatal errors that force quit the pipeline. However, if re-initiate, these jobs runs smoothly without error. Current suspection of the failure are related to the latency of creating files on hpc. Trying out different --latency-wait value to see if it helps (with example data, latency = 600 finished successfully, =5 (default) failed). The error messages suggested that it is within each job that failed, however, executing these jobs separately in terminal gives no errors. Log files are empty for the failed jobs. Comtributions are welcome!"**


 
