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
	snakemake   # this automatically search for snakefile and execute it

**run on cluster: qsub**

(This is also recorded in qsub_snake.sh in the repo)

	snakemake --jobs 100 --latency-wait 600 --use-conda --cluster "qsub -l h_rt=02:00:00 -l h_vmem=30G  -l tmpspace=100G -o [your path]/concall/log/$(date "+%Y%m%d-%H%M%S")snake_stdout.txt -e [your path]/concall/log/$(date "+%Y%m%d-%H%M%S")snake_stderr.txt -pe threaded {threads}"`


# Development
**Currently the qsub functionality is still facing fatal errors that force quit the pipeline. However, if re-initiate, these jobs runs smoothly without error. Current suspection of the failure are related to the latency of creating files on hpc. Trying out different --latency-wait value to see if it helps. The error messages suggested that it is within each job that failed, however, executing these jobs separately in terminal gives no errors. Log files are empty for the failed jobs. 
- 1. tried without specifying latency-wait: smolecule_ins/ smolecule_bb 2 out of 6 jobs failed
- 2. latency-wait =600s, successfully finished.
- 3. latency-wait =60s, failed at one of the smolecule_ins/smolecule_bb rule.
- 4. retry latency-wait = 600s, failed at group job bowtie_split (all jobs failed).
- 5? No idea what to suspect for the errors now...

Comtributions are welcome!"**


 
