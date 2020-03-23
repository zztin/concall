TITLE=$1
snakemake --cluster sge_wrapper.py --jobs 50 --latency-wait 240 --use-conda --use-singularity --rerun-incomplete --keep-going --configfile config-$TITLE.yaml --snakefile tidehunter.smk
snakemake --snakefile tidehunter.smk --configfile config-$TITLE.yaml --report $(date +"%Y%m%d-%H%M%S")-$TITLE.html
