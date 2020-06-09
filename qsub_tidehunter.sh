TITLE=$1
snakemake --cluster "qsub -V -b y -S /bin/bash" --jobs 50 --latency-wait 240 --use-conda --rerun-incomplete --keep-going --configfile ./configfiles/config-$TITLE.yaml --snakefile ./tidehunter_only.smk
#snakemake --snakefile ./tidehunter_only.smk --configfile ./configfiles/config-$TITLE.yaml --report ./reports/$(date +"%Y%m%d-%H%M%S")-$TITLE.html
