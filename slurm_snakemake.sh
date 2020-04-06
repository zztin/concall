TITLE=$1
snakemake --cluster slurm_wrapper.py --jobs 50 --latency-wait 240 --use-conda --use-singularity --rerun-incomplete --keep-going --restart-times 3 --quiet --configfile ./configfiles/config-$TITLE.yaml
snakemake --report ./reports/$(date +"%Y%m%d-%H%M%S")-$TITLE.html --configfile ./configfiles/config-$TITLE.yaml
echo "report generated at ./datetime+title.html"
