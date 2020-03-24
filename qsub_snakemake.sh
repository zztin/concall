TITLE=$1
snakemake --cluster sge_wrapper.py --jobs 50 --latency-wait 240 --use-conda --use-singularity --rerun-incomplete --keep-going --restart-times 3 --configfile config-$TITLE.yaml
snakemake --report $(date +"%Y%m%d-%H%M%S")-$TITLE.html
echo "report generated at ./datetime+title.html"
