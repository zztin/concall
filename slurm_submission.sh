snakemake --cluster slurm_wrapper.py --jobs 50 --configfile $1 --latency-wait 240 --keep-going --restart-times 3 --quiet
# snakemake --configfile $1 --report $1.report.html
