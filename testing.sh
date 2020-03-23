#/bin/bash
snakemake --restart-times 3 --configfile config_test2.yaml  --use-singularity
#snakemake --configfile config_test2.yaml --report
