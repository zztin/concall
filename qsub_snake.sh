# 20200103
# configfile: config_test.yaml
snakemake --jobs 4 --latency-wait 60 --use-conda --cluster "qsub -o /hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/concall/log/$(date +"%Y%m%d-%H%M%S")snake_stdout.txt -e /hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/concall/log/$(date +"%Y%m%d-%H%M%S")snake_stderr.txt -pe threaded {threads}"
