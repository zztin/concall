# 20200103
# configfile: config_test.yaml
snakemake --jobs 8 --latency-wait 120 --use-conda --verbose --cluster "qsub -l h_rt=00:30:00 -l h_vmem=30G -l tmpspace=40G -o /hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/concall/log/$(date +"%Y%m%d-%H%M%S")snake_stdout.txt -e /hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/concall/log/$(date +"%Y%m%d-%H%M%S")snake_stderr.txt -pe threaded {threads}"
