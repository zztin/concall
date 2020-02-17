# 20200103
# configfile: config_test.yaml

TITLE=$1
snakemake --jobs 50 --latency-wait 180 --use-conda --rerun-incomplete --keep-going --cluster "qsub -l h_rt=6:00:00 -l h_vmem=80G -l tmpspace=200G -o /hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/concall/log/$(date +"%Y%m%d-%H%M%S")snake_stdout.txt -e /hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/concall/log/$(date +"%Y%m%d-%H%M%S")snake_stderr.txt -pe threaded {threads}"
#snakemake --report $(date +"%Y%m%d-%H%M%S")-$TITLE.html
#snakemake --configfile $CONFIG --jobs 50 --latency-wait 120 --use-conda --rerun-incomplete --keep-going --report report-$(date +"%Y%m%d-%H%M%S")-$TITLE.html --cluster "qsub -l h_rt=4:00:00 -l h_vmem=60G -l tmpspace=150G -o /hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/concall/log/$(date +"%Y%m%d-%H%M%S")snake_stdout.txt -e /hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/concall/log/$(date +"%Y%m%d-%H%M%S")snake_stderr.txt -pe threaded {threads}"
