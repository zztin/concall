#!/bash/bin

snakemake --configfile ./configfiles/config-$1.yaml --snakefile tidehunter_only.smk --cores 1 --use-conda
#snakemake --configfile ./configfiles/config-$1.yaml --snakefile tidehunter_only.smk --report ./reports/$1_report.html
