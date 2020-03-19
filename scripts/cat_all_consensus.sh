#!bin/bash
SAMPLE_NAME=$1
OUT_NAME=$2

mkdir -p /hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/conbow2/results_all/$SAMPLE_NAME


cat /hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/conbow2/results/$SAMPLE_NAME/*/*.csv > /hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/conbow2/results_all/$SAMPLE_NAME/$OUT_NAME.csv
cat /hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/conbow2/results/$SAMPLE_NAME/*/*_consensus_bb.fasta > /hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/conbow2/results_all/$SAMPLE_NAME/$OUT_NAME.fasta
cd /hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/conbow2/results_all/$SAMPLE_NAME/
