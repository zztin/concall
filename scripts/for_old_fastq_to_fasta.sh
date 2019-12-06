#!/bin/bash
IN=$2
OUT=$1
for i in /hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/conbow2/results/${IN}/*/*_ins.fastq; do
    BASENAME=`basename $i .fastq`    
    sed -n '1~4s/^@/>/p;2~4p' $i > ${OUT}/${BASENAME}.fasta
done

for i in /hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/conbow2/results/${IN}/*/*_bb.fastq; do
    BASENAME=`basename $i .fastq`
    sed -n '1~4s/^@/>/p;2~4p' $i > ${OUT}/${BASENAME}.fasta
done
