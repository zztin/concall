#!/bin/bash

IN_DIR=$1
OUT_DIR=$2

mkdir -p $OUT_DIR
for FILE in $IN_DIR/*.fastq
do
	filename=$(basename $FILE .fastq)
	cat $FILE | split -l 160000 - ${OUT_DIR}/${filename}x --additional-suffix=.fastq --numeric-suffixes
done
