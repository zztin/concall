#!/bin/bash/
# this is for old medaka function which build one consensus out of 1 file containing repeats within 1 read.
DIR_OUTPUT=`realpath $1` # Directory to store all output in
BASE_NAME=`basename  $DIR_OUTPUT`

awk 'function basename(file, a, n) {n = split(file, a, "/"); return a[n-1]}((FNR) % 2 == 1){print $0"_"basename(FILENAME)}((FNR) % 2 != 1){print $0}' $DIR_OUTPUT/medaka_bb/consensus/*/consensus.fasta > $DIR_OUTPUT/${BASE_NAME}_consensus_bb.fasta

