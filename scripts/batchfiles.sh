#!/usr/bin/env bash

X=$1
Y=$2
batch=$3
filepath=$4




for i in $batch
do
  tail -n +$X batch$i $filepath | head -n $((Y-X+1))
done

echo 'finished.'

