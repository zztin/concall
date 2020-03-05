FILE_FASTA=$1
FILE_SAM=$2
FILE_REF=$3
NAME=$4
THREADS=$5
FILE_BAM=$6
FILE_SORTED_BAM=$7
bwa mem -t $5 -c 100 -M -R"@RG\\tID:$NAME\\tSM:$NAME\\tPL:NANOPORE\\tLB:$NAME" $3 $1 > $2;
samtools view -h -F 256 $2 > $6;
samtools sort -l 9 $6 -o $7;
samtools index $7;
# samtools view -F 256 <-- only primary reads
