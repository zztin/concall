#$ -l h_rt=2:00:00
#$ -l h_vmem=8G
#$ -t 1-40:10
#$ -pe threaded 2
#$ -cwd
#$ -M l.t.chen-4@umcutrecht.nl
#$ -e ./std_err_p2
#$ -o ./std_out_p2

# bash /hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/concall/scripts/postprocessing.sh 02_split 03_bwa bb
FILE_ALL_FASTA=$1
DIR_TXT=`realpath $2`
DIR_FASTA=`realpath $3`
DIR_OUTPUT=`realpath $4` # without bb or ins
INS_OR_BB=$5
REF=$6

mkdir -p $DIR_OUTPUT

OUT=$DIR_OUTPUT

python3 /hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/concall/scripts/split_fasta.py ${DIR_FASTA} ${DIR_TXT} ${FILE_ALL_FASTA}

for NUM in {1..40}; do
bwa bwasw \
        -b 5 \
        -q 2 \
        -r 1 \
        -z 10 \
        -T 15 \
        -t 4 \
        -f ${DIR_OUTPUT}/$NUM.sam \
        ${FILE_REF} \
        ${DIR_FASTA}/consensus_${INS_OR_BB}_${NUM}.fasta


# convert sam file into bam file and sort the bam files produce bam and bai files.

/hpc/local/CentOS7/cog_bioinf/samtools-1.9/samtools view -h -F 256 $DIR_OUTPUT/$NUM.sam > $DIR_OUTPUT/$NUM.bam
/hpc/local/CentOS7/cog_bioinf/samtools-1.9/samtools sort -l 9 $DIR_OUTPUT/$NUM.bam -o $DIR_OUTPUT/$NUM.sorted.bam
/hpc/local/CentOS7/cog_bioinf/samtools-1.9/samtools index $DIR_OUTPUT/$NUM.sorted.bam
done
