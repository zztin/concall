#$ -l h_rt=2:00:00
#$ -l h_vmem=8G
#$ -t 1-40:10
#$ -pe threaded 2
#$ -cwd
#$ -M l.t.chen-4@umcutrecht.nl
#$ -e ./std_err_p2
#$ -o ./std_out_p2

# bash /hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/concall/scripts/postprocessing.sh 02_split 03_bwa bb

DIR_INPUT=`realpath $1`
DIR_OUTPUT=`realpath $2` # Directory to store all output in
INS_OR_BB=$3


mkdir -p $DIR_OUTPUT_${INS_OR_BB}

OUT=$DIR_OUTPUT_$INS_OR_BB

for NUM in {1..40}; do
/hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/concall/scripts/bwa.sh \
        ${DIR_INPUT}_${INS_OR_BB}/consensus_${INS_OR_BB}_${NUM}.fasta \
        ${DIR_OUTPUT}_${INS_OR_BB}/$NUM.sam \
        /hpc/cog_bioinf/GENOMES/Cyclomics_reference_genome/version9/Homo_sapiens.GRCh37.GATK.illumina_cyclomics_backbone.fasta

# convert sam file into bam file and sort the bam files produce bam and bai files.

/hpc/local/CentOS7/cog_bioinf/samtools-1.9/samtools view -h $DIR_OUTPUT_${INS_OR_BB}/$NUM.sam > $DIR_OUTPUT_${INS_OR_BB}/$NUM.bam
/hpc/local/CentOS7/cog_bioinf/samtools-1.9/samtools sort -l 9 $DIR_OUTPUT_${INS_OR_BB}/$NUM.bam -o $DIR_OUTPUT_${INS_OR_BB}/$NUM.sorted.bam
/hpc/local/CentOS7/cog_bioinf/samtools-1.9/samtools index $DIR_OUTPUT_${INS_OR_BB}/$NUM.sorted.bam
done
