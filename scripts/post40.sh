#$ -l h_rt=2:00:00
#$ -l h_vmem=8G
#$ -t 1-40:10
#$ -pe threaded 2
#$ -cwd
#$ -M l.t.chen-4@umcutrecht.nl
#$ -e ./std_err_p2
#$ -o ./std_out_p2
 
DIR_INPUT=`realpath $1`
DIR_OUTPUT=`realpath $2` # Directory to store all output in
INS_OR_BB=$3


mkdir -p $DIR_OUTPUT

NUM=40more

./bwa.sh \
        $DIR_INPUT/consensus_${INS_OR_BB}_40more.fasta \
        $DIR_OUTPUT/consensus_${INS_OR_BB}_40more.sam \
        /hpc/cog_bioinf/GENOMES/Cyclomics_reference_genome/version9/Homo_sapiens.GRCh37.GATK.illumina_cyclomics_backbone.fasta

# convert sam file into bam file and sort the bam files produce bam and bai files.

/hpc/local/CentOS7/cog_bioinf/samtools-1.9/samtools view -h $DIR_OUTPUT/consensus_${INS_OR_BB}_40more.sam > $DIR_OUTPUT/consensus_${INS_OR_BB}_40more.bam
/hpc/local/CentOS7/cog_bioinf/samtools-1.9/samtools sort -l 9 $DIR_OUTPUT/consensus_${INS_OR_BB}_40more.bam -o $DIR_OUTPUT/consensus_${INS_OR_BB}_40more.sorted.bam
/hpc/local/CentOS7/cog_bioinf/samtools-1.9/samtools index $DIR_OUTPUT/consensus_${INS_OR_BB}_40more.sorted.bam
