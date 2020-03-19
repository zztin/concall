#$ -l h_rt=08:00:00
#$ -l h_vmem=16G
#$ -pe threaded 4

FILE_IN=$1
FILE_OUT=$2
FILE_REF=$3 #/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta 

# Based on Amins MC4C work, using his toolset as well
DIR_AMIN=/hpc/cog_bioinf/ridder/users/aallahyar/
BWA=$DIR_AMIN/My_Works/Useful_Sample_Codes/BWA/bwa/bwa


# Map split reads to reference genome
#echo "Mapping split reads to reference genome"
$BWA bwasw \
	-b 5 \
	-q 2 \
	-r 1 \
	-z 10 \
	-T 15 \
	-t 4 \
 	-f ${FILE_OUT} \
	${FILE_REF} \
	${FILE_IN}
	
