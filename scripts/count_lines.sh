PATH=$1
for N in {1..40}
do 
samtools view $PATH/$N.sorted.bam | wc -l 
done

