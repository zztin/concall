# first qlogin 
conda activate snake-mdk
# 0c51ad27f9ea9eae49f9b4140524acd33f5f0ee2 -- works.
snakemake --use-conda --use-singularity --configfiles ./configfiles/config-test40unzip.yaml
echo $?
if [ "$?" == 0 ];
    echo "local test40unzip  done. successful"
#2a20ebe7810a2595816ef2c823465e83273de5d4 -- works
sh ./qsub_snakemake.sh test40unzip
echo $?
if [ "$?" == 0 ];
    echo "sh qsub_snakemake.sh test40unzip done. successful"
#2a20ebe7810a2595816ef2c823465e83273de5d4 -- works
sh ./qsub_snakemake.sh test40000

echo $?
if [ "$?" == 0 ];
    echo "sh qsub_snakemake.sh test40000zip done. successful"
#5d55c11c10cd48aeba184be1301af2dfcb49ede4 --  test2 always work with this version
sh ./qsub_snakemake.sh test2
echo $?
if [ "$?" == 0 ];
    echo "sh qsub_snakemake.sh test40gz done. successful"
