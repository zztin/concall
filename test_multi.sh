# first qlogin 
conda activate snake-mdk
snakemake --use-conda --use-singularity --configfiles ./configfiles/config-test40unzip.yaml
echo $?
if [ "$?" == 0 ];
    echo "local test40unzip  done. successful"
sh ./qsub_snakemake.sh test40unzip
echo $?
if [ "$?" == 0 ];
    echo "sh qsub_snakemake.sh test40unzip done. successful"

#5d55c11c10cd48aeba184be1301af2dfcb49ede4 --  test2 always work with this version
sh ./qsub_snakemake.sh test2
echo $?
if [ "$?" == 0 ];
    echo "sh qsub_snakemake.sh test40gz done. successful"
