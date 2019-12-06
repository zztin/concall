#$ -l h_rt=01:00:00
#$ -l h_vmem=4G
#$ -pe threaded 1
#$ -cwd
#$ -o /hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/concall/log/
#$ -e /hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/concall/log/
#$ -M lchen4@umcutrecht.nl
#$ -m a

#echo "test qsub"
conda activate snake-mdk

/hpc/cog_bioinf/ridder/users/lchen/miniconda37_2/envs/snake-mdk/bin/medaka smolecule --length 50 --depth 20 --threads 1 /hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/concall/data/samples/03_split_fasta_formatted/MUT_refWT_b00/FAK80297_b08ac56b5a71e0628cfd2168a44680a365dc559f_301_bb.fasta /hpc/cog_bioinf/ridder/users/lchen/Projects/Medaka_t/concall/data/samples/04_medaka_smolecule/MUT_refWT_b00/FAK80297_b08ac56b5a71e0628cfd2168a44680a365dc559f_301_bb_d20_t1_medaka_full_path/
echo "bash script finished."
conda deactivate
