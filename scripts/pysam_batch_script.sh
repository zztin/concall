#!/bin/bash

for i in {1..39}; do 
pysamstats --type variation --fasta ~/Projects/references/Homo_sapiens.GRCh37.GATK.illumina.fasta --chromosome 17 --start 7573944 --end 7574081 ./exact/bin_consensus_max40/consensus_$i.bam > ./pysamstats_results/amp5_con$i.txt
pysamstats --type variation --fasta ~/Projects/references/Homo_sapiens.GRCh37.GATK.illumina.fasta --chromosome 17 --start 7576800 --end 7576965 ./exact/bin_consensus_max40/consensus_$i.bam > ./pysamstats_results/amp9_con$i.txt
pysamstats --type variation --fasta ~/Projects/references/Homo_sapiens.GRCh37.GATK.illumina.fasta --chromosome 17 --start 7577000 --end 7577160 ./exact/bin_consensus_max40/consensus_$i.bam > ./pysamstats_results/amp12_con$i.txt; done

pysamstats --type variation --fasta ~/Projects/references/Homo_sapiens.GRCh37.GATK.illumina.fasta --chromosome 17 --start 7573944 --end 7574081 ./exact/bin_consensus_max40/consensus_40more.bam > ./pysamstats_results/amp5_con40more.txt

pysamstats --type variation --fasta ~/Projects/references/Homo_sapiens.GRCh37.GATK.illumina.fasta --chromosome 17 --start 7574000 --end 7576965 ./exact/bin_consensus_max40/consensus_40more.bam > ./pysamstats_results/amp9_con40more.txt

pysamstats --type variation --fasta ~/Projects/references/Homo_sapiens.GRCh37.GATK.illumina.fasta --chromosome 17 --start 7577000 --end 7577160 ./exact/bin_consensus_max40/consensus_40more.bam > ./pysamstats_results/amp12_con40more.txt


