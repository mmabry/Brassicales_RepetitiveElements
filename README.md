# Brassicales_RepetitiveElements
Scripts used for Repetitive element content not correlated with whole-genome duplication or reflect phylogeny in the Brassicales (In prep)

#1. Filter adaptors from raw reads
'''bash
#! /bin/bash

#SBATCH -J TrimAdapt
#SBATCH -o TrimAdapt.o%J
#SBATCH -e TrimAdapt.e%J
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -p BioCompute,hpc5
#SBATCH --account=biosci
#SBATCH -t 2-00:00:00
#SBATCH --mem=4000

module load python/python-2.7.13
export PATH=/home/mmabry/scratch/ncbi-blast-2.5.0+/bin/:$PATH


python /home/mmabry/yangya-phylogenomic_dataset_construction-489685700c2a/filter_fastq.py Cleomella_serrulata_JHall_R1.fastq Cleomella_serrulata_JHall_R2.fastq /home/mmabry/yangya-phylogenomic_dataset_construction-489685700c2a/data/UniVec-TruSeq_adapters 2
'''
