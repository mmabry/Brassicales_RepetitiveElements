## Brassicales_RepetitiveElements
Scripts used for Repetitive element content not correlated with whole-genome duplication or reflect phylogeny in the Brassicales (In prep)

# 1. Filter adaptors from raw reads

```bash
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
```

# 2. Run Trinity *note: this also trims to remove poor quality reads
#### download and un tar the latest version of trinity
```bash
#! /bin/bash

#SBATCH -J Trinity
#SBATCH -o Trinity.o%J
#SBATCH -e Trinity.e%J
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -p BioCompute,hpc5
#SBATCH --account=biosci
#SBATCH -t 2-00:00:00
#SBATCH --mem=80000

export PATH=/home/mmabry/software/trinityrnaseq/:$PATH
module load bowtie/bowtie-1.1.2
module load java/openjdk
module load trimmomatic/trimmomatic-0.35
module load tbb/tbb-4.4.4

Trinity --seqType fq --trimmomatic --quality_trimming_params 'SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25' --max_memory 80G --CPU 12 --full_cleanup --output Cleomella_serrulata_JHall.trinity --left Brassicales_RNA_RawReads/Cleomella_serrulata_JHall_R1.fastq.filtered --right Brassicales_RNA_RawReads/Cleomella_serrulata_JHall_R2.fastq.filtered
```

# 3. BUSCO check for Transcriptome completeness
#### install busco https://gitlab.com/ezlab/busco/blob/master/BUSCO_v3_userguide.pdf
#### make sure paths are correct in config file, use command below to get correct path
```bash
module show <name of program> 
#copy PATH and place it in config file
```
