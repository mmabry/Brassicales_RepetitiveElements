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
```bash
#! /bin/bash

#SBATCH -J BUSCO
#SBATCH -o BUSCO.o%J
#SBATCH -e BUSCO.e%J
#SBATCH -N 1
#SBATCH -n 14
#SBATCH -p BioCompute,Lewis
#SBATCH --account=biosci
#SBATCH -t 2-00:00:00
#SBATCH --mem=8000

module load python/python-2.7.13

python /home/mmabry/software/busco/scripts/run_BUSCO.py -i Cleomella_serrulata_JHall.trinity.Trinity.fasta -o Cleomella_serrulata -l /home/mmabry/software/busco/embryophyta_odb9 -m tran
```

#### Plot results
```bash
#create folder with all BUSCO_short_summaries
cp run_05103/short_summary_05103.txt BUSCO_summaries/.
```
```bash
#! /bin/bash

#SBATCH -J plotBUSCO
#SBATCH -o plotBUSCO.o%J
#SBATCH -e plotBUSCO.e%J
#SBATCH -N 1
#SBATCH -n 14
#SBATCH -p BioCompute,Lewis
#SBATCH --account=biosci
#SBATCH -t 2-00:00:00
#SBATCH --mem=8000

module load python/python-2.7.13

python /home/mmabry/scratch/busco/scripts/generate_plot.py -wd /home/mmabry/scratch/BUSCO_summaries
```
# 4. Translate assembled transcriptomes to amino acids for downstream analyses
#### download and unzip TransDecoder, build it by typing make in the base installation directiory
```bash
wget https://github.com/TransDecoder/TransDecoder/archive/v3.0.1.tar.gz
tar -xvzf <filename>
```
#### first script use LongORFs
```bash
#! /bin/bash

#SBATCH -J LongOrfs
#SBATCH -o LongOrfs.o%J
#SBATCH -e LongOrfs.e%J
#SBATCH -N 1
#SBATCH -n 14
#SBATCH -p BioCompute,Lewis
#SBATCH --account=biosci
#SBATCH -t 2-00:00:00
#SBATCH --mem=80000

export PATH=/home/mmabry/software/TransDecoder-3.0.1/:$PATH

for file in *.fasta; do /home/mmabry/software/TransDecoder-3.0.1/TransDecoder.LongOrfs -t ${file}; done
```
#### second script uses Predict
```bash
#! /bin/bash

#SBATCH -J Predict
#SBATCH -o Predict.o%J
#SBATCH -e Predict.e%J
#SBATCH -N 1
#SBATCH -n 14
#SBATCH -p BioCompute,Lewis
#SBATCH --account=biosci
#SBATCH -t 2-00:00:00
#SBATCH --mem=8000


export PATH=/home/mmabry/software/TransDecoder-3.0.1/:$PATH

for file in *.fasta; do /home/mmabry/software/TransDecoder-3.0.1/TransDecoder.Predict -t ${file}; done
```
# 5. Rename files
#### remane all transcriptomes from transdecoder to species_name.faa and then run script renameSeq.py to rename all transcripts species_name_1, species_name_2, species_name_3...etc
```python
#!/usr/bin/env python

from __future__ import print_function
from sys import argv

fn = argv[1]
SpeciesName = fn.split(".")[0]
counter = 0
with open(fn) as f:
        for line in f:
                if line[0] == '>':
                        counter += 1
                        print(">", SpeciesName, "_", counter, sep='')
                else:
                     	print(line.strip())
```
```bash
for file in *.faa; do python RenameSeqs.py $file > Fixed_names/$file; done
```
