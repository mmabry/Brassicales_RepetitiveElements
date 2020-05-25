## Brassicales_RepetitiveElements
Scripts used for Repetitive element content not correlated with whole-genome duplication or reflect phylogeny in the Brassicales (In prep)

# Table of contents
- [1. Phylogeny](#1-phylogeny)
  - [A. Filter adaptors from raw reads](#a-Filter-adaptors-from-raw-reads)
  - [B. Trinity](#b-run-trinity-note-this-also-trims-to-remove-poor-quality-reads)
  - [C. BUSCO check for Transcriptome completeness](#c-BUSCO-check-for-Transcriptome-completeness)
  - [D. Translate assembled transcriptomes](#d-Translate-assembled-transcriptomes-to-amino-acids-for-downstream-analyses)
  - [E. Rename files](#e-Rename-files)
  - [F. Run OrthoFinder to determine OrthoGroups](#f-Run-OrthoFinder-to-determine-OrthoGroups)
  - [G. Filter alignments by taxa occupancy](#g-filter-alignments-by-taxa-occupancy-httpsgithubcommu-ircffilter_by_ortho_group)
  - [H. Filter alignments based on quality/gaps](#h-filter-alignments-based-on-qualitygaps-httpsgithubcommu-ircffilter_by_gap_fraction)
  - [I. Rename headers for PhylotreePruner](#i-Rename-headers-for-PhylotreePruner-(in-later-in-the-steps)-before-running-first-set-of-trees)
  - [J. Run RAxML for tree inference](#j-Run-RAxML-for-tree-inference)
  - [K. PhyloTreePruner](#k-PhyloTreePruner)
  - [L. Final gene tree estimation](#l-Finally,-run-RAxML-one-more-time-to-get-final-gene-trees)
  - [M. Species tree inference with ASTRAL](#m-Species-tree-inference-with-ASTRAL)
- [2. Repetitive Element Clustering](#2-Repetitive-Element-Clustering)
- [3. Regression Analyses](#3-Regression-Analyses)
- [4. Hierarchical Clustering](#4-Hierarchical-Clustering)
- [5. Ultrametric Tree](#5-Ultrametric-Tree) 
  - [A. Concatenate alignments](#a-concatenate-alignments-using-the-concatenate_matricespy-script-from-httpsbitbucketorgwashjaketranscriptome_phylogeny_toolssrcmaster)
  - [B. Run RAxML to optimize Branch lengths](#b-Run-RAxML-to-optimize-Branch-lengths-and-model-parameters-using-the-concatenated-alignment-and-ASTRAL-phylogeny-as-a-fixed-input-tree.)
  - [C. Use TreePL to time calibrate phylogeny](#c-use-treepl-to-time-calibrate-phylogeny-httpsgithubcomblackrimtreeplwikiquick-run)
- [6. Bayou](#6-bayou-httpsgithubcomuyedajbayoublobmastertutorialmd)
- [7. Owie](#7-owie--example-from-from-httpwwwphytoolsorgcordoba2017ex10multi-regimehtml)



# 1. Phylogeny 
## A. Filter adaptors from raw reads
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

## B. Run Trinity *note: this also trims to remove poor quality reads
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

## C. BUSCO check for Transcriptome completeness
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
## D. Translate assembled transcriptomes to amino acids for downstream analyses
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
## E. Rename files
#### remane all transcriptomes from transdecoder (*.pep file) to species_name.faa and then run script renameSeq.py to rename all transcripts species_name_1, species_name_2, species_name_3...etc
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
## F. Run OrthoFinder to determine OrthoGroups
#### first run it using diamond to get groups
```bash
#! /bin/bash

#SBATCH -J Ortho_Cleo
#SBATCH -o Ortho_Cleo.o%J
#SBATCH -e Ortho_Cleo.e%J
#SBATCH -N 1
#SBATCH -n 56
#SBATCH -p BioCompute,Lewis
#SBATCH --qos=long
#SBATCH -t 7-00:00:00
#SBATCH --mem=400G


module load python/python-2.7.13
module load py-scipy/py-scipy-0.18.1
module load py-numpy/py-numpy-1.11.2
module load mafft/mafft-7.299
module load ircf/ircf-modules 
module load diamond/diamond-0.9.22
export PATH="/home/mmabry/software/mcl-14-137/bin/":$PATH
export PATH=/home/mmabry/software/FastME/tarball/fastme-2.1.5.1/binaries/:$PATH
export PATH=/home/mmabry/software/:$PATH #this is for running diamond and FastTree


/home/mmabry/software/OrthoFinder-2.2.6/OrthoFinder/orthofinder/orthofinder.py -f /storage/htc/pireslab/mmabry_htc/Brassicales_Transcriptomes_aa/FixedNames/renameSeq/Family_level/Cleomaceae/ -t 56 -S diamond
```
#### run it a second time to get alignments needed for filtering in the next step

```bash
#! /bin/bash

#SBATCH -J Ortho2_Capp
#SBATCH -o Ortho2_Capp.o%J
#SBATCH -e Ortho2_Capp.e%J
#SBATCH -N 1
#SBATCH -n 56
#SBATCH -p BioCompute,Lewis
#SBATCH -t 2-00:00:00
#SBATCH --mem=400G

module load python/python-2.7.13
module load py-scipy/py-scipy-0.18.1
module load py-numpy/py-numpy-1.11.2
module load mafft/mafft-7.299
module load ircf/ircf-modules 
module load fasttree/fasttree-2.1.9 
export PATH="/home/mmabry/software/mcl-14-137/bin/":$PATH
export PATH=/home/mmabry/software/FastME/tarball/fastme-2.1.5.1/binaries/:$PATH
export PATH=/home/mmabry/software/:$PATH #this is for running diamond and FastTree


/home/mmabry/software/OrthoFinder-2.2.6/OrthoFinder/orthofinder/orthofinder.py -fg /storage/htc/pireslab/mmabry_htc/Brassicales_Transcriptomes_aa/FixedNames/renameSeq/Family_level/Capparidaceae/Results_Feb13/ -M msa -t 56 -ot
```
## G. Filter alignments by taxa occupancy https://github.com/MU-IRCF/filter_by_ortho_group
```perl
#!/usr/bin/env perl
use v5.10;
use strict;
use warnings;
use autodie;

use File::Copy qw(move);
#use Data::Show;

my $DEBUG = 0;

MAIN(@ARGV);

sub MAIN {

    my $in_file_name = shift // die 'Filename and number of taxa required';
    my $min_num_taxa = shift // die 'Number of taxa required';
    my $out_dir      = shift // 'FilteredTaxaAlignments';

    open(my $fh, '<', $in_file_name);

    my %ortho_groups;

    system("mkdir -p $out_dir");

    while (my $line = readline $fh) {

        #if ($line =~ m{\A [>] ( (?:\w|[-])+ ) _* }xms ) {
        if ($line =~ m{\A [>] ([a-zA-Z0-9]+ _ [a-zA-Z]+ ) }xms ) {
                say "MATCH" if $DEBUG;
                $ortho_groups{$1}++;
        }
	else {
            say "$line did not match" if $DEBUG;
        }
    }

    close($fh);

    #show %ortho_groups if $DEBUG;

    my $num_taxa = scalar keys %ortho_groups;

    #show $num_taxa if $DEBUG;

    if ($num_taxa >= $min_num_taxa ) {
        warn "Moving '$in_file_name' to '$out_dir' because it has $num_taxa (i.e. >= minimum of $min_num_taxa)\n";
        move($in_file_name, $out_dir);
    }
    else {
	warn "Not moving '$in_file_name' because it has $num_taxa (i.e. <= minimum of $min_num_taxa)\n";
    }
}

sub length_of_seq {
    my $string = shift;
    my @nts = grep { $_ ne '-'} split //,$string;
}
```
#### use script below to actually run perl script above on files, change # (i.e. 59) to number corresponding to percent occupancy wanted
```bash
srun --pty -p Interactive -c2 --mem 4G /bin/bash
module load perl/perl-5.26.2-intel
for file in *.fa; do perl filter_by_ortho_group.pl $file 56; done &> movelog_filter_by_ortho_group &
```

## H. Filter alignments based on quality/gaps https://github.com/MU-IRCF/filter_by_gap_fraction
```perl
#!/usr/bin/env perl
use v5.10;
use strict;
use warnings;
use autodie;

use File::Copy qw(move);
#use Data::Show;

my $DEBUG = 0;

MAIN(@ARGV);

sub MAIN {

    my $in_file_name       = shift // die 'Filename and number of taxa required';
    my $max_fraction_gaps  = shift // die 'Maximum fraction gaps';
    my $out_dir            = shift // 'AlignmentsFilteredByQuality';
    
    open(my $fh, '<', $in_file_name);
    
    system("mkdir -p $out_dir");
    
    my $total_gaps;
    my $total_not_gaps;

    while (my $line = readline $fh) {

        next if substr($line,0,1) eq '>';

        # remove newline
        chomp $line;

        my $num_gaps = count_gaps_in($line);

        my $length   = length($line);

        my $num_not_gaps = $length - $num_gaps;
   
        $total_gaps     += $num_gaps;
        $total_not_gaps += $num_not_gaps;
    }

    close($fh);

    my $fraction_gaps = $total_gaps/($total_gaps + $total_not_gaps);

    if ($fraction_gaps >= $max_fraction_gaps ) {
        warn "Not moving '$in_file_name' to '$out_dir' because it is $fraction_gaps gaps (i.e. >= maximum of $max_fraction_gaps)\n";
    }
    else {
        warn "Moving '$in_file_name' to '$out_dir' because it is $fraction_gaps gaps (i.e. <= maximum of $max_fraction_gaps)\n";
        move($in_file_name, $out_dir); 
    }

}

sub count_gaps_in {
    my $string = shift;

    my $num_gaps = $string =~ tr/-//;

    return $num_gaps;
}
```
#### use script below to actually run perl script above on files, change # (i.e. 0.4) to number corresponding to percent gaps user wants to allow
```bash
srun --pty -p Interactive -c2 --mem 4G /bin/bash
module load perl/perl-5.26.2-intel
for file in *.fa; do perl /storage/htc/biocompute/scratch/mmabry/filter_by_alignment_quality/filter_by_alignment_quality $file 0.4 ; done &> movelog_quality0.4_final &
```

## I. Rename headers for PhylotreePruner (in later in the steps) before running first set of trees
```bash
srun --pty -p Interactive -c2 --mem 4G /bin/bash  #this works if I do not want to be on the head node, specifically on Lewis
for file in *.fa; do sed -e "s/\(.*\)_/\1@/g" $file > PTP_fixedNames/$file ; done
```

## J. Run RAxML for tree inference
```bash
#! /bin/bash

#SBATCH -J RAxML
#SBATCH -o RAxML.o%J
#SBATCH -e RAxML.e%J
#SBATCH -N 1
#SBATCH -n 14
#SBATCH -p hpc5,BioCompute,Lewis
#SBATCH -t 1-00:00:00
#SBATCH --mem=8000

export PATH=/home/mmabry/software/standard-RAxML/:$PATH

file=${1}
file2=$(basename $file | sed 's/\.fa$//g')

raxmlHPC-PTHREADS -T 14 -f a -N 100 -p 12345 -x 12345 -n "${file2}.tre" -s "${file}" -m PROTCATWAG
```
#### call files in htc, write in hpc
```bash
for file in /storage/htc/pireslab/mmabry_htc/Brassicales_Transcriptomes_aa/FixedNames/renameSeq/Family_level/Brassicaceae/Results_Feb13/Orthologues_Feb15/Alignments/FilteredTaxaAlignments/AlignmentsFilteredByQuality/PTP_fixedNames/*.fa; do sbatch RAxML.sh ${file}; done
```

## K. PhyloTreePruner
```bash
#! /bin/bash

#SBATCH -J PTP
#SBATCH -o PTP.o%J
#SBATCH -e PTP.e%J
#SBATCH -N 1
#SBATCH -n 14
#SBATCH -p hpc5,BioCompute,Lewis
#SBATCH -t 2-00:00:00
#SBATCH --mem=8000

module load java/openjdk/java-1.8.0-openjdk

file=${1}
file2=$(basename $file | sed 's/\.fa$//g')

java -cp /home/mmabry/software/PhyloTreePruner/src_and_wrapper_scripts/ PhyloTreePruner /group/pireslab/mmabry_hpc/Brassicales/Brassicaceae/RAxML_bipartitions.${file2}.tre 10 ${file} 0.5 u
```
```bash
for file in /storage/htc/pireslab/mmabry_htc/Brassicales_Transcriptomes_aa/FixedNames/renameSeq/Family_level/Brassicaceae/Results_Feb13/Orthologues_Feb15/Alignments/FilteredTaxaAlignments/AlignmentsFilteredByQuality/PTP_fixedNames/*.fa; do sbatch PhyloTreePruner.sh ${file}; done
```
#### before running the next step remove the @ symbols
```bash
for file in *.fa; do cat $file | cut -f1 -d '@' > ${file}_cut.fa; done
```

## L. Finally, run RAxML one more time to get final gene trees
```bash
#! /bin/bash

#SBATCH -J RAxML
#SBATCH -o RAxML.o%J
#SBATCH -e RAxML.e%J
#SBATCH -N 1
#SBATCH -n 14
#SBATCH -p hpc5,BioCompute,Lewis
#SBATCH -t 1-00:00:00
#SBATCH --mem=8000

export PATH=/home/mmabry/software/standard-RAxML/:$PATH

file=${1}
file2=$(basename $file | sed 's/\.fa$//g')

raxmlHPC-PTHREADS -T 14 -f a -N 100 -p 12345 -x 12345 -n "${file2}.tre" -s "${file}" -m PROTCATWAG
```
#### call files in htc, write in hpc
```bash
for file in /storage/htc/pireslab/mmabry_htc/Brassicales_Transcriptomes_aa/FixedNames/renameSeq/Family_level/Res_Bat_Mor_Car/Results_Feb13/Orthologues_Feb13_1/Alignments/FilteredTaxaAlignments/AlignmentsFilteredByQuality/PTP_fixedNames/FinalGeneTreeAlignments/*.fa_cut.fa; do sbatch RAxML.sh ${file}; done
```

## M. Species tree inference with ASTRAL

#### in the folder with the besttree from RAxML run this to make a file with all the trees in one
```bash
cat *bestTree* > Capparidaceae.tre
```
```bash
#! /bin/bash

#SBATCH -J ASTRAL
#SBATCH -o ASTRAL.o%J
#SBATCH -e ASTRAL.e%J
#SBATCH -N 1
#SBATCH -n 14
#SBATCH --account=biosci
#SBATCH -p hpc5,BioCompute,Lewis
#SBATCH -t 01:00:00
#SBATCH --mem=80G

module load java/openjdk/java-1.8.0-openjdk

java -jar /home/mmabry/software/ASTRAL_III/Astral/astral.5.6.1.jar -i Brassicales.tre -t 2 -o Brassicales_Discordance_ASTRAL.tre
```


# 2. Repetitive Element Clustering

# 3. Regression Analyses

# 4. Hierarchical Clustering

# 5. Ultrametric Tree 
## A. Concatenate alignments using the 'concatenate_matrices.py' script from https://bitbucket.org/washjake/transcriptome_phylogeny_tools/src/master/ 
```bash
#! /bin/bash

#SBATCH -p BioCompute,hpc5  # partition
#SBATCH -J ConcatMatrix  # give the job a custom name
#SBATCH -o results-%j.out  # give the job output a custom name
#SBATCH -t 2-00:00:00  # two day time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 14  # number of cores (AKA tasks)
#SBATCH --mem=80G #memory

export PYTHONPATH="/home/mmabry/software/biopython-1.68/":$PYTHONPATH
export PATH="/home/mmabry/software/phyutility/":$PATH

python /home/mmabry/software/washjake-transcriptome_phylogeny_tools-1d8e6439be95/concatenate_matrices.py 1404Alignments 5 2 aa BrassicalesREconcatMatrix
```
## B. Run RAxML to optimize Branch lengths and model parameters using the concatenated alignment and ASTRAL phylogeny as a fixed input tree.
```bash
#! /bin/bash

#SBATCH -p BioCompute,hpc5,Lewis  # partition
#SBATCH -J RAxML  # give the job a custom name
#SBATCH -o RAxML-%j.out  # give the job output a custom name
#SBATCH -t 2-00:00:00  # two day time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 24  # number of cores (AKA tasks)
#SBATCH --mem=80G #memory

export PATH=/home/mmabry/software/standard-RAxML/:$PATH

raxmlHPC-PTHREADS -T 24 -f e -p 12345 -m PROTCATWAG -q BrassicalesREconcatMatrix.model -s BrassicalesREconcatMatrix.phy -t ../RepElem_Brassicales_ASTRAL.tre -n ConcatBrassicalesRE -o Moringa_sp_Moringa_sp,Carica_sp_Carica_sp 
```
## C. Use TreePL to time calibrate phylogeny https://github.com/blackrim/treePL/wiki/Quick-run
#### First run the script below with thorough and prime, then run a second time with the prime options
```txt
treefile = RAxML_result.ConcatBrassicalesRE.tre
smooth = 1000
numsites = 765087
mrca = PALAEOCLEOME Aethionema_arabicum_Aethionema_arabicum Cleome_violaceae_Cleome_violaceae
min = PALAEOCLEOME 47.8
max = PALAEOCLEOME 52.58
mrca = DRESSIANTHA Carica_sp_Carica_sp Batis_maritima_Batis_maritima
min = DRESSIANTHA 89.9
max = DRESSIANTHA 98.78
outfile = BrassicalesRE.dated.tre
thorough
#prime
opt = 2
moredetail
optad = 2
moredetailad
optcvad = 2
moredetailcvad
```

# 6. Bayou https://github.com/uyedaj/bayou/blob/master/tutorial.md
```R
install_github("uyedaj/bayou")
library(bayou)
library(viridis)
```
```R
##set working directory
setwd("/Users/mem2c2/OneDrive - University of Missouri/Computer/Projects/BrassicalesRE/")

####Simulating a multi-optima OU process on a phylogeny

##get Brassicales RE tree from treePL
tree <- ladderize(read.tree("treePL/BrassicalesRE.dated.tre"))

tree <- drop.tip(tree, c("Cleomella_serrulata_Cleomella_serrulata")) ## we do not have RE data for this one :(

##for Brassicaceae family tree
tree <- extract.clade(tree, 77)

##for Cleomaceae family tree
tree <- extract.clade(tree, 122)

##forjust BMAP taxa
tree <- keep.tip(tree, c("Aethionema_arabicum_Aethionema_arabicum",
"Cakile_maritima_Cakile_maritima", 
"Caulanthus_amplexicaulis_Caulanthus_amplexicaulis", 
"Crambe_hispanica_Crambe_hispanica",
"Descurania_pinnata_Descurania_pinnata", 
"Descurainia_sophioides_Descurainia_sophioides",
"Diptychocarpus_strictus_Diptychocarpus_strictus",
"Eruca_vesicaria_Eruca_vesicaria",
"Euclidium_syriacum_Euclidium_syriacum",
"Iberis_amara_Iberis_amara",
"Isatis_tinctoria_Isatis_tinctoria",
"Lepidium_sativum_Lepidium_sativum",
"Lunaria_annua_Lunaria",
"Malcomia_maritima_Malcomia_maritima",
"Myagrum_perfoliatum_Myagrum_perfoliatum",
"Rorippa_islandica_Rorippa_islandica",
"Sinapis_alba_Sinapis_alba",
"Stanleya_pinnata_Stanleya_pinnata",
"Thlaspi_arvense_Thlaspi_arvense"))


##for any tree do the following
tree <- reorder(tree, "postorder")

plot(tree, cex = 0.5)

##get RE data
df <- read.csv("data.summary2.csv")

df <- df[df$Species %in% tree$tip.label, ]

## change Total_RE to Copia or Gypsy when using those.
dat <- as.vector(df$Total_RE)
names(dat) <- df$Species

hist(dat)

###make priors
priorOU <- make.prior(tree, dists=list(dalpha="dlnorm", 
                                       dsig2="dlnorm",
                                       dk="cdpois", 
                                       dtheta="dnorm", 
                                       dsb = "dsb"), 
                    param=list(dalpha=list(), 
                               dsig2=list(), 
                               dk=list(lambda=1,kmax=2*length(tree$tip.label)-2),
                               dtheta=list(mean=mean(dat), sd=2*sd(dat))))

###initiate the MCMC chain with some starting values
startpars <- priorSim(priorOU, tree, plot=TRUE)$pars[[1]]
priorOU(startpars)

mcmcOU <- bayou.makeMCMC(tree, dat, prior=priorOU, 
                         new.dir=TRUE, outname="modelOU_RE", plot.freq=NULL) # Set up the MCMC
mcmcOU$run(1000000) # Run the MCMC

###full MCMC results are written to a set of files
chainOU <- mcmcOU$load()

###We can set a "burnin" parameter that tells the package coda to discard the first bit of the chain
chainOU <- set.burnin(chainOU, 0.3)
summary(chainOU)
plot(chainOU, auto.layout=FALSE)
dev.off()

###check out results vs truth
#par(mfrow=c(3,1))
plotSimmap.mcmc(chainOU, burnin = 0.3, lwd = 2, pp.cutoff = 0.3, pal = viridis, circle.col = "#33638DFF")
plotBranchHeatMap(tree, chainOU, "theta", burnin = 0.3, pal = viridis, edge.width = 2)
phenogram.density(tree, dat, burnin = 0.3, chainOU, pp.cutoff = 0.3)
```

# 7. Owie : Example from from http://www.phytools.org/Cordoba2017/ex/10/Multi-regime.html
```R
library(phytools)
```
```R
##set working directory
setwd("/Users/mem2c2/OneDrive - University of Missouri/Computer/Projects/BrassicalesRE/")

tree <- ladderize(read.tree("treePL/BrassicalesRE.dated.tre"))

tree <- drop.tip(tree, c("Cleomella_serrulata_Cleomella_serrulata"))

##for Brassicaceae family tree
tree <- extract.clade(tree, 78)

##for Cleomaceae family tree
tree <- extract.clade(tree, 122)

##for BMAP taxa
tree <- keep.tip(tree, c("Aethionema_arabicum_Aethionema_arabicum", "Cakile_maritima_Cakile_maritima", "Caulanthus_amplexicaulis_Caulanthus_amplexicaulis", 
"Crambe_hispanica_Crambe_hispanica" , "Descurania_pinnata_Descurania_pinnata" , "Descurainia_sophioides_Descurainia_sophioides", 
"Diptychocarpus_strictus_Diptychocarpus_strictus", "Eruca_vesicaria_Eruca_vesicaria" , "Euclidium_syriacum_Euclidium_syriacum", 
"Iberis_amara_Iberis_amara", "Isatis_tinctoria_Isatis_tinctoria", "Lepidium_sativum_Lepidium_sativum", "Lunaria_annua_Lunaria", 
"Malcomia_maritima_Malcomia_maritima", "Myagrum_perfoliatum_Myagrum_perfoliatum", "Rorippa_islandica_Rorippa_islandica", "Sinapis_alba_Sinapis_alba", 
"Stanleya_pinnata_Stanleya_pinnata", "Thlaspi_arvense_Thlaspi_arvense"))

##for data
df <- read.csv("data.summary2.csv")
df <- df[df$Species %in% tree$tip.label, ]

plot(tree, cex = 0.5)

##Seletion regime for testing atAlpha vs all other Brassicales
atAlpha <- as.vector(df$Atalpha)
names(atAlpha) <- df$Species

atAlphaTree <- make.simmap(tree, atAlpha, model="ER")
plot(atAlphaTree)

##Seletion regime for testing thAlpha vs all other Cleomeaceae
Thalpha <-as.vector(df$Thalpha)
names(Thalpha) <- df$Species

ThalphaTree <- make.simmap(tree, Thalpha, model="ER")
plot(ThalphaTree)

##Seletion regime for testing Brassiceae vs all other Brassicaceae
Brassiceae <-as.vector(df$Brassiceae)
names(Brassiceae) <- df$Species

BrassiceaeTree <- make.simmap(tree, Brassiceae, model="ER")
plot(BrassiceaeTree)

##Seletion regime for testing only known ploidy of BMAP taxa
Neo <- as.vector(c("Diploid",
                   "NeoPolyploid",
                   "NeoPolyploid",
                   "NeoPolyploid",
                   "NeoPolyploid",
                   "Diploid",
                   "Diploid",
                   "NeoPolyploid",
                   "Diploid", 
                   "NeoPolyploid",
                   "NeoPolyploid",
                   "NeoPolyploid",
                   "NeoPolyploid", 
                   "Diploid",
                   "Diploid",
                   "Diploid",
                   "NeoPolyploid", 
                   "Diploid",
                   "Diploid"))
names(Neo) <- c("Aethionema_arabicum_Aethionema_arabicum",
                "Cakile_maritima_Cakile_maritima",
                "Caulanthus_amplexicaulis_Caulanthus_amplexicaulis",
                "Crambe_hispanica_Crambe_hispanica",
                "Descurania_pinnata_Descurania_pinnata",
                "Descurainia_sophioides_Descurainia_sophioides",
                "Diptychocarpus_strictus_Diptychocarpus_strictus",
                "Eruca_vesicaria_Eruca_vesicaria",
                "Euclidium_syriacum_Euclidium_syriacum", 
                "Iberis_amara_Iberis_amara",
                "Isatis_tinctoria_Isatis_tinctoria",
                "Lepidium_sativum_Lepidium_sativum",
                "Lunaria_annua_Lunaria", 
                "Malcomia_maritima_Malcomia_maritima",
                "Myagrum_perfoliatum_Myagrum_perfoliatum",
                "Rorippa_islandica_Rorippa_islandica",
                "Sinapis_alba_Sinapis_alba", 
                "Stanleya_pinnata_Stanleya_pinnata",
                "Thlaspi_arvense_Thlaspi_arvense")

NeoTree <- make.simmap(tree, Neo, model="ER")
plot(NeoTree)


##Setting data
TEs <- as.vector(df$Total_RE)
names(TEs) <- df$Species

Gypsy <- as.vector(df$Gypsy)
names(Gypsy) <- df$Species

Copia <- as.vector(df$Copia)
names(Copia) <- df$Species
```
```R
install.packages("OUwie")
library(OUwie)

##change depending on which selection regime is being used
OU_df<-data.frame(df$Species, atAlpha, TEs)
OU_df<-data.frame(df$Species, atAlpha, Gypsy)
OU_df<-data.frame(df$Species, atAlpha, Copia)

##change tree depending on which selection regime is being used

#Ornstein-Uhlenbeck model with different state means and a single alpha and sigma^2 acting all selective regimes
fitOUM<-OUwie(atAlphaTree, OU_df,model="OUM",simmap.tree=TRUE) 
fitOUM

#Ornstein-Uhlenbeck model with a single optimum for all species 
fitOU1<-OUwie(atAlphaTree, OU_df,model="OU1",simmap.tree=TRUE) 
fitOU1

#single-rate Brownian motion
fitBM1<-OUwie(atAlphaTree, OU_df,model="BM1",simmap.tree=TRUE) 
fitBM1

#Brownian motion with different rate parameters for each state on a tree
fitBMS<-OUwie(atAlphaTree, OU_df,model="BMS",simmap.tree=TRUE, root.station=FALSE) 
fitBMS

#new Ornstein-Uhlenbeck models that assume different state means as well as multiple sigma^2
fitOUMV<-OUwie(atAlphaTree, OU_df,model="OUMV",simmap.tree=TRUE) 
fitOUMV

#new Ornstein-Uhlenbeck models that assume different state means as well as multiple alpha
fitOUMA<-OUwie(atAlphaTree, OU_df,model="OUMA",simmap.tree=TRUE) 
fitOUMA

#new Ornstein-Uhlenbeck models that assume different state means as well as multiple alpha and sigma^2 per selective regime 
fitOUMVA<-OUwie(atAlphaTree, OU_df,model="OUMVA",simmap.tree=TRUE) 
fitOUMVA

#gather AICc scores
ouwie_aicc<-c(fitOUM$AICc, fitOU1$AICc, fitBM1$AICc, fitBMS$AICc, fitOUMV$AICc, fitOUMA$AICc, fitOUMVA$AICc)
names(ouwie_aicc)<-c("fitOUM","fitOU1","fitBM1", "fitBMS",  "fitOUMV", "fitOUMA", "fitOUMVA")

#determine which one is weighted highest for all tests run
aic.w(ouwie_aicc)
```


