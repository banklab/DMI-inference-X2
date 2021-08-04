#!/bin/sh

###### This shell script contains most key simulation and analyses scripts. To avoid confusion, small scripts are only available upon request.

wdir=`pwd`
code=$wdir/code
data=$wdir/data/
output=$wdir/output/

####################################################
### Toy model: haploid population with recombination.
#### R/Rscript version 4.
:<<A
Two-locus DMI: Assume parental hapotype are Ab and aB.
DMI happens between A and B. Direct selection coefficient at A is aplha, at B is beta.
Epistatic coefficient is -gamma. Recombination between two loci is r.
Admixture proportion at generation 0 is f(aB)
Migration rate can be changed within the code or add as a command argument
Run: Rscript $code/haploidWrecombining.R aplha beta gamma r f output_dir
A

Rscript $code/haploidWrecombining.R 0.001 0.002 0.5 0.1 0.5 $output/1/

########################################################
### SLiM simulations.
### SLiM version version 3.5.
:<<B
Why did we use the “recessive” scenario?
Modeling DMIs in a diploid population is more complicated than suggested by the toy model, which assumes a haploid population with recombination. In a diploid population, the epistatic interaction does depend on the dominance of each incompatible allele. Turelli and Orr cleverly avoided the continuous distribution of dominance in their paper by focusing on backcross and Haldane’s rule inferred from sibling species deference mapping or introgression experiments (Turelli and Orr 2000). We can see that DMIs are caused by interactions that include at least one recessive allele. The incompatibilities may be caused by two homozygous loci (homozygous-homozygous genotype, h0 as the selection coefficient) like Haldane’s rule, which is that (h0 is much smaller than the selection coefficient of homozygous-heterozygous genotype (h1 as the selection coefficient, h0<<h1). Or it may be caused by heterozygous-homozygous genotypes, as in the backcross, in which h1 is much smaller than the selection coefficient of heterozygous-heterozygous genotypes (h2 as the selection coefficient, h1<<h2, defined as recessive scenario in Blanckaert and Bank 2018). The codominance DMI scenario is proposed (h0<h1<h2, Bank et al. 2012, Blanckaert and Bank 2018). All the three scenarios are listed in Table S2. Recessive can protect disadvantageous alleles from quick elimination. Since the fitness of F1 individuals affects the formation of hybrid populations, we use the h1<<h2 scenario to demonstrate statistical power and demographic effects.
B
##### 1. Wright-Fisher model, non-overlapping generations;
:<<C
Here we provide our script on two DMIs randomly distributed on four chromosomes, and with two neutral chromosome in the same gemone. See method section "False and true positive rates with one and two DMI pairs". All the other senario can be adapt from this script in our paper. 
Two DMI (four DMI loci) are randomly distributed on four chromosome. Each chromosome carries 100 makers, the recombination rate betweeen two neighboring markers is 1%. Namely, the genetic map is one marker per centimorgan. In addtion, two neutral crhomosomes are simulated. Therefore, 600 markers in total.
Chromosomes cannot be simulated independently. They are seperated by recombination rate =0.5. Here, 6 chromosomes are simulated in one genome. Marker 0-399 are four chomorosomes with random two DMIs. 400-599 are two neutral chromosomes. In "1. Parameters", DMI1_loc1(randomNum) is draw from a uniform distribution [0,1], which is the relative position at four chromosomes. DMI1_loc1(randomNum)X400 is DMI1_loc1(genomeCoordinate). SLiM assumes a chromosome starts from 0. For convience, in our script, we assume a chromosome starts from 1.
In this script, genomes of 500 individuals from generation 30 and generation 50 are recorded.
Main SLiM script: code/dmi.WF.twoDMI.6chr.h0h1.slim
Check the script for paramters.
Output include
1. Parameters
DMI1_Epi(gamma)	DMI1_loc1(randomNum)	DMI1_loc2	DMI2_Epi	DMI2_loc1	DMI2_loc2	DMI1_loc1(genomeCoordinate)	DMI1_loc2	DMI2_loc1	DMI2_loc2	DMI1_alpha	DMI1_beta	DMI2_alpha	DMI2_beta
-0.327676	0.046949	0.988007	-0.353637	0.698941	0.0342439	18	395	279	13	0.00118887	0.00277652	0.00223266	0.00277652
2.DMI frequencies
3. Genomes of 500 individuals.
C
# for example: migration rate is 0.01, admixture proportion is 0.3 and 0.7.
m=0.01
startingFreq=0.3
slim -d "outPath='$output/2/sim.demo'" -d mP1=$m -d freqP1=$startingFreq $code/dmi.WF.twoDMI.6chr.h0h1.slim

prefix=$output/2/sim.demo.F30
numLoci=600
## write indiviudal genomes.
perl $code/genome-01.pl $prefix $prefix $numLoci
## use 200 genomes for DMI statistics.
dRegion=400
head -200 $prefix.01 | perl -nae 'for($i=0; $i<'$dRegion'; $i+=1){print " $F[$i]";}print "\n";' >$prefix.200
Rscript $code/genome2stat.R $prefix.200 $prefix.ms.stats


##### 2. Moran model, overlapping generations
### One DMI. Add age structure.
### Parameters are very similar to others. See scripts.
m=0.01
startingFreq=0.3
slim -d "outPath='$output/3/'" -d mP1=$m -d freqP1=$startingFreq $code/dmi.Moran.oneDMI.slim


#######################################################
### Custom DMI simulation script for sensitivity and specificity
### python version 3
##### each run output 100 simulations with the same parameters
##### Check parameters: python3 simDMI.oneDMI.WF.py -h

python3 $code/simDMI.oneDMI.WF.py --gamma -0.5 --output $output/4/sim100
Rscript $code/cal.basic.sensSpec.R $output/4/sim100 $output/4/sim100.stat

#######################################################
### Demo: statistics for fish data.
### genotype 0 is allele from Xiphophorus birchmanni and 1 from X. malinche
### Input is the plink .ped and .map files. Heterozygous genotypes have been changed into missing data. Here, we extracted hundrands markers on chromosome 20 from TLMC as an example.

prefix=$output/5/TLMC_chr20
Rscript $code/fishStat.permutation.R $prefix.map $prefix.ped $prefix.out

#######################################################
