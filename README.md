# Background
DeepBSA is a novel bulked segregant analysis (BSA) software for the dissection of complex traits. Two brand-new algorithms are developed in DeepBSA named deep learning (DL) and k-value (K), which can be applied on different number (at least 2) of bulked pools. DeepBSA also integrates five widely used algorithms - ED4, delta SNP_index, G', Ridit and SmoothLOD, and DL performs better than them with absolute bias and signal-to-noise ratio in our simulation. Overall, DeepBSA provides a user-friendly, OS-compatible, and all-in-one pipeline, which do not need sophisticated bioinformatics skills for BSA.

# Installation
DeepBSA is available for both Windows and Linux, and the download link is: http://zeasystemsbio.hzau.edu.cn/tools.html. The alternate cloud download link is:
链接：https://pan.baidu.com/s/1PbqOu5fDXK2RU5Hi3G4p6A?pwd=c71e 
提取码：c71e

# Update history
2022.11.15 version1.4: Improving the function of Simulator and offering the software for Linux.

2022.08.30 version1.3: Adding PDF file of mapping result and CSV file of algorithm value.

2022.08.16 version1.2

2022.07.25 version1.1

# Input
The input file for DeepBSA is the VCF file, which contains genomic variants for all bulked pools. For the genomic variant calling, we'd love to recommendate using GATK using the guided bioinformatic pipeline as follows:

```
***Taking two mixed pools as examples***
##building reference index
samtools faidx Referencegenome.fa
bwa index Referencegenome.fa

##mapping
bwa mem -t 8 -M -P Referencegenome.fa High_Forward.fastq High_Reverse.fastq >bsa_H.sam
bwa mem -t 8 -M -P Referencegenome.fa Low_Forward.fastq Low_Reverse.fastq >bsa_L.sam

##pretreatment for GATK SNP calling for hight pool
java -jar ${EBROOTPICARD}/picard.jar CleanSam INPUT=bsa_H.sam OUTPUT=bsa_H_cleaned.sam
java -jar ${EBROOTPICARD}/picard.jar FixMateInformation INPUT=bsa_H_cleaned.sam OUTPUT=bsa_H_cleaned_fixed.sam SO=coordinate
java -jar ${EBROOTPICARD}/picard.jar AddOrReplaceReadGroups INPUT=bsa_H_cleaned_fixed.sam OUTPUT=bsa_H_cleaned_fixed_group.bam LB=bsa_H SO=coordinate RGPL=illumina PU=barcode SM=bsa_H
samtools index bsa_H_cleaned_fixed_group.bam
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicatesWithMateCigar INPUT=bsa_H_cleaned_fixed_group.bam OUTPUT=bsa_H_cleaned_fixed_group_DEDUP.bam M=bsa_H_cleaned_fixed_group_DEDUP.mx AS=true REMOVE_DUPLICATES=true MINIMUM_DISTANCE=500
samtools index bsa_H_cleaned_fixed_group_DEDUP.bam

##same pretreatment for GATK SNP calling for low pool

##genomic variant calling
java -Xmx64g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R Referencegenome.fa -nct 8 -I bsa_H_cleaned_fixed_group_DEDUP.bam -I bsa_L_cleaned_fixed_group_DEDUP.bam -o bsa_H_L_snps_indels.vcf
```

# Usage
## For windows

The “Instruction or Manual” file can be download in github and it is also packed into the DeepBSA_windows.zip.
## For linux

### Requirment
R and Python 3.7(or greater) should be installed. Other require python packages can be quickly installed by running "./requirment.txt" in main dictory as follows.
```
#Install
wget http://zeasystemsbio.hzau.edu.cn/Tools/DeepBSA_linux_v1.4.tar.gz
tar -xvzf DeepBSA_linux_v1.4.tar.gz
cd DeepBSA_linux_v1.4/
./requirment.txt

#QTL mapping 
cd bin/
python3 main.py -h

#usage: main.py [-h] --i I [--m M] [--p P] [--p1 P1] [--p2 P2] [--p3 P3] [--s S] [--w W] [--t T]
 optional arguments:
  -h, --help  show this help message and exit
  --i I       The input file path(vcf/csv).
  --m M       The algorithm(DL/K/ED4/SNP/SmoothG/SmoothLOD/Ridit) used. Default is DL.
  --p P       Whether to pretreatment data(1[True] or 0[False]). Default is True.
  --p1 P1     Pretreatment step 1: Number of read thread, the SNP whose number lower than it will be filtered. Default is 0.
  --p2 P2     Pretreatment step 2: Chi-square test(1[True] or 0[False]). Default is 1[True].
  --p3 P3     Pretreatment step 3: Continuity test(1[True] or 0[False]). Default is 1[True].
  --s S       The function to smooth the result(Tri-kernel-smooth\LOWESS\Moving Average), Defalut is LOWESS
  --w W       Windows size of LOESS. The number is range from 0-1. 0 presents the best size for minimum AICc. Default is 0(auto).
  --t T       The threshold to find peaks(float). Default is 0(auto)

#Data simulation
cd DeepBSA_linux_v1.4/bin/
python3 simulate_progress.py -h

#usage: simulate_progress.py [-h] --i I --p P --r R --e E --s S
 optional arguments:
  -h, --help  show this help message and exit
  --i I       individual
  --p P       pools
  --r R       ratio
  --e E       effective points
  --s S       save path

```
More details for parameters can be got in the “Instruction or Manual” file.

# Cite

Li Z., Chen X., Shi S., Zhang H., Wang X., Chen H., Li W., and Li L. (2022). DeepBSA: A deep-learning algorithm improves bulked segregant analysis for dissecting complex traits. Mol. Plant. doi: https://doi.org/10.1016/j.molp.2022.08.004.

