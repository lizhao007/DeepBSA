# DeepBSA
DeepBSA is a novel bulked segregant analysis (BSA) software for the dissection of complex traits. Two brand-new algorithms are developed in DeepBSA named deep learning (DL) and k-value (K), which can be applied on different number (at least 2) of bulked pools. DeepBSA also integrates five widely used algorithms - ED4, delta SNP_index, G', Ridit and SmoothLOD, and DL performs better than them with absolute bias and signal-to-noise ratio in our simulation. Overall, DeepBSA provides a user-friendly, OS-compatible, and all-in-one pipeline, which do not need sophisticated bioinformatics skills for BSA.

Installation
The “Instruction or Manual” file of software was packed into the DeepBSA.zip.

Input
The input file for DeepBSA is the VCF file, which contains genomic variants for all bulked pools. For the genomic variant calling, we'd love to recommendate using GATK using the guided bioinformatic pipeline as follows:

#Taking two mixed pools as examples
#mapping

bwa mem -t 8 -M -P Referencegenome.fa High_Forward.fastq High_Reverse.fastq >bsa_H.sam
bwa mem -t 8 -M -P Referencegenome.fa Low_Forward.fastq Low_Reverse.fastq >bsa_L.sam

#pretreatment for GATK SNP calling for hight pool

java -jar ${EBROOTPICARD}/picard.jar CleanSam INPUT=bsa_H.sam OUTPUT=bsa_H_cleaned.sam
java -jar ${EBROOTPICARD}/picard.jar FixMateInformation INPUT=bsa_H_cleaned.sam OUTPUT=bsa_H_cleaned_fixed.sam SO=coordinate

java -jar ${EBROOTPICARD}/picard.jar AddOrReplaceReadGroups INPUT=bsa_H_cleaned_fixed.sam OUTPUT=bsa_H_cleaned_fixed_group.bam LB=bsa_H SO=coordinate RGPL=illumina PU=barcode SM=bsa_H

samtools index bsa_H_cleaned_fixed_group.bam

java -jar ${EBROOTPICARD}/picard.jar MarkDuplicatesWithMateCigar INPUT=bsa_H_cleaned_fixed_group.bam OUTPUT=bsa_H_cleaned_fixed_group_DEDUP.bam M=bsa_H_cleaned_fixed_group_DEDUP.mx AS=true REMOVE_DUPLICATES=true MINIMUM_DISTANCE=500

samtools index bsa_H_cleaned_fixed_group_DEDUP.bam

#pretreatment for GATK SNP calling for low pool

Similar to high pool
#genomic variant calling

java -Xmx64g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R Referencegenome.fa -nct 8 -I bsa_H_cleaned_fixed_group_DEDUP.bam -I bsa_L_cleaned_fixed_group_DEDUP.bam -o bsa_H_L_snps_indels.vcf

Cite
Li Z., Chen X., Shi S., Zhang H., Wang X., Chen H., Li W., and Li L. (2022). DeepBSA: A deep-learning algorithm improves bulked segregant analysis for dissecting complex traits. Mol. Plant. doi: https://doi.org/10.1016/j.molp.2022.08.004.
https://www.cell.com/molecular-plant/abstract/S1674-2052(22)00267-2
