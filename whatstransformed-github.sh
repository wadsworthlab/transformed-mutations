#####################################################################
#GC Transformant Analysis: What's Transformed
#This script can be used to identify DNA tracts in recipient cell lines that have been inherited from donor strains.
#See: Wadsworth CB, Sater MRA, Bhattacharyya R, Grad Y. 2019. Impact of population structure in the design of RNA-based diagnostics for antibiotic resistance in Neisseria gonorrhoeae. Antimicrobial Agents and Chemotherapy 63(8):e00549-19. https://doi.org/10.1128/AAC.00549-19
##############################################################################


############################################
#de novo assembly via Spades
############################################

#Constructing a 28-BL Reference Genome using PE reads in Spades de novo
#.yaml file for 2 libraries 
[
{
orientation: "fr",
type: "paired-end",
right reads: [ 
"28-BL_S16_L001_R1_001.fastq"
],
left reads: [
"28-BL_S16_L001_R2_001.fastq"
]
},
{
orientation: "fr",
type: "paired-end",
right reads: [ 
"28-BL_S14_L001_R1_001.fastq"
],
left reads: [
"28-BL_S14_L001_R2_001.fastq"
]
}
]

module load dev/python/2.7.6 
SPAdes-3.7.0-Linux/bin/spades.py -t 4 --dataset 28-BL.yaml -o assembly_28-BL_reference

#Generate library stats with screed
module load dev/python/2.7.6 
python screed/assemstats.py 20 assembly_28-BL_reference/scaffolds.fasta > assembly_28-BL_reference/28-BL_lib_stats.txt

############################################
#Improve 28-BL reference with Pilon
############################################

module load seq/bowtie/2.2.4
bowtie2-build assembly_28-BL_reference/scaffolds.fasta assembly_28-BL_reference/28-BL_Index

bowtie2 --end-to-end --very-sensitive -x assembly_28-BL_reference/28-BL_Index -1 28-BL_S16_L001_R1_001.fastq -2 28-BL_S16_L001_R2_001.fastq -S 28-BL_improved.sam

module load seq/samtools/1.3
samtools view -bS -o 28-BL_improved.bam 28-BL_improved.sam
samtools sort -o 28-BL_improved.sorted.bam 28-BL_improved.bam
samtools index 28-BL_improved.sorted.bam

module load dev/java/jdk1.8
java -Xmx16G -jar pilon-1.16.jar --genome assembly_28-BL_reference/scaffolds.fasta --frags 28-BL_improved.sorted.bam --output 28-BL_improved --vcf --fix all --tracks --minqual 15


############################################
#Pilon to call variants from improved 28-BL reference
#complete for transformants and donor strains
############################################

#Build reference
bowtie2-build 28-BL_improved.fasta assembly_28-BL_reference/28-BL_Index_Improved

#!/bin/bash
module load seq/bowtie/2.2.4
for i in $(ls *_001.fastq | rev | cut -c 14- | rev | uniq | cut -c 73-)
do
bsub -q long -W 99:00 -R "rusage[mem=12000]" bowtie2 --end-to-end --very-sensitive -x assembly_28-BL_reference/28-BL_Index_Improved -1 transformants/${i}_R1_001.fastq -2 transformants/${i}_R2_001.fastq -S ${i}_improved.sam
done

chmod +x align28-BL.sh
./align28-BL.sh

#!/bin/bash
module load seq/samtools/1.3
for i in $(ls *.sam | rev | cut -c 5- | rev)
do
samtools view -bS -o ${i}.bam ${i}.sam
done

chmod +x bs.sh
bsub -q priority -W 4:00 -R "rusage[mem=20000]" ./bs.sh

#!/bin/bash
module load seq/samtools/1.3
for i in $(ls *.sam | rev | cut -c 5- | rev)
do
samtools sort -o ${i}.sorted.bam ${i}.bam
done

chmod +x sort.sh
bsub -q priority -W 4:00 -R "rusage[mem=20000]" ./sort.sh

#!/bin/bash
less 
for i in $(ls *.sorted.bam)
do
samtools index ${i}
done 

chmod +x bamindex.sh
bsub -q priority -W 4:00 ./bamindex.sh

#!/bin/bash
module load dev/java/jdk1.8
for i in $(ls *.sorted.bam | rev | cut -c 12- | rev)
do
java -Xmx16G -jar pilon-1.16.jar --genome assembly_28-BL_reference/28-BL_improved.fasta --frags ${i}.sorted.bam --output ${i} --vcf --tracks --minqual 15
done

chmod +x pilon.sh
./pilon.sh

#Select the sites that are different from the reference

#!/bin/bash
for i in $(ls *.vcf)
do
perl -anle 'if ($F[4]!~/\./){print;}' ${i} | grep -E '(PASS)' > ${i}.filtered
done

chmod +x compile_pilon.sh
bsub -q priority -W 20:00 ./compile_pilon.sh

############################################
#Identify  transformed regions that have been inherited in all transformant lines (replicates), to ID mutations presumably causal to a phenotype
############################################

#Run in R 
#Remove variant calls for 28Bl mapped back to 28Bl (presumably biased sequencing errors)
BL_14 <- read.csv("28-BL_S14_L001_improved_improved.vcf.filtered", sep="\t", head=F)
BL_16 <- read.csv("28-BL_S16_L001_improved_improved.vcf.filtered", sep="\t", head=F)
BL_14 <- BL_14[,-3:-10]
BL_16 <- BL_16[,-3:-10]
BL_14$node_position <- with(BL_14, paste(V1, V2, sep="_"))
BL_16$node_position <- with(BL_16, paste(V1, V2, sep="_"))
incorrect <- merge(BL_16, BL_14, by="node_position", all.y=F, all.x=F)
incorrect$x <- rep("x",3)
incorrect <- incorrect[,-2:-5]

p1_r1_p1 <- read.csv("p1-r1-p1_S1_L001_improved_improved.vcf.filtered", sep="\t", head=F)
p1_r1_p1$node_position <- with(p1_r1_p1, paste(V1, V2, sep="_"))
p1_r1_p1_x <- merge(p1_r1_p1, incorrect, by="node_position", all.x=T, all.y=T)
wordList <- c("x")
p1_r1_p1_final <- subset(p1_r1_p1_x, !(x %in% wordList))
write.csv(p1_r1_p1_final, "~/Desktop/pilon_improved_cleaned_28-BL/p1_r1_p1_improved_cleaned.csv", row.names = FALSE)

p1_r4_p1 <- read.csv("p1_r4_p1_improved_improved.csv", sep="\t", head=F)
p1_r4_p1$node_position <- with(p1_r4_p1, paste(V1, V2, sep="_"))
p1_r4_p1_x <- merge(p1_r4_p1, incorrect, by="node_position", all.x=T, all.y=T)
wordList <- c("x")
p1_r4_p1_final <- subset(p1_r4_p1_x, !(x %in% wordList))
write.csv(p1_r4_p1_final, "~/Desktop/pilon_improved_cleaned_28-BL/p1_r4_p1_improved_cleaned.csv", row.names = FALSE)

#Repeat for all transforming replicates, then merge to find shared sites

shared_sites <- merge(p1_r4_p1_final, p1_r4_p2_final, by="node_position", all.x=F, all.y=F)
shared_sites <- merge(shared_sites, p1_r4_p3_final, by="node_position", all.x=F, all.y=F)
shared_sites <- merge(shared_sites, p1_r4_p4_final, by="node_position", all.x=F, all.y=F)
shared_sites <- merge(shared_sites, p2_r1_p1_final, by="node_position", all.x=F, all.y=F)
shared_sites <- merge(shared_sites, p2_r2_p1_final, by="node_position", all.x=F, all.y=F)
shared_sites <- merge(shared_sites, p2_r2_p2_final, by="node_position", all.x=F, all.y=F)
shared_sites <- merge(shared_sites, p2_r2_p3_final, by="node_position", all.x=F, all.y=F)
shared_sites <- merge(shared_sites, p1_r4_p1_final, by="node_position", all.x=F, all.y=F)

write.csv(shared_sites, "~/Desktop/pilon_improved_cleaned_28-BL/LAX/shared_sites_LAX.csv",row.names = FALSE)

