###################################################################################################################
# genemap (RNAseq)
###################################################################################################################
cd /media/pericles/CLASH/RNAseq/;
mkdir Discas/;
mkdir fil/;
for X in SRR609666 SRR609665 ; do
cutadapt -j 10 -m 38 -a AGATCGGAAGAG -o Discas/${X}_trimed2.fastq.gz ${X}.fastq.gz;
cutadapt -j 10 -m 30 -u 4 -u -4 -o Discas/${X}_trimed.fastq.gz Discas/${X}_trimed2.fastq.gz;
cutadapt -j 10 -q 20 -m 25 -o Discas/${X}_trim.fastq.gz Discas/${X}_trimed.fastq.gz;
done;

for X in SRR609666 SRR609665 ; do
seqkit rmdup -s Discas/${X}_trim.fastq.gz -o Discas/${X}_trim2.fastq.gz;
bowtie2 -p 20 -x /media/pericles/CLASH/database/Bowtie2Index/dm6_genome -U Discas/${X}_trim2.fastq.gz -S Discas/${X}.sam;
samtools view -@ 20 -bS Discas/${X}.sam > Discas/${X}.bam;
samtools sort -@ 20 -m 4G Discas/${X}.bam > genome/${X}_collapse_sort.bam;
samtools index genome/${X}_collapse_sort.bam;
done;

cd /media/pericles/CLASH/RNAseq/;
mkdir GROcount/;
# count gene
for X in SRR609666 SRR609665 ; do
bedtools intersect -v -abam genome/${X}_collapse_sort.bam \
-b /media/pericles/TEfind/rmblast/dm6_rmblast_mask.bed > \
genome/${X}_remove.bam;
samtools sort -@ 20 -m 4G genome/${X}_remove.bam > genome/${X}_remove_sort.bam;
samtools index genome/${X}_remove_sort.bam;
source activate meme
featureCounts -T 24 -s 1 -t gene -g gene_id --fraction -M -O -d 20 -D 2000 \
-a /media/pericles/CLASH/database/gtf/dm6_use.gtf -o Discas/${X}_featurecount.txt \
genome/${X}_remove_sort.bam;
cut -f1,6,7 Discas/${X}_featurecount.txt \
> GROcount/${X}_featurecount.txt;
conda deactivate;
done;

# TE count
mkdir TE/;
for X in SRR609666 SRR609665 ; do
bowtie2 -p 20 -a -x /media/pericles/CLASH/database/Bowtie2Index/TEref \
-U Discas/${X}_trim2.fastq.gz -S Discas/${X}_TE.sam;
samtools view -@ 20 -bS Discas/${X}_TE.sam > Discas/${X}_TE.bam;
samtools sort -@ 20 -m 4G Discas/${X}_TE.bam > TE/${X}_TE_sort.bam;
samtools index TE/${X}_TE_sort.bam;
samtools idxstats TE/${X}_TE_sort.bam > GROcount/${X}_TE.txt;
done;




