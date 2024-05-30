###################################################################################################################
# genemap (CLIP)
###################################################################################################################
cd /media/pericles/CLASH/CLIPseq/;
mkdir Discas/;
fastqc -t 20 *.fastq.gz -o Discas/;
mkdir gene/;
mkdir bedgraph/;
mkdir bedgraph/;
mkdir fil/;
for X in PiwiCLIP_IN PiwiCLIP_IP_rep2 PiwiCLIP_IP_rep1 ; do
cutadapt -j 10 -q 20 -m 100 -o Discas/${X}_fill.fastq.gz ${X}.fastq.gz;
cutadapt -j 10 -m 27 --max-n 0 -a TGGAATTCTCGGGTGCCAAGG -o Discas/${X}_trimed.fastq.gz Discas/${X}_fill.fastq.gz;
seqkit rmdup -s Discas/${X}_trimed.fastq.gz -o Discas/${X}_trimed2.fastq.gz;
cutadapt -j 10 -u 9 -m 18 -o Discas/${X}_f.fastq.gz Discas/${X}_trimed2.fastq.gz;
bowtie2 -p 20 --very-sensitive -N 1 -a -x /media/pericles/CLASH/database/Bowtie2Index/dm6_gene -U Discas/${X}_f.fastq.gz -S Discas/${X}.sam;
samtools view -@ 20 -F 4 -F 16 -bS Discas/${X}.sam > Discas/${X}.bam;
samtools sort -@ 20 -m 4G Discas/${X}.bam > Discas/${X}_gene_sort.bam;
samtools index Discas/${X}_gene_sort.bam;
bedtools intersect -v -abam Discas/${X}_gene_sort.bam \
-b /media/pericles/CLASH/database/anno/dm6_gene_ano2_pre.bed > fil/${X}_gene_fil.bam;
samtools index fil/${X}_gene_fil.bam;
samtools idxstats fil/${X}_gene_fil.bam > gene/${X}_count.txt;
source activate macs;
bamCoverage -p 20 -of bedgraph -bs 1 --normalizeUsing None -b fil/${X}_gene_fil.bam -o bedgraph/${X}.bedgraph;
conda deactivate;
done;

seqkit stats -j 20 Discas/*_f.fastq.gz -T > read_info.tsv;
cd /media/pericles/CLASH/CLIPseq/;
mkdir fastqc/;
#fastqc --quiet -t 20 *.fastq.gz -o fastqc/;
