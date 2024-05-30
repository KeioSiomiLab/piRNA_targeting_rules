

###################################################################################################################
# testis Aub
###################################################################################################################

#testis Aub
cd /media/hermione/otherdata_compare/Aub_testis;

mkdir fastqc/;
fastqc --quiet -t 20 *.fastq.gz -o fastqc/;

mkdir Discas/;
mkdir trim/;
mkdir figures/;
mkdir otherdata/;
mkdir rnaplex/;
mkdir vienna/;
mkdir blastout/;

for X in SRR12213370 SRR12213371  ; do
cutadapt -j 4 --max-n 0 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -m 23 -M 30 \
-o Discas/${X}_trim.fastq.gz  ${X}.fastq.gz;
bowtie2 -p 8 -x /media/pericles/CLASH/database/Bowtie2Index/dm6_rmstRNA --un Discas/${X}_removed.fastq \
-U Discas/${X}_trim.fastq.gz -S Discas/${X}_rmstRNA.sam;
bowtie2 -p 8 -x /media/pericles/CLASH/database/Bowtie2Index/miRNA --un Discas/${X}_removed2.fastq \
-U Discas/${X}_removed.fastq -S Discas/${X}_rmstRNA.sam;
seqkit fq2fa Discas/${X}_removed2.fastq -o ${X}.fa;
done;

/usr/bin/Rscript  Aub_testis_make_piRNA_db.R;

blastn -word_size 15 -evalue 1e-3 -outfmt 6 -db /media/hermione/piRNAmodel/blastnIndex/TE \
 -query piRNA_all.fasta -out blastout/piRNA.txt;
blastn -word_size 15 -evalue 1e-3 -outfmt 6 -db /media/hermione/piRNAmodel/blastnIndex/gene \
 -query piRNA_all.fasta -out blastout/gene.txt;

/usr/bin/Rscript  /media/hermione/otherdata_compare/fly_piRNA_plex_prepare.R /media/hermione/otherdata_compare/Aub_testis/;

parallel --max-procs=16 \
 'RNAplex --temp=26 < rnaplex/{/.}.fasta > vienna/{/.}.rnaplex 2> /dev/null;
 /usr/bin/Rscript /media/hermione/otherdata_compare/piRNA_vienna_prepare.R \
 {/.} /media/hermione/otherdata_compare/Aub_testis/;' ::: model_RNAplex model_RNAplex2;

bedtools intersect -wa -wb -a rnaplex/gene_remove1.bed \
-b /media/pericles/CLASH/database/anno/dm6_gene_ano.bed > rnaplex/gene_remove1_overlap.tsv;

/usr/bin/Rscript /media/hermione/otherdata_compare/piRNA_vienna_table.R /media/hermione/otherdata_compare/Aub_testis/;

bedtools merge -i /media/hermione/otherdata_compare/Aub_testis/vienna/All_model_vienna_pre.bed \
> /media/hermione/otherdata_compare/Aub_testis/vienna/All_model_vienna_merge.bed;

cd /media/hermione/otherdata_compare/Aub_testis;
mkdir RNA/;
mkdir featurecount/;
mkdir TE/;
for X in  SRR12216401 SRR12216402 SRR12216403 ;do
cutadapt -j 10 -q 20 -m 20 -o Discas/${X}_trim.fastq.gz ${X}.fastq.gz;
STAR runMode alignReads --runThreadN 16 --genomeDir /media/pericles/CLASH/database/star_genome_dm6/ \
--readFilesIn Discas/${X}_trim.fastq.gz  \
--outSAMattributes All --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM \
--outFilterScoreMin 10  --readFilesCommand zcat --alignEndsType EndToEnd --outFileNamePrefix RNA/${X}_STAR;
samtools index RNA/${X}_STARAligned.sortedByCoord.out.bam;
source activate meme
featureCounts -T 24 -t exon -g gene_id --fraction -M -O -d 50 -D 2000 \
-a /media/pericles/CLASH/database/gtf/dm6_use.gtf -o Discas/${X}_featurecount.txt \
RNA/${X}_STARAligned.sortedByCoord.out.bam;
cut -f1,6,7 Discas/${X}_featurecount.txt \
> featurecount/${X}_featurecount.txt;
conda deactivate;
bowtie2 -p 20 -a -x /media/pericles/CLASH/database/Bowtie2Index/TEref -U Discas/${X}_trim.fastq.gz -S Discas/${X}_TE.sam;
samtools view -@ 20 -bS Discas/${X}_TE.sam > Discas/${X}_TE.bam;
samtools sort -@ 20 -m 4G Discas/${X}_TE.bam > Discas/${X}_TE_sort.bam;
samtools index Discas/${X}_TE_sort.bam;
samtools idxstats Discas/${X}_TE_sort.bam > TE/${X}_TEcount.txt;
done;
seqkit stats -j 20  Discas/*_trim.fastq.gz -T  > read_info.tsv;

/usr/bin/Rscript /media/hermione/otherdata_compare/Aub_testis/Aub_testis_analysis.R


###################################################################################################################
# embryo Aub
###################################################################################################################

#embryo Aub

cd /media/hermione/otherdata_compare/Aub_embryo/;

mkdir fastqc/;
fastqc --quiet -t 20 *.fastq.gz -o fastqc/;

mkdir Discas/;
mkdir trim/;
mkdir figures/;
mkdir otherdata/;
mkdir rnaplex/;
mkdir vienna/;
mkdir blastout/;

for X in SRR3051363  ; do
cutadapt -j 4 --max-n 0 -m 23 -M 30 \
-o Discas/${X}_trim.fastq.gz  ${X}.fastq.gz;

bowtie2 -p 8 -x /media/pericles/CLASH/database/Bowtie2Index/dm6_rmstRNA --un Discas/${X}_removed.fastq \
-U Discas/${X}_trim.fastq.gz -S Discas/${X}_rmstRNA.sam;
bowtie2 -p 8 -x /media/pericles/CLASH/database/Bowtie2Index/miRNA --un Discas/${X}_removed2.fastq \
-U Discas/${X}_removed.fastq -S Discas/${X}_rmstRNA.sam;

seqkit fq2fa Discas/${X}_removed2.fastq  -o ${X}.fa;
done;

/usr/bin/Rscript  Aub_embryo_make_piRNA_db.R;

blastn -word_size 15 -evalue 1e-3 -outfmt 6 -db /media/hermione/piRNAmodel/blastnIndex/TE \
 -query piRNA_all.fasta -out blastout/piRNA.txt;
blastn -word_size 15 -evalue 1e-3 -outfmt 6 -db /media/hermione/piRNAmodel/blastnIndex/gene \
 -query piRNA_all.fasta -out blastout/gene.txt;

/usr/bin/Rscript  /media/hermione/otherdata_compare/fly_piRNA_plex_prepare.R /media/hermione/otherdata_compare/Aub_embryo/;

parallel --max-procs=16 \
 'RNAplex --temp=26 < rnaplex/{/.}.fasta > vienna/{/.}.rnaplex 2> /dev/null;
 /usr/bin/Rscript /media/hermione/otherdata_compare/piRNA_vienna_prepare.R \
 {/.} /media/hermione/otherdata_compare/Aub_embryo/;' ::: model_RNAplex model_RNAplex2;

bedtools intersect -wa -wb -a rnaplex/gene_remove1.bed \
-b /media/pericles/CLASH/database/anno/dm6_gene_ano.bed > rnaplex/gene_remove1_overlap.tsv;

/usr/bin/Rscript /media/hermione/otherdata_compare/piRNA_vienna_table.R /media/hermione/otherdata_compare/Aub_embryo/;

bedtools merge -i /media/hermione/otherdata_compare/Aub_embryo/vienna/All_model_vienna_pre.bed \
> /media/hermione/otherdata_compare/Aub_embryo/vienna/All_model_vienna_merge.bed;

mkdir RNA/;
mkdir featurecount/;
mkdir TE/;
for X in  SRR3051356 SRR3051357 SRR3051361 SRR3051362 ;do
cutadapt -j 10 -q 20 -m 20 -o Discas/${X}_trim.fastq.gz ${X}.fastq.gz;
STAR runMode alignReads --runThreadN 16 --genomeDir /media/pericles/CLASH/database/star_genome_dm6/ \
--readFilesIn Discas/${X}_trim.fastq.gz  \
--outSAMattributes All --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM \
--outFilterScoreMin 10  --readFilesCommand zcat --alignEndsType EndToEnd --outFileNamePrefix RNA/${X}_STAR;
samtools index RNA/${X}_STARAligned.sortedByCoord.out.bam;
source activate meme
featureCounts -T 24 -t exon -g gene_id --fraction -M -O -d 50 -D 2000 \
-a /media/pericles/CLASH/database/gtf/dm6_use.gtf -o Discas/${X}_featurecount.txt \
RNA/${X}_STARAligned.sortedByCoord.out.bam;
cut -f1,6,7 Discas/${X}_featurecount.txt \
> featurecount/${X}_featurecount.txt;
conda deactivate;
bowtie2 -p 20 -a -x /media/pericles/CLASH/database/Bowtie2Index/TEref -U Discas/${X}_trim.fastq.gz -S Discas/${X}_TE.sam;
samtools view -@ 20 -bS Discas/${X}_TE.sam > Discas/${X}_TE.bam;
samtools sort -@ 20 -m 4G Discas/${X}_TE.bam > Discas/${X}_TE_sort.bam;
samtools index Discas/${X}_TE_sort.bam;
samtools idxstats Discas/${X}_TE_sort.bam > TE/${X}_TEcount.txt;
done;
seqkit stats -j 20  Discas/*_trim.fastq.gz -T  > read_info.tsv;
/usr/bin/Rscript /media/hermione/otherdata_compare/Aub_embryo/Aub_embryo_analysis.R
###################################################################################################################
# Miwi testis
###################################################################################################################

# Miwi testis
fasterq-dump -e 8 SRR363958;

# spike-in
#UGCUAGUCUGUUAUCGACCUGACCUCAUAG
#UGCUAGUCUGUUCGAUACCUGACCUCAUAG
#UGCUAGUCUGUUGUCACGAAGACCUCAUAG
#UGCUAGUCUUAUCGACCUCCUCAUAG
#UGCUAGUCUUCGAUACCUCCUCAUAG
#UGCUAGUCUUGUCACGAACCUCAUAG
#UGCUAGUUAUCGACCUUCAUAG
#UGCUAGUUCGAUACCUUCAUAG
#UGCUAGUUGUCACGAAUCAUAG


#NNNCGANNNUACNNN
#NNNAUCNNNAGUNNN


#SRR21528491 SRR21528492 SRR21528493 SRR21528494 SRR21528495 SRR21528497 SRR21528498 SRR21528499 \
#SRR21528500 SRR21528501 SRR21528502 SRR21528503
cd /media/hermione/otherdata_compare/Miwi/;
mkdir fastqc/;
fastqc --quiet -t 20 *.fastq.gz -o fastqc/;

mkdir Discas/;
mkdir trim/;
mkdir figures/;
mkdir otherdata/;
mkdir rnaplex/;
mkdir vienna/;
mkdir blastout/;

#miRNA
#https://mirbase.org/download/


cd /media/hermione/otherdata_compare/Miwi/;
seqkit seq --rna2dna database/miRNA_pre.fasta -o database/miRNA.fasta

mkdir database/Bowtie2Index/;
bowtie2-build --quiet -f database/miRNA.fasta database/Bowtie2Index/miRNA;
bowtie2-build --quiet -f database/rRNA.fasta database/Bowtie2Index/rRNA;


cd /media/hermione/otherdata_compare/Miwi/;
mkdir database/blastnIndex/;
makeblastdb -in /media/qqprosperodd/miranda/Refee/database/dfam_TE.fasta -dbtype nucl \
-out database/blastnIndex/TE -parse_seqids;
makeblastdb -in /media/qqprosperodd/miranda/Refee/database/mm10.fa -dbtype nucl \
-out database/blastnIndex/genome -parse_seqids;
seqkit fx2tab /media/qqprosperodd/miranda/Refee/database/dfam_TE.fasta | cut -f 1-2 > TE.tsv;
seqkit fx2tab /media/qqprosperodd/miranda/Refee/database/mm10.fa | cut -f 1-2 > genome.tsv;

mkdir database/RepMas/;
source activate repeatmasker;
RepeatMasker -s -pa 18 -e rmblast -gff -norna -nolow -no_is \
-lib /media/qqprosperodd/miranda/Refee/database/dfam_TE.fasta \
-dir database/RepMas/ \
/media/qqprosperodd/miranda/Refee/database/mm10.fa \
> Discas/Repmas_of_log.txt;
conda deactivate;

/usr/bin/Rscript  Mm_db.R;




fastqc --quiet -t 20 SRR21528*.fastq.gz -o fastqc/;

# 5' cut seq
#NNNCGANNNUACNNN
#NNNAUCNNNAGUNNN
# 3' cut seq
#NNNGTCNNNTAGNNN

cd /media/hermione/otherdata_compare/Miwi/;
for X in SRR21528491 SRR21528492 SRR21528493 SRR21528494 SRR21528495 SRR21528497 SRR21528498 SRR21528499 \
SRR21528500 SRR21528501 SRR21528502 SRR21528503  ; do
cutadapt -j 10 -q 30 -m 50 -o Discas/${X}_c.fastq.gz ${X}.fastq.gz;
cutadapt -j 10 -a TGGAATTCTCGGGTGCCAAGG -m 26  -o Discas/${X}_c2.fastq.gz Discas/${X}_c.fastq.gz;
seqkit rmdup -s Discas/${X}_c2.fastq.gz -o Discas/${X}_c3.fastq.gz;
seqkit grep -j 10 -sirp '^[ATGC][ATGC][ATGC]CGA[ATGC][ATGC][ATGC]TAC' Discas/${X}_c3.fastq.gz -o Discas/${X}_c4_adaptor1.fastq.gz;
seqkit grep -j 10 -sirp '^[ATGC][ATGC][ATGC]ATC[ATGC][ATGC][ATGC]AGT' Discas/${X}_c3.fastq.gz -o Discas/${X}_c4_adaptor2.fastq.gz;
zcat Discas/${X}_c4_adaptor*.fastq.gz  | pigz -p 6 > Discas/${X}_c4.fastq.gz;
cutadapt -u 15 -u -15 -j 10 --max-n 0 -m 22 -M 34 -o Discas/${X}_trim.fastq.gz  Discas/${X}_c4.fastq.gz;

bowtie2 -p 8 -x database/Bowtie2Index/rRNA --un Discas/${X}_removed.fastq \
-U Discas/${X}_trim.fastq.gz -S Discas/${X}_rmstRNA.sam;
bowtie2 -p 8 -x database/Bowtie2Index/miRNA --un Discas/${X}_removed2.fastq \
-U Discas/${X}_removed.fastq -S Discas/${X}_rmstRNA.sam;
seqkit fq2fa Discas/${X}_removed.fastq  -o ${X}.fa;
done;

/usr/bin/Rscript  Miwi_make_piRNA_db.R;

blastn -word_size 15 -evalue 1e-3 -outfmt 6 -db database/blastnIndex/TE \
 -query piRNA_all.fasta -out blastout/piRNA.txt;
blastn -word_size 15 -evalue 1e-3 -outfmt 6 -db database/blastnIndex/genome \
 -query piRNA_all.fasta -out blastout/gene.txt;

/usr/bin/Rscript  /media/hermione/otherdata_compare/Mm_piRNA_plex_prepare.R /media/hermione/otherdata_compare/Miwi/;

parallel --max-procs=16 \
 'RNAplex --temp=26 < rnaplex/{/.}.fasta > vienna/{/.}.rnaplex 2> /dev/null;
 /usr/bin/Rscript /media/hermione/otherdata_compare/piRNA_vienna_prepare.R \
 {/.} /media/hermione/otherdata_compare/Miwi/;' ::: model_RNAplex model_RNAplex2;

bedtools intersect -wa -wb -a rnaplex/gene_remove1.bed \
-b /media/hermione/otherdata_compare/Miwi/database/Mm_TE_ano.bed > rnaplex/gene_remove1_overlap.tsv;

bedtools intersect -wa -wb -a rnaplex/gene_remove1.bed \
-b /media/hermione/otherdata_compare/Miwi/database/Mm_gene_ano.bed > rnaplex/exon_overlap.tsv;

bedtools intersect -wa -wb -a rnaplex/gene_remove1.bed \
-b /media/hermione/otherdata_compare/Miwi/database/gene.bed > rnaplex/gene_overlap.tsv;


/usr/bin/Rscript /media/hermione/otherdata_compare/piRNA_vienna_table_Miwi.R /media/hermione/otherdata_compare/Miwi/;

bedtools merge -i /media/hermione/otherdata_compare/Miwi/vienna/All_model_vienna_pre.bed \
> /media/hermione/otherdata_compare/Miwi/vienna/All_model_vienna_merge.bed;



#RNA
mkdir RNA/;
mkdir featurecount/;
mkdir TE/;
for X in  SRR363965 SRR363967 ;do
cutadapt -j 10 -q 20 -m 20 -o Discas/${X}_trim.fastq.gz ${X}.fastq.gz;
STAR runMode alignReads --runThreadN 16 --genomeDir /media/qqprosperodd/miranda/Refee/database/STAR_index_TE/ \
--readFilesIn Discas/${X}_trim.fastq.gz --outReadsUnmapped Fastx \
--outFileNamePrefix TE/${X}_STAR --outSAMattributes All --outSAMtype BAM SortedByCoordinate  \
--outFilterScoreMin 10  --readFilesCommand zcat --alignEndsType EndToEnd ;
pigz -f -p 16 TE/${X}_STARUnmapped.out*;
samtools index TE/${X}_STARAligned.sortedByCoord.out.bam;

STAR runMode alignReads --runThreadN 16 --genomeDir /media/qqprosperodd/miranda/Refee/database/STAR_index/ \
--readFilesIn TE/${X}_STARUnmapped.out.mate1.gz \
 --outFileNamePrefix RNA/${X}_STAR --outSAMattributes All --outSAMtype BAM SortedByCoordinate  \
--outFilterScoreMin 10  --readFilesCommand zcat --alignEndsType EndToEnd;
samtools index RNA/${X}_STARAligned.sortedByCoord.out.bam;
source activate meme
featureCounts -T 24 -t exon -g gene_id --fraction -M -O -d 50 -D 2000 \
-a /media/qqprosperodd/miranda/Refee/database/gencode.vM25.annotation.gtf -o Discas/${X}_featurecount.txt \
RNA/${X}_STARAligned.sortedByCoord.out.bam;
cut -f1,6,7 Discas/${X}_featurecount.txt \
> featurecount/${X}_featurecount.txt;
conda deactivate;
samtools idxstats TE/${X}_STARAligned.sortedByCoord.out.bam > TE/${X}_TEcount.txt;
done;
seqkit stats -j 20  Discas/*_trim.fastq.gz -T  > read_info.tsv;



/usr/bin/Rscript /media/hermione/otherdata_compare/Miwi/Miwi_analysis.R



