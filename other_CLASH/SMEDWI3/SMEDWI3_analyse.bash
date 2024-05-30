#######################################################################################################
# SMEDWI3 CLASH
#######################################################################################################
#Planaria CLASH
#genome
#https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_022537955.1/
# annotation
# https://zenodo.org/records/5149023#.Ybk9jn3MK3I

#unused genome
#https://planosphere.stowers.org/pub/fasta/smed_chr_ref_v1.fa

# rRNA annotation
#https://www.sciencedirect.com/science/article/pii/S2213596017300223




# planaria piRNA
cd ~/Desktop/srrfiles/;
for X in SRR8162647 SRR8162648 SRR8162649;do
fasterq-dump -e 8 ${X};
done;
#CLASH
cd ~/Desktop/srrfiles/;
for X in SRR8842976 SRR8842977 SRR8842978 SRR8842979 SRR8842980 SRR8842981;do
fasterq-dump -e 8 ${X};
done;
cd ~/Desktop/srrfiles/;
for X in SRR8842982 SRR8842983 SRR8842984 SRR8842985 SRR8842986 SRR8842987;do
fasterq-dump -e 8 ${X};
done;



cd /media/hermione/SMEDWI3/;
mkdir Discas/;
mkdir Process/;
mkdir endtoend/;
mkdir fastqc/;
fastqc -t 12 *.fastq.gz -o fastqc/;



cd /media/hermione/SMEDWI3/;
mkdir trimed/;
for X in SRR8162647 SRR8162648 SRR8162649;do
cutadapt -j 10 -q 30 -m 50 -o Discas/${X}_c.fastq.gz ${X}.fastq.gz;
cutadapt -j 10 -a TGGAATTCTCGGGTGCCAAGG -m 26  -o Discas/${X}_c2.fastq.gz Discas/${X}_c.fastq.gz;
zcat Discas/${X}_c2.fastq.gz | awk -f solexa2fasta.awk | \
awk -f fasta2tab.awk > Discas/${X}_clipped_qf.tab;
perl make_comp_fasta.pl Discas/${X}_clipped_qf.tab > Discas/${X}_comp.fasta;
seqkit fx2tab Discas/${X}_comp.fasta  | cut -f 1-2  > Process/${X}_comp.tsv;
done;


cd /media/hermione/SMEDWI3/;
mkdir ref/;
mkdir Bowtie2Index/;
mkdir Process2/;
seqkit fx2tab database/smed_genome.fa | cut -f 1-2 > database/smed_genome.tsv;
bowtie2-build --quiet -f database/smed_genome.fa ref/smed;

# rRNA annotation
mkdir blastnIndex/;
makeblastdb -in database/smed_genome.fa -dbtype nucl \
-out blastnIndex/genome -parse_seqids;
blastn -word_size 20 -evalue 1e-3 -outfmt 6 -db blastnIndex/genome \
 -query database/raw/rRNA_others_raw.fasta -out database/raw/rRNA_align_out.txt;

/usr/bin/Rscript  rRNA_blast.R

bedtools getfasta -s -fi database/smed_genome.fa -bed database/raw/rRNA_align.bed -fo database/raw/rRNA_all.fasta;
bowtie2-build --quiet -f database/raw/rRNA_all.fasta ref/rRNA;


samtools faidx database/raw/dd_Smed_v6.pcf.contigs.fasta;
seqkit seq --rna2dna database/raw/arb-silva.de_2024-01-22_id1299249_tax_silva_trunc.fasta -o database/raw/rRNA_silva.fasta



for X in SRR8842976 SRR8842977 SRR8842978 SRR8842979 SRR8842980 SRR8842981 \
SRR8842982 SRR8842983 SRR8842984 SRR8842985 SRR8842986 SRR8842987 ;do
cutadapt -j 10 -q 30  -a TGGAATTCTCGGGTGCCAAGG -m 45 \
-o Discas/${X}_trim_pre.fastq.gz ${X}.fastq.gz;
done;
cd /media/hermione/SMEDWI3/;
for X in SRR8842976 SRR8842977 SRR8842978 SRR8842979 SRR8842980 SRR8842981 \
SRR8842982 SRR8842983 SRR8842984 SRR8842985 SRR8842986 SRR8842987 ;do
bowtie2 -p 20 -x ref/smed --un Discas/${X}_trim.fastq \
-U Discas/${X}_trim_pre.fastq.gz -S Discas/${X}_rmstRNA.sam;
pigz -p 10 Discas/${X}_trim.fastq;
done;

cd /media/hermione/SMEDWI3/;
for X in SRR8842976 SRR8842977 SRR8842978 SRR8842979 SRR8842980 SRR8842981 \
SRR8842982 SRR8842983 SRR8842984 SRR8842985 SRR8842986 SRR8842987;do
zcat Discas/${X}_trim.fastq.gz | awk -f solexa2fasta.awk | \
awk -f fasta2tab.awk > Discas/${X}_clipped_qf.tab;
perl make_comp_fasta.pl Discas/${X}_clipped_qf.tab > Discas/${X}_comp.fasta;
seqkit fx2tab Discas/${X}_comp.fasta | cut -f 1-2  > Discas/${X}_del.tsv;
bowtie2 -p 20 --very-sensitive -f -U Discas/${X}_comp.fasta --reorder \
 -x ref/smed -S Discas/${X}.sam;
samtools view -@ 20 -F 4 Discas/${X}.sam | cut -f 1-6 > Discas/${X}_map2.sam;
/usr/bin/Rscript  filter_endtoend.R Discas/${X}_map2.sam Discas/${X}_del.tsv Discas/${X}_comp2.fasta;
seqkit fx2tab Discas/${X}_comp2.fasta  | cut -f 1-2 > Discas/${X}_comp2.tsv;
done;

#go to pi_database_make
/usr/bin/Rscript  pi_database_make.R;

mkdir BowtieIndex/;
find ./Discas -name "*_comp2.fasta" | sed 's/_comp2.fasta//' | parallel --max-procs=16 \
 'bowtie-build --quiet -f Discas/{/.}_comp2.fasta BowtieIndex/{/.}_CLASH_ref';
for Xdef in $(ls Discas/. | grep _comp2.fasta$) ; do
X=`echo ${Xdef} | sed 's/_comp2.fasta//'`;
bowtie -p 20 -f -a -v 0 -S BowtieIndex/${X}_CLASH_ref database/piRNA_all.fasta Discas/${X}.sam;
samtools view -@ 20 -F 4 -F 16 Discas/${X}.sam | cut -f 1-6 > Process2/${X}_map.sam;
done;
find ./Discas -name "*_comp2.fasta" | sed 's/_comp2.fasta//' | parallel --max-procs=16 \
  'seqkit fx2tab Discas/{/.}_comp2.fasta | cut -f 1-2 > Process2/{/.}_comp.tsv; \
  /usr/bin/Rscript  piclash_pra200308.R Process2/{/.}_map.sam Process2/{/.}_comp.tsv database/piRNA_all.tsv Process2/{/.}_comp_tf.fasta';



###################
cd /media/hermione/SMEDWI3/;
mkdir Process3/;
mkdir vienna/;
for Xdef in $(ls Process2/. | grep _comp_tf.fasta$) ; do
X=`echo ${Xdef} | sed 's/_comp_tf.fasta//'`;
bowtie2 -p 20 -f -x ref/smed -U Process2/${X}_comp_tf.fasta -S Discas/${X}_tf.sam;
samtools view -@ 20 -F 4 Discas/${X}_tf.sam | cut -f 1-6 > Process3/${X}_tf_map.sam;
done;

find ./Process2 -name "*_comp_tf.fasta" | sed 's/_comp_tf.fasta//' | parallel --max-procs=16 \
 'seqkit fx2tab Process2/{/.}_comp_tf.fasta | cut -f 1-2 > Process3/{/.}_comp_tf.tsv;
 /usr/bin/Rscript  piclash_tf_sep3.R Process3/{/.}_tf_map.sam Process3/{/.}_comp_tf.tsv database/piRNA_all.tsv Process3/{/.}_chimera;
  RNAplex --temp=26 < Process3/{/.}_chimera_rnaplex.fasta > vienna/{/.}.rnaplex 2> /dev/null;
  /usr/bin/Rscript  piclash_make_vienna2.R Process3/{/.}_chimera.tsv vienna/{/.}.rnaplex vienna/{/.};';


mkdir vienna1fig/;

# piRNA mapping to rRNA
mkdir rRNA_filter/;
bowtie2 -p 20 -f -x ref/rRNA -U database/piRNA_all.fasta -S Discas/rRNA_piRNA.sam;
samtools view -@ 20 -F 4 Discas/rRNA_piRNA.sam | cut -f 1-6 > rRNA_filter/rRNA_piRNA.sam;

makeblastdb -in database/raw/rRNA_all.fasta -dbtype nucl \
-out blastnIndex/rRNA -parse_seqids;
blastn -word_size 10 -evalue 1e-3 -outfmt 6 -db blastnIndex/rRNA \
 -query database/piRNA_all.fasta -out database/raw/piRNA_rRNA.txt;

makeblastdb -in database/raw/rRNA_silva.fasta -dbtype nucl \
-out blastnIndex/rRNA_silva -parse_seqids;
blastn -word_size 10 -evalue 1e-3 -outfmt 6 -db blastnIndex/rRNA_silva \
 -query database/piRNA_all.fasta -out database/raw/piRNA_rRNA_silva.txt;

bowtie2 -p 20 -f -x ref/smed -U database/piRNA_all.fasta -S Discas/piRNA_genome.sam;
samtools view -@ 20 -f 4 Discas/piRNA_genome.sam | cut -f 1-6 > rRNA_filter/piRNA_genome_remove.sam;

/usr/bin/Rscript  pivienna_analyse.R;

