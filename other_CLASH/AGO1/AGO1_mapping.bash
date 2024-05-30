
cd /media/hermione/AGO1_CLASH/;
mkdir Discas/;
mkdir Process/;
mkdir endtoend/;
mkdir fastqc/;
fastqc -t 6 *.fastq -o fastqc/;
for X in SRR959751 SRR959752 SRR959753 SRR959754 SRR959755 SRR959756 SRR959757 SRR959758 SRR959759;do
fasterq-dump -e 6 ${X};
done;

for X in  SRR959760 SRR959761 SRR959762 SRR959763 ;do
fasterq-dump -e 6 ${X};
done;
for X in SRR959751 SRR959756;do
fasterq-dump -e 6 ${X};
done;

cd /media/hermione/AGO1_CLASH/;
mkdir trimed/;
for X in SRR959751;do
cutadapt -j 10 -q 30 -a TGGAATTCTCGGGTGCCAAGGC -m 21 -o trimed/${X}_trim.fastq.gz ${X}.fastq.gz;
done;
for X in SRR959752 SRR959753;do
seqkit grep -j 10 -sirp ^AC ${X}.fastq.gz -o Discas/${X}fill.fastq.gz;
cutadapt -j 10 -u 2 -o Discas/${X}_trim.fastq.gz Discas/${X}fill.fastq.gz;
cutadapt -j 10 -q 30 -a TGGAATTCTCGGGTGCCAAGGC -m 21 -o trimed/${X}_trim.fastq.gz Discas/${X}_trim.fastq.gz;
done;
for X in SRR959754 SRR959755;do
seqkit grep -j 10 -sirp ^GA ${X}.fastq.gz -o Discas/${X}fill.fastq.gz;
cutadapt -j 10 -u 2 -o Discas/${X}_trim.fastq.gz Discas/${X}fill.fastq.gz;
cutadapt -j 10 -q 30 -a TGGAATTCTCGGGTGCCAAGGC -m 21 -o trimed/${X}_trim.fastq.gz Discas/${X}_trim.fastq.gz;
done;
for X in SRR959756 SRR959757;do
seqkit grep -j 10 -sirp ^CACAGC ${X}.fastq.gz -o Discas/${X}fill.fastq.gz;
cutadapt -j 10 -u 6 -o Discas/${X}_trim.fastq.gz Discas/${X}fill.fastq.gz;
cutadapt -j 10 -q 30 -a TGGAATTCTCGGGTGCCAAGGC -m 21 -o trimed/${X}_trim.fastq.gz Discas/${X}_trim.fastq.gz;
done;
for X in SRR959758 ;do
cutadapt -j 10 -u 3 ${X}.fastq.gz -o Discas/${X}fill1.fastq.gz;
seqkit grep -j 10 -sirp ^TAAGC Discas/${X}fill1.fastq.gz -o Discas/${X}fill.fastq.gz;
cutadapt -j 10 -u 5 -o Discas/${X}_trim.fastq.gz Discas/${X}fill.fastq.gz;
cutadapt -j 10 -q 30 -a TGGAATTCTCGGGTGCCAAGGC -m 21 -o trimed/${X}_trim.fastq.gz Discas/${X}_trim.fastq.gz;
done;
for X in SRR959759 ;do
cutadapt -j 10 -u 3 ${X}.fastq.gz -o Discas/${X}fill1.fastq.gz;
seqkit grep -j 10 -sirp ^ATTAGC Discas/${X}fill1.fastq.gz -o Discas/${X}fill.fastq.gz;
cutadapt -j 10 -u 6 -o Discas/${X}_trim.fastq.gz Discas/${X}fill.fastq.gz;
cutadapt -j 10 -q 30 -a TGGAATTCTCGGGTGCCAAGGC -m 21 -o trimed/${X}_trim.fastq.gz Discas/${X}_trim.fastq.gz;
done;

cd /media/hermione/AGO1_CLASH/;
bowtie2-build --quiet -f ref/hOH7.fasta ref/bowtie2/hOH7;
bowtie2-build --quiet -f ref/hOH7-microRNA.fasta ref/bowtie2/miRNA;
for X in hOH7 hOH7-microRNA; do
seqkit fx2tab ref/${X}.fasta | cut -f 1-2 > ref/${X}.tsv
done;
for X in SRR959751 SRR959752 SRR959753 SRR959754 SRR959755 SRR959756 SRR959757 SRR959758 SRR959759;do
zcat trimed/${X}_trim.fastq.gz | awk -f solexa2fasta.awk | \
awk -f fasta2tab.awk > Process/${X}_clipped_qf.tab;
perl make_comp_fasta.pl Process/${X}_clipped_qf.tab > Process/${X}_comp.fasta;
done;
cd /media/hermione/AGO1_CLASH/;
mkdir Bowtie2Index/;
mkdir Process2/;
for Xdef in $(ls Process/. | grep _comp.fasta$) ; do
X=`echo ${Xdef} | sed 's/_comp.fasta//'`;
bowtie2-build --quiet -f Process/${X}_comp.fasta Bowtie2Index/${X}_CLASH_ref;
bowtie2 -p 20 -f -a -x Bowtie2Index/${X}_CLASH_ref -U ref/hOH7-microRNA.fasta -S Discas/${X}.sam;
samtools view -@ 20 -F 4 -F 16 Discas/${X}.sam  | cut -f 1-6 > Process2/${X}_map.sam;
done;
find ./Process -name "*_comp.fasta" | sed 's/_comp.fasta//' | parallel --max-procs=16 \
'seqkit fx2tab Process/{/.}_comp.fasta | cut -f 1-2  > Process2/{/.}_comp.tsv;
/usr/bin/Rscript  miclash_pra200308.R Process2/{/.}_map.sam Process2/{/.}_comp.tsv ref/hOH7-microRNA.tsv Process2/{/.}_comp_tf.fasta;';


###################
cd /media/hermione/AGO1_CLASH/;
mkdir Process3/;
mkdir vienna/;
for Xdef in $(ls Process2/. | grep _comp_tf.fasta$) ; do
X=`echo ${Xdef} | sed 's/_comp_tf.fasta//'`;
seqkit fx2tab Process2/${X}_comp_tf.fasta | cut -f 1-2 > Process3/${X}_comp_tf.tsv;
bowtie2 -p 20 -f -a -x ref/bowtie2/hOH7 -U Process2/${X}_comp_tf.fasta -S Discas/${X}_tf.sam;
samtools view -@ 20 -F 4 -F 16 Discas/${X}_tf.sam  | cut -f 1-6 > Process3/${X}_tf_map.sam;
done;

find ./Process2 -name "*_comp_tf.fasta" | sed 's/_comp_tf.fasta//' | parallel --max-procs=16 \
 '/usr/bin/Rscript  miclash_tf_sep3.R Process3/{/.}_tf_map.sam Process3/{/.}_comp_tf.tsv ref/hOH7-microRNA.tsv Process3/{/.}_chimera;
 RNAplex --temp=26 < Process3/{/.}_chimera_rnaplex.fasta > vienna/{/.}.rnaplex 2> /dev/null;
 /usr/bin/Rscript --silent --slave --vanilla miclash_make_vienna2.R Process3/{/.}_chimera.tsv vienna/{/.}.rnaplex vienna/{/.};';

/usr/bin/Rscript  mivienna_analyse.R;
