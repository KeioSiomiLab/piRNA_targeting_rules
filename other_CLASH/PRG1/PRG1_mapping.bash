cd /media/hermione/PRG1_CLASH/;
mkdir Discas/;
mkdir Process/;
mkdir endtoend/;
mkdir fastqc/;

for X in  SRR6512652 SRR6512653 SRR6512654 SRR6512655 SRR538357 SRR2140770 ;do
fasterq-dump -e 6 ${X};
done;
fastqc -t 6 *.fastq -o fastqc/;

pigz -p 6 *.fastq;
for X in  SRR6512653 SRR6512654 SRR538357;do
fasterq-dump -e 6 ${X};
done;
fastqc -t 6 *.fastq -o fastqc/;
pigz -p 6 *.fastq;



cd /media/hermione/PRG1_CLASH/;
mkdir trimed/;
for X in SRR538357;do
cutadapt -j 10 -q 30 -m 36 -o Discas/${X}_t.fastq.gz ${X}.fastq.gz;
cutadapt -j 10 -u -4 -u 11 -m 21 -o trimed/${X}_trim.fastq.gz Discas/${X}_t.fastq.gz;
zcat trimed/${X}_trim.fastq.gz | awk -f solexa2fasta.awk | \
awk -f fasta2tab.awk > Process/${X}_clipped_qf.tab;
perl make_comp_fasta.pl Process/${X}_clipped_qf.tab > Process/${X}_comp.fasta;
seqkit fx2tab Process/${X}_comp.fasta  | cut -f 1-2  > Process/${X}_comp.tsv;
done;
for X in SRR2140770;do
cutadapt -j 10 -q 30 -m 50 -o Discas/${X}_t.fastq.gz ${X}.fastq.gz;
cutadapt -j 10 -u -4 -u 25 -m 21 -o trimed/${X}_trim.fastq.gz Discas/${X}_t.fastq.gz;
zcat trimed/${X}_trim.fastq.gz | awk -f solexa2fasta.awk | \
awk -f fasta2tab.awk > Process/${X}_clipped_qf.tab;
perl make_comp_fasta.pl Process/${X}_clipped_qf.tab > Process/${X}_comp.fasta;
seqkit fx2tab Process/${X}_comp.fasta  | cut -f 1-2  > Process/${X}_comp.tsv;
done;

for X in SRR6512652 SRR6512653 SRR6512654 SRR6512655;do
cutadapt -j 10 -q 30  -a AGATCCTCGGCCGCGACC -m 22 -o trimed/${X}_trim.fastq.gz ${X}.fastq.gz;
done;
mkdir stats/;
fastqc  --quiet -t 20 trimed/*_trim.fastq.gz -o stats/;

seqkit fx2tab database/WS276_genome.fa | cut -f 1-2 > database/WS276_genome.tsv;

cd /media/hermione/PRG1_CLASH/;
mkdir ref/;
mkdir Bowtie2Index/;
mkdir Process2/;
bowtie2-build --quiet -f database/WS276_genome.fa ref/WS276;
for X in SRR6512652 SRR6512653 SRR6512654 SRR6512655;do
zcat trimed/${X}_trim.fastq.gz | awk -f solexa2fasta.awk | \
awk -f fasta2tab.awk > Process/${X}_clipped_qf.tab;
perl make_comp_fasta.pl Process/${X}_clipped_qf.tab > Process/${X}_comp.fasta;
seqkit fx2tab Process/${X}_comp.fasta | cut -f 1-2  > Process/${X}_del.tsv;
bowtie2 -p 20 --very-sensitive -f -U Process/${X}_comp.fasta --reorder \
 -x ref/WS276 -S endtoend/${X}.sam;
samtools view -@ 20 -F 4 endtoend/${X}.sam | cut -f 1-6 > endtoend/${X}_map2.sam;
/usr/bin/Rscript  filter_endtoend.R endtoend/${X}_map2.sam Process/${X}_del.tsv Process/${X}_comp2.fasta;
seqkit fx2tab Process/${X}_comp2.fasta  | cut -f 1-2 > Process2/${X}_comp2.tsv;
done;

#go to pi_database_make
/usr/bin/Rscript  pi_database_make.R;
###################
cd /media/hermione/PRG1_CLASH/;
mkdir Process3/;
mkdir vienna/;
for Xdef in $(ls Process2/. | grep _comp_tf.fasta$) ; do
X=`echo ${Xdef} | sed 's/_comp_tf.fasta//'`;
bowtie2 -p 20 -f -a -x ref/WS276 -U Process2/${X}_comp_tf.fasta -S Discas/${X}_tf.sam;
samtools view -@ 20 -F 4 Discas/${X}_tf.sam | cut -f 1-6 > Process3/${X}_tf_map.sam;
done;

find ./Process2 -name "*_comp_tf.fasta" | sed 's/_comp_tf.fasta//' | parallel --max-procs=16 \
 '/usr/bin/Rscript  piclash_tf_sep3.R Process3/{/.}_tf_map.sam Process2/{/.}_pimap.tsv Process3/{/.}_chimera;
  RNAplex --temp=26 < Process3/{/.}_chimera_rnaplex.fasta > vienna/{/.}.rnaplex 2> /dev/null;
  /usr/bin/Rscript  piclash_make_vienna2.R Process3/{/.}_chimera.tsv vienna/{/.}.rnaplex vienna/{/.};';


/usr/bin/Rscript  pivienna_analyse.R;

