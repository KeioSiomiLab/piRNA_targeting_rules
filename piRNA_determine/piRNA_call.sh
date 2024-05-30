
###################################################################################################################
# database (piRNA)
###################################################################################################################
cd /media/pericles/CLASH/piRNA/;
mkdir Discas/;
mkdir trim/;
mkdir figures/;
mkdir otherdata/;
for X in SRR2749802 SRR2749801 ; do
cutadapt -j 10 -q 30 -m 50 -o Discas/${X}_c.fastq.gz fastq/${X}.fastq.gz;
cutadapt -j 10 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -m 26  -o Discas/${X}_c2.fastq.gz Discas/${X}_c.fastq.gz;
seqkit rmdup -s Discas/${X}_c2.fastq.gz -o Discas/${X}_c3.fastq.gz;
cutadapt -u 4 -u -4 -j 10 --max-n 0 -m 23 -M 30 -o trim/${X}_trim.fastq.gz  Discas/${X}_c3.fastq.gz;
done;
for X in SRR9158321 ; do
cutadapt -j 10 -q 30 -m 50 -o Discas/${X}_c.fastq.gz fastq/${X}.fastq.gz;
cutadapt -j 10 -a TGGAATTCTCGGGTGCCAAGG -m 26  -o Discas/${X}_c2.fastq.gz Discas/${X}_c.fastq.gz;
seqkit rmdup -s Discas/${X}_c2.fastq.gz -o Discas/${X}_c3.fastq.gz;
cutadapt -u 4 -u -4 -j 10 --max-n 0 -m 23 -M 30 -o trim/${X}_trim.fastq.gz  Discas/${X}_c3.fastq.gz;
done;
for X in piRNA1 piRNA2 piRNA3 ; do
cutadapt -j 10 -q 30 -m 50 -o Discas/${X}_c.fastq.gz fastq/${X}.fastq.gz;
cutadapt -j 10 -a TGGAATTCTCGGGTGCCAAGG -m 26  -o Discas/${X}_c2.fastq.gz Discas/${X}_c.fastq.gz;
seqkit rmdup -s Discas/${X}_c2.fastq.gz -o Discas/${X}_c3.fastq.gz;
cutadapt -u 4 -u -4 -j 10 --max-n 0 -m 23 -M 30 -o trim/${X}_trim.fastq.gz  Discas/${X}_c3.fastq.gz;
done;
find ./trim -name "*_trim.fastq.gz" | sed 's/_trim.fastq.gz//' | \
 parallel --max-procs=16 'seqkit fq2fa trim/{/.}_trim.fastq.gz -o {/.}.fa';

/usr/bin/Rscript  /media/pericles/CLASH/make_piRNA_db.R;
