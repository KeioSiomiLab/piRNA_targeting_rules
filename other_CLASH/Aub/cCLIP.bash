
###################################################################################################################
# cCLIP analysis
###################################################################################################################
cd /media/hermione/cCLIP/ovary/;
mkdir fastqc/;
fastqc --quiet -t 20 ovary/*.fastq.gz -o fastqc/;

cd /media/hermione/cCLIP/;
mkdir Discas/;
mkdir trim/;
mkdir figures/;
mkdir otherdata/;
for X in SRR3051366 SRR3051369  ; do
cutadapt -j 10 --max-n 0 -m 23 -M 30 -o trim/${X}_trim.fastq.gz  ovary/${X}.fastq.gz;
seqkit fq2fa trim/${X}_trim.fastq.gz -o ${X}.fa;
done;
/usr/bin/Rscript  cclip_make_piRNA_db.R;


# TE  mapping
cd /media/hermione/cCLIP/;
mkdir TE;
for X in SRR3051366 SRR3051369 piRNA_all; do
bowtie2 -p 20 -a -f -x /media/pericles/CLASH/database/Bowtie2Index/TEref -U ${X}.f{a,asta} -S Discas/${X}_TE.sam;
samtools view -@ 20 -bS Discas/${X}_TE.sam > Discas/${X}_TE.bam;
samtools sort -@ 20 -m 4G Discas/${X}_TE.bam > TE/${X}_TE_sort.bam;
samtools index TE/${X}_TE_sort.bam;
source activate macs;
for bin in 1 10;do
for A in forward reverse;do
bamCoverage -p 20 -of bigwig -bs ${bin} --filterRNAstrand ${A} --normalizeUsing CPM \
-b TE/${X}_TE_sort.bam \
-o TE/${X}_${bin}_${A}.bigwig;
done;
done;
conda deactivate;
done;



cd /media/hermione/cCLIP/;
mkdir Discas/;
mkdir Process/;
mkdir endtoend/;
mkdir stats/;
mkdir shuffle/;
mkdir test/;

for X in SRR3051373 SRR3051374 SRR3051375  ; do
cutadapt -j 10 -m 40 -q 30 -o Process/${X}_2.fastq.gz ovary/${X}.fastq.gz;
cutadapt -j 10 -m 22 --max-n 0 -a GTGTCAGTCACTTCCAGCGGTCGTATGCCGTCTTCTGCTTG --discard-untrimmed -o Process/${X}_3.fastq.gz Process/${X}_2.fastq.gz;
zcat Process/${X}_3.fastq.gz | awk -f cclip_solexa2fasta.awk | \
   awk -f cclip_fasta2tab.awk > Process/${X}_5.tab;
perl cclip_make_comp_fasta.pl Process/${X}_5.tab > Process/${X}_comp.fasta;
seqkit fx2tab Process/${X}_comp.fasta | cut -f 1-2 > Process/${X}_del.tsv;
bowtie2 -p 20 --very-sensitive -f -U Process/${X}_comp.fasta --reorder \
 -x /media/pericles/CLASH/database/Bowtie2Index/dm6_gene2 -S endtoend/${X}.sam;
samtools view -@ 20 -F 4 -F 16 endtoend/${X}.sam | cut -f 1-6 > endtoend/${X}_map2.sam;
done;


find ./endtoend -name "*_map2.sam" | sed 's/_map2.sam//' | parallel --max-procs=16 \
 '/usr/bin/Rscript  cclip_filter_endtoend.R endtoend/{/.}_map2.sam Process/{/.}_del.tsv Process/{/.}_comp2.fasta';


cd /media/hermione/cCLIP/;
mkdir Bowtie2Index/;
mkdir BowtieIndex/;
mkdir Process2/;
seqkit fx2tab piRNA_all.fasta | cut -f 1-2 > piRNA_all.tsv;
find ./Process -name "*_comp2.fasta" | sed 's/_comp2.fasta//' | parallel --max-procs=16 \
 'bowtie-build --quiet -f Process/{/.}_comp2.fasta BowtieIndex/{/.}_CLASH_ref';
for Xdef in $(ls Process/. | grep _comp.fasta$) ; do
X=`echo ${Xdef} | sed 's/_comp.fasta//'`;
bowtie -p 20 -f -a -v 0 -S -x BowtieIndex/${X}_CLASH_ref piRNA_all.fasta Discas/${X}.sam;
samtools view -@ 20 -F 4 -F 16 Discas/${X}.sam | cut -f 1-6 > Process2/${X}_map.sam;
done;
find ./Process -name "*_comp2.fasta" | sed 's/_comp2.fasta//' | parallel --max-procs=16 \
  'seqkit fx2tab Process/{/.}_comp2.fasta | cut -f 1-2 > Process2/{/.}_comp.tsv; \
  /usr/bin/Rscript  cclip_clash_pra200308.R Process2/{/.}_map.sam Process2/{/.}_comp.tsv piRNA_all.tsv Process2/{/.}_comp_tf.fasta';


cd /media/hermione/cCLIP/;
mkdir Process3/;
mkdir vienna/;
mkdir vienna1fig/;
mkdir typeanno/;
mkdir vienna2fig/;
mkdir vienna4fig/;
for Xdef in $(ls Process2/. | grep _comp_tf.fasta$) ; do
X=`echo ${Xdef} | sed 's/_comp_tf.fasta//'`;
bowtie2 -p 20 -f -a -x /media/pericles/CLASH/database/Bowtie2Index/dm6_gene2 -U Process2/${X}_comp_tf.fasta -S Discas/${X}_tf.sam;
samtools view -@ 20 -F 4 -F 16 Discas/${X}_tf.sam | cut -f 1-6 > Process3/${X}_tf_map.sam;
done;


find ./Process2 -name "*_comp_tf.fasta" | sed 's/_comp_tf.fasta//' | parallel --max-procs=16 \
  'seqkit fx2tab Process2/{/.}_comp_tf.fasta | cut -f 1-2 > Process3/{/.}_comp_tf.tsv;
  /usr/bin/Rscript  cclip_clash_tf_sep3.R Process3/{/.}_tf_map.sam Process3/{/.}_comp_tf.tsv piRNA_all.tsv Process3/{/.}_chimera;
  bedtools intersect -wa -wb -a Process3/{/.}_chimera_target.bed -b /media/pericles/CLASH/database/anno/dm6_gene_ano.bed > typeanno/{/.}_overlap.tsv;
  RNAplex --temp=26 < Process3/{/.}_chimera_rnaplex.fasta > vienna/{/.}.rnaplex 2> /dev/null;
  /usr/bin/Rscript --silent --slave --vanilla cclip_clash_make_vienna2.R Process3/{/.}_chimera.tsv vienna/{/.}.rnaplex vienna/{/.};';

cd /media/hermione/cCLIP/;
/usr/bin/Rscript  cclip_vienna_table_process.R;


