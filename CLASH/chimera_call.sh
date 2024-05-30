

###################################################################################################################
# CLASH
###################################################################################################################
cd /media/pericles/CLASH/;
mkdir Discas/;
mkdir Process/;
mkdir endtoend/;
mkdir stats/;
mkdir shuffle/;
mkdir test/;

for Y in $(ls fastq/. | grep _R1.fastq.gz$) ; do
X=`echo ${Y} | sed 's/_R1.fastq.gz//'`;
cutadapt -j 10 -m 40 -q 30 -o Process/${X}_2.fastq.gz fastq/${X}_R1.fastq.gz;
cutadapt -j 10 -m 30 -a TGGAATTCTCGGGTGCCAAGGC --discard-untrimmed -o Process/${X}_3.fastq.gz Process/${X}_2.fastq.gz;
cutadapt -j 10 --max-n 0 -o Process/${X}_4.fastq.gz Process/${X}_3.fastq.gz;
zcat Process/${X}_4.fastq.gz | seqkit rmdup -s -o Process/${X}_5.fastq.gz;
cutadapt -j 10 -m 22 -u 4 -u -4 -o Process/${X}_6.fastq.gz Process/${X}_5.fastq.gz;
zcat Process/${X}_6.fastq.gz | awk -f solexa2fasta.awk | \
   awk -f fasta2tab.awk > Process/${X}_5.tab;
perl make_comp_fasta.pl Process/${X}_5.tab > Process/${X}_comp.fasta;
seqkit fx2tab Process/${X}_comp.fasta | cut -f 1-2 > Process/${X}_del.tsv;
hisat2 -p 20 --known-splicesite-infile /media/pericles/CLASH/database/anno/dm6_gene_splicesite.txt --rna-strandness F -k 2000 \
 -f -x /media/pericles/CLASH/database/hisat2Index/dm6_gene -U Process/${X}_comp.fasta -S endtoend/${X}.sam;
samtools view -@ 20 -F 4 -F 16 endtoend/${X}.sam | cut -f 1-6 > endtoend/${X}_map2.sam;
done;
find ./endtoend -name "*_map2.sam" | sed 's/_map2.sam//' | parallel --max-procs=16 \
 '/usr/bin/Rscript  filter_endtoend.R endtoend/{/.}_map2.sam Process/{/.}_del.tsv Process/{/.}_comp2.fasta';



cd /media/pericles/CLASH/;
mkdir Bowtie2Index/;
mkdir BowtieIndex/;
mkdir Process2/;
seqkit fx2tab piRNA/anno/piRNA_all2.fasta | cut -f 1-2 > piRNA/anno/piRNA_all2.tsv;
find ./Process -name "*_comp2.fasta" | sed 's/_comp2.fasta//' | parallel --max-procs=16 \
 'bowtie-build --quiet -f Process/{/.}_comp2.fasta BowtieIndex/{/.}_CLASH_ref';
for Xdef in $(ls Process/. | grep _comp.fasta$) ; do
X=`echo ${Xdef} | sed 's/_comp.fasta//'`;
bowtie -p 20 -f -a -v 0 -S BowtieIndex/${X}_CLASH_ref piRNA/anno/piRNA_all2.fasta Discas/${X}.sam;
samtools view -@ 20 -F 4 -F 16 Discas/${X}.sam | cut -f 1-6 > Process2/${X}_map.sam;
done;
find ./Process -name "*_comp2.fasta" | sed 's/_comp2.fasta//' | parallel --max-procs=16 \
  'seqkit fx2tab Process/{/.}_comp2.fasta | cut -f 1-2 > Process2/{/.}_comp.tsv; \
  /usr/bin/Rscript  /media/pericles/CLASH/clash_pra200308.R Process2/{/.}_map.sam Process2/{/.}_comp.tsv piRNA/anno/piRNA_all2.tsv Process2/{/.}_comp_tf.fasta';


###################
cd /media/pericles/CLASH/;
mkdir Process3/;
mkdir vienna/;
mkdir vienna1fig/;
mkdir typeanno/;
mkdir vienna2fig/;
mkdir vienna4fig/;
for Xdef in $(ls Process2/. | grep _comp_tf.fasta$) ; do
X=`echo ${Xdef} | sed 's/_comp_tf.fasta//'`;
bowtie2 -p 20 -f -a -x database/Bowtie2Index/dm6_gene -U Process2/${X}_comp_tf.fasta -S Discas/${X}_tf.sam;
samtools view -@ 20 -F 4 -F 16 Discas/${X}_tf.sam | cut -f 1-6 > Process3/${X}_tf_map.sam;
done;

find ./Process2 -name "*_comp_tf.fasta" | sed 's/_comp_tf.fasta//' | parallel --max-procs=16 \
  'seqkit fx2tab Process2/{/.}_comp_tf.fasta | cut -f 1-2 > Process3/{/.}_comp_tf.tsv;
  /usr/bin/Rscript --silent --slave --vanilla clash_tf_sep3.R Process3/{/.}_tf_map.sam Process3/{/.}_comp_tf.tsv piRNA/anno/piRNA_all2.tsv Process3/{/.}_chimera;
  bedtools intersect -wa -wb -a Process3/{/.}_chimera_target.bed -b database/anno/dm6_gene_ano.bed > typeanno/{/.}_overlap.tsv;
  bedtools intersect -wa -wb -a Process3/{/.}_chimera_target.bed \
  -b /media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_ins_mapped2.bed > typeanno/{/.}_overlap2.tsv;
  RNAplex --temp=26 < Process3/{/.}_chimera_rnaplex.fasta > vienna/{/.}.rnaplex 2> /dev/null;
  /usr/bin/Rscript --silent --slave --vanilla clash_make_vienna2.R Process3/{/.}_chimera.tsv vienna/{/.}.rnaplex vienna/{/.};';

for Xdef in $(ls Process2/. | grep _comp_tf.fasta$) ; do
X=`echo ${Xdef} | sed 's/_comp_tf.fasta//'`;
/usr/bin/Rscript --silent --slave --vanilla rnaseT1_extract.R vienna/${X}
done;
rm vienna/*_vienna.tsv;


#######################
#splice
cd /media/pericles/CLASH/;
mkdir splice/;

find ./Process2 -name "*_comp_tf.fasta" | sed 's/_comp_tf.fasta//' | parallel --max-procs=16 \
'/usr/bin/Rscript  filter_forsplice.R Process3/{/.}_chimera.tsv Process3/{/.}_comp_tf.tsv splice/{/.}_comp_tf_splice.tsv;
seqkit tab2fx splice/{/.}_comp_tf_splice.tsv > splice/{/.}_comp_tf_splice.fasta;';
for Xdef in $(ls Process2/. | grep _comp_tf.fasta$) ; do
X=`echo ${Xdef} | sed 's/_comp_tf.fasta//'`;
hisat2 -p 2 --known-splicesite-infile /media/pericles/CLASH/database/anno/dm6_gene_splicesite.txt --rna-strandness F -k 2000 \
 -f -x /media/pericles/CLASH/database/hisat2Index/dm6_gene -U splice/${X}_comp_tf_splice.fasta -S Discas/${X}_tf_splice.sam;
samtools view -@ 20 -F 4 -F 16 Discas/${X}_tf_splice.sam | cut -f 1-6 > splice/${X}_tf_splice_map.sam;
done;
find ./Process2 -name "*_comp_tf.fasta" | sed 's/_comp_tf.fasta//' | parallel --max-procs=16 \
'/usr/bin/Rscript  clash_tf_sep3_splice.R splice/{/.}_tf_splice_map.sam splice/{/.}_comp_tf_splice.tsv piRNA/anno/piRNA_all2.tsv splice/{/.}_chimera;
RNAplex --temp=26 < splice/{/.}_chimera_rnaplex.fasta > splice/{/.}.rnaplex 2> /dev/null;
/usr/bin/Rscript --silent --slave --vanilla clash_make_vienna2.R splice/{/.}_chimera.tsv splice/{/.}.rnaplex splice/{/.};';

for Xdef in $(ls Process2/. | grep _comp_tf.fasta$) ; do
X=`echo ${Xdef} | sed 's/_comp_tf.fasta//'`;
/usr/bin/Rscript --silent --slave --vanilla rnaseT1_extract.R splice/${X}
done;
rm splice/*_vienna.tsv;

###################################################################################################################
# vienna-analysis & shuffle
###################################################################################################################
# vienna-analysis
cd /media/pericles/CLASH/;
/usr/bin/Rscript  vienna_table_process.R;
# shuffle
cd /media/pericles/CLASH/;
find ./shuffle -name "*_rnaplex.fasta" | sed 's/_rnaplex.fasta//' | parallel --max-procs=16 \
 'RNAplex --temp=26 < shuffle/{/.}_rnaplex.fasta > shuffle/{/.}.rnaplex 2> /dev/null;'
for X in Vienna_sle ; do
/usr/bin/Rscript --silent --slave --vanilla clash_make_vienna_shuf2.R shuffle/${X}_chimera.tsv shuffle/ shuffle/${X}_vienna_shuffle.tsv;
done;

###################################################################################################################
# stats 
###################################################################################################################
cd /media/pericles/CLASH/;
mkdir fastqc/;

# fastqc --quiet -t 20 Process/*.fastq -o fastqc/;
# fastqc --quiet -t 20 fastq/*.fastq.gz -o fastqc/;
cd /media/pericles/CLASH/;
mkdir fasta_tsv/;
seqkit stats -j 20 */*.fast{a,q.gz} -T > fasta_tsv/comp_stats.tsv;


