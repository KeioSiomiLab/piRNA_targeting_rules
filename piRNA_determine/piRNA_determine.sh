
###################################################################################################################
# database (hybdata)
###################################################################################################################
cd /media/pericles/CLASH/database/;
mkdir reffasta/;
mkdir anno/;
find ./refseq -name "*-r6.36.fasta.gz" | sed 's/-r6.36.fasta.gz//' | sed 's/dmel-all-//' | \
 parallel --max-procs=16 'gunzip -c refseq/dmel-all-{/.}-r6.36.fasta.gz > reffasta/{/.}.fasta';
cp -u /media/pericles/TEfind/database/TE_dm.fasta /media/pericles/CLASH/database/reffasta/ensembl.fasta;
cp -u /media/pericles/CLASH/piRNA/piRNA_all.fasta /media/pericles/CLASH/database/reffasta/piRNA_all.fasta;
cp -u /media/pericles/CLASH/piRNA/piRNA_all_tes.fasta /media/pericles/CLASH/database/reffasta/piRNA_all_tes.fasta;
cp -u /media/pericles/TEfind/database/TE_dm.fasta /media/pericles/CLASH/database/ensembl.fasta;
cp -u /media/pericles/CLASH/piRNA/piRNA_all.fasta /media/pericles/CLASH/database/piRNA_all.fasta;
cp -u /media/pericles/CLASH/piRNA/piRNA_all_tes.fasta /media/pericles/CLASH/database/piRNA_all_tes.fasta;

find ./reffasta -name "*.fasta" | sed 's/.fasta//' | \
 parallel --max-procs=16 'seqkit fx2tab reffasta/{/.}.fasta | cut -f 1-2 > reffasta/{/.}.tsv';
/usr/bin/Rscript  /media/pericles/CLASH/construct_gene_database.R;

cd /media/pericles/CLASH/database/;
mkdir Bowtie2Index/;
bowtie2-build --quiet -f dm6_gene.fasta Bowtie2Index/dm6_gene;
bowtie2-build --quiet -f reffasta/gene.fasta Bowtie2Index/genes;
bowtie2-build --quiet -f reffasta/ensembl.fasta Bowtie2Index/TEref;
bowtie2-build --quiet -f reffasta/miRNA.fasta Bowtie2Index/miRNA;
bowtie-build --quiet -f reffasta/ensembl.fasta BowtieIndex/TEref;
bowtie2-build --quiet -f dm6_rmstRNA.fasta Bowtie2Index/dm6_rmstRNA;
samtools faidx dm6_gene.fasta;
seqkit fx2tab dm6_gene.fasta | cut -f 1-2 > dm6_gene.tsv;

cp -u /media/pericles/TEfind/database/dm6_chr.fasta /media/pericles/CLASH/database/dm6_chr.fasta;
bowtie2-build --quiet -f dm6_chr.fasta Bowtie2Index/dm6_genome;

mkdir BowtieIndex/;
bowtie-build --quiet -f dm6_gene.fasta BowtieIndex/dm6_gene;
bowtie-build --quiet -f dm6_gene2.fasta BowtieIndex/dm6_gene2;
bowtie2-build --quiet -f dm6_gene2.fasta Bowtie2Index/dm6_gene2;
mkdir hisat2Index/;
hisat2-build --quiet -f /media/pericles/CLASH/database/ensembl.fasta hisat2Index/TE;
hisat2-build --quiet -f /media/pericles/CLASH/database/dm6_gene.fasta hisat2Index/dm6_gene;


cd /media/pericles/CLASH/database/anno/;
mkdir /media/pericles/TEfind/vcfprocess/telomere/;
bedtools intersect -wa -wb -a /media/pericles/CLASH/database/anno/gene.bed  \
 -b /media/pericles/TEfind/rmblast/dm6_rmblast_mask.bed > /media/pericles/CLASH/database/anno/repeat_gene.txt;
bedtools intersect -wa -wb -a /media/pericles/CLASH/database/anno/gene.bed \
 -b /media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_ins.bed > /media/pericles/CLASH/database/anno/insertion_gene.txt;
bedtools intersect -wa -wb -a /media/pericles/TEfind/vcfprocess/pbsv_vcf_ins.bed \
 -b /media/pericles/TEfind/database/centromere_telomere.txt > \
 /media/pericles/TEfind/vcfprocess/telomere/pbsv_vcf_ins_telomere.txt;
bedtools sort -i /media/pericles/CLASH/database/anno/gene.bed  \
> /media/pericles/CLASH/database/anno/gene_sort.bed;
# closest gene of insertion!
bedtools closest -io -k 2 -a /media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_ins_ref.bed \
-b /media/pericles/CLASH/database/anno/gene_sort.bed \
> /media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_ins_closest.bed;
bedtools intersect -wa -wb -a /media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_ins_ref.bed \
-b /media/pericles/CLASH/database/anno/gene_sort.bed \
> /media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_ins_inside.bed;
bedtools intersect -wa -wb -a /media/pericles/TEfind/newTE/OSCrepeat.bed \
-b /media/pericles/CLASH/database/anno/gene_sort.bed \
> /media/pericles/TEfind/newTE/OSC_inside.bed;

cd /media/pericles/CLASH/database/anno/;
bedtools intersect -wa -wb -a gene.bed -b /media/pericles/TEfind/vcfprocess/embl_pbsv_del_repeat.bed > repeat_gene2.txt;
gunzip -c /media/pericles/CLASH/database/gtf/dmel-all-r6.36.gtf.gz > /media/pericles/CLASH/database/gtf/dm6.gtf;
/usr/bin/Rscript  /media/pericles/CLASH/GTF_process.R;
/usr/bin/Rscript  /media/pericles/CLASH/get_gene_database_annotation.R;
gffread /media/pericles/CLASH/database/gtf/dm6_igv.gtf \
 -g /media/pericles/TEfind/database/dm6_chr.fasta -E -o /media/pericles/CLASH/database/gtf/dm6_igv.gff3;

python ~/anaconda3/bin/hisat2_extract_splice_sites.py /media/pericles/CLASH/database/gtf/dm6.gtf > \
  /media/pericles/CLASH/database/gtf/splicesites.txt;


cd /media/pericles/CLASH/RNAseq/;
mkdir genome;
mkdir rsem/;
mkdir /media/pericles/CLASH/database/star_genome_dm6/;
STAR --runThreadN 16 --runMode genomeGenerate --genomeSAindexNbases 12 \
--genomeDir /media/pericles/CLASH/database/star_genome_dm6/ \
--genomeFastaFiles /media/pericles/CLASH/database/dm6_chr.fasta \
--sjdbGTFfile /media/pericles/CLASH/database/gtf/dm6_use.gtf;
STAR --runThreadN 16 --runMode genomeGenerate --genomeSAindexNbases 8 \
--genomeDir /media/pericles/CLASH/database/star_TE_dm6/ \
--genomeFastaFiles /media/pericles/CLASH/database/reffasta/ensembl.fasta ;
source activate meme;
rsem-prepare-reference --gtf /media/pericles/CLASH/database/gtf/dm6_use.gtf -p 10 \
/media/pericles/CLASH/database/dm6_chr.fasta \
/media/pericles/CLASH/RNAseq/rsem;
conda deactivate;

###################################################################################################################
# genemap (piRNA)
###################################################################################################################

cd /media/pericles/CLASH/piRNA/;
mkdir Discas/;
mkdir gene/;
mkdir anno/;
mkdir figout/;
mkdir gene_reg/;
find ./trim -name "*_trim.fastq.gz" | sed 's/_trim.fastq.gz//' | \
 parallel --max-procs=16 'seqkit fq2fa trim/{/.}_trim.fastq.gz -o {/.}.fa;
 seqkit fx2tab {/.}.fa | cut -f 1-2 > {/.}.tsv;';
for X in SRR2749802 SRR2749801 SRR9158321 piRNA1 piRNA2 piRNA3; do
bowtie2 -f -p 20 -a -x /media/pericles/CLASH/database/Bowtie2Index/dm6_gene \
 -U piRNA_${X}.fa -S Discas/${X}.sam;
samtools view -@ 20 -F 4 -F 16 Discas/${X}.sam | cut -f 1-6 > gene/${X}_map.sam;
done;
find ./trim -name "*_trim.fastq.gz" | sed 's/_trim.fastq.gz//' | \
 parallel --max-procs=16 '/usr/bin/Rscript  /media/pericles/CLASH/piRNA_selection.R \
 gene/{/.}_map.sam gene/{/.}_map2;
 /usr/bin/Rscript  /media/pericles/CLASH/piRNAQC.R \
 gene/{/.}_map2.tsv /media/pericles/TEfind/newTE/allTE.tsv piRNA_{/.}.tsv figures/ {/.};';


cd /media/pericles/CLASH/piRNA/;
find ./trim -name "*_trim.fastq.gz" | sed 's/_trim.fastq.gz//' | \
 parallel --max-procs=16 'bedtools getfasta -fi /media/pericles/CLASH/database/dm6_gene.fasta \
 -bed gene_reg/{/.}_gene.bed -fo gene_reg/{/.}_gene.fasta;
 seqkit fx2tab gene_reg/{/.}_gene.fasta | cut -f 1-2 > gene_reg/{/.}_gene.tsv;
 /usr/bin/Rscript  /media/pericles/CLASH/piRNAQC2.R piRNA_{/.}.tsv gene_reg/{/.}_gene.tsv gene_reg/ {/.};';


for X in piRNA_all ; do
seqkit fx2tab ${X}.fasta | cut -f 1-2 > ${X}.tsv;
bowtie2 -f -p 20 -a -x /media/pericles/CLASH/database/Bowtie2Index/dm6_gene \
 -U ${X}.fasta -S Discas/${X}.sam;
samtools view -@ 20 -F 4 -F 16 Discas/${X}.sam | cut -f 1-6 > gene/${X}_map.sam;
done;




parallel --max-procs=16 \
 '/usr/bin/Rscript  /media/pericles/CLASH/piRNA_selection.R gene/{/.}_map.sam gene/{/.}_map2;
 bedtools intersect -wa -wb -a gene/{/.}_map2.bed -b /media/pericles/CLASH/database/anno/dm6_gene_ano.bed > gene/{/.}_overlap.tsv;
 bedtools intersect -wa -wb -a gene/{/.}_map2.bed -b /media/pericles/CLASH/database/anno/dm6_gene_ano2.bed > gene/{/.}_overlap2.tsv;
 /usr/bin/Rscript  /media/pericles/CLASH/piRNA_extract_target.R gene/{/.}_map2 \
  /media/pericles/TEfind/newTE/allTE.tsv {/.}.tsv anno/{/.} figures/{/.};
  bedtools getfasta -fi /media/pericles/CLASH/database/dm6_gene.fasta \
  -bed gene_reg/{/.}_gene.bed -fo gene_reg/{/.}_gene.fasta;
  seqkit fx2tab gene_reg/{/.}_gene.fasta | cut -f 1-2 > gene_reg/{/.}_gene.tsv;' ::: piRNA_all;
/usr/bin/Rscript  /media/pericles/CLASH/calculate_piRNA_CPM.R;



###################################################################################################################
# piRNA expression in TEs
###################################################################################################################
cd /media/pericles/CLASH/piRNA;
mkdir TE;
mkdir gene;
mkdir TE/mismatch/;
for X in piRNA_all; do
bowtie -a --best --strata --threads 20 -l 20 -f -n 3 -S /media/pericles/CLASH/database/BowtieIndex/TEref \
${X}.fasta Discas/${X}_TE.sam;
samtools view -@ 20 -f 16 -bS Discas/${X}_TE.sam > Discas/${X}_TE.bam;
samtools sort -@ 20 -m 4G Discas/${X}_TE.bam > TE/${X}_TE_sort.bam;
samtools index TE/${X}_TE_sort.bam;
bedtools bamtobed -i TE/${X}_TE_sort.bam > TE/${X}_TE.bed;
bowtie -a --best --strata --threads 20 -l 20 -f -n 3 -S /media/pericles/CLASH/database/BowtieIndex/dm6_gene \
${X}.fasta Discas/${X}_gene.sam;
samtools view -@ 20 -f 16 -bS Discas/${X}_gene.sam > Discas/${X}_gene.bam;
samtools sort -@ 20 -m 4G Discas/${X}_gene.bam > gene/${X}_gene_sort.bam;
samtools index gene/${X}_gene_sort.bam;
samtools idxstats TE/${X}_TE_sort.bam > TE/${X}_count.txt;
source activate macs;
for bin in 1 10;do
bamCoverage --filterRNAstrand forward --binSize ${bin} -p 20 --normalizeUsing None --scaleFactor -1 \
 -b TE/${X}_TE_sort.bam -o TE/${X}_reverse_${bin}.bw ;
bamCoverage --filterRNAstrand reverse --binSize ${bin} -p 20 --normalizeUsing None  \
 -b TE/${X}_TE_sort.bam -o TE/${X}_forward_${bin}.bw  ;
bigWigToBedGraph TE/${X}_reverse_${bin}.bw TE/${X}_reverse_${bin}.bedgraph;
bigWigToBedGraph TE/${X}_forward_${bin}.bw TE/${X}_forward_${bin}.bedgraph;
done;
conda deactivate;
for Y in 0 1 2 3;do
bowtie -a --best --strata --threads 20 -l 20 -f -v ${Y} --best -S /media/pericles/CLASH/database/BowtieIndex/TEref \
${X}.fasta Discas/${X}_${Y}_TE.sam;
samtools view -@ 20 -f 16 -bS Discas/${X}_${Y}_TE.sam > Discas/${X}_${Y}_TE.bam;
samtools sort -@ 20 -m 4G Discas/${X}_${Y}_TE.bam > TE/mismatch/${X}_${Y}_TE_sort.bam;
samtools index TE/mismatch/${X}_${Y}_TE_sort.bam;
bedtools bamtobed -i TE/mismatch/${X}_${Y}_TE_sort.bam > TE/mismatch/${X}_${Y}_TE.bed;
samtools idxstats TE/mismatch/${X}_${Y}_TE_sort.bam > TE/mismatch/${X}_${Y}_count.txt;
source activate macs;
for bin in 1 10;do
bamCoverage --filterRNAstrand forward --binSize ${bin} -p 20 --normalizeUsing None --scaleFactor -1 \
 -b TE/mismatch/${X}_${Y}_TE_sort.bam -o TE/mismatch/${X}_${Y}_reverse_${bin}.bw ;
bamCoverage --filterRNAstrand reverse --binSize ${bin} -p 20 --normalizeUsing None  \
 -b TE/mismatch/${X}_${Y}_TE_sort.bam -o TE/mismatch/${X}_${Y}_forward_${bin}.bw  ;
bigWigToBedGraph TE/mismatch/${X}_${Y}_reverse_${bin}.bw TE/mismatch/${X}_${Y}_reverse_${bin}.bedgraph;
bigWigToBedGraph TE/mismatch/${X}_${Y}_forward_${bin}.bw TE/mismatch/${X}_${Y}_forward_${bin}.bedgraph;
done;
conda deactivate;
done
done;

cd /media/pericles/CLASH/piRNA;
for X in SRR2749802 SRR2749801 SRR9158321 piRNA1 piRNA2 piRNA3; do
bowtie -a --best --strata --threads 20 -l 20 -f -n 3 -S /media/pericles/CLASH/database/BowtieIndex/TEref \
${X}.fa  Discas/${X}_TE.sam;
samtools view -@ 20 -f 16 -bS Discas/${X}_TE.sam > Discas/${X}_TE.bam;
samtools sort -@ 20 -m 4G Discas/${X}_TE.bam > TE/${X}_TE_sort.bam;
samtools index TE/${X}_TE_sort.bam;
samtools idxstats TE/${X}_TE_sort.bam > TE/${X}_count.txt;
source activate macs;
for bin in 1 10;do
bamCoverage --filterRNAstrand forward --binSize ${bin} -p 20 --normalizeUsing None --scaleFactor -1 \
 -b TE/${X}_TE_sort.bam -o TE/${X}_reverse_${bin}.bw ;
bamCoverage --filterRNAstrand reverse --binSize ${bin} -p 20 --normalizeUsing None  \
 -b TE/${X}_TE_sort.bam -o TE/${X}_forward_${bin}.bw  ;
bigWigToBedGraph TE/${X}_reverse_${bin}.bw TE/${X}_reverse_${bin}.bedgraph;
bigWigToBedGraph TE/${X}_forward_${bin}.bw TE/${X}_forward_${bin}.bedgraph;
done;
conda deactivate;
for Y in 0 1 2 3;do
bowtie -a --best --strata --threads 20 -l 20 -f -v ${Y} --best -S /media/pericles/CLASH/database/BowtieIndex/TEref \
${X}.fa Discas/${X}_${Y}_TE.sam;
samtools view -@ 20 -f 16 -bS Discas/${X}_${Y}_TE.sam > Discas/${X}_${Y}_TE.bam;
samtools sort -@ 20 -m 4G Discas/${X}_${Y}_TE.bam > TE/mismatch/${X}_${Y}_TE_sort.bam;
samtools index TE/mismatch/${X}_${Y}_TE_sort.bam;
samtools idxstats TE/mismatch/${X}_${Y}_TE_sort.bam > TE/mismatch/${X}_${Y}_count.txt;
source activate macs;
for bin in 1 10;do
bamCoverage --filterRNAstrand forward --binSize ${bin} -p 20 --normalizeUsing None --scaleFactor -1 \
 -b TE/mismatch/${X}_${Y}_TE_sort.bam -o TE/mismatch/${X}_${Y}_reverse_${bin}.bw ;
bamCoverage --filterRNAstrand reverse --binSize ${bin} -p 20 --normalizeUsing None  \
 -b TE/mismatch/${X}_${Y}_TE_sort.bam -o TE/mismatch/${X}_${Y}_forward_${bin}.bw  ;
bigWigToBedGraph TE/mismatch/${X}_${Y}_reverse_${bin}.bw TE/mismatch/${X}_${Y}_reverse_${bin}.bedgraph;
bigWigToBedGraph TE/mismatch/${X}_${Y}_forward_${bin}.bw TE/mismatch/${X}_${Y}_forward_${bin}.bedgraph;
done;
conda deactivate;
done
done;
seqkit stats -j 20  *.fa -T  > read_info.tsv;



