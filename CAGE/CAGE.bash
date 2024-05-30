
# CAGE
cd /media/hermione/CAGE/;
for X in chr2L chr2R chr3L chr3R chrX chr4 chrY;do
cp -u /media/pericles/TEfind/database/allchromosome/${X}.fasta \
~/R_libs/BSgenome/extdata/${X}.fasta;
pigz -f -p 20 ~/R_libs/BSgenome/extdata/${X}.fasta;
done;
/usr/bin/Rscript  /media/hermione/CAGE/CAGEBSgenome_build.R;

cd /media/hermione/CAGE/;
mkdir Discas/;
mkdir fastqc/;
mkdir trim/;
mkdir figures/;
fastqc --quiet -t 24 *.fastq.gz -o fastqc/;

cd /media/hermione/CAGE/;
for Y in $(ls . | grep F.fastq.gz$) ; do
X=`echo ${Y} | sed 's/F.fastq.gz//'`;
cutadapt -j 10 --pair-filter=any --max-n 0 -m 100 -a GAAGAGCACACGTCTGAAC -A GAAGAGCGTCGTGTAGGGA \
-o trim/${X}_1_trim.fastq.gz -p trim/${X}_2_trim.fastq.gz ${X}F.fastq.gz ${X}R.fastq.gz;
done;


cd /media/hermione/CAGE/;
mkdir flash/;
mkdir flashlog/;
source activate macs;
for Y in $(ls . | grep F.fastq.gz$) ; do
X=`echo ${Y} | sed 's/F.fastq.gz//'`;
flash2 -t 20 -d /media/hermione/CAGE/flash/ -o ${X} \
${X}F.fastq.gz ${X}R.fastq.gz 2>&1 | tee flashlog/${X}_flash.log;
pigz -f -p 20 flash/*.fastq;
done;
conda deactivate;
rm flash/*;




cd /media/hermione/CAGE/;
mkdir genome;
mkdir TE;
ulimit -n 10000;
for Y in $(ls . | grep F.fastq.gz$) ; do
X=`echo ${Y} | sed 's/F.fastq.gz//'`;
STAR runMode alignReads --runThreadN 24 --genomeDir /media/pericles/CLASH/database/star_genome_dm6/ \
--readFilesIn trim/${X}_1_trim.fastq.gz trim/${X}_2_trim.fastq.gz \
--outReadsUnmapped None  \
--outSAMattributes All --outSAMtype BAM SortedByCoordinate \
 --readFilesCommand zcat --outSAMunmapped None --outFilterMultimapNmax 20 --alignIntronMax 1000000 \
--outFileNamePrefix /media/hermione/CAGE/genome/${X}_STAR;
samtools index /media/hermione/CAGE/genome/${X}_STARAligned.sortedByCoord.out.bam;
samtools view -@ 20 -f 0x2 -bS -uq 20 /media/hermione/CAGE/genome/${X}_STARAligned.sortedByCoord.out.bam \
| samtools sort -@ 20 -m 4G -n - \
| bedtools bamtobed -mate1 -bedpe -i stdin > /media/hermione/CAGE/genome/${X}.bed;
source activate macs;
for bin in 1 10;do
for A in forward reverse;do
bamCoverage -p 20 -of bigwig -bs ${bin} --filterRNAstrand ${A} --normalizeUsing CPM \
-b genome/${X}_STARAligned.sortedByCoord.out.bam \
-o genome/${X}_${bin}_${A}.bigwig;
done;
done;
samtools flagstat /media/hermione/CAGE/genome/${X}_STARAligned.sortedByCoord.out.bam;
STAR runMode alignReads --runThreadN 24 --genomeDir /media/pericles/CLASH/database/star_TE_dm6/ \
--readFilesIn trim/${X}_1_trim.fastq.gz trim/${X}_2_trim.fastq.gz \
--outReadsUnmapped None --limitBAMsortRAM 31000000000 \
--outSAMattributes All --outSAMtype BAM SortedByCoordinate \
 --readFilesCommand zcat --outSAMunmapped None --outFilterMultimapNmax 20 --alignIntronMax 1000000 \
--outFileNamePrefix /media/hermione/CAGE/TE/${X}_TE;
samtools index /media/hermione/CAGE/TE/${X}_TEAligned.sortedByCoord.out.bam;
samtools view -@ 20 -f 0x2 -bS -uq 20 /media/hermione/CAGE/TE/${X}_TEAligned.sortedByCoord.out.bam \
| samtools sort -@ 20 -m 4G -n - \
| bedtools bamtobed -mate1 -bedpe -i stdin > /media/hermione/CAGE/TE/${X}_TE.bed;
done;




rm Discas/*;
# https://github.com/alexdobin/STAR/issues/528
cd /media/hermione/CAGE/;
mkdir fastqc2;
fastqc --quiet -t 24 Discas/*_1_trim.fastq.gz Discas/*_2_trim.fastq.gz -o fastqc2/;


cd /media/hermione/CAGE/;
mkdir stat/;
for Y in $(ls . | grep F.fastq.gz$) ; do
X=`echo ${Y} | sed 's/F.fastq.gz//'`;
samtools stats -@ 20 -i 8000 /media/hermione/CAGE/genome/${X}_STARAligned.sortedByCoord.out.bam > stat/${X}_stat.txt;
plot-bamstats -l -p stat/${X}/${X} stat/${X}_stat.txt
grep ^IS stat/${X}_stat.txt | cut -f 2- > stat/${X}_len.tsv;
done;



mkdir ctss/;
/usr/bin/Rscript  /media/hermione/CAGE/CAGE_process.R;

/usr/bin/Rscript  /media/hermione/CAGE/CAGE_normalize.R;


for X in figures ;do
mkdir /media/hermione/CAGE/${X}/ctss/;
cat /media/hermione/CAGE/${X}/siEGFP.tagClusters.bed \
/media/hermione/CAGE/${X}/siPiwi.tagClusters.bed > /media/hermione/CAGE/${X}/premerge.bed;
bedtools sort -i /media/hermione/CAGE/${X}/premerge.bed \
> /media/hermione/CAGE/${X}/premerge_sort.bed;
bedtools merge -c 6 -o distinct -s -i /media/hermione/CAGE/${X}/premerge_sort.bed \
> /media/hermione/CAGE/${X}/mergeTSS.bed;
bedtools intersect -v -a /media/hermione/CAGE/${X}/mergeTSS.bed \
-b /media/pericles/TEfind/rmblast/dm6_rmblast_mask.bed > \
/media/hermione/CAGE/${X}/CAGE_remove.bed;
for A in siEGFP_rep1 siEGFP_rep2 siEGFP_rep3 siPIWI_rep1 siPIWI_rep2 siPIWI_rep3 ;do
bedtools intersect -wa -wb -a /media/hermione/CAGE/${X}/CAGE_remove.bed \
-b /media/hermione/CAGE/ctss/cage${A}.bed > /media/hermione/CAGE/${X}/ctss/cage${A}.cov;
done;
awk '{print $1 "\t" $2 "\t" $3 "\t" "cage" "\t" "0" "\t" $4}' \
/media/hermione/CAGE/${X}/CAGE_remove.bed > /media/hermione/CAGE/${X}/CAGE_remove2.bed;
bedtools closest -s -d -k 1 -a /media/hermione/CAGE/${X}/CAGE_remove2.bed \
-b /media/pericles/CLASH/database/gtf/dm6_tss.bed > /media/hermione/CAGE/${X}/CAGE_annotate.bed;
bedtools intersect -s -wa -wb -a /media/hermione/CAGE/${X}/CAGE_remove2.bed \
-b /media/pericles/CLASH/database/gtf/dm6_5UTR.bed > \
/media/hermione/CAGE/${X}/CAGE_inter.bed;
done;
bedtools sort -i /media/hermione/CAGE/figures/All.samples.tagClusters.qLow0.1_qUp0.9.bed | \
bedtools merge -c 6 -o collapse -s -i - > /media/hermione/CAGE/figures/all_peaks_merged.bed;

/usr/bin/Rscript  /media/hermione/CAGE/CAGE_annotate.R;

# fusion analysis of CAGE
cd /media/hermione/CAGE/;
mkdir fusion/;
ulimit -n 10000;
STAR --runThreadN 16 --runMode genomeGenerate --genomeSAindexNbases 12 --genomeDir /media/ariel/HN00139860/fusion/normal_star_in/ \
--genomeFastaFiles /media/pericles/TEfind/database/fusion_genome_TE.fasta \
--sjdbGTFfile /media/ariel/HN00139860/fusion/siPIWI_isoforms_fusion.gtf;
for Y in $(ls . | grep F.fastq.gz$) ; do
X=`echo ${Y} | sed 's/F.fastq.gz//'`;
STAR runMode alignReads --runThreadN 24 --genomeDir /media/ariel/HN00139860/fusion/normal_star_in/ \
--readFilesIn trim/${X}_1_trim.fastq.gz trim/${X}_2_trim.fastq.gz --chimNonchimScoreDropMin 10 --outSAMstrandField intronMotif \
--outReadsUnmapped None --twopassMode Basic --chimSegmentMin 12 --chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10 \
--chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4 --chimMultimapNmax 20 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 \
--alignInsertionFlush Right --alignSplicedMateMapLminOverLmate 0 --alignSplicedMateMapLmin 30 \
--outSAMattributes All --outSAMtype BAM SortedByCoordinate --alignMatesGapMax 200000 --chimSegmentReadGapMax 3 \
--readFilesCommand zcat --outSAMunmapped None --outFilterMultimapNmax 20 --alignIntronMax 200000 --alignSJstitchMismatchNmax 5 -1 5 5 \
--outFileNamePrefix /media/hermione/CAGE/fusion/${X}_STAR;
samtools index /media/hermione/CAGE/fusion/${X}_STARAligned.sortedByCoord.out.bam;
done;

cd /media/hermione/CAGE/;
mkdir overlap_res/;
mkdir chimera_res/;
mkdir overlap/;
mkdir fusion_figure/;
find ./fusion -name "*_STARChimeric.out.junction" | sed 's/_STARChimeric.out.junction//' | parallel --max-procs=6 \
'/usr/bin/Rscript chimera_process1.R fusion/{/.}_STARChimeric.out.junction fusion/{/.}_STARSJ.out.tab overlap/{/.};
bedtools intersect -s -wa -wb -a /media/pericles/CLASH/database/anno/gene_sort.bed \
-b overlap/{/.}_junction.bed > overlap_res/{/.}_gene_overlap.bed;
/usr/bin/Rscript chimera_process2.R overlap_res/{/.}_gene_overlap.bed \
overlap/{/.}_junction.result chimera_res/{/.}_chimera.tsv;';

seqkit stats -j 20 trim/*_trim.fastq.gz -T > comp_stats.tsv;


