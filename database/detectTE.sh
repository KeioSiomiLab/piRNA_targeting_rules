###################################################################################################################
# tebreak
###################################################################################################################
cd /media/pericles/TEfind/;
mkdir IlluminaProcess/;
mkdir bwaIndex/;
mkdir tebreakres/;
bwa index database/dm6_chr.fasta;
samtools faidx database/dm6_chr.fasta;

for X in Input-siE-1 Input-siE-2 ; do
cutadapt -m 50 -j 10 -a AGATCGGAAGAG -A AGATCGGAAGAG \
 -o IlluminaProcess/${X}_1_trim.fastq.gz -p IlluminaProcess/${X}_2_trim.fastq.gz \
  illuminaDNA/${X}_R1.fastq.gz illuminaDNA/${X}_R2.fastq.gz;
bwa mem -M -Y -t 20 bwaIndex/dm6_chr IlluminaProcess/${X}_1_trim.fastq.gz IlluminaProcess/${X}_2_trim.fastq.gz | \
samtools view -b - > IlluminaProcess/${X}.bam;
samtools sort -@ 20 -m 4G IlluminaProcess/${X}.bam > IlluminaProcess/${X}_sort.bam;
samtools index IlluminaProcess/${X}_sort.bam;
picard MarkDuplicates -I IlluminaProcess/${X}_sort.bam -O IlluminaProcess/${X}_rmdup.bam -M metrics.out;
samtools index IlluminaProcess/${X}_rmdup.bam;
done;
rm IlluminaProcess/*;

tebreak -b IlluminaProcess/Input-siE-1_rmdup.bam,IlluminaProcess/Input-siE-2_rmdup.bam \
 -r database/dm6_chr.fasta \
 -p 20 -d database/annotateTE/TEannotate_tebreak.txt \
 -m database/centromere_telomere.txt \
 --max_ins_reads 500 \
 -i database/TE_dm.fasta \
 -o tebreakres/test.tsv;
/usr/bin/Rscript  TEbreak_res_process.R;

###################################################################################################################
# svcall, Sniffle and pbsv
###################################################################################################################
cd /media/pericles/TEfind/;
mkdir svcall/;
mkdir vcfprocess/;
seqkit fq2fa pacbioDNA/OSC_genome_HiFi.fastq.gz -o pacbioDNA/OSC_genome_HiFi.fasta;

mkdir qual/;
mkdir qual2/;
source activate canu;
NanoStat --fastq pacbioDNA/OSC_genome_HiFi.fastq.gz -t 16 > pacbioDNA/nanostat.txt;
NanoPlot --fastq pacbioDNA/OSC_genome_HiFi.fastq.gz --loglength -t 16 -o qual;
conda deactivate;

longdir=svcall
ref=database/dm6_chr.fasta
for X in pacbio ; do
source activate canu;
pbmm2 align -j 16 ${ref} pacbioDNA/OSC_genome_HiFi.bam ${longdir}/${X}_pbsv.bam --sort --median-filter --sample sample1;
conda deactivate;
samtools index ${longdir}/${X}_pbsv.bam;
source activate canu;
pbindex ${longdir}/${X}_pbsv.bam;
pbsv discover ${longdir}/${X}_pbsv.bam ${longdir}/${X}.svsig.gz;
pbsv call --ccs -x 40 --max-ins-length 15K ${ref} ${longdir}/${X}.svsig.gz ${longdir}/${X}_pbsv.vcf;
conda deactivate;
done;
/usr/bin/Rscript  vcf_process.R;
for X in pbsv_vcf_del_homo ; do
Y=`echo ${X} | sed 's/_vcf_del_homo//' | sed 's/_vcf_del_zeroone//'`;
bedtools intersect -f 0.90 -v -a rmblast/dm6_rmblast_mask.bed -b vcfprocess/${X}.bed > vcfprocess/embl_${Y}_del_repeat.bed;
done;

# pacbio data coverage
samtools depth /media/pericles/TEfind/svcall/pacbio_pbsv.bam | \
 awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' > /media/pericles/TEfind/OSC_coverage.txt;
# https://www.biostars.org/p/5165/

# get fasta seq
bedtools getfasta -s -fi database/dm6_chr.fasta \
 -bed vcfprocess/embl_pbsv_del_repeat.bed \
 > rmblast/TE_genome.fasta;
seqkit fx2tab rmblast/TE_genome.fasta | cut -f 1-2 > rmblast/TE_genome.tsv;

# pbsv insertion
cd /media/pericles/TEfind/;
mkdir vcfprocess/rmblast/;
mkdir vcfprocess/fig/;
source activate repeatmasker
for X in pbsv ;do
RepeatMasker -s -pa 18 -e rmblast -gff -norna -nolow -no_is \
-lib database/TE_dm.fasta \
-dir vcfprocess/rmblast/ \
vcfprocess/${X}_vcf_ins_read.fasta \
> vcfprocess/rmblast/${X}_of_log.txt;
done;
conda deactivate;
/usr/bin/Rscript  vcf_repmask.R;
bedtools intersect -wo -a /media/pericles/TEfind/database/annotateTE/TE_bed6_full.bed \
 -b /media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_ins_locus.bed > \
/media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_ins_annotate.bed;
bedtools intersect -wo -a /media/pericles/TEfind/database/annotateTE/TE_bed6_full.bed \
 -b /media/pericles/TEfind/vcfprocess/rmblast/dm6_rmblast_locus.bed > \
/media/pericles/TEfind/vcfprocess/rmblast/dm6_rmblast_annotate.bed;
/usr/bin/Rscript  vcf_annotate.R;

# fullness consideration
bedtools intersect -f 0.80 -wo -a /media/pericles/TEfind/database/annotateTE/TE_bed6_full.bed \
 -b /media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_ins_locus.bed > \
/media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_ins_annotate2.bed;

cd /media/pericles/TEfind/;
mkdir minimap2/;
for A in rmblast/TE_genome vcfprocess/pbsv_vcf_ins_read ;do
X=`echo ${A} | sed 's/.*\///'`;
minimap2 -t 20 -ax splice:hq \
 /media/pericles/TEfind/database/TE_dm.fasta \
 /media/pericles/TEfind/${A}.fasta > \
 /media/pericles/TEfind/minimap2/${X}.sam;
samtools view -@ 20 -bS minimap2/${X}.sam > minimap2/${X}.bam;
samtools sort -@ 20 -m 4G minimap2/${X}.bam > minimap2/${X}_sort.bam;
samtools index minimap2/${X}_sort.bam;
done;

cd /media/pericles/TEfind/;
mkdir compare/;
cat vcfprocess/rmblast/dm6_pbsv_ins_venn.bed tebreakres/TEinssite.bed > compare/ins3.bed;
bedtools sort -i compare/ins3.bed > compare/ins4.bed;
bedtools cluster -d 50 -i compare/ins4.bed > compare/cluster2.bed;
bedtools intersect -v -a compare/cluster2.bed -b database/centromere_telomere.txt > compare/cluster2_remove.bed;



###################################################################################################################
# merge vcf and tldr
###################################################################################################################

cd /media/pericles/CLASH/database/anno/;
bedtools intersect -wa -wb -a /media/pericles/CLASH/database/anno/gene.bed  \
 -b /media/pericles/TEfind/rmblast/dm6_rmblast_mask.bed > /media/pericles/CLASH/database/anno/repeat_gene.txt;
bedtools intersect -wa -wb -a /media/pericles/CLASH/database/anno/gene.bed \
 -b /media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_ins.bed > /media/pericles/CLASH/database/anno/insertion_gene.txt;

mkdir /media/pericles/TEfind/vcfprocess/telomere/;
bedtools intersect -wa -wb -a /media/pericles/TEfind/vcfprocess/pbsv_vcf_ins.bed \
 -b /media/pericles/TEfind/database/centromere_telomere.txt > \
 /media/pericles/TEfind/vcfprocess/telomere/pbsv_vcf_ins_telomere.txt;


###################################################################################################################
# genome assembly
###################################################################################################################
cd /media/pericles/TEfind/;
mkdir canu;
mkdir kmc/;
cd /media/pericles/TEfind/kmc/;
mkdir kmc/;
mkdir kmc2/;
mkdir kmc3/;
source activate busco
kmc -k32 -t16 -m64 -cs10000 -fq /media/pericles/TEfind/pacbioDNA/OSC_genome_HiFi.fastq.gz database /media/pericles/TEfind/kmc/;
kmc_tools transform database histogram reads.histo -cx10000;
conda deactivate;

/usr/bin/Rscript  ~/genomescope2.0/genomescope.R -i reads.histo -k 32 -p 1 -o kmc > kmc/log.txt;
/usr/bin/Rscript  ~/genomescope2.0/genomescope.R -i reads.histo -k 32 -p 2 -o kmc2 > kmc2/log.txt;
/usr/bin/Rscript  ~/genomescope2.0/genomescope.R -i reads.histo -k 32 -p 4 -o kmc3 > kmc3/log.txt;

# hifiasm
cd /media/pericles/TEfind/;
mkdir hifiasm;
cd hifiasm/;
source activate hifiasm
hifiasm -t 28 --primary -D 10 -o OSC_hifiasm.asm \
--h1 /media/hermione/HiC_OSC/siEGFP1/siEGFP1_1.fastq.gz,/media/hermione/HiC_OSC/siEGFP2/siEGFP2_1.fastq.gz \
--h2 /media/hermione/HiC_OSC/siEGFP1/siEGFP1_2.fastq.gz,/media/hermione/HiC_OSC/siEGFP2/siEGFP2_2.fastq.gz \
/media/pericles/TEfind/pacbioDNA/OSC_genome_HiFi.fastq.gz 2> log.txt;
conda deactivate;
awk '/^S/{print ">"$2;print $3}' OSC_hifiasm.asm.hic.p_ctg.gfa > OSC_hifiasm.asm.hic.p_ctg.fa;
awk '/^S/{print ">"$2;print $3}' OSC_hifiasm.asm.hic.hap1.p_ctg.gfa > OSC_hifiasm.asm.hic.hap1.p_ctg.fa;
awk '/^S/{print ">"$2;print $3}' OSC_hifiasm.asm.hic.hap2.p_ctg.gfa > OSC_hifiasm.asm.hic.hap2.p_ctg.fa;
# Real time: 94666.542 sec; CPU: 2200364.883 sec; Peak RSS: 78.777 GB
# use this for flam haploid analysis and other datas.
cd /media/pericles/TEfind/hifiasm/;
mkdir quast/;

# comparison of other datas.
cd /media/pericles/TEfind/hifiasm/;
source activate canu;
quast.py /media/pericles/TEfind/database/dmel-all-chromosome-r6.36.fasta \
-R /media/pericles/TEfind/database/dmel-all-chromosome-r6.36.fasta \
-g /media/pericles/CLASH/database/gtf/dm6.gtf -t 16 --fragmented -o quast/quast2;
conda deactivate;

source activate canu;
quast.py /media/pericles/TEfind/hifiasm/OSC_hifiasm.asm.hic.p_ctg.fa \
/media/pericles/TEfind/hifiasm/OSC_hifiasm.asm.hic.hap1.p_ctg.fa \
/media/pericles/TEfind/hifiasm/OSC_hifiasm.asm.hic.hap2.p_ctg.fa \
-R /media/pericles/TEfind/database/dmel-all-chromosome-r6.36.fasta \
-g /media/pericles/CLASH/database/gtf/dm6.gtf -t 24 --fragmented -o quast/quast7;
conda deactivate;


cd /media/pericles/TEfind/hifiasm/;
source activate canu;
busco -f -i /media/pericles/TEfind/database/dmel-all-chromosome-r6.36.fasta -o busco_dm6_ref \
-l diptera_odb10 -m genome;
busco -f -i /media/pericles/TEfind/hifiasm/OSC_hifiasm.asm.hic.p_ctg.fa -o busco_hifiasm \
-l diptera_odb10 -m genome;
busco -f -i /media/pericles/TEfind/hifiasm/OSC_hifiasm.asm.hic.hap1.p_ctg.fa -o busco_hifiasm_hap1 \
-l diptera_odb10 -m genome;
busco -f -i /media/pericles/TEfind/hifiasm/OSC_hifiasm.asm.hic.hap2.p_ctg.fa -o busco_hifiasm_hap2 \
-l diptera_odb10 -m genome;
conda deactivate;



cd /media/pericles/TEfind/hifiasm/;
mkdir dot/;
minimap2 -t 20 -x asm5 \
/media/pericles/TEfind/database/dm6_chr.fasta \
/media/pericles/TEfind/hifiasm/OSC_hifiasm.asm.hic.p_ctg.fa > \
/media/pericles/TEfind/hifiasm/dot/hifiasm.paf;
minimap2 -t 20 -x asm5 \
/media/pericles/TEfind/database/dm6_chr.fasta \
/media/pericles/TEfind/hifiasm/OSC_hifiasm.asm.hic.hap1.p_ctg.fa > \
/media/pericles/TEfind/hifiasm/dot/hap1.paf;
minimap2 -t 20 -x asm5 \
/media/pericles/TEfind/database/dm6_chr.fasta \
/media/pericles/TEfind/hifiasm/OSC_hifiasm.asm.hic.hap2.p_ctg.fa > \
/media/pericles/TEfind/hifiasm/dot/hap2.paf;

cd /media/pericles/TEfind/hifiasm/;
cd dot/;
for X in hifiasm hap1 hap2 ;do
/usr/bin/Rscript ~/dotPlotly/pafCoordsDotPlotly.R \
-i /media/pericles/TEfind/hifiasm/dot/${X}.paf -o /media/pericles/TEfind/hifiasm/dot/${X}_large.plot \
-m 20000 -q 500000 -s -t -l -a 5 -w 4.5;
done;



###################################################################################################################
# update flam and 20A sequence.
# flam(dm6) is far different from OSC flam sequences.
###################################################################################################################

cd /media/pericles/TEfind/;
mkdir flam/;
mkdir Discas/;
cat /media/pericles/CLASH/database/reffasta/gene.fasta | seqkit grep -r -p FBgn0267704 > flam/flam.fasta;
cat /media/pericles/CLASH/database/reffasta/gene.fasta | seqkit grep -r -p FBgn0287603 > flam/20A.fasta;
mkdir flam/genes/;
cat /media/pericles/CLASH/database/reffasta/gene.fasta | seqkit grep -r -p FBgn0024807 > flam/genes/DIP1.fasta;
cat /media/pericles/CLASH/database/reffasta/gene.fasta | seqkit grep -r -p FBgn0285970 > flam/genes/CG32500.fasta;
cat /media/pericles/CLASH/database/reffasta/gene.fasta | seqkit grep -r -p FBgn0085521 > flam/genes/CG40813.fasta;
cat /media/pericles/CLASH/database/reffasta/gene.fasta | seqkit grep -r -p FBgn0040028 > flam/genes/CG17450.fasta;
cat /media/pericles/CLASH/database/reffasta/gene.fasta | seqkit grep -r -p FBgn0052820 > flam/genes/CG32820.fasta;
cat /media/pericles/CLASH/database/reffasta/gene.fasta | seqkit grep -r -p FBgn0052857 > flam/genes/CG32857.fasta;
cat /media/pericles/CLASH/database/reffasta/gene.fasta | seqkit grep -r -p FBgn0052819 > flam/genes/CG32819.fasta;
cat /media/pericles/CLASH/database/reffasta/gene.fasta | seqkit grep -r -p FBgn0027588 > flam/genes/CG14476.fasta;
cat /media/pericles/CLASH/database/reffasta/gene.fasta | seqkit grep -r -p FBgn0031182 > flam/genes/CG1644.fasta;
cat /media/pericles/CLASH/database/reffasta/gene.fasta | seqkit grep -r -p FBgn0031181 > flam/genes/CG14584.fasta;
cat flam/genes/*.fasta > flam/gene_all.fasta;

cd /media/pericles/TEfind/;
mkdir flam/hap/;
for Y in flam 20A;do
for X in p;do
for Z in hap1 hap2;do
minimap2 -ax map-pb  flam/${Y}.fasta \
/media/pericles/TEfind/hifiasm/OSC_hifiasm.asm.hic.${Z}.p_ctg.fa > flam/hap/overlaps_${Y}_${Z}.sam;
samtools view -@ 20 -F 4 flam/hap/overlaps_${Y}_${Z}.sam | cut -f 1-6 > flam/hap/overlaps_${Y}_${Z}_hit.sam;
samtools view -@ 20 -bS flam/hap/overlaps_${Y}_${Z}.sam > Discas/overlaps_${Y}_${Z}.bam;
samtools sort -@ 20 -m 4G Discas/overlaps_${Y}_${Z}.bam > Discas/overlaps_${Y}_${Z}_sort.bam;
samtools index Discas/overlaps_${Y}_${Z}_sort.bam;
bedtools bamtobed -i Discas/overlaps_${Y}_${Z}_sort.bam > flam/hap/overlaps_${Y}_${Z}_sort;
done;
done;
done;
seqkit grep -r -p h1tg000004l /media/pericles/TEfind/hifiasm/OSC_hifiasm.asm.hic.hap1.p_ctg.fa | seqkit seq -w 0 > flam/hap/flam_hap1.fasta;
samtools faidx flam/hap/flam_hap1.fasta;
seqkit grep -r -p h2tg000007l /media/pericles/TEfind/hifiasm/OSC_hifiasm.asm.hic.hap2.p_ctg.fa | seqkit seq -w 0 > flam/hap/flam_hap2.fasta;
samtools faidx flam/hap/flam_hap2.fasta;
for X in hap1 hap2 ;do
minimap2 -ax map-pb   flam/hap/flam_${X}.fasta flam/flam.fasta > flam/hap/flam1_${X}.sam;
samtools view -@ 20 -bS flam/hap/flam1_${X}.sam > flam/hap/flam1_${X}.bam;
samtools sort -@ 20 -m 4G flam/hap/flam1_${X}.bam > flam/hap/flam1_${X}_sort.bam;
samtools index flam/hap/flam1_${X}_sort.bam;
bedtools bamtobed -i flam/hap/flam1_${X}_sort.bam > flam/hap/flam1_${X}_sort;
done;

cd /media/pericles/TEfind/;
mkdir flam/hap/rmblast2/;
mkdir /media/pericles/TEfind/flam/hap/blastnIndex/;
mkdir /media/pericles/TEfind/flam/hap/Discas/;
for X in hap1 hap2 ;do
source activate repeatmasker;
RepeatMasker -s -pa 12 -e rmblast -gff -norna -nolow -no_is \
-lib database/TE_dm.fasta \
-dir flam/hap/rmblast2/ \
flam/hap/flam_${X}.fasta \
> flam/hap/rmblast2/of_log.txt;
conda deactivate;
makeblastdb -in flam/hap/flam_${X}.fasta -dbtype nucl \
-out flam/hap/blastnIndex/${X} -parse_seqids;
blastn -word_size 100 -evalue 1e-40 -outfmt 6 -db flam/hap/blastnIndex/${X} \
 -query flam/gene_all.fasta -out flam/hap/gene_res_${X}.txt;
done;

for X in hap1 hap2;do
  blastn -word_size 100 -evalue 1e-40 -outfmt 6 -db flam/hap/blastnIndex/${X} \
   -query flam/gene_all.fasta -out flam/hap/gene_res_${X}.txt;
done;

/usr/bin/Rscript /media/pericles/TEfind/repmask2.R;

cd /media/pericles/TEfind/;
mkdir flam/hap/piRNA;
mkdir flam/hap/BowtieIndex/;
for Y in hap1 hap2 ;do
bowtie-build --quiet -f flam/hap/flam_${Y}.fasta flam/hap/BowtieIndex/flam_${Y};
for X in piRNA_all; do
bowtie -a --best --strata --threads 20 -l 20 -f -n 3 -S -x flam/hap/BowtieIndex/flam_${Y} \
/media/pericles/CLASH/piRNA/${X}.fasta > Discas/${X}_${Y}_ctg.sam;
samtools view -@ 20 -bS Discas/${X}_${Y}_ctg.sam > Discas/${X}_${Y}_ctg.bam;
samtools sort -@ 20 -m 4G Discas/${X}_${Y}_ctg.bam > flam/hap/piRNA/${X}_${Y}_ctg_sort.bam;
samtools index flam/hap/piRNA/${X}_${Y}_ctg_sort.bam;
grep -v "XS:" Discas/${X}_${Y}_ctg.sam > Discas/${X}_${Y}uniq.sam;
samtools view -@ 20 -bS Discas/${X}_${Y}uniq.sam > Discas/${X}_${Y}uniq.bam;
samtools sort -@ 20 -m 4G Discas/${X}_${Y}uniq.bam > flam/hap/piRNA/${X}_${Y}uniq_sort.bam;
samtools index flam/hap/piRNA/${X}_${Y}uniq_sort.bam;
source activate macs;
for bin in 1 10;do
bamCoverage --filterRNAstrand forward --binSize ${bin} -p 20 --normalizeUsing None --scaleFactor -1 \
 -b flam/hap/piRNA/${X}_${Y}_ctg_sort.bam -o flam/hap/piRNA/${X}_${Y}_reverse_${bin}.bw ;
bamCoverage --filterRNAstrand reverse --binSize ${bin} -p 20 --normalizeUsing None  \
 -b flam/hap/piRNA/${X}_${Y}_ctg_sort.bam -o flam/hap/piRNA/${X}_${Y}_forward_${bin}.bw  ;
bamCoverage --filterRNAstrand forward --binSize ${bin} -p 20 --normalizeUsing None --scaleFactor -1 \
 -b flam/hap/piRNA/${X}_${Y}uniq_sort.bam -o flam/hap/piRNA/${X}_${Y}uniq_reverse_${bin}.bw ;
bamCoverage --filterRNAstrand reverse --binSize ${bin} -p 20 --normalizeUsing None  \
 -b flam/hap/piRNA/${X}_${Y}uniq_sort.bam -o flam/hap/piRNA/${X}_${Y}uniq_forward_${bin}.bw  ;
done;
conda deactivate;
done;
done;





cd /media/pericles/TEfind/;
for Y in hap1 hap2 ;do
for X in OSC_genome_HiFi; do
minimap2 -t 20 -ax map-hifi --sam-hit-only --secondary=no -k 25 --MD /media/pericles/TEfind/flam/hap/flam_${Y}.fasta \
 /media/pericles/TEfind/pacbioDNA/OSC_genome_HiFi.fastq.gz > \
 /media/pericles/TEfind/flam/Discas/${X}.sam;
samtools view -@ 20 -bS flam/Discas/${X}.sam > flam/Discas/${X}.bam;
samtools sort -@ 20 -m 4G flam/Discas/${X}.bam > flam/${X}_${Y}_sort.bam;
samtools index flam/${X}_${Y}_sort.bam;
done;
done;
rm flam/Discas/*;

cd /media/pericles/TEfind/;
mkdir flam/new_flam/;
cat flam/hap/flam_hap1.fasta | seqkit subseq -r 22482068:22871420 | seqkit seq -t DNA | sed 's/>h1tg000004l/>flam_hap1/g' > flam/new_flam/flam_hap1_region.fasta;
cat flam/hap/flam_hap1.fasta | seqkit subseq -r 22367007:22404696 | seqkit seq -t DNA | sed 's/>h1tg000004l/>20A_hap1/g'> flam/new_flam/20A_hap1_region.fasta;
cat flam/hap/flam_hap2.fasta | seqkit subseq -r 22424732:23147602 | seqkit seq -t DNA | sed 's/>h2tg000007l/>flam_hap2/g' > flam/new_flam/flam_hap2_region.fasta;
cat flam/hap/flam_hap2.fasta | seqkit subseq -r 22302607:22340362 | seqkit seq -t DNA | sed 's/>h2tg000007l/>20A_hap2/g'> flam/new_flam/20A_hap2_region.fasta;
for X in flam_hap1 flam_hap2 20A_hap1 20A_hap2 ;do
seqkit fx2tab flam/new_flam/${X}_region.fasta | cut -f 1-2 > flam/new_flam/${X}_region.tsv;
samtools faidx flam/new_flam/${X}_region.fasta;
done;

cd /media/pericles/TEfind/;
mkdir flam/rmblast/;
source activate repeatmasker;
for X in flam_hap1 flam_hap2 20A_hap1 20A_hap2;do
RepeatMasker -s -pa 18 -e rmblast -gff -norna -nolow -no_is \
-lib database/TE_dm.fasta \
-dir flam/rmblast/ \
flam/new_flam/${X}_region.fasta \
> flam/rmblast/${X}_log.txt;
done;
conda deactivate;

/usr/bin/Rscript  repmask.R;

# piRNA
cd /media/pericles/TEfind/;
mkdir flam/new_piRNA;
mkdir flam/new_piRNA/Bowtie2Index/;
mkdir flam/new_piRNA/Discas/;
for Y in flam_hap1 flam_hap2 20A_hap1 20A_hap2;do
bowtie2-build --quiet -f flam/new_flam/${Y}_region.fasta flam/new_piRNA/Bowtie2Index/${Y};
for X in piRNA_all; do
bowtie2 -f -p 20 -a -x flam/new_piRNA/Bowtie2Index/${Y} \
-U /media/pericles/CLASH/piRNA/${X}.fasta -S flam/new_piRNA/Discas/${X}_${Y}_ctg.sam;
samtools view -@ 20 -bS flam/new_piRNA/Discas/${X}_${Y}_ctg.sam > flam/new_piRNA/Discas/${X}_${Y}_ctg.bam;
samtools sort -@ 20 -m 4G flam/new_piRNA/Discas/${X}_${Y}_ctg.bam > flam/new_piRNA/Discas/${X}_${Y}_ctg_sort.bam;
samtools index flam/new_piRNA/Discas/${X}_${Y}_ctg_sort.bam;
grep -v "XS:" flam/new_piRNA/Discas/${X}_${Y}_ctg.sam > flam/new_piRNA/Discas/${X}_${Y}uniq.sam;
samtools view -@ 20 -bS flam/new_piRNA/Discas/${X}_${Y}uniq.sam > flam/new_piRNA/Discas/${X}_${Y}uniq.bam;
samtools sort -@ 20 -m 4G flam/new_piRNA/Discas/${X}_${Y}uniq.bam > flam/new_piRNA/Discas/${X}_${Y}uniq_sort.bam;
samtools index flam/new_piRNA/Discas/${X}_${Y}uniq_sort.bam;
bedtools bamtobed -i flam/new_piRNA/Discas/${X}_${Y}_ctg_sort.bam > flam/new_piRNA/${X}_${Y}_ctg.bed;
bedtools bamtobed -i flam/new_piRNA/Discas/${X}_${Y}uniq_sort.bam > flam/new_piRNA/${X}_${Y}uniq.bed;
done;
done;

cd /media/pericles/TEfind/;
for Y in flam_hap1 flam_hap2 20A_hap1 20A_hap2 ;do
for X in OSC_genome_HiFi; do
minimap2 -t 20 -ax map-hifi --sam-hit-only --secondary=no -k 25 --MD flam/new_flam/${Y}_region.fasta \
 /media/pericles/TEfind/pacbioDNA/OSC_genome_HiFi.fastq.gz > \
 /media/pericles/TEfind/flam/Discas/${X}.sam;
 grep -v -E '[0-9]{4,}S' flam/Discas/${X}.sam > flam/Discas/${X}_soft.sam;
samtools view -@ 20 -bS flam/Discas/${X}_soft.sam > flam/Discas/${X}.bam;
samtools sort -@ 20 -m 4G flam/Discas/${X}.bam > flam/${X}_${Y}_sort.bam;
samtools index flam/${X}_${Y}_sort.bam;
done;
done;
rm flam/Discas/*;



###################################################################################################################
# Oat, ex, CG11459, Adgf-A
###################################################################################################################

cd /media/pericles/TEfind/;
mkdir chimera/;
cat /media/pericles/CLASH/database/reffasta/gene.fasta | seqkit grep -r -p FBgn0022774 > chimera/Oat.fasta;
cat /media/pericles/CLASH/database/reffasta/gene.fasta | seqkit grep -r -p FBgn0005564 > chimera/Shal.fasta;
cat /media/pericles/CLASH/database/reffasta/gene.fasta | seqkit grep -r -p FBgn0036749 > chimera/CG7460.fasta;
cat /media/pericles/CLASH/database/reffasta/gene.fasta | seqkit grep -r -p FBgn0285958 > chimera/Fuca.fasta;

cat chimera/Oat.fasta chimera/Shal.fasta chimera/CG7460.fasta chimera/Fuca.fasta > chimera/some_gene.fasta;

cd /media/pericles/TEfind/;
mkdir chimera/blastnIndex/;
for X in hap1 hap2 ;do
makeblastdb -in /media/pericles/TEfind/hifiasm/OSC_hifiasm.asm.hic.${X}.p_ctg.fa -dbtype nucl \
-out chimera/blastnIndex/${X} -parse_seqids;
blastn -word_size 100 -evalue 1e-40 -outfmt 6 -db chimera/blastnIndex/${X} \
 -query chimera/some_gene.fasta -out chimera/gene_res_${X}.txt;
done;

mkdir chimera/chimera_ins/;
seqkit grep -r -p h2tg000030l /media/pericles/TEfind/hifiasm/OSC_hifiasm.asm.hic.hap2.p_ctg.fa | seqkit seq -w 0 | \
seqkit subseq -r 4099750:4108350 | seqkit seq -t DNA | sed 's/>h2tg000030l/>Oat/g' > chimera/chimera_ins/Oat_hap2.fasta;
samtools faidx chimera/chimera_ins/Oat_hap2.fasta;
seqkit grep -r -p h1tg000008l /media/pericles/TEfind/hifiasm/OSC_hifiasm.asm.hic.hap1.p_ctg.fa | seqkit seq -w 0 | \
seqkit subseq -r 4320000:4357000 | seqkit seq -t DNA | sed 's/>h1tg000008l/>Shal/g' > chimera/chimera_ins/Shal_hap1.fasta;
samtools faidx chimera/chimera_ins/Shal_hap1.fasta;
seqkit grep -r -p h1tg000008l /media/pericles/TEfind/hifiasm/OSC_hifiasm.asm.hic.hap1.p_ctg.fa | seqkit seq -w 0 | \
seqkit subseq -r 6310000:6332000 | seqkit seq -t DNA | sed 's/>h1tg000008l/>CG7460/g' > chimera/chimera_ins/CG7460_hap1.fasta;
samtools faidx chimera/chimera_ins/CG7460_hap1.fasta;
seqkit grep -r -p h1tg000005l /media/pericles/TEfind/hifiasm/OSC_hifiasm.asm.hic.hap1.p_ctg.fa | seqkit seq -w 0 | \
seqkit subseq -r 4783000:4796000 | seqkit seq -t DNA | sed 's/>h1tg000005l/>Fuca/g' > chimera/chimera_ins/Fuca_hap1.fasta;
samtools faidx chimera/chimera_ins/Fuca_hap1.fasta;


