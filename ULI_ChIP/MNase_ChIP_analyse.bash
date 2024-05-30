

cd /media/hermione/MNaseChIP/
mkdir map/;
mkdir Discas/;
mkdir genome/;
for X in sense anti ;do
for Y in ni input H3K4me3 H3K9me3;do
cutadapt -j 10 -q 20 -m 50 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
-o map/${X}_${Y}_1.fastq.gz ${X}_${Y}_1.fastq.gz;

seqkit stats -j 20 map/*_1.fastq.gz | awk ' BEGIN { OFS = ";" } { print($1, $4) }' > read_info.csv;

bowtie2 -p 20 -x /media/pericles/CLASH/database/Bowtie2Index/dm6_genome \
map/${X}_${Y}_1.fastq.gz -S Discas/${X}_${Y}.sam;
samtools view -@ 20 -q 15 -bS Discas/${X}_${Y}.sam > Discas/${X}_${Y}.bam;
samtools sort -@ 20 -m 4G Discas/${X}_${Y}.bam > genome/${X}_${Y}_sort.bam;
samtools index genome/${X}_${Y}_sort.bam;
done;
done;

for X in sense anti ;do
for Y in ni input H3K4me3 H3K9me3;do
bowtie2 -p 20 -a -x /media/hermione/211027_M02097_0280_000000000-K54Y5/Bowtie2Index/luc_ChIP_${X} \
 map/${X}_${Y}_1.fastq.gz -S Discas/${X}_${Y}_luc.sam;
samtools view -@ 20 -bS Discas/${X}_${Y}_luc.sam > Discas/${X}_${Y}_luc.bam;
samtools sort -@ 20 -m 4G Discas/${X}_${Y}_luc.bam > genome/${X}_${Y}_luc_sort.bam;
samtools index genome/${X}_${Y}_luc_sort.bam;
done;
done;


mkdir bigwig/;
source activate macs;
for X in sense anti ;do
for Y in ni input H3K4me3 H3K9me3;do
for bin in 1 10;do
bamCoverage --binSize ${bin} -p 20 --normalizeUsing None  \
 -b genome/${X}_${Y}_luc_sort.bam -o bigwig/${X}_${Y}_luc_${bin}.bw  ;
bigWigToBedGraph bigwig/${X}_${Y}_luc_${bin}.bw bigwig/${X}_${Y}_luc_${bin}.bedgraph;
bamCoverage --binSize ${bin} -p 20 --normalizeUsing None  \
 -b genome/${X}_${Y}_sort.bam -o bigwig/${X}_${Y}_${bin}.bw  ;
bigWigToBedGraph bigwig/${X}_${Y}_${bin}.bw bigwig/${X}_${Y}_${bin}.bedgraph;
done;
done;
done;
for X in sense anti ;do
for Y in H3K4me3 H3K9me3;do
for bin in 100 10 50 5;do
bamCompare -b1 genome/${X}_${Y}_sort.bam -b2 genome/${X}_input_sort.bam  \
-o bigwig/${X}_${Y}_compare_${bin}.bigwig -of bigwig --binSize=${bin} --minMappingQuality 15 --normalizeUsing CPM --scaleFactorsMethod None -p 20 --operation log2;
done;
done;
done;
conda deactivate


