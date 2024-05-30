###################################################################################################################
# Isoseq analysis
###################################################################################################################
cd /media/ariel/HN00139860/;
unzip siEGFP_OSC.zip;
unzip siPIWI_OSC.zip;

# move isoseq folder

# https://www.pacb.com/wp-content/uploads/Single-Cell-Iso-Seq-Library-Preparation-Using-SMRTbell-Express-Template-Prep-Kit-2.0-Customer-Training.pdf

###################################################################################################################
# Processing Iso-seq data
###################################################################################################################
cd /media/ariel/Isoseq/;
mkdir process/;
mkdir dataset/;
mkdir sam/;
source activate canu
for X in siPIWI siEGFP; do
  if [[ ${X} = siEGFP ]]; then
    A=005007
  elif [[ ${X} = siPIWI ]]; then
    A=211937
  fi
ccs ${X}_OSC/m54229_201205_${A}.subreads.bam process/${X}.ccs.bam --all --num-threads 16;
bam2fastq -o process/${X}_ccs process/${X}.ccs.bam;
lima process/${X}.ccs.bam primers.fasta process/${X}_demux_ccs.bam --isoseq --dump-clips --num-threads 16;
isoseq3 refine process/${X}_demux_ccs.primer_5p--primer_3p.bam primers.fasta process/${X}_flnc.bam --require-polya -j 16;
bam2fastq -o process/${X}_flnc process/${X}_flnc.bam;
gunzip -k -f /media/ariel/Isoseq/pocess/${X}_flnc.fastq.gz;
samtools index ${X}_OSC/m54229_201205_${A}.subreads.bam;
done;

for X in siPIWI2 siEGFP2; do
  if [[ ${X} = siEGFP2 ]]; then
    A=221109_230936
  elif [[ ${X} = siPIWI2 ]]; then
    A=221110_193755
  fi
ccs ${X}_OSC/m54209_${A}.subreads.bam process/${X}.ccs.bam --all  --num-threads 16;
bam2fastq -o process/${X}_ccs process/${X}.ccs.bam;
lima process/${X}.ccs.bam primers.fasta process/${X}_demux_ccs.bam --isoseq --dump-clips  --num-threads 16;
isoseq3 refine process/${X}_demux_ccs.primer_5p--primer_3p.bam primers.fasta process/${X}_flnc.bam --require-polya -j 16;
bam2fastq -o process/${X}_flnc process/${X}_flnc.bam;
gunzip -k -f /media/ariel/Isoseq/process/${X}_flnc.fastq.gz;
samtools index ${X}_OSC/m54209_${A}.subreads.bam;
done;

ls process/*_flnc.bam > master_flnc.fofn;
isoseq3 cluster master_flnc.fofn process/polish.bam --verbose --use-qvs --log-file process/cluster.log -j 16;
dataset create --type SubreadSet --force process/merged.subreadset.xml \
siEGFP_OSC/*.subreadset.xml siPIWI_OSC/*.subreadset.xml siEGFP2_OSC/*.subreadset.xml siPIWI2_OSC/*.subreadset.xml;
isoseq3 summarize process/polish.bam process/summary.csv;
pbmm2 align process/polish.bam /media/pericles/TEfind/database/dm6_chr.fasta sam/gene_sort.bam --preset ISOSEQ --sort -j 16;
pbmm2 align process/polish.bam /media/pericles/TEfind/database/TE_dm.fasta sam/TE_sort.bam --preset ISOSEQ --sort -j 16;
isoseq3 collapse sam/gene_sort.bam process/merged.subreadset.xml sam/gene.gff \
--max-fuzzy-junction 0 --min-aln-coverage 0.99 --min-aln-identity 0.95;
isoseq3 collapse sam/TE_sort.bam process/merged.subreadset.xml sam/TE.gff \
--max-fuzzy-junction 0 --min-aln-coverage 0.99 --min-aln-identity 0.95;
conda deactivate;

# for QC
cd /media/ariel/Isoseq/;
mkdir ccs_sam/;
for X in siPIWI siEGFP siPIWI2 siEGFP2; do
minimap2 -t 20 -ax splice:hq -uf \
 /media/pericles/TEfind/database/dm6_chr.fasta \
 process/${X}_flnc.fastq.gz > \
 Discas/${X}_genome.sam;
minimap2 -t 20 -ax splice:hq -uf \
 /media/pericles/TEfind/database/TE_dm.fasta \
 process/${X}_flnc.fastq.gz > \
 Discas/${X}_TE.sam;
for Y in TE genome;do
samtools view -@ 20 -bS Discas/${X}_${Y}.sam > Discas/${X}_${Y}.bam;
samtools sort -@ 20 -m 4G Discas/${X}_${Y}.bam > ccs_sam/${X}_${Y}_sort.bam;
samtools index ccs_sam/${X}_${Y}_sort.bam;
done;
done;
# get coverage
for X in siPIWI siEGFP siPIWI2 siEGFP2; do
bedtools genomecov -bga -split -strand + \
-ibam ccs_sam/${X}_TE_sort.bam > ccs_sam/${X}_TE.bedgraph;
done;

source activate macs;
for X in siPIWI siEGFP siPIWI2 siEGFP2; do
bamCoverage --samFlagExclude 16 --binSize 10 -p 20 --normalizeUsing None --outFileFormat bigwig \
 -b ccs_sam/${X}_TE_sort.bam -o ccs_sam/${X}_TE_10.bigwig ;
done;
for X in siPIWI siEGFP; do
bigwigCompare -b1 ccs_sam/${X}_TE_10.bigwig -b2 ccs_sam/${X}2_TE_10.bigwig \
--operation add -p 20 -o ccs_sam/${X}_TE_10.bedgraph -of bedgraph ;
done;
conda deactivate;


# polish cal

for X in polish ;do
minimap2 -t 20 -ax splice -uf --secondary=no -C5 \
 /media/pericles/TEfind/database/dm6_chr.fasta \
 process/${X}.hq.fasta.gz > \
 sam/gene_cupcake.sam;
minimap2 -t 20 -ax splice -uf --secondary=no -C5 \
 /media/pericles/TEfind/database/TE_dm.fasta \
 process/${X}.hq.fasta.gz > \
 sam/TE_cupcake.sam;
gunzip -c process/${X}.hq.fasta.gz > process/${X}.hq.fasta;
for Y in gene TE;do
samtools view -@ 20 -bS sam/${Y}_cupcake.sam > sam/${Y}_cupcake.bam;
samtools sort -@ 20 -m 4G sam/${Y}_cupcake.bam > sam/${Y}_cupcake_sort.bam;
samtools index sam/${Y}_cupcake_sort.bam;
sort -k 3,3 -k 4,4n sam/${Y}_cupcake.sam > sam/${Y}_cupcake_sort.sam;
source activate SQANTI3.env;
export PYTHONPATH=$PYTHONPATH:~/cDNA_Cupcake/sequence/;
mkdir sam/${Y}/;
collapse_isoforms_by_sam.py --input process/${X}.hq.fasta --gen_mol_count \
-s sam/${Y}_cupcake_sort.sam --dun-merge-5-shorter -o sam/${Y}/${Y};
get_abundance_post_collapse.py sam/${Y}/${Y}.collapsed  process/${X}.cluster_report.csv;
conda deactivate;
done;
done;



cd /media/ariel/Isoseq/;
/usr/bin/Rscript  /media/ariel/Isoseq/cupcake_make_flnc_report_normal.R;
source activate SQANTI3.env;
export PYTHONPATH=$PYTHONPATH:~/cDNA_Cupcake/sequence/;
for X in gene TE;do
filter_away_subset.py sam/${X}/${X}.collapsed;
~/cDNA_Cupcake/post_isoseq_cluster/demux_isoseq_with_genome.py \
--mapped_fafq sam/${X}/${X}.collapsed.filtered.rep.fa \
--read_stat sam/${X}/${X}.collapsed.read_stat.txt \
--classify_csv process/flnc_report.csv \
-o sam/${X}/${X}.collapsed.filtered.fl_count.txt;
done;
# mapped_id_rex = re.compile('(PB\S*.\d+[.\d+]?)')   (demux_isoseq_with_genome.py)
# this is fault!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# change demux_isoseq_with_genome.py like this....
# mapped_id_rex = re.compile('(PB.\d+.\d+|PBfusion.\d+)')
conda deactivate;
for X in gene;do
sed -E 's/("([^"]*)")?,/\2\t/g' sam/${X}/${X}.collapsed.filtered.fl_count.txt >\
 sam/${X}/${X}.collapsed.filtered.fl_count.tsv
done;



# SQANTI library modification needed!
# insert this into SQANTI3_report.R & SQANTI3_filter_report.R & compare_ML_variante.R & SQANTI3_rules_filter.R
# .libPaths("/home/qqprosperodd/anaconda3/envs/SQANTI3.env/lib/R/library")

# remove this from  SQANTI3_filter_report.R & SQANTI3_MLfilter.R & SQANTI3_rules_filter.R
#!/bin/bash Rscript

cd /media/ariel/Isoseq/;
mkdir merge/;
cd /media/ariel/Isoseq/merge/;
source activate SQANTI3.env;
export PYTHONPATH=$PYTHONPATH:~/cDNA_Cupcake/sequence/;
for X in gene;do
python ~/SQANTI3/sqanti3_qc.py /media/ariel/Isoseq/sam/${X}/${X}.collapsed.filtered.gff \
/media/pericles/CLASH/database/gtf/dm6_use.gtf \
/media/pericles/TEfind/database/dm6_chr.fasta \
--fl_count /media/ariel/Isoseq/sam/${X}/${X}.collapsed.filtered.fl_count.txt \
--CAGE_peak /media/hermione/CAGE/figures/All.samples.tagClusters_SQANTI.bed \
--polyA_motif /media/pericles/CLASH/database/polyA.txt -t 8 \
--report skip;

mkdir /media/ariel/Isoseq/merge/rules/;
python ~/SQANTI3/sqanti3_filter.py rules \
/media/ariel/Isoseq/merge/${X}.collapsed.filtered_classification.txt \
--gtf /media/ariel/Isoseq/merge/${X}.collapsed.filtered_corrected.gtf \
--faa /media/ariel/Isoseq/merge/${X}.collapsed.filtered_corrected.faa \
--dir /media/ariel/Isoseq/merge/rules/  --skip_report -v;
done;
conda deactivate;

cd /media/ariel/Isoseq/;
mkdir mergeTE/;
cd /media/ariel/Isoseq/mergeTE/;
source activate SQANTI3.env;
export PYTHONPATH=$PYTHONPATH:~/cDNA_Cupcake/sequence/;
for X in TE;do
python ~/SQANTI3/sqanti3_qc.py /media/ariel/Isoseq/sam/${X}/${X}.collapsed.filtered.gff \
/media/pericles/CLASH/database/gtf/dm6_TE.gtf \
/media/pericles/TEfind/database/TE_dm.fasta  \
--fl_count /media/ariel/Isoseq/sam/${X}/${X}.collapsed.filtered.fl_count.txt \
--polyA_motif /media/pericles/CLASH/database/polyA.txt -t 8 \
--report skip;
done;
conda deactivate;



# https://github.com/PacificBiosciences/pbbioconda/issues/164
cd /media/ariel/Isoseq/;
mkdir qual/;
mkdir nanostat/;
source activate nanostat;

for X in siPIWI siEGFP; do
NanoStat --fasta ${X}_OSC/${X}_OSC_subreads.fasta.gz -t 16 > nanostat/${X}_subread_nanostat.txt;
NanoPlot --fasta ${X}_OSC/${X}_OSC_subreads.fasta.gz --loglength -t 16 -o qual/${X}_subread;
done;
for X in siPIWI siEGFP; do
NanoStat --fasta ${X}2_OSC/${X}r2_subreads.fasta.gz -t 16 > nanostat/${X}2_subread_nanostat.txt;
NanoPlot --fasta ${X}2_OSC/${X}r2_subreads.fasta.gz --loglength -t 16 -o qual/${X}2_subread;
done;
for X in siPIWI siEGFP siPIWI2 siEGFP2; do
for Y in ccs flnc ;do
NanoStat --fastq process/${X}_${Y}.fastq.gz -t 16 > nanostat/${X}_${Y}_nanostat.txt;
NanoPlot --fastq process/${X}_${Y}.fastq.gz --loglength -t 16 -o qual/${X}_${Y};
done;
done;
NanoStat --fasta process/polish.hq.fasta.gz -t 16 > nanostat/polish_nanostat.txt;
NanoPlot --fasta process/polish.hq.fasta.gz --loglength -t 16 -o qual/polish;
conda deactivate;

###################################################################################################################
# detect fusion Transcritome
###################################################################################################################

cd /media/ariel/Isoseq/;
mkdir fusion/;
for X in polish ;do
minimap2 -t 20 -ax splice -uf --secondary=no -C5 \
 /media/pericles/TEfind/database/fusion_genome_TE.fasta \
 process/${X}.hq.fasta.gz > \
 /media/ariel/Isoseq/fusion/${X}_isoforms.sam;
 sort -k 3,3 -k 4,4n /media/ariel/Isoseq/fusion/${X}_isoforms.sam > \
 /media/ariel/Isoseq/fusion/${X}_isoforms_sort.sam;
source activate SQANTI3.env;
export PYTHONPATH=$PYTHONPATH:~/cDNA_Cupcake/sequence/;
fusion_finder.py --input /media/ariel/Isoseq/process/${X}.hq.fasta \
 -s /media/ariel/Isoseq/fusion/${X}_isoforms_sort.sam \
-o /media/ariel/Isoseq/fusion/${X}_isoforms_fusion \
--min_locus_coverage_bp 20 -d 1000 --is_flnc;
get_abundance_post_collapse.py /media/ariel/Isoseq/fusion/${X}_isoforms_fusion \
 /media/ariel/Isoseq/process/${X}.cluster_report.csv;
conda deactivate;
done;

cd /media/ariel/Isoseq/;
for X in polish ;do
source activate SQANTI3.env;
export PYTHONPATH=$PYTHONPATH:~/cDNA_Cupcake/sequence/;
~/cDNA_Cupcake/post_isoseq_cluster/demux_isoseq_with_genome.py \
--mapped_fafq fusion/${X}_isoforms_fusion.rep.fa \
--read_stat fusion/${X}_isoforms_fusion.read_stat.txt \
--classify_csv /media/ariel/Isoseq/process/flnc_report.csv  \
-o /media/ariel/Isoseq/fusion/${X}_isoforms_fusion.fl_count.txt;
conda deactivate;
sed -E 's/("([^"]*)")?,/\2\t/g' fusion/${X}_isoforms_fusion.fl_count.txt >\
 fusion/${X}_isoforms_fusion.fl_count.tsv;
mkdir /media/ariel/Isoseq/fusion/chimera/;
cd /media/ariel/Isoseq/fusion/chimera/;
source activate SQANTI3.env;                 
export PYTHONPATH=$PYTHONPATH:~/cDNA_Cupcake/sequence/;
python ~/SQANTI3/sqanti3_qc.py /media/ariel/Isoseq/fusion/${X}_isoforms_fusion.gff \
/media/pericles/CLASH/database/gtf/fusion_genome_TE.gtf \
/media/pericles/TEfind/database/fusion_genome_TE.fasta \
--is_fusion --fl_count /media/ariel/Isoseq/fusion/${X}_isoforms_fusion.fl_count.txt \
--polyA_motif /media/pericles/CLASH/database/polyA.txt --report skip;
conda deactivate;
done;

mkdir /media/pericles/TEfind/overlap/;
for X in polish; do
bedtools window -a /media/pericles/TEfind/newTE/OSCrepeat_pbsv.bed \
-b /media/ariel/Isoseq/fusion/chimera/${X}_isoforms_fusion_corrected.gtf.cds.gff -w 10000 \
> /media/pericles/TEfind/overlap/${X}_overlap.txt;
done;


###################################################################################################################
# chimera view
###################################################################################################################

cd /media/ariel/Isoseq/;

# add reference
for Z in Shal_hap1 CG7460_hap1 Fuca_hap1 Oat_hap2 ;do
minimap2 -t 20 -ax splice -uf --sam-hit-only -C5 \
 /media/pericles/TEfind/chimera/chimera_ins/${Z}.fasta \
 chimera_view/${Z}_selected_gene1.fasta > \
 Discas/${Z}_selected_gene1.sam;
samtools view -bS Discas/${Z}_selected_gene1.sam | \
samtools sort - > Discas/${Z}_selected_gene1.bam;
samtools index Discas/${Z}_selected_gene1.bam;
bedtools bamtobed -bed12 -color 0,153,73 -split \
-i Discas/${Z}_selected_gene1.bam > chimera_view/${Z}_selected_gene1.bed;
done;

# define remove region
for X in polish ;do
for Z in Shal_hap1 CG7460_hap1 Fuca_hap1 Oat_hap2 ;do
minimap2 -t 20 -ax splice -uf --sam-hit-only -C5 \
 /media/pericles/TEfind/chimera/chimera_ins/${Z}.fasta \
 process/${X}.hq.fasta.gz > \
 Discas/${X}_${Z}.sam;
samtools view -@ 20 -bS Discas/${X}_${Z}.sam > Discas/${X}_${Z}.bam;
bedtools intersect -v -f 0.9 -abam Discas/${X}_${Z}.bam \
-b chimera_view/${Z}_remove.bed > Discas/${X}_${Z}_remove.bam;
samtools sort -@ 20 -m 4G Discas/${X}_${Z}_remove.bam > chimera_view/${X}_${Z}_sort.bam;
samtools index chimera_view/${X}_${Z}_sort.bam;
samtools view -h chimera_view/${X}_${Z}_sort.bam > Discas/${X}_${Z}_sort.sam;
grep -v -E '[0-9]{2,}S' Discas/${X}_${Z}_sort.sam > chimera_view/${X}_${Z}_sort.sam;
done;
done;


# edit manually
for Z in Shal_hap1 CG7460_hap1 Fuca_hap1 Oat_hap2 ;do
samtools view -bS chimera_view/polish_${Z}_selected.sam | samtools sort - > Discas/polish_${Z}_selected.bam;
samtools index Discas/polish_${Z}_selected.bam;
bedtools bamtobed -bed12 -color 51,102,255 -split \
-i Discas/polish_${Z}_selected.bam > chimera_view/polish_${Z}_selected.bed;
done;

# using bed is 
# Shal_hap1 CG7460_hap1 Fuca_hap1 Oat_hap2 


