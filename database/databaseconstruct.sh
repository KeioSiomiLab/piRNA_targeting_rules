#########################################
# data base for TE ins detect
#########################################
# https://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/cytoBand.txt.gz
cd /media/pericles/TEfind/;
mkdir database/eachTE/;
mkdir database/annotateTE/;
mkdir Discas/;
mkdir database/allchromosome/;
/usr/bin/Rscript  TEembl_construct.R database/transposon_sequence_set.embl.txt;
genomedata=database/dmel-all-chromosome-r6.36.fasta
genome=dm6
seqkit fx2tab ${genomedata} | cut -f 1-2 > database/${genome}_process.tsv;
/usr/bin/Rscript  chromosome_construct.R database/${genome}_process.tsv;
/usr/bin/Rscript  dfamtobed.R database/dm6_dfam.nrph.hits database/dm6_dfam_nrph;
cat /media/pericles/TEfind/database/cytoBand_igv.txt | sed 's/NA//g'  > /media/pericles/TEfind/database/cytoBand_igv2.txt;
samtools faidx /media/pericles/TEfind/database/dm6_chr.fasta;
awk -v OFS='\t' {'print $1,$2'} /media/pericles/TEfind/database/dm6_chr.fasta.fai \
> /media/pericles/TEfind/database/dm6_chr.genome;

#########################################
# Repeatmasker for embl data
#########################################
cd /media/pericles/TEfind/;
mkdir embl_mask/;
for engine in rmblast; do
mkdir ${engine}/;
source activate repeatmasker;
RepeatMasker -s -pa 18 -e ${engine} -gff -norna -nolow -no_is \
-lib database/TE_dm.fasta \
-dir ${engine}/ \
database/dm6_chr.fasta \
> ${engine}/of_log.txt;
conda deactivate;
done;
/usr/bin/Rscript  Mask_dm6_genome_process.R;
bedtools merge -i /media/pericles/TEfind/rmblast/dm6_rmblast_mask.bed > \
/media/pericles/TEfind/rmblast/dm6_rmblast_merge.bed;

source activate macs;
findRestSite --fasta /media/pericles/TEfind/database/dm6_chr.fasta \
--searchPattern 'AAGCTT' -o /media/pericles/TEfind/database/HindIIIrestsites.bed;
findRestSite --fasta /media/pericles/TEfind/database/dm6_chr.fasta \
--searchPattern 'GATC' -o /media/pericles/TEfind/database/DpnIIrestsites.bed;
conda deactivate;

