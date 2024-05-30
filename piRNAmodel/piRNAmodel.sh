

###################################################################################################################
# model calculation
###################################################################################################################

# first, piRNA is mapped to TEs.

cd /media/hermione/piRNAmodel;
mkdir bowtieIndex/;
mkdir blastnIndex/;
mkdir blastout/;
makeblastdb -in /media/pericles/CLASH/database/reffasta/ensembl.fasta -dbtype nucl \
-out blastnIndex/TE -parse_seqids;
makeblastdb -in /media/pericles/CLASH/database/reffasta/gene.fasta -dbtype nucl \
-out blastnIndex/gene -parse_seqids;
seqkit fx2tab /media/pericles/CLASH/database/reffasta/ensembl.fasta | cut -f 1-2 > TE.tsv;
seqkit fx2tab /media/pericles/CLASH/database/reffasta/gene.fasta | cut -f 1-2 > gene.tsv;

blastn -word_size 15 -evalue 1e-3 -outfmt 6 -db blastnIndex/TE \
 -query /media/pericles/CLASH/piRNA/piRNA_all.fasta -out blastout/piRNA.txt;
blastn -word_size 15 -evalue 1e-3 -outfmt 6 -db blastnIndex/gene \
 -query /media/pericles/CLASH/piRNA/piRNA_all.fasta -out blastout/gene.txt;


mkdir rnaplex/;

/usr/bin/Rscript  piRNA_plex_prepare.R;
mkdir vienna/;


parallel --max-procs=16 \
 'RNAplex --temp=26 < rnaplex/{/.}.fasta > vienna/{/.}.rnaplex 2> /dev/null;
 /usr/bin/Rscript piRNA_vienna_prepare.R {/.};' ::: model_RNAplex model_RNAplex2;

bedtools intersect -wa -wb -a /media/hermione/piRNAmodel/rnaplex/gene_remove1.bed \
-b /media/pericles/CLASH/database/anno/dm6_gene_ano.bed > /media/hermione/piRNAmodel/rnaplex/gene_remove1_overlap.tsv;

/usr/bin/Rscript piRNA_vienna_table.R


