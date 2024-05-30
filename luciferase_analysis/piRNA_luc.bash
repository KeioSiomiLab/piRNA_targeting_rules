


cd /media/hermione/piRNA_luc;
mkdir TE;
mkdir TE/mismatch/;
mkdir bowtieIndex/;
mkdir blastnIndex/;
mkdir Discas/;
mkdir rnaup/;
mkdir rnaup_temp/;
mkdir rnaplex/;
mkdir rnaplex_temp/;
mkdir hyb_temp/;
mkdir binding/;
mkdir figure/;
mkdir target_sequence/;
mkdir process/;
mkdir flam/;
mkdir plasmid/;

bowtie-build --quiet -f /media/pericles/TEfind/flam/new_flam/flam_hap1_region.fasta bowtieIndex/flam_region_hap1;
bowtie-build --quiet -f /media/pericles/TEfind/flam/new_flam/flam_hap2_region.fasta bowtieIndex/flam_region_hap2;
bowtie -a --best --strata --threads 4 -l 20 -f -S -x bowtieIndex/flam_region_hap1 \
piRNA_query.fasta > process/piRNA_query1.sam;
bowtie -a --best --strata --threads 4 -l 20 -f -S -x bowtieIndex/flam_region_hap2 \
piRNA_query.fasta > process/piRNA_query2.sam;

/usr/bin/Rscript  flam_piRNA_processs.R

cd /media/hermione/piRNA_luc;
mkdir flam_mapping/;
cat /media/pericles/TEfind/flam/new_flam/flam_hap1_region.fasta | seqkit subseq -r 7108:7195 | seqkit seq -t DNA | \
sed 's/>flam_hap1/>flam_hap1_piRNA134303/g' > flam_mapping/flam_hap1_piRNA134303.fasta;
cat /media/pericles/TEfind/flam/new_flam/flam_hap1_region.fasta | seqkit subseq -r 79590:79675 | seqkit seq -t DNA | \
sed 's/>flam_hap1/>flam_hap1_piRNA268036/g' > flam_mapping/flam_hap1_piRNA268036.fasta;
cat flam_mapping/flam_*.fasta > flam_mapping/all_target.fasta;
makeblastdb -in flam_mapping/all_target.fasta -dbtype nucl \
-out blastnIndex/flam_target -parse_seqids;

seqkit fx2tab flam_mapping/all_target.fasta | cut -f 1-2 > flam_mapping/flam_mapping_target.tsv;
blastn -word_size 14 -evalue 1e-3 -outfmt 6 -db blastnIndex/flam_target \
 -query /media/pericles/CLASH/piRNA/piRNA_all.fasta -out flam_mapping/flam_mapping.txt;


# plasmid target (piRNA134303, piRNA268036, mdg1_anti)
makeblastdb -in target_sequence/piRNA_plasmid_target.fasta -dbtype nucl \
-out blastnIndex/piRNA_plasmid_target -parse_seqids;
seqkit fx2tab target_sequence/piRNA_plasmid_target.fasta | cut -f 1-2 > plasmid/piRNA_plasmid_target.tsv;
for X in piRNA_all; do
blastn -word_size 14 -evalue 1e-3 -outfmt 6 -db blastnIndex/piRNA_plasmid_target \
 -query /media/pericles/CLASH/piRNA/${X}.fasta -out plasmid/plasmid.txt;
done;


# mdg1 mutation calculation
cd /media/hermione/piRNA_luc;
seqkit fx2tab target_sequence/mdg1_anti_mut.fasta | cut -f 1-2 > mdg1_anti_mut.tsv;
/usr/bin/Rscript  piRNA_make_rnaplex.R
find ./rnaplex_temp -name "*.fasta" | sed 's/.fasta//' | parallel --max-procs=16 \
  'RNAplex --temp=26 < rnaplex_temp/{/.}.fasta > rnaplex/{/.}.rnaplex 2> /dev/null;
  /usr/bin/Rscript --silent --slave --vanilla piRNA_make_binding.R hyb_temp/{/.}.tsv rnaplex/{/.}.rnaplex binding/{/.};';

find ./rnaplex_temp2 -name "*.fasta" | sed 's/.fasta//' | parallel --max-procs=16 \
  'RNAplex --temp=26 < rnaplex_temp2/{/.}.fasta > rnaplex2/{/.}.rnaplex 2> /dev/null;
  /usr/bin/Rscript --silent --slave --vanilla piRNA_make_binding.R hyb_temp2/{/.}.tsv rnaplex2/{/.}.rnaplex binding2/{/.};';


cd /media/hermione/piRNA_luc;
/usr/bin/Rscript  /media/hermione/piRNA_luc/Vienna.R
/usr/bin/Rscript  /media/hermione/piRNA_luc/Vienna2.R

/usr/bin/Rscript  flam_piRNA.R;

