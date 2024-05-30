

bash databasecostruct.sh 
bash piRNA_call.sh 
bash TEfind.bash

# database making and piRNA decision
bash piRNA_determine.sh 



# CLIP-seq analysis
bash CLIP_mapping.bash 
# GRO-seq analysis
bash GRO_mapping.bash 
# Iso-seq analysis
bash Pacbio_RNA.bash
# CAGE analysis
bash CAGE.bash 
# ChIP-seq analysis
bash MNase_ChIP_analyse.bash 
# ChIP-qPCR
/usr/bin/Rscript  ChIP_qPCR.R 


bash piRNA_luc.bash 
# flam_piRNA.R make many figures
# luciferase result analysis
/usr/bin/Rscript  sample_analyse.R

# CLASH
bash chimera_call.bash 

#other Piwi data analysis
bash AGO1_mapping.bash
bash PRG1_mapping.bash 
bash cCLIP_mapping.bash 
bash SMEDWI3_analyse.bash 
/usr/bin/Rscript  compare.R 

# piRNA_model
bash piRNA_model.sh 


/usr/bin/Rscript  RNA.R;

# get coverage!!!
/usr/bin/Rscript chimera_polish_count.R


/usr/bin/Rscript  fig1.R;
/usr/bin/Rscript  fig2.R;
/usr/bin/Rscript  fig3.R;
/usr/bin/Rscript  fig4.R;

/usr/bin/Rscript  fig5.R;
/usr/bin/Rscript  fig6.R;

# other Piwi comparison
bash comparison_analysis.sh
/usr/bin/Rscript  vector_simulation.R;
