# Slideseq_SSC_Niche

### This repository contains custom code for analyzing data reported in the following paper:
     
     Shreya, R., et al. (2023). Dissecting the Spermatogonial Stem Cell Niche Using Spatial Transcriptomics.

For the previously generated, processed wild type and diabetic mouse Slide-seq datasets, please go to https://www.dropbox.com/s/ygzpj0d0oh67br0/Testis_Slideseq_Data.zip?dl=0.

For the newly generated, proccessed wild type mouse Slide-seq dataset, please go to https://cloud.biohpc.swmed.edu/index.php/s/bb6Jsj5JjM4pezW

For the processed normal human Slide-seq datasets, please go to https://www.dropbox.com/s/q5djhy006dq1yhw/Human.7z?dl=0.   
   
    
    1. MappedDGE is the digital gene expression matrix;
    
    2. Beadlocations is the bead location matrix;
    

### Seminiferous Tubule Assignment Workflow.ipynb 
   
    Convert Slide-seq data to image data and segment individual seminiferous tubules.
   
### RCTD.R

    Calculate cell type weights for each Slide-seq bead.
    
### Post-RCTD Analyses.ipynb

    Assign cell type based on the cell type weights calculated by RCTD, plot data, and calculate the neighborhood enrichment score for the SPG subtypes.

### NICHES.R

    Identify stage-depdent and cell type-specific LR pairs.

### Post-RCTD Analyses.ipynb

    Calculate the ES purity score for wild type and diabetic seminiferous tubules. 
    
### Pairwise Spatial Contact Frequency Analysis.m 

    Calculate the pairwise spatial contact frequency for wild type and diabetic seminiferous tubules. 


