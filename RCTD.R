#Load packages
library(spacexr)
library(Matrix)
library(doParallel)
library(ggplot2)

#Set home directory
home="your/home/directory"

#Read in the scRNA-seq reference data
counts <- read.csv(paste0(home,'scRNAseq_Testis_Ref.csv'))
counts[1:5,1:2]
rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
counts[1:5,1:2]

#Read in the cell type info for the scRNA-seq data
metaData <- read.csv(paste0(home,'Cell_Ident.csv'))
metaData [1:5,1:2]

cell_types <- metaData$CellType; names(cell_types) <- metaData$barcode # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- colSums(counts)

#Generate reference object
reference <- Reference(counts, cell_types, nUMI)
print(dim(reference@counts)) #observe Digital Gene Expression matrix
table(reference@cell_types) #number of occurences for each cell type
saveRDS(reference, paste0(home,"SCRef.RDS"))

#Read in the gene expression matrix and the bead coordinate matrix of the Slide-seq data
counts <- read.csv(paste0(home,"MappedDGEForR.csv")) # load in counts matrix
coords <- read.csv(paste0(home,"BeadLocationsForR.csv")) # load in coordinate matrix

counts[1:5,1:2]
rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
counts[1:5,1:2]
dim(counts)

coords[1:5,1:2]
rownames(coords) <- coords$barcodes; coords$barcodes <- NULL # Move barcodes to rownames
coords[1:5,1:2]
dim(coords)

nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI

### Create SpatialRNA object
puck <- SpatialRNA(coords, counts, nUMI)

## Examine SpatialRNA object (optional)
print(dim(puck@counts)) # observe Digital Gene Expression matrix
hist(log(puck@nUMI,2)) # histogram of log_2 nUMI
print(head(puck@coords)) # start of coordinate data.frame
barcodes <- colnames(puck@counts) # pixels to be used (a list of barcode names). 

saveRDS(puck, file = paste0(home,'Puck.rds'))

myRCTD <- create.RCTD(puck, reference, max_cores = 8)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')

saveRDS(myRCTD,paste0(home,'myRCTD_testes.rds'))

results <- myRCTD@results

# normalize the cell type proportions to sum to 1.
norm_weights = normalize_weights(results$weights) 
cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- myRCTD@spatialRNA
resultsdir <- 'Your/result/directory/' 
dir.create(resultsdir)

cell_type_df <-results$results_df
write.csv(cell_type_df, paste0(resultsdir,"RCTD_Cell_Type.csv"))
write.csv(norm_weights, paste0(resultsdir,"RCTD_Cell_Type_Weights.csv"))

# make the plots 
# Plot the confident weights for each cell type as in full_mode (saved as 
# 'results/cell_type_weights.pdf')
plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights) 

# Plot all weights for each cell type as in full_mode. (saved as 
# 'results/cell_type_weights_unthreshold.pdf')
plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 

# Plot the weights for each cell type as in doublet_mode. (saved as 
# 'results/cell_type_weights_doublets.pdf')
plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, 
                     results$results_df) 

# Plot the number of confident pixels of each cell type in 'full_mode'. (saved as 
# 'results/cell_type_occur.pdf')
plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)

# make a map of all cell types, (saved as 
# 'results/all_cell_types.pdf')
plot_all_cell_types(results$results_df, spatialRNA@coords, cell_type_names, resultsdir) 

# doublets
#obtain a dataframe of only doublets
doublets <- results$results_df[results$results_df$spot_class == "doublet_certain",] 

# Plot all doublets in space (saved as 
# 'results/all_doublets.pdf')
plot_doublets(spatialRNA, doublets, resultsdir, cell_type_names) 

# Plot all doublets in space for each cell type (saved as 
# 'results/all_doublets_type.pdf')
plot_doublets_type(spatialRNA, doublets, resultsdir, cell_type_names) 

# a table of frequency of doublet pairs 
doub_occur <- table(doublets$second_type, doublets$first_type) 

# Plots a stacked bar plot of doublet ocurrences (saved as 
# 'results/doublet_stacked_bar.pdf')
plot_doub_occur_stack(doub_occur, resultsdir, cell_type_names) 
