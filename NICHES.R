### load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(SeuratWrappers)
library(NICHES)
library(viridis)
library(cowplot)
library(patchwork)

### path
home="Your/home/directory/"

#Load in the Slide-seq gene expression matrix
counts <- read.csv(paste0(home,'MappedDGEForR_RCTD_Weights.csv'))
counts[1:5,1:2]
rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames

counts_t = t(counts)
dim(counts_t) 
counts_t[1:5,1:2]

#Load in the Slide-seq bead coordinate matrix
metaData <- read.csv(paste0(home,'Puck_RCTD_Cell_Type_and_Coords_Weights.csv'))
metaData [1:5,1:2]
rownames(metaData) <- metaData [,1]; metaData [,1] <- NULL # Move first column to rownames
metaData [1:5,1:2]
dim(metaData)

#Create Seurat object
dge <- CreateSeuratObject(counts=counts_t, meta.data = metaData, min.cells = 5, min.features = 0)

rm(counts)
rm(counts_t)

dge <- SetIdent(dge, value = dge@meta.data$CellType)
Idents(dge)

dge <- NormalizeData(dge)

#Impute the data
dge <- SeuratWrappers::RunALRA(dge)
saveRDS(dge, file = paste0(home,"Seurat_Puck_imputed.rds"))

#Run NICHES
NICHES_output <- RunNICHES(object = dge,
                           LR.database = "omnipath",
                           species = "mouse",
                           assay = "alra",
                           position.x = 'x',
                           position.y = 'y',
                           meta.data.to.map = c('orig.ident','CellType'),
                           CellToCell = F,CellToSystem = F,SystemToCell = F,
                           CellToCellSpatial = F,CellToNeighborhood = F,NeighborhoodToCell = T)

#NICHES outputs a list of objects. Each object contains a certain style of cell-system signaling atlas. Above, we have 
#only calculated a single one of interest, namely, individual cellular microenvironment. We next isolate this output and 
#embed using UMAP to visualize the microenvironemnt of each cell.

niche <- NICHES_output[['NeighborhoodToCell']]
Idents(niche) <- niche[['ReceivingType']]

# Add Niches output as an assay
niches.data <- GetAssayData(object =  niche[['NeighborhoodToCell']], slot = 'data')
colnames(niches.data) <- niche[['ReceivingCell']]$ReceivingCell
dge[["NeighborhoodToCell"]] <- CreateAssayObject(data = niches.data )
DefaultAssay(dge) <- "NeighborhoodToCell"

#sub <- ScaleData(sub)
niche_sig = dge$NeighborhoodToCell@data
dim(niche_sig)
niche_sig[1:5,1:5]

write.csv(niche_sig, paste0(home,'puck_niche_signaling_raw_data.csv'))

dge <- ScaleData(dge)
niche_scale = dge$NeighborhoodToCell@scale.data
dim(niche_scale)
niche_scale[1:5,1:5]

write.csv(niche_scale, paste0(home,'puck_niche_signaling_scaled_data.csv'))


# Scale and visualize
niche <- ScaleData(niche)
niche <- FindVariableFeatures(niche,selection.method = "disp")
niche <- RunPCA(niche)
ElbowPlot(niche,ndims = 50)

my_cols <- c('SPG1'='#aeadb3','SPG2'='#E6C122','SPG3'='#1FA195','SPG4'='#B95FBB','Leydig'='#D4D915',
             'Endothelial'='#ff9a36','Sertoli'='#CCB1F1','ES'='#faf4cf',
             'RS'='#A4DFF2','SPC'='#F68282', 'T'='#2FF18B')

niche <- RunUMAP(niche,dims = 1:10)
DimPlot(niche,reduction = 'umap',pt.size = 0.1,shuffle = T, label = F, cols = my_cols) +NoAxes()+
  ggtitle('Cellular Microenvironment')+guides(colour=guide_legend(ncol=1,override.aes = list(size=6)))

# Find cell type markers
mark <- FindAllMarkers(niche,min.pct = 0.25,only.pos = T,test.use = "roc")
GOI_niche <- mark %>% group_by(cluster) %>% top_n(5,myAUC)
DoHeatmap(niche,features = unique(GOI_niche$gene))+ 
  scale_fill_viridis(option = "E")

write.csv(mark, paste0(home,'niche_signaling_markers.csv'))

# Find stage-dependent markers
stage <- read.csv(paste0(home,'Puck_Cluster_and_Stage_Weights.csv'))
stage[1:5,1:2]
rownames(stage) <- stage[,1]; stage[,1] <- NULL # Move first column to rownames
stage[1:5,]

cells.use <- stage$barcode
subset_dge <- subset(dge, cells = cells.use)
subset_dge 

subset_dge<-AddMetaData(subset_dge , metadata=stage, col.name = 'Stage')
subset_dge@meta.data[1:5, 1:7]

subset_dge <- SetIdent(subset_dge, value = subset_dge@meta.data$Stage)
Idents(subset_dge)

stage_mark <- FindAllMarkers(subset_dge,min.pct = 0.25,only.pos = T,test.use = "roc")
GOI_niche <- stage_mark %>% group_by(cluster)
DoHeatmap(subset_dge,features = unique(GOI_niche$gene))+ 
  scale_fill_viridis(option = "E")

write.csv(stage_mark, paste0(home,'niche_signaling_stage_markers.csv'))

#Let's only look at a subset of cell types with sufficient number of cells
sub <- subset(dge,idents = c('SPG1','SPG2','SPG3','SPG4','Leydig',
                             'Endothelial','Sertoli'))

NICHES_output <- RunNICHES(object = sub,
                           LR.database = "omnipath",
                           species = "mouse",
                           assay = "alra",
                           position.x = 'x',
                           position.y = 'y',
                           meta.data.to.map = c('orig.ident','CellType'),
                           CellToCell = FALSE,CellToSystem = F,SystemToCell = F,
                           CellToCellSpatial = F,CellToNeighborhood = F,NeighborhoodToCell = TRUE)

#NICHES outputs a list of objects. Each object contains a certain style of cell-system signaling atlas. Above, we have 
#only calculated a single one of interest, namely, individual cellular microenvironment. We next isolate this output and 
#embed using UMAP to visualize the microenvironemnt of each cell.

niche <- NICHES_output[['NeighborhoodToCell']]
Idents(niche) <- niche[['ReceivingType']]

# Scale and visualize
niche <- ScaleData(niche)
niche <- FindVariableFeatures(niche,selection.method = "disp")
niche <- RunPCA(niche)
ElbowPlot(niche,ndims = 50)

my_cols <- c('SPG1'='#F68282','SPG2'='#E6C122','SPG3'='#1FA195','SPG4'='#B95FBB','Leydig'='#D4D915',
             'Endothelial'='#ff9a36','Sertoli'='#CCB1F1') 

niche <- RunUMAP(niche,dims = 1:10)
DimPlot(niche,reduction = 'umap',pt.size = 0.5,shuffle = T, label = F, cols = my_cols) +NoAxes()+
  ggtitle('Cellular Microenvironment')+guides(colour=guide_legend(ncol=1,override.aes = list(size=6)))

# Find markers
mark <- FindAllMarkers(niche,min.pct = 0.25,only.pos = T,test.use = "roc")
write.csv(mark, paste0(home,'niche_signaling_puck_markers.csv'))

GOI_niche <- mark %>% group_by(cluster) %>% top_n(5,myAUC)
DoHeatmap(niche,features = unique(GOI_niche$gene))+ 
  scale_fill_viridis(option = "E")


#This gives us a sense of the character of each celltype niche. However, it doesn’t allow us to see 
#where the signal is coming from. To ask these kinds of questions, we need to characterize cell-cell 
#relationships, which we do as follows:

scc.imputed <- RunNICHES(object = sub,
                         LR.database = "omnipath",
                         species = "mouse",
                         assay = "alra",
                         position.x = 'x',
                         position.y = 'y',
                         meta.data.to.map = c('orig.ident','CellType'),
                         CellToCell = FALSE,CellToSystem = F,SystemToCell = F,
                         CellToCellSpatial = TRUE,CellToNeighborhood = F,
                         NeighborhoodToCell = F)

demo.2 <- scc.imputed$CellToCell
demo.2 <- ScaleData(demo.2)
demo.2 <- FindVariableFeatures(demo.2)
demo.2 <- RunPCA(demo.2)
ElbowPlot(demo.2,ndims=50)

PCHeatmap(demo.2,dims = 1:10,balanced = T,cells = 100)

demo.2 <- RunUMAP(demo.2,dims = 1:10)
DimPlot(demo.2,reduction = 'umap',group.by = 'VectorType',label = F, label.size = 6) + 
  NoAxes()+ ggtitle('SPG Niche Cell-Cell Signaling')+
  guides(colour=guide_legend(ncol=3,override.aes = list(size=6)))

#We are now prepared to fully dissect a given cellular niche. 
#Let’s see how different cells influence the SPG1 population within this system:

Idents(demo.2) <- demo.2[['ReceivingType']]
SPG1.network <- subset(demo.2,idents ='SPG1')
Idents(SPG1.network) <- SPG1.network[['VectorType']]
mark.SPG1 <- FindAllMarkers(SPG1.network,
                            logfc.threshold = 0.5, 
                            min.pct = 0.5,
                            only.pos = T,
                            test.use = 'roc')
# Pull markers of interest to plot
mark.SPG1$ratio <- mark.SPG1$pct.1/mark.SPG1$pct.2
marker.list.SPG1 <- mark.SPG1 %>% group_by(cluster) %>% top_n(5,avg_log2FC)
write.csv(mark.SPG1, paste0(home,"Niche_Cell_to_SPG1_Signaling_Markers.csv"))

#Plot in Heatmap form
DoHeatmap(SPG1.network,features = marker.list.SPG1$gene,cells = WhichCells(SPG1.network,downsample = 100))+ 
  scale_fill_viridis(option = "E")

