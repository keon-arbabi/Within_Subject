library(tidyverse)
library(Seurat)

# set directory 
setwd("/")

# get reference seurat object 
seu_ref = readRDS("/external/rprshnas01/netdata_kcni/stlab/Public/Seurat_objects/Seu_AIBS_obj_update_07JUN21.rds")
# get seurat object you want to map to
load("/external/rprshnas01/kcni/karbabi/R Projects/Within_Subject/Revision_Analysis/integration_scRNAseq/OUTPUTS/SME_Data_Integrated.RData")
seu_map = SME_Data_Integrated
rm(SME_Data_Integrated)

# find variable features for transferring 
seu_ref = FindVariableFeatures(seu_ref, selection.method = "vst", nfeatures = 2000, verbose = TRUE)
seu_map = FindVariableFeatures(seu_map, selection.method = "vst", nfeatures = 2000, verbose = TRUE)

# transfer
tanchors = FindTransferAnchors(reference = seu_ref, query = seu_map, dims = 1:30) # find anchors for mapping/transferring/integrating data
predictions = TransferData(anchorset = tanchors, 
                           refdata = seu_ref$subclass_label_expanded_L35IT, #this is the key part of specifying what class/subclass/cluster labels you want to map onto
                           dims = 1:30)
# get relevant variables 
metada_to_add = predictions[, c("predicted.id", "prediction.score.max")]
# add to metadata 
seu_map = AddMetaData(seu_map, metadata = metada_to_add)
# save 
save(seu_map, file = "/external/rprshnas01/kcni/karbabi/R Projects/Within_Subject/keon_analyses/output/SME_Data_Integrated_Our_Labels.rds")

# compare 
seu_map = readRDS("/external/rprshnas01/kcni/karbabi/R Projects/Within_Subject/keon_analyses/output/SME_Data_Integrated_Our_Labels.rds")
DimPlot(object = seu_map, reduction = "umap", label = TRUE, pt.size = 0.5) +
  theme(legend.position="none")