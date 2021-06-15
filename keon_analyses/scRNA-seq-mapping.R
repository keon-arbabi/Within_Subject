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
new_meta = predictions[, c("predicted.id", "prediction.score.max")]
# add to metadata and set to active ident
seu_map = AddMetaData(seu_map, metadata = new_meta)
Idents(seu_map) = "predicted.id"
# save 
saveRDS(seu_map, file = "/external/rprshnas01/kcni/karbabi/R Projects/Within_Subject/keon_analyses/output/SME_Data_Integrated_Our_Labels.rds")

# compare 
#seu_map = readRDS("/external/rprshnas01/kcni/karbabi/R Projects/Within_Subject/keon_analyses/output/SME_Data_Integrated_Our_Labels.rds")
p1 = DimPlot(object = seu_map, reduction = "umap", label = TRUE, pt.size = 0.5) +
  theme(legend.position="none")

load("/external/rprshnas01/kcni/karbabi/R Projects/Within_Subject/Revision_Analysis/integration_scRNAseq/OUTPUTS/SME_Data_Integrated_Labels.RData")
p2 = DimPlot(object = SME_Data_Integrated, reduction = "umap", label = TRUE, pt.size = 0.5) +
  theme(legend.position="none")

p1 + p2

seu_ref = ScaleData(seu_ref, verbose = FALSE)
seu_ref  = RunPCA(seu_ref, npcs = 30, verbose = FALSE)
seu_ref  = RunUMAP(seu_ref, reduction = "pca", dims = 1:30, verbose = FALSE)
DimPlot(object = seu_ref, reduction = "umap", label = TRUE, pt.size = 0.5) +
  theme(legend.position="none")


seu_ref = RunUMAP(seu_ref, dims = 1:30, reduction = "pca", return.model = TRUE)
tanchors = FindTransferAnchors(reference = seu_ref, query = seu_map, dims = 1:30, reference.reduction = "pca") 
seu_map = MapQuery(anchorset = tanchors, reference = seu_ref, query = seu_map, 
                    refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap")


p1 <- DimPlot(seu_ref, reduction = "umap", group.by = "subclass_label_expanded_L35IT", label = TRUE, label.size = 3, 
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(seu_map, reduction = "ref.umap", group.by = "predicted.id", label = TRUE, 
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

















