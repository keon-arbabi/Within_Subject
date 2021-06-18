library(tidyverse)
library(Seurat)

# set directory 
setwd("/")

# get reference seurat object 
seu_ref = readRDS("/external/rprshnas01/netdata_kcni/stlab/Public/Seurat_objects/Seu_AIBS_obj_update_07JUN21.rds")
#seu_ref = subset(seu_ref, region_label == "MTG")

# get seurat object you want to map to
load("/external/rprshnas01/kcni/karbabi/R Projects/Within_Subject/Revision_Analysis/integration_scRNAseq/OUTPUTS/SME_Data_Integrated.RData")
seu_query = SME_Data_Integrated
rm(SME_Data_Integrated)

# find variable features for transferring 
seu_ref = FindVariableFeatures(seu_ref, selection.method = "vst", nfeatures = 5000, verbose = TRUE)
seu_query = FindVariableFeatures(seu_query, selection.method = "vst", nfeatures = 5000, verbose = TRUE)

# standard work-flow
seu_ref = ScaleData(seu_ref, verbose = TRUE)
seu_ref = RunPCA(seu_ref, npcs = 30, verbose = TRUE)
seu_ref = RunUMAP(seu_ref, reduction = "pca", dims = 1:30, verbose = TRUE)
# view reference clusters 
DimPlot(object = seu_ref, reduction = "umap", label = TRUE, pt.size = 0.5) +
  theme(legend.position="none")

# transfer
tanchors = FindTransferAnchors(reference = seu_ref, 
                               query = seu_query, 
                               normalization.method = "LogNormalize", # Name of normalization method used: LogNormalize or SCT
                               recompute.residuals = TRUE, # If using SCT as a normalization method, compute query Pearson residuals using the reference SCT model parameters
                               reduction = "pcaproject", # Project the PCA from the reference onto the query. We recommend using PCA when reference and query datasets are from scRNA-seq
                               reference.reduction = NULL, # Name of dimensional reduction to use from the reference if running the pcaproject workflow (default is PCA)
                               project.query = FALSE, # Project the PCA from the query dataset onto the reference. Use only in rare cases where the query dataset has a much larger cell number, but the reference dataset has a unique assay for transfer
                               features = NULL, #Features to use for dimensional reduction. If not specified, set as variable features of the reference object which are also present in the query
                               scale = TRUE, # Scale query data
                               npcs = 30, # Number of PCs to compute on reference if reference.reduction is not provided
                               l2.norm = TRUE, # Perform L2 normalization on the cell embeddings after dimensional reduction
                               dims = 1:30, # Which dimensions to use from the reduction to specify the neighbor search space
                               k.anchor = 5, # How many neighbors (k) to use when finding anchors
                               k.filter = 200, # How many neighbors (k) to use when filtering anchors. Set to NA to turn off filtering
                               k.score = 30, # How many neighbors (k) to use when scoring anchors
                               max.features = 200, # The maximum number of features to use when specifying the neighborhood search space in the anchor filtering
                               nn.method = "annoy", # Method for nearest neighbor finding. Options include: rann, annoy
                               n.trees = 50, # More trees gives higher precision when using annoy approximate nearest neighbor search
                               eps = 0, # Error bound on the neighbor finding algorithm (from RANN or RcppAnnoy)
                               approx.pca = TRUE, # Use truncated singular value decomposition to approximate PCA
                               mapping.score.k = NULL, # Compute and store nearest k query neighbors in the AnchorSet object that is returned.
                               verbose = TRUE
                               )
                               
predictions = TransferData(anchorset = tanchors, 
                           refdata = seu_ref$subclass_label_expanded_L35IT, 
                           reference = NULL,
                           query = NULL,
                           weight.reduction = "pcaproject",
                           l2.norm = FALSE,
                           dims = NULL,
                           k.weight = 50,
                           sd.weight = 1,
                           eps = 0,
                           n.trees = 50,
                           verbose = TRUE,
                           slot = "data",
                           prediction.assay = FALSE,
                           store.weights = TRUE
                           )

# get relevant variables 
new_meta = predictions[, c("predicted.id", "prediction.score.max")]
# add to metadata and set to active ident
seu_query = AddMetaData(seu_query, metadata = new_meta)
Idents(seu_query) = "predicted.id"

# read in cell-type identifiers from Berto et al. 2021
load("/external/rprshnas01/kcni/karbabi/R Projects/Within_Subject/Revision_Analysis/integration_scRNAseq/OUTPUTS/SME_Data_Integrated_Labels.RData")

# confusion matrix
conf_table = round(prop.table(table(SME_Data_Integrated$Cell, seu_query$predicted.id),2)*100, 1)
conf_table[conf_table == 0] = "-"
conf_table
# save 
write.csv(conf_table, file = "/external/rprshnas01/kcni/karbabi/R Projects/Within_Subject/keon_analyses/output/confusion_table.csv")

# plot prediction.score.max 
plot(density(predictions$prediction.score.max),
     main = "prediction.score.max")
# filter by prediction.score.max
conf_cells = predictions %>%
  rownames_to_column(var = "id") %>%
  filter(predicted.id %in% c("L2/3 IT","L3/5 IT","VIP") & prediction.score.max < 0.6) %>%
  pull(id)
seu_query_filt = seu_query[,!colnames(seu_query) %in% conf_cells]

# compare
p1 = DimPlot(object = seu_query_filt, reduction = "umap", label = TRUE, pt.size = 0.5) 
p2 = DimPlot(object = SME_Data_Integrated, reduction = "umap", label = TRUE, pt.size = 0.5) 

plot_list = p1 + p2
ggsave(plot = plot_list, file = "/external/rprshnas01/kcni/karbabi/R Projects/Within_Subject/keon_analyses/output/mapping_comparison.png",
       device = "png", width = 18, height = 9, units = "in", dpi = 600)

# save final seurat objects
saveRDS(seu_query, file = "/external/rprshnas01/kcni/karbabi/R Projects/Within_Subject/keon_analyses/output/seu_query.rds")
saveRDS(seu_query_filt, file = "/external/rprshnas01/kcni/karbabi/R Projects/Within_Subject/keon_analyses/output/seu_query_filt.rds")


