library(tidyverse)
library(broom)
library(here)
library(cowplot)
library(ggpubr)
library(markerGeneProfile)
library(gplots)

#==== load data ====
# set directory 
setwd("/external/rprshnas01/kcni/karbabi/R Projects/Within_Subject")
# get some functions from Berto et al., 2021
source("Utils.R")

# Load adjusted data 
load("rawdata/expAdj.RData")
# format expression data
data_working = boot_data_all %>% 
  as.data.frame() %>%
  rownames_to_column(var = "Gene.Symbol")
# format meta data
data_meta = phenoAll %>%
  rownames_to_column(var = "sample_id")

# Load processed SME data
load("rawdata/SME_Values.RData")
data_sme = Lateral_resected %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var = "UT.code")

# Load in marker genes 
marker_data = read.csv("https://raw.githubusercontent.com/sonnyc247/MarkerSelection/master/Data/Outputs/CSVs_and_Tables/Markers/All_hodge_regions/new_ALLReg_results_ITexpand_WL35IT_lfct15_minpct25_dup.csv")[,-1]

# Manually add inhibitory/excitatory suffix to subclass labels 
marker_data %<>%
  mutate(class = case_when(
    subclass == "LAMP5" ~ "Inh",
    subclass == "L5 ET" ~ "Exc",
    subclass == "PAX6" ~ "Inh",
    subclass == "L5/6 IT Car3" ~ "Exc",
    subclass == "VIP" ~ "Inh",
    subclass == "L6 CT" ~ "Exc",
    subclass == "L6 IT"  ~ "Exc",
    subclass == "L6b" ~ "Exc",
    subclass == "PVALB" ~ "Inh",
    subclass == "L5/6 NP" ~ "Exc",
    subclass == "SST" ~ "Inh",
    subclass == "L4 IT" ~ "Exc",
    subclass == "L3/5 IT" ~ "Exc",
    subclass == "L2/3 IT" ~ "Exc")
    ) %>% 
  relocate(class, .before = "subclass") %>%
  unite(subclass, c(class, subclass), sep = "_", remove = F, na.rm = T)
marker_data$subclass = gsub(" ", "_", marker_data$subclass)
marker_data$class[is.na(marker_data$class)] = "NonN"
# fix labels 
marker_data$subclass = gsub("\\.", "/", marker_data$subclass)

paste("marker matches in data: ", length(intersect(unlist(data_working$Gene.Symbol), unlist(marker_data$gene))), "/",
      nrow(marker_data))

# Get vector of unique cell types 
(cell_types = marker_data$subclass %>% unique())
# Organize markers into a list 
marker_list = lapply(cell_types, function(cell_type){
  return(marker_data %>% filter(subclass == cell_type) %>% pull(gene) %>% unlist())
})
names(marker_list) = cell_types

#==== marker gene profile ====
# Run MGP analysis
estimations =  mgpEstimate(
  exprData = data_working,
  genes = marker_list,
  geneColName = 'Gene.Symbol',
  outlierSampleRemove = FALSE, # should outlier samples removed. This is done using boxplot stats
  geneTransform = NULL, # this is the default option for geneTransform
  groups = NULL, # if there are experimental groups provide them here. if not desired set to NULL
  seekConsensus = FALSE, # ensures gene rotations are positive in both of the groups
  removeMinority = TRUE)

# merge cell type proportions with sample metadata
mgp_estimates = as.data.frame(estimations$estimates) %>%
  rownames_to_column(var = "sample_id")
names(mgp_estimates) = gsub("\\.", "/", names(mgp_estimates))
mgp_df = data_meta %>%
  inner_join(data_sme, by = "UT.code") %>%
  inner_join(mgp_estimates, by = "sample_id") %>%
  pivot_longer(-c(colnames(data_meta), colnames(data_sme)),
               names_to = "cell_type",
               values_to = "cell_proportion") %>%
  pivot_longer(-c(colnames(data_meta), cell_type, cell_proportion),
               names_to = "oscillation",
               values_to = "oscillation_value")

# make objects for plotting
plot_celltypes = names(mgp_estimates)[-1]
plot_oscillation = unique(mgp_df$oscillation)
plot_list = list()

# scatter plots for oscillation vs. mgps 
for(i in 1:length(plot_celltypes)){
  plot_list[[i]] = ggplot(
    mgp_df %>% filter(cell_type == plot_celltypes[i]),
    aes(x = oscillation_value, y = cell_proportion)) +
    geom_smooth(method = "lm", se = F) + 
    geom_point() +
    ylab(paste(plot_celltypes[i], " MGP")) + xlab ("Oscillation") +
    theme_bw() +
    facet_wrap(~oscillation, scales = "fixed")
}

#pdf("keon_analyses/output/scatterplots-subclass.pdf",width=25,height=15)
ggarrange(plotlist = plot_list)
#dev.off()            

#==== correlation analysis ====
# reformat dataframe
mgp_estimates_sme = mgp_estimates %>%
  filter(sample_id %in% rownames(pheno.ws)) %>%
  column_to_rownames(var = "sample_id") %>%
  t() 
# quick spearman correlation  
res = list()
for(i in 1:nrow(Lateral_resected)) {
  res[[i]] = FastCor(mgp_estimates_sme,Lateral_resected[i,],method="spearman",alternative="two.sided",cores=1,override=TRUE)%>%
    as.data.frame() %>%
    rownames_to_column('Gene') 
}
# add column for wave 
SME_cor_LR = data.frame()
for(i in 1:length(res)){
  res[[i]]$Waves = rep(paste(rownames(Lateral_resected)[i]),nrow(res[[i]]))
}
# combine list into dataframe 
SME_cor_LR = do.call(rbind,res)
SME_cor_LR = SME_cor_LR[c(1,4,2,3)] %>%
  arrange(Rho)

# format for heatmap
SME_cor_LR_wide = SME_cor_LR %>%
  select(-Pval) %>%
  pivot_wider(names_from = Waves, values_from = Rho) %>%
  column_to_rownames(var = "Gene") %>%
  mutate_at(plot_oscillation, ~as.numeric(as.character(.))) %>%
  as.matrix()

#pdf("keon_analyses/output/correlation-heatmap-subclass.pdf",width=10,height=8)
heatmap.2(SME_cor_LR_wide, scale = "none",
          col = bluered(100), trace = "none", density.info = "none",
          margins = c(10,10),
          key.xlab = "Rho")
#dev.off()    

#==== linear models ====
plot_list = list()
for(i in 1:length(plot_oscillation)){
  
  lm_df = mgp_df %>%
    filter(oscillation == plot_oscillation[i]) %>%
    group_by(cell_type) %>%
    do(print(tidy(lm(scale(cell_proportion) ~ scale(oscillation_value),  data = .)))) %>%
    ungroup() %>%
    mutate(padj = p.adjust(`p.value`, method = 'BH')) %>%
    mutate(term = recode(term, 
                         `scale(oscillation_value)` = "oscillation")) %>%
    mutate(class = case_when(
      str_detect(cell_type, "Inh") ~ "Inhibitory",
      str_detect(cell_type, "Exc") ~ "Excitatory",
      TRUE ~ "Non-Neuronal"
    ))
  
  # Beta coeffs per cell type for phenotype and age effects
  beta_plot = lm_df %>% 
    filter(term %in% 'oscillation') %>% 
    mutate(cell_type = fct_reorder(cell_type, estimate))
  
  plot_list[[i]] = ggplot(beta_plot, aes(x = cell_type, y = estimate)) + 
    geom_hline(yintercept = 0) + 
    geom_bar(stat = "identity") + 
    geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error)) + 
    ylab('') + 
    xlab('') + 
    ggtitle(paste(plot_oscillation[i])) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_wrap(~class, drop = T, scales = "free_x")
}

#pdf("keon_analyses/output/linear-models-subclass.pdf",width=15,height=9)
p = ggarrange(plotlist = plot_list)
annotate_figure(p, left = "Std. Beta coeff.", bottom = "Cell type proportions")
#dev.off()            

#==== bulk vs. sc ====

# load-in processed seurat object 
seu_obj = readRDS(file = "keon_analyses/output/seu_query_filt.rds")
# get cell-type proportions by counting cells of each type within subjects 
ctp_df = data.frame(sample_id = seu_obj@meta.data$orig.ident,
                    cell_type = seu_obj@meta.data$predicted.id) %>%
  group_by(sample_id) %>%
  count(cell_type) %>%
  na.omit() %>%
  mutate(sum = sum(n)) %>%
  ungroup() %>%
  mutate(percent = n/sum*100) 
# add class suffix
ctp_df = ctp_df %>%
  mutate(class = case_when(
    cell_type == "LAMP5" ~ "Inh",
    cell_type == "L5 ET" ~ "Exc",
    cell_type == "PAX6" ~ "Inh",
    cell_type == "L5/6 IT Car3" ~ "Exc",
    cell_type == "VIP" ~ "Inh",
    cell_type == "L6 CT" ~ "Exc",
    cell_type == "L6 IT"  ~ "Exc",
    cell_type == "L6b" ~ "Exc",
    cell_type == "PVALB" ~ "Inh",
    cell_type == "L5/6 NP" ~ "Exc",
    cell_type == "SST" ~ "Inh",
    cell_type == "L4 IT" ~ "Exc",
    cell_type == "L3/5 IT" ~ "Exc",
  cell_type == "L2/3 IT" ~ "Exc")
  ) %>% 
  relocate(class, .before = "cell_type") %>%
  unite(cell_type, c(class, cell_type), sep = "_", remove = T, na.rm = T)
ctp_df$cell_type = gsub(" ", "_", ctp_df$cell_type)

# reformat mgp dataframe 
mgp_estimates = as.data.frame(estimations$estimates) %>%
  rownames_to_column(var = "sample_id")
names(mgp_estimates) = gsub("\\.", "/", names(mgp_estimates))

mgp_df = data_meta %>%
  inner_join(mgp_estimates, by = "sample_id") %>%
  pivot_longer(-c(colnames(data_meta)),
               names_to = "cell_type",
               values_to = "cell_proportion")
mgp_df$sample_id = gsub("_","", mgp_df$sample_id)

# bring together mgp and sc-ctps 
combined_df = inner_join(mgp_df, ctp_df, by = c("sample_id","cell_type"))

pdf("keon_analyses/output/bulk-vs-sc-comparison.pdf",width=15,height=9)
  ggplot(combined_df, aes(x = percent, y = cell_proportion)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ylab(paste(plot_celltypes[i], "bulk-MGP")) + xlab ("sc-CTP") +
  stat_cor(method = "spearman") +
  theme_bw() +
  facet_wrap(~cell_type, scale = "free")
dev.off() 

