library(Matrix)
library(data.table)
library(Seurat)
library(dplyr)
library(gt)
library(SPOTlight)
library(igraph)
library(RColorBrewer)

#single cell data
cortex_sc = Read10X("scRNAseq/filtered_feature_bc_matrix/")
cortex_sc = CreateSeuratObject(counts = cortex_sc)
cortex_sc

#spatial data
anterior = Load10X_Spatial("spatialData/", filename = "filtered_feature_bc_matrix.h5")
anterior

#extract count data from seurat object
counts = as.matrix(cortex_sc@assays$RNA@counts)

library(scMAGIC)
source("D:/SNU/scMAGIC.R")
source("D:/SNU/scMAGIC_atlas.R")

data = data("HCL_ref")
pdf("scMAGIC.pdf", width = 100, height = 100)

#function for celltype identification
  output.scMAGIC = scMAGIC_atlas(counts, HCL_ref, atlas = 'HCL',
                               type_ref = 'sum-counts', use_RUVseq = F,
                               min_cell = 5, num_threads = 8)
dev.off()
conf = output.scMAGIC$ConfidenceScore
celltype = output.scMAGIC$scMAGIC.tag
dfMAGIC = cbind.data.frame(colnames(counts), celltype)
dfMAGIC = cbind.data.frame(dfMAGIC, conf)
rownames(dfMAGIC) = dfMAGIC[,1]
colnames(dfMAGIC) = c("ID", "scMAGIC", "conf")

x = sapply(dfMAGIC[,2], function(x) strsplit(x, "_")[[1]][1], USE.NAMES=FALSE)
x = sapply(x, function(x) strsplit(x, "[(]")[[1]][1], USE.NAMES=FALSE)
dfMAGIC = cbind.data.frame(dfMAGIC, x)
colnames(dfMAGIC)[4] = "scMAGIC.l1"

freq = as.data.frame(table(dfMAGIC$scMAGIC.l1))
num = which(freq[,2] < 10)
replaceT = freq[num, 1]
replaceT = paste(replaceT, collapse = "|")
if(replaceT == ""){
  dfMAGIC$scMAGIC.l1 = gsub("Unknown", "Unassigned", dfMAGIC$scMAGIC.l1)
} else {
  dfMAGIC$scMAGIC.l1 = gsub(replaceT, "Unassigned", dfMAGIC$scMAGIC.l1)
  dfMAGIC$scMAGIC.l1 = gsub("Unknown", "Unassigned", dfMAGIC$scMAGIC.l1)
}

identical(rownames(dfMAGIC), rownames(cortex_sc@meta.data))
cortex_sc@meta.data = cbind.data.frame(cortex_sc@meta.data, dfMAGIC)

library(randomcoloR)
n = length(unique(cortex_sc@meta.data$scMAGIC.l1))
palette = distinctColorPalette(n)

#process single-cell data
cortex_sc = Seurat::SCTransform(cortex_sc, verbose = FALSE) %>%
  Seurat::RunPCA(., verbose = FALSE) %>%
  Seurat::RunUMAP(., dims = 1:30, verbose = FALSE)

DimPlot(cortex_sc, group.by = 'scMAGIC.l1', cols = palette, label = T)

# computing marker genes
Seurat::Idents(object = cortex_sc) <- cortex_sc@meta.data$scMAGIC.l1
cluster_markers_all <- Seurat::FindAllMarkers(object = cortex_sc, 
                                              assay = "SCT",
                                              slot = "data",
                                              verbose = TRUE, 
                                              only.pos = TRUE)

saveRDS(object = cluster_markers_all,
        file = here::here("markers_sc.RDS"))

#SPOTlight decomposition

spotlight_ls <- SPOTlight(
  x = cortex_sc,
  y = anterior@assays$Spatial@counts,
  mgs = cluster_markers_all,
  n_top = NULL,
  gene_id = "gene",
  group_id = "cluster",
  weight_id = "p_val_adj",
  min_prop = 0.01,
  verbose = TRUE,
)

saveRDS(object = spotlight_ls, file = here::here("spotlight_ls.rds"))

res = readRDS("spotlight_ls.rds")
mod <- res$NMF

plotTopicProfiles(
  x = mod,
  y = cortex_sc$scMAGIC.l1,
  facet = FALSE,
  min_prop = 0.01,
  ncol = 1) +
  theme(aspect.ratio = 1)

plotTopicProfiles(
  x = mod,
  y = cortex_sc$scMAGIC.l1,
  facet = TRUE,
  min_prop = 0.01,
  ncol = 6)

library(NMF)
sign <- basis(mod)
colnames(sign) <- paste0("Topic", seq_len(ncol(sign)))

mat <- res$mat
plotCorrelationMatrix(mat)

plotInteractions(mat, "heatmap")
plotInteractions(mat, "network")

ct <- colnames(mat)
mat[mat < 0.1] <- 0
# Define color palette
# (here we use 'paletteMartin' from the 'colorBlindness' package)
paletteMartin <- c(
  "#000000", "#004949", "#009292", "#ff6db6", "#ffb6db", 
  "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", 
  "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d")
pal <- colorRampPalette(paletteMartin)(length(ct))
names(pal) <- ct
plotSpatialScatterpie(
  x = anterior,
  y = mat,
  cell_types = colnames(y),
  img = FALSE,
  scatterpie_alpha = 1,
  pie_scale = 0.4) +
  scale_fill_manual(
    values = pal,
    breaks = names(pal))
