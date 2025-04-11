library(ArchR)
library(PCAtools)
library(umap)
proj_used <- readRDS("~/polyATAC/multiomi/merge_wolf_archr_subset/Save-ArchR-Project.rds")
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj_used,
  useMatrix = "PeakMatrix",
  groupBy = "celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
proj_used <- addMotifAnnotations(ArchRProj = proj_used, motifSet = "cisbp", name = "Motif")
proj_used <- addBgdPeaks(proj_used)
proj_used <- addDeviationsMatrix(
  ArchRProj = proj_used,
  peakAnnotation = "Motif",
  force = TRUE
)
chromvar_mt = getMatrixFromProject(proj_used,"MotifMatrix",threads = 10)
saveRDS(chromvar_mt,'/data/jiangjunyao/polyATAC/multiomi/merge_wolf/chromvar_dev.rds')
saveRDS(proj_used,"~/polyATAC/multiomi/merge_wolf_archr_subset/Save-ArchR-Project.rds")
obj@assays[['TFActivity_all_peak']] = CreateAssayObject(chromvar_mt@assays@data$z)
### pca
chromvar_pca = pca(t(chromvar_mt@assays@data$z))
plot(chromvar_pca$sdev[1:50])
chromvar_umap = umap(chromvar_pca$loadings[1:10])
chromvar_umap = chromvar_umap$layout
chromvar_umap = as.data.frame(chromvar_umap)
obj@reductions['chromvar_pca_umap'] = obj@reductions$umap
obj@reductions['chromvar_pca_umap']$chromvar_pca_umap@cell.embeddings[,1] = chromvar_umap[colnames(obj),1]
obj@reductions['chromvar_pca_umap']$chromvar_pca_umap@cell.embeddings[,2] = chromvar_umap[colnames(obj),2]
obj@reductions['chromvar_pca'] = obj@reductions$pca
for (i in 1:50) {
  obj@reductions['chromvar_pca']$chromvar_pca@cell.embeddings[,i] = chromvar_pca$loadings[colnames(obj),i]
}
saveRDS(obj,"~/polyATAC/polyatac_V2.rds")
### chormvar pc
list_pc = list()
for (i in unique(obj$ct_leiden)) {
  cell_use = rownames(obj@meta.data[obj$ct_leiden==i,])
  pc = colMeans(obj@reductions$chromvar_pca@cell.embeddings[cell_use,1:10])
  list_pc[[i]] = pc
}
df_pc = as.data.frame(list_pc)
pheatmap::pheatmap(cor(df_pc),border_color = NA)
