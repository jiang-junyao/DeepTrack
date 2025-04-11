library(spacexr)
sc = subset(brain_sub,celltype=='Spinal cord')
sc = RunUMAP(sc,reduction = 'scvi',dims = 1:10,min.dist = 1)
UMAPPlot(sc,group.by='celltype2')+scale_color_manual(values = jdb_palette('corona'))+
  theme_void()+theme(text = element_text(size=16))
ggsave('/data/jiangjunyao/polyATAC/polyatac_figure/single cell fig/subcluster lineage analysis/spinal cord AP.pdf',width = 10,height = 8)
###
scvi_latetn = brain_sub@reductions$scvi@cell.embeddings[,1:10]
scvi_average = list()
for (i in unique(brain_sub$celltype2)) {
  scvi_average[[i]] = colMeans(scvi_latetn[rownames(brain_sub@meta.data[brain_sub@meta.data$celltype2==i,]),])
}
latent_average = as.data.frame(scvi_average)
colnames(latent_average) = unique(brain_sub$celltype2)
lineage_simi = pheatmap(cor(latent_average),clustering_method = 'mcquitty')$tree_row
label_use <- data.frame(label =lineage_simi$labels )
label_use$celltype = label_use[,1]
rownames(label_use) = label_use[,1]
p1 = ggtree(lineage_simi,ladderize = 1,size=0.7)
p1%<+%label_use+geom_tiplab(offset=0.005)+xlim(NA,0.15)
rctd_decon <- function(ref,querry,nCore = 10){
  # extract information to pass to the RCTD Reference function
  counts <- ref@assays$RNA@counts
  cluster <- as.factor(ref$celltype2)
  names(cluster) <- colnames(ref)
  nUMI <- ref$nCount_RNA
  names(nUMI) <- colnames(ref)
  reference <- Reference(counts, cluster, nUMI)

  # set up query with the RCTD function SpatialRNA
  counts <- querry@assays$RNA@counts
  coords <- data.frame(querry$X,querry$Y)
  colnames(coords) <- c("x", "y")
  coords[is.na(colnames(coords))] <- NULL
  query <- SpatialRNA(coords, counts, colSums(counts))

  RCTD <- create.RCTD(query, reference, max_cores = nCore)
  RCTD <- run.RCTD(RCTD, doublet_mode = "multi")
  return(RCTD)
}
e95_10x = readRDS('/data/jiangjunyao/polyATAC/E9_10 visium/e95_sc_sub.rds')
rctd_result_full = rctd_decon(sc,e95_10x)
rctd_pred = rctd_result_full@results
all_weight = list()
cell_type = list()
first_type = c()
for (J in 1:length(rctd_pred)) {
  all_weight[[J]] = rctd_pred[[J]]$all_weights
  cell_type[[J]] = rctd_pred[[J]]$cell_type_list
  first_type = c(first_type,rctd_pred[[J]]$cell_type_list[1])
}
weight_df = as.data.frame(t(as.data.frame(all_weight)))
rownames(weight_df) = colnames(e95_10x)
e95_10x@meta.data = cbind(e95_10x@meta.data[,1:7],weight_df)
e95_10x$first_type = first_type
FeaturePlot(e95_10x,colnames(weight_df),reduction = 'spatial',pt.size=2.5,min.cutoff = 0.25)
