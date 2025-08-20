library(scales)
library(ArchR)
library(Hmisc)
library(Seurat)
library(diffusionMap)
library(destiny)
library(BuenColors)
library(pheatmap)
library(PseudotimeDE)
library(corrplot)
library(org.Mm.eg.db)
library(monocle)
library(factoextra)
source('/data/jiangjunyao/polyATAC/script/smooth_by_bin.R')
gene_enrich <- function(gene,org.db,enrich.db,use_internal_data = TRUE,
                        organism = NULL,pvalueCutoff = 0.05){
  library(clusterProfiler)
  gene1 <- clusterProfiler::bitr(gene, fromType = "SYMBOL",
                                 toType = c("SYMBOL", "ENTREZID"),
                                 OrgDb = org.db)
  if (enrich.db =='KEGG') {
    k1 <- clusterProfiler::enrichKEGG(gene = gene1$ENTREZID,
                                      pvalueCutoff = pvalueCutoff
                                      ,use_internal_data = use_internal_data,
                                      minGSSize = 2)
  }else if(enrich.db =='GO'){
    k1 = clusterProfiler::enrichGO(gene = gene1$ENTREZID,
                                   OrgDb = org.db,
                                   keyType = "ENTREZID",
                                   ont = "BP",
                                   pvalueCutoff = pvalueCutoff,minGSSize = 2)
  }
  return(k1)
}

extract_expressed_TFs <- function(seurat_object,TFs,cells_quantile = 0.05){
  TFs <- TFs[TFs%in%rownames(seurat_object)]
  matrix_tf <- seurat_object@assays$RNA@counts[TFs,]
  if (cells_quantile==0) {
    TfExp <- matrix_tf[rowSums(as.matrix(matrix_tf))>0,]
  }else{
    quantile_exp <- ncol(as.matrix(matrix_tf))/(1/cells_quantile)
    TfExp <- matrix_tf[ncol(as.matrix(matrix_tf))-rowSums(as.matrix(matrix_tf==0))>quantile_exp,]}
  return(TfExp)
}
diffgenetest_pseudotime <- function(monocle_object){
  mo <- detectGenes(monocle_object, min_expr = 1)
  #mo <- estimateSizeFactors(mo)
  diff1 <- monocle::differentialGeneTest(mo,
                                         fullModelFormulaStr = "~Pseudotime",
                                         relative_expr = TRUE,cores=1
  )
  ed <- c()
  for (i in diff1$num_cells_expressed) {
    a <- log(i) * 0.95 - log(i) * 0.05
    ed <- c(ed, a)
  }
  diff1$expression_difference <- ed
  return(diff1)
}
sparse.cor <- function(x){
  n <- nrow(x)
  cMeans <- colMeans(x)
  cSums <- colSums(x)
  # Calculate the population covariance matrix.
  # There's no need to divide by (n-1) as the std. dev is also calculated the same way.
  # The code is optimized to minize use of memory and expensive operations
  covmat <- tcrossprod(cMeans, (-2*cSums+n*cMeans))
  crossp <- as.matrix(crossprod(x))
  covmat <- covmat+crossp
  sdvec <- sqrt(diag(covmat)) # standard deviations of columns
  covmat/crossprod(t(sdvec)) # correlation matrix
}
min_max_scale <- function(predicted_values){
  v1=(predicted_values - min(predicted_values)) / (max(predicted_values) - min(predicted_values))
  return(v1)
}
get_pseudotime_de <- function(data, reverse = FALSE,gene.use = NULL) {
  pd <- data.frame(Pseudotime=1:ncol(data),sample='meta',row.names = colnames(data))
  pd <- new("AnnotatedDataFrame", data = pd)
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new("AnnotatedDataFrame", data = fData)
  data1=as(as.matrix(data), "sparseMatrix")
  monocle_cds <- monocle::newCellDataSet(data1,
                                         phenoData = pd,
                                         featureData = fd,
                                         lowerDetectionLimit = 0.1,
                                         expressionFamily = VGAM::negbinomial.size()
  )
  cds <- BiocGenerics::estimateSizeFactors(monocle_cds)
  deg=diffgenetest_pseudotime(cds)
  return(deg)
}
get_pseudotime_de_seurat <- function(data, reverse = FALSE,gene.use = NULL) {
  pd <- data.frame(Pseudotime=data$scale_pseudotime,
                   sample='meta',row.names = colnames(data))
  pd <- new("AnnotatedDataFrame", data = pd)
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new("AnnotatedDataFrame", data = fData)
  #data1=as(data@assays$RNA@counts, "sparseMatrix")
  monocle_cds <- monocle::newCellDataSet(data@assays$RNA@data,
                                         phenoData = pd,
                                         featureData = fd,
                                         lowerDetectionLimit = 0.1,
                                         expressionFamily = VGAM::negbinomial.size()
  )
  cds <- BiocGenerics::estimateSizeFactors(monocle_cds)
  deg=diffgenetest_pseudotime(cds)
  return(deg)
}
gene_cor = function(candidate_gene,obj,col1){
  library(foreach)
  library(doParallel)
  registerDoParallel(cores=10)
  pval = c()
  corv = c()
  results <- foreach(i = candidate_gene, .combine=rbind) %dopar% {
    cor_test <- cor.test(obj@assays$RNA@data[i,], obj@meta.data[,col1])
    cbind(gene=i, pval=cor_test$p.value, corv=cor_test$estimate)
  }
  stopImplicitCluster()
  cor_df <- as.data.frame(results)
  colnames(cor_df) <- c("candidate_gene", "pval", "corv")
  return(cor_df)
}
gene_cor_df = function(mt){
  library(foreach)
  library(doParallel)
  registerDoParallel(cores=10)
  pval = c()
  corv = c()
  results <- foreach(i = rownames(mt), .combine=rbind) %dopar% {
    cor_test <- cor.test(mt[i,], 1:ncol(mt))
    cbind(gene=i, pval=cor_test$p.value, corv=cor_test$estimate)
  }
  stopImplicitCluster()
  cor_df <- as.data.frame(results)
  colnames(cor_df) <- c("candidate_gene", "pval", "corv")
  cor_df[,3] = as.numeric(cor_df[,3])
  cor_df[,2] = as.numeric(cor_df[,2])
  return(cor_df)
}
Color1 <- c(rgb(72 / 255, 85 / 255, 167 / 255), rgb(255 / 255, 255 / 255, 255 / 255),
            rgb(239 / 255, 58 / 255, 37 / 255))
fate_col = c('#4DB6AC',
             '#0D47A1','#1B5E20','#9ECAE1','#C5E1A5')
names(fate_col) = c('Balance','Paraxial mesoderm','Spinal cord',"Paraxial_mesoderm_bias","Spinal_cord_bias")
ct_col = c('#0D47A1','#1B5E20','#4DB6AC')
names(ct_col) = c('Paraxial mesoderm','Spinal cord','Neuromesodermal progenitors')

obj = readRDS('/data/jiangjunyao/polyATAC/nmp_final/nmp_sp_pm_only_rna.rds')
### cal diffusion map
dm_dim = DiffusionMap(obj@reductions$scvi@cell.embeddings,distance='cosine',k=50,n_eigs = 5,rotate = T)
dm_dim2 = as.data.frame(dm_dim@eigenvectors[,c(1,2)])
colnames(dm_dim2) = c('embed1','embed2')
dm_dim2$celltype = obj$celltype
ggplot(dm_dim2,aes(x=embed1,y=embed2,color=celltype))+geom_point(size=0.8)+theme_classic()
###
obj@reductions$diffmap = obj@reductions$umap
obj@reductions$diffmap@cell.embeddings[,1] = ddetree_dim2[,1]
obj@reductions$diffmap@cell.embeddings[,2] = ddetree_dim2[,2]
pm_pseudotime <- read.csv("~/polyATAC/multiomi/nmp_fate_pseudotime/pm_pseudotime.csv", row.names=1)
sp_pseudotime <- read.csv("~/polyATAC/multiomi/nmp_fate_pseudotime/sp_pseudotime.csv", row.names=1)
pm_multi_cell = intersect(rownames(pm_pseudotime),colnames(obj))
sp_multi_cell = colnames(obj)[!colnames(obj) %in% pm_multi_cell]
sp_obj = subset(obj,cells=sp_multi_cell)
pm_obj = subset(obj,cells=pm_multi_cell)
sp_obj@meta.data[sp_multi_cell,'Pseudotime'] = sp_pseudotime[sp_multi_cell,'pseudotime']
pm_obj@meta.data[pm_multi_cell,'Pseudotime'] = pm_pseudotime[pm_multi_cell,'pseudotime']
sp_obj@meta.data[sp_multi_cell,'sp_score'] = sp_pseudotime[sp_multi_cell,'fate_map_transition_map_Spinal.cord']
pm_obj@meta.data[pm_multi_cell,'sp_score'] = pm_pseudotime[pm_multi_cell,'fate_map_transition_map_Spinal.cord']
sp_obj@meta.data[sp_multi_cell,'pm_score'] = sp_pseudotime[sp_multi_cell,'fate_map_transition_map_Paraxial.mesoderm']
pm_obj@meta.data[pm_multi_cell,'pm_score'] = pm_pseudotime[pm_multi_cell,'fate_map_transition_map_Paraxial.mesoderm']
pm_obj$Pseudotime = -pm_obj$Pseudotime
sp_obj$Pseudotime = -sp_obj$Pseudotime
pm_obj$scale_pseudotime = min_max_scale(pm_obj$Pseudotime)
sp_obj$scale_pseudotime = min_max_scale(sp_obj$Pseudotime)
FeaturePlot(pm_obj,'scale_pseudotime')+scale_color_gradientn(colors = jdb_palette("flame_light"))
FeaturePlot(sp_obj,'scale_pseudotime')+scale_color_gradientn(colors = jdb_palette("flame_light"))
UMAPPlot(pm_obj,group.by='celltype')+scale_color_manual(values = fate_col)
UMAPPlot(sp_obj,group.by='celltype')+scale_color_manual(values = fate_col)
obj$Pseudotime = obj$scale_pseudotime
obj@meta.data[colnames(pm_obj),'scale_pseudotime'] = pm_obj$scale_pseudotime
obj@meta.data[colnames(sp_obj),'scale_pseudotime'] = sp_obj$scale_pseudotime
###
obj$celltype_ori = obj$celltype
obj$sp_score = NA
obj@meta.data[obj$celltype %in% c('Balance',
                                  'Paraxial_mesoderm_bias',
                                  'Spinal_cord_bias'),'celltype_ori'] = 'Neuromesodermal progenitors'
obj@meta.data[colnames(sp_obj),'sp_score'] = sp_obj$sp_score
obj@meta.data[colnames(pm_obj),'sp_score'] = pm_obj$sp_score
obj@meta.data[colnames(sp_obj),'pm_score'] = sp_obj$pm_score
obj@meta.data[colnames(pm_obj),'pm_score'] = pm_obj$pm_score
saveRDS(obj,'/data/jiangjunyao/polyATAC/nmp_final/nmp_sp_pm_only_rna_2025_0709.rds')
### plot umap
UMAPPlot(obj,group.by='celltype')+scale_color_manual(values = fate_col)
ggsave('/data/jiangjunyao/polyATAC/figure_after_science/single cell fig/NMP_fate_analysis/embedding/NMP umap embedding celltype.pdf',
       width = 10,height = 8)
UMAPPlot(obj,group.by='celltype_ori')+scale_color_manual(values = ct_col)
ggsave('/data/jiangjunyao/polyATAC/figure_after_science/single cell fig/NMP_fate_analysis/embedding/NMP umap embedding celltype_original.pdf',
       width = 10,height = 8)
FeaturePlot(obj,'scale_pseudotime')+
  scale_color_gradientn(colors = jdb_palette("flame_light"))
ggsave('/data/jiangjunyao/polyATAC/figure_after_science/single cell fig/NMP_fate_analysis/embedding/NMP umap embedding pseudotime.pdf',
       width = 10,height = 8)
UMAPPlot(obj,group.by='time')
ggsave('/data/jiangjunyao/polyATAC/figure_after_science/single cell fig/NMP_fate_analysis/embedding/NMP umap embedding time.pdf',
       width = 10,height = 8)
FeaturePlot(obj,'sp_score')+scale_colour_viridis_c()+theme_void()
ggsave('/data/jiangjunyao/polyATAC/figure_after_science/single cell fig/NMP_fate_analysis/embedding/sp_bias_score.pdf',width = 10,height = 8)
FeaturePlot(obj,'pm_score')+scale_colour_viridis_c()+theme_void()
ggsave('/data/jiangjunyao/polyATAC/figure_after_science/single cell fig/NMP_fate_analysis/embedding/pm_bias_score.pdf',width = 10,height = 8)
table(obj$celltype)
obj_only_nmp = subset(obj, celltype %in% c('Balance',
                                          'Paraxial_mesoderm_bias',
                                          'Spinal_cord_bias'))
FeaturePlot(obj_only_nmp,'sp_score')+scale_color_gradientn(colors = jdb_palette("brewer_green"))+theme_void()
ggsave('/data/jiangjunyao/polyATAC/figure_after_science/single cell fig/NMP_fate_analysis/embedding/sp_bias_score_only_nmp.pdf',width = 10,height = 8)
FeaturePlot(obj_only_nmp,'pm_score')+scale_color_gradientn(colors = jdb_palette("brewer_purple"))+theme_void()
ggsave('/data/jiangjunyao/polyATAC/figure_after_science/single cell fig/NMP_fate_analysis/embedding/pm_bias_score_only_nmp.pdf',width = 10,height = 8)
### plot diffusion map
DimPlot(obj,group.by='celltype',reduction = 'diffmap')+scale_color_manual(values = fate_col)+NoLegend()
ggsave('/data/jiangjunyao/polyATAC/figure_after_science/single cell fig/NMP_fate_analysis/embedding/NMP diffusion map embedding celltype.pdf',
       width = 10,height = 8)
DimPlot(obj,group.by='celltype_ori',reduction = 'diffmap')+scale_color_manual(values = ct_col)
ggsave('/data/jiangjunyao/polyATAC/figure_after_science/single cell fig/NMP_fate_analysis/embedding/NMP diffusion map embedding celltype_original.pdf',
       width = 10,height = 8)
DimPlot(obj,group.by='time',reduction = 'diffmap')+
  scale_color_manual(values = jdb_palette('brewer_spectra'))
ggsave('/data/jiangjunyao/polyATAC/figure_after_science/single cell fig/NMP_fate_analysis/embedding/NMP diffusion map embedding real time.pdf',
       width = 10,height = 8)


plot_df = cbind(obj@reductions$diffmap@cell.embeddings,obj$scale_pseudotime,obj$sp_score,obj$pm_score)
colnames(plot_df)[3:5] = c('Pseudotime','sp_score','pm_score')
ggplot(plot_df,aes(x=DiffusionMap_1,y=DiffusionMap_2,color=Pseudotime))+
  geom_point(size=0.9)+scale_color_gradientn(colors = rev(jdb_palette("brewer_marine")))+
  theme_classic()+theme(text = element_text(size=16))
ggsave('/data/jiangjunyao/polyATAC/figure_after_science/single cell fig/NMP_fate_analysis/embedding/NMP diffusion map embedding pseudotime.pdf',
       width = 10,height = 8)
ggplot(plot_df,aes(x=DiffusionMap_1,y=DiffusionMap_2,color=sp_score))+
  geom_point(size=0.9)+scale_color_viridis_c()+
  theme_classic()+theme(text = element_text(size=16))
ggsave('/data/jiangjunyao/polyATAC/figure_after_science/single cell fig/NMP_fate_analysis/embedding/NMP diffusion map embedding spinal cord bias score.pdf',
       width = 10,height = 8)
ggplot(plot_df,aes(x=DiffusionMap_1,y=DiffusionMap_2,color=pm_score))+
  geom_point(size=0.9)+scale_color_viridis_c()+
  theme_classic()+theme(text = element_text(size=16))
ggsave('/data/jiangjunyao/polyATAC/figure_after_science/single cell fig/NMP_fate_analysis/embedding/NMP diffusion map embedding paraxial mesoderm bias score.pdf',
       width = 10,height = 8)
### no legend
###
ggplot(plot_df,aes(x=DiffusionMap_1,y=DiffusionMap_2,color=Pseudotime))+
  geom_point(size=0.9)+scale_color_gradientn(colors = rev(jdb_palette("brewer_marine")))+
  theme_classic()+theme(text = element_text(size=16))+guides(color = "none")
ggsave('/data/jiangjunyao/polyATAC/figure_after_science/single cell fig/NMP_fate_analysis/embedding/NMP diffusion map embedding pseudotime_nolegend.pdf',
       width = 10,height = 8)
ggplot(plot_df,aes(x=DiffusionMap_1,y=DiffusionMap_2,color=sp_score))+
  geom_point(size=0.9)+scale_color_viridis_c()+
  theme_classic()+theme(text = element_text(size=16))+guides(color = "none")
ggsave('/data/jiangjunyao/polyATAC/figure_after_science/single cell fig/NMP_fate_analysis/embedding/NMP diffusion map embedding spinal cord bias score_nolegend.pdf',
       width = 10,height = 8)
ggplot(plot_df,aes(x=DiffusionMap_1,y=DiffusionMap_2,color=pm_score))+
  geom_point(size=0.9)+scale_color_viridis_c()+
  theme_classic()+theme(text = element_text(size=16))+guides(color = "none")
ggsave('/data/jiangjunyao/polyATAC/figure_after_science/single cell fig/NMP_fate_analysis/embedding/NMP diffusion map embedding paraxial mesoderm bias score_nolegend.pdf',
       width = 10,height = 8)
#########################
### gene module analysis
#########################
pm_gene = rownames(extract_expressed_TFs(pm_obj,rownames(pm_obj),0.01))
sp_gene = rownames(extract_expressed_TFs(sp_obj,rownames(sp_obj),0.01))

pm_sp_gene = unique(c(pm_gene,sp_gene))
pm_obj$celltype2 = pm_obj$celltype
pm_obj$celltype2[pm_obj$celltype2!='Paraxial mesoderm'] = 'NMP'
pm_smooth = smooth_pseudotime(as.matrix(pm_obj@assays$RNA@data),pm_obj@meta.data,
                              pm_sp_gene,100,NULL,c('NMP','Paraxial mesoderm'))
colnames(pm_smooth[[1]]) = paste0('metacell',1:ncol(pm_smooth[[1]]))

pm_de = get_pseudotime_de(pm_smooth[[1]])
pm_de = pm_de[pm_de$num_cells_expressed>2,]
pm_de = pm_de[pm_de$pval<=0.05,]
pm_de = pm_de[pm_de$gene_short_name %in% names(pm_smooth[[3]])[pm_smooth[[3]]>0.5],]


sp_obj$celltype2 = sp_obj$celltype
sp_obj$celltype2[sp_obj$celltype2!='Spinal cord'] = 'NMP'
sp_smooth = smooth_pseudotime(as.matrix(sp_obj@assays$RNA@data),sp_obj@meta.data,
                              pm_sp_gene,100,NULL,c('NMP','Spinal cord'))
colnames(sp_smooth[[1]]) = paste0('metacell',1:ncol(sp_smooth[[1]]))

sp_de = get_pseudotime_de(sp_smooth[[1]])
sp_de = sp_de[sp_de$num_cells_expressed>2,]
sp_de = sp_de[sp_de$pval<=0.05,]
sp_de = sp_de[sp_de$gene_short_name %in% names(sp_smooth[[3]])[sp_smooth[[3]]>0.5],]


pheatmap(sp_smooth[[2]][sp_de$gene_short_name,],cluster_cols = F,border_color = NA,show_colnames = F)
pheatmap(pm_smooth[[2]][pm_de$gene_short_name,],cluster_cols = F,border_color = NA,show_colnames = F)

final_fea = unique(c(sp_de$gene_short_name,pm_de$gene_short_name))
final_fea = final_fea[final_fea %in% rownames(pm_smooth[[1]])]
final_fea = final_fea[final_fea %in% rownames(sp_smooth[[1]])]

dual_lineage = cbind(pm_smooth[[2]][final_fea,100:1],
                     sp_smooth[[2]][final_fea,])
dual_lineage2 = cbind(pm_smooth[[1]][final_fea,100:1],
                     sp_smooth[[1]][final_fea,])
#dual_lineage = scale_mt(dual_lineage)
pheatmap(dual_lineage,cluster_cols = F,border_color = NA,show_colnames = F)
colnames(dual_lineage) = 1:200
saveRDS(dual_lineage,'/data/jiangjunyao/polyATAC/nmp_final/nmp_pseudotime_module/dual_lineage_mt_0709.rds')
write.csv(pm_smooth,'/data/jiangjunyao/polyATAC/nmp_final/nmp_pseudotime_module/only_rna/pm_smooth_0709.csv')
write.csv(sp_smooth,'/data/jiangjunyao/polyATAC/nmp_final/nmp_pseudotime_module/only_rna/sp_smooth_0709.csv')
write.csv(sp_de,'/data/jiangjunyao/polyATAC/nmp_final/nmp_pseudotime_module/only_rna/sp_deg_0709.csv')
write.csv(pm_de,'/data/jiangjunyao/polyATAC/nmp_final/nmp_pseudotime_module/only_rna/pm_deg_0709.csv')

### modulize gene
fviz_nbclust(dual_lineage, kmeans, method = "wss",k.max = 15)
group = kmeans(dual_lineage,centers = 8)$cluster
names(group) = rownames(dual_lineage)
group = sort(group,decreasing = T)
pheatmap(dual_lineage[names(group),],
         cluster_cols = F,border_color = NA,
         show_colnames = F,cluster_rows = F,gaps_col = 100)

re_cluster_group = group[group%in%c(3,2,1)]
group2 = kmeans(dual_lineage[names(re_cluster_group),],centers = 5)$cluster
names(group2) = names(re_cluster_group)
group2 = sort(group2,decreasing = T)
pheatmap(dual_lineage[names(group2),],
         cluster_cols = F,border_color = NA,
         show_colnames = F,cluster_rows = F,gaps_col = 100)


group_df = data.frame(group)
group_df$module = '114514'
group_df$module[group_df[,1]==8] = 'early mesoderm regulate'
group_df$module[group_df[,1]==7] = 'late neuron regulate'
group_df$module[group_df[,1]==6] = 'late mesoderm regulate'
group_df$module[group_df[,1]==5] = 'late neuron regulate'
group_df$module[group_df[,1]==4] = 'late mesoderm regulate'
group_df[names(group2[group2==5]),2] = 'early neuron regulate'
group_df[names(group2[group2==4]),2] = 'early mesoderm regulate'
group_df[names(group2[group2==3]),2] = 'Balance regulate'
group_df[names(group2[group2==2]),2] = 'early neuron regulate'
group_df[names(group2[group2==1]),2] = 'early mesoderm regulate'

group_df = group_df[order(group_df$module),]
table(group_df$module)
write.csv(group_df,'/data/jiangjunyao/polyATAC/nmp_final/nmp_pseudotime_module/only_rna/group_df_20241115.csv')

### do heatmap
dual_lineage = readRDS('/data/jiangjunyao/polyATAC/nmp_final/nmp_pseudotime_module/dual_lineage_mt_0709.rds')
group_df_20241115 = read.csv('/data/jiangjunyao/polyATAC/nmp_final/nmp_pseudotime_module/only_rna/group_df_20250709.csv')
legend_colo = list(module=c('Balance regulate'='#78cdd1',
                              'late mesoderm regulate'='#33a3dc',
                              'late neuron regulate'='#77ac98',
                              'early mesoderm regulate'='#90d7ec',
                              "early neuron regulate"='#cde6c7'),
                   celltype=c('Balance NMP'='#4DB6AC',
                              'Paraxial mesoderm biased NMP'='#9ECAE1',
                              'Spinal cord biased NMP'='#C5E1A5',
                              'Paraxial mesoderm'='#0D47A1',
                              'Spinal cord'='#1B5E20'))
col <- colorRampPalette(c('#2166AC', '#4393C3','#92C5DE', '#D1E5F0', '#F7F7F7', '#FDDBC7', '#F4A582', '#D6604D', '#B2182B'))(100)
anno_row = data.frame(group_df_20241115$module)
rownames(anno_row) = rownames(dual_lineage)
colnames(anno_row) = 'module'
anno_col = data.frame(c(rep('Paraxial mesoderm',67),
                              rep('Paraxial mesoderm biased NMP',23),
                              rep('Balance NMP',16),
                              rep('Spinal cord biased NMP',7),
                              rep('Spinal cord',87)),row.names = 1:200)
colnames(anno_col) = 'celltype'
pdf('/data/jiangjunyao/polyATAC/figure_after_science/single cell fig/NMP_fate_analysis/gene_module_heatmap2.pdf',width = 10,height = 8)
colnames(dual_lineage) = 1:200
pheatmap(dual_lineage[group_df_20241115$X,],
         show_rownames = F,cluster_cols = F,cluster_rows = F,gaps_col = 100,
         show_colnames = F,annotation_row = anno_row,annotation_col = anno_col,
         color = col,annotation_colors = legend_colo,gaps_row = c(20,67,87,136))
dev.off()
sort_genes_by_expression_pattern <- function(dataframe) {

  gene_center <- apply(dataframe, 1, function(x) {
    sum(x * seq_along(x)) / sum(x)  # 加权平均位置

  })
  sorted_genes <- names(sort(gene_center))

  dataframe_sorted <- dataframe[sorted_genes, ]

  return(dataframe_sorted)
}
# list1 = list()
# for (i in unique(group_df_20241115$module)) {
#   print(i)
#   group_use = group_df_20241115[group_df_20241115$module==i,]
#   list1[[i]] = as.data.frame(sort_genes_by_expression_pattern(dual_lineage[group_use$X,]))
# }
lm = group_use = group_df_20241115[group_df_20241115$module=="late mesoderm regulate",]
ln = group_use = group_df_20241115[group_df_20241115$module=="late neuron regulate",]
em = group_use = group_df_20241115[group_df_20241115$module=="early mesoderm regulate",]
en = group_use = group_df_20241115[group_df_20241115$module=="early neuron regulate",]
em = as.data.frame(sort_genes_by_expression_pattern(dual_lineage[em$X,1:100]))
en = as.data.frame(sort_genes_by_expression_pattern(dual_lineage[en$X,200:101]))
lm = as.data.frame(sort_genes_by_expression_pattern(dual_lineage[lm$X,1:100]))
ln = as.data.frame(sort_genes_by_expression_pattern(dual_lineage[ln$X,200:101]))

final_idx = c(group_df_20241115[group_df_20241115$module=="Balance regulate",2],
              rownames(em),rownames(en),rownames(lm),rownames(ln))
#rownames(group_df_20241115)=group_df_20241115$X
#group_df_20241115 = group_df_20241115[rownames(final_df),]
#rownames(group_df_20241115) == rownames(final_df)
anno_row = data.frame(group_df_20241115$module)
rownames(anno_row) = final_idx
colnames(anno_row) = 'module'
pheatmap(dual_lineage[final_idx,],
         show_rownames = F,cluster_cols = F,cluster_rows = F,gaps_col = 100,
         show_colnames = F,annotation_row = anno_row,annotation_col = anno_col,
         color = col,annotation_colors = legend_colo,gaps_row = c(20,66,87,133))


### do line plot
gene_plot = c('Cdx2','Cdx4','Wnt5a','Wnt3a','Hoxc5','Hoxa9','Hoxa7','Hoxb9','Hoxb8',
              'Hoxc4','Hoxc9','Hes1','Sema6a','Zeb2','Crabp2','Ptn','Efna5','Nin',
              'Id2','Sema3a','Sema3e','Ptk7','Il17rd','Fn1','Cyp26a1','Hoxa3','Ifitm1','T','Fgf8')
for (i in group_df_20241115$X ) {
  selected_gene = i
  module_name = group_df_20241023[group_df_20241023$X==i,3]
  mt1 = as.matrix(obj@assays$RNA@data[selected_gene,])
  embed = as.data.frame(obj@reductions$diffmap@cell.embeddings)
  embed$expression = mt1
  ggplot(embed,aes(x=DiffusionMap_1,y=DiffusionMap_2,color=expression))+geom_point()+
    scale_color_gradientn(colors = jdb_palette("brewer_celsius"))+theme_classic()+
    theme(text = element_text(size=16))+ggtitle(i)
  ggsave(paste0('/data/jiangjunyao/polyATAC/figure_after_science/single cell fig/NMP_fate_analysis/gene_pseudotime_dynamics/',module_name,'_',
                i,'_UMAP.pdf'),width = 10,height = 8)
  df_exp = data.frame(c(pm_smooth[[2]][selected_gene,],
                        sp_smooth[[2]][selected_gene,]),
                      c(1:100,1:100),
                      c(rep('Mesoderm lineage',100),rep('Neuron lineage',100)))
  colnames(df_exp) = c('expression','pseudotime','lineage')
  ggplot(df_exp,aes(x=pseudotime,y=expression,color=lineage))+geom_line(size=1)+
    theme_classic()+theme(text = element_text(size=16),plot.title = element_text(hjust = 0.5))+
    scale_color_manual(values = c('#6950a1','#84bf96'))+ggtitle(selected_gene)+
    ylab('Scaled expression')+scale_y_continuous(limits = c(0, 2))
  ggsave(paste0('/data/jiangjunyao/polyATAC/figure_after_science/single cell fig/NMP_fate_analysis/gene_pseudotime_dynamics/',module_name,'_',
                i,'_line.pdf'),width = 10,height = 8)
}
### average module
list_aver = list()
for (i in unique(group_df_20241115$module)){
  gene_use = group_df_20241115[group_df_20241115$module==i,1]
  mt1 = as.matrix(dual_lineage2[gene_use,])
  list_aver[[i]] = colMeans(mt1)
}
for (i in unique(group_df_20241115$module)) {
  gene_use = group_df_20241115[group_df_20241115$module==i,1]
  mt1 = as.matrix(obj@assays$RNA@data[gene_use,])
  mt1 = colMeans(mt1)
  embed = as.data.frame(obj@reductions$diffmap@cell.embeddings)
  embed$expression = mt1
  ggplot(embed,aes(x=DiffusionMap_1,y=DiffusionMap_2,color=expression))+geom_point()+
    scale_color_gradientn(colors = jdb_palette("brewer_celsius"))+theme_classic()+
    theme(text = element_text(size=16))+ggtitle(i)
  ggsave(paste0('/data/jiangjunyao/polyATAC/figure_after_science/single cell fig/NMP_fate_analysis/module_average_expression/diffusionmap embedding_',i,'.pdf'),
         width = 10,height = 8)
  ggplot(embed,aes(x=DiffusionMap_1,y=DiffusionMap_2,color=expression))+geom_point()+
    scale_color_gradientn(colors = jdb_palette("brewer_celsius"))+theme_classic()+
    theme(text = element_text(size=16))+ggtitle(i)+guides(color = "none")
  ggsave(paste0('/data/jiangjunyao/polyATAC/figure_after_science/single cell fig/NMP_fate_analysis/module_average_expression/diffusionmap embedding_',i,'_nolenged.pdf'),
         width = 10,height = 8)
  df_exp = data.frame(c(colMeans(pm_smooth[[2]][gene_use,]),
                        colMeans(sp_smooth[[2]][gene_use,])),
                      c(1:100,1:100),
                      c(rep('Mesoderm lineage',100),rep('Neuron lineage',100)))
  colnames(df_exp) = c('expression','pseudotime','lineage')
  ggplot(df_exp,aes(x=pseudotime,y=expression,color=lineage))+geom_line(size=1)+
    theme_classic()+theme(text = element_text(size=16),plot.title = element_text(hjust = 0.5))+
    scale_color_manual(values = c('#6950a1','#84bf96'))+ggtitle(i)+
    ylab('Scaled expression')+scale_y_continuous(limits = c(0, 2))
  ggsave(paste0('/data/jiangjunyao/polyATAC/figure_after_science/single cell fig/NMP_fate_analysis/module_average_expression/lineplot_',i,'.pdf'),
         width = 10,height = 8)
}
### co-expression
no_ba = group_df_20241115[group_df_20241115[,3] != 'Balance regulate',1]
dual_lineage2 = dual_lineage2[no_ba,]
anno_row2 = data.frame(group_df_20241115[group_df_20241115[,3] != 'Balance regulate',3])
colnames(anno_row2) = 'module'
rownames(anno_row2) = rownames(dual_lineage2)
pdf('/data/jiangjunyao/polyATAC/figure_after_science/single cell fig/NMP_fate_analysis/inter module regulation_single gene.pdf',width = 10,height = 8)
pheatmap(cor(t(dual_lineage2)),cluster_cols = F,cluster_rows = F,
         show_rownames = F,annotation_row =anno_row2,annotation_col = anno_row ,
         annotation_colors = legend_colo,show_colnames = F,color = col)
dev.off()
### inter-module level co-expression
module_expression = do.call(bind_cols,list_aver[-1])
result <- rcorr(as.matrix(module_expression))
correlation_matrix <- result$r
adjusted_pvalues <- p.adjust(result$P, method = "fdr")

adjusted_pvalue_matrix <- matrix(adjusted_pvalues, nrow = nrow(correlation_matrix), ncol = ncol(correlation_matrix))
colnames(adjusted_pvalue_matrix) = colnames(correlation_matrix)
rownames(adjusted_pvalue_matrix) = rownames(correlation_matrix)
adjusted_pvalue_matrix[is.na(adjusted_pvalue_matrix)] = 0
pdf('/data/jiangjunyao/polyATAC/figure_after_science/single cell fig/NMP_fate_analysis/inter module regulation.pdf',width = 10,height = 8)
corrplot(correlation_matrix,
         method = "circle",
         type = "upper",
         p.mat = adjusted_pvalue_matrix,
         sig.level = 0.01,insig = "blank",
         pch.cex = 0.8,tl.col = "black",col = col)
dev.off()
### enrichment for each module### enrichment for each module### enrichment for each module
lm_go = gene_enrich(group_df_20241115[group_df_20241115$module=='late mesoderm regulate',1],org.Mm.eg.db,'GO')
lm_go = setReadable(lm_go,OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

ln_go = gene_enrich(group_df_20241115[group_df_20241115$module=='late neuron regulate',1],org.Mm.eg.db,'GO')
ln_go = setReadable(ln_go,OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

em_go = gene_enrich(group_df_20241115[group_df_20241115$module=='early mesoderm regulate',1],org.Mm.eg.db,'GO')
em_go = setReadable(em_go,OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

en_go = gene_enrich(group_df_20241115[group_df_20241115$module=='early neuron regulate',1],org.Mm.eg.db,'GO')
en_go = setReadable(en_go,OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

bi_go = gene_enrich(group_df_20241115[group_df_20241115$module=='Balance regulate',1],org.Mm.eg.db,'GO')
bi_go = setReadable(bi_go,OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

selected_func = rbind(
                      em_go@result[c(4,26,30,67,73,111),],
                      en_go@result[c(1,4,16,26),])

selected_func$logq = -log10(selected_func$qvalue)
selected_func$rank = 10:1
selected_func$group = c(c(
                              rep('early mesoderm regulate',6),
                          rep('early neuron regulate',4)))
write.csv(selected_func,'/data/jiangjunyao/polyATAC/nmp_final/nmp_pseudotime_module/only_rna/selected_module_GO_terms_20250709.csv')
ggplot(selected_func,aes(y=reorder(Description,rank),x=logq,fill=group))+
  geom_col()+theme_classic()+scale_x_continuous(expand = c(0,0))+
  xlab('-log10 qvalue')+
  ylab('GO terms')+theme(text = element_text(size=16))+scale_fill_manual(values = legend_colo$module)
ggsave('/data/jiangjunyao/polyATAC/figure_after_science/single cell fig/NMP_fate_analysis/module_GO.pdf',
       width = 10,height = 8)
write.csv(bi_go@result,
          '/data/jiangjunyao/polyATAC/nmp_final/nmp_pseudotime_module/only_rna/Bipotent_go_20250709.csv')
write.csv(lm_go@result,
          '/data/jiangjunyao/polyATAC/nmp_final/nmp_pseudotime_module/only_rna/lm_go_20250709.csv')
write.csv(ln_go@result,
          '/data/jiangjunyao/polyATAC/nmp_final/nmp_pseudotime_module/only_rna/ln_go_20250709.csv')
write.csv(em_go@result,
          '/data/jiangjunyao/polyATAC/nmp_final/nmp_pseudotime_module/only_rna/em_go_20250709.csv')
write.csv(en_go@result,
          '/data/jiangjunyao/polyATAC/nmp_final/nmp_pseudotime_module/only_rna/en_go_20250709.csv')


