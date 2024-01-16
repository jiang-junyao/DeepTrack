multi_relationship = function(obj){
  library(ggplot2)
  library(FateMapper)
  ### seurat object with cell_fate & cell_type
  obj = NormalizeData(obj,assay = 'activity',normalization.method = 'LogNormalize')
  obj = NormalizeData(obj,assay = 'RNA',normalization.method = 'LogNormalize')
  cell_use = rownames(obj@meta.data)[!is.na(obj@meta.data$cell_fate)]
  obj_plot = subset(obj,cells=cell_use)
  obj_plot@active.ident = as.factor(obj_plot$cell_type)
  ave_exp = AverageExpression(obj_plot,slot = 'data')$RNA
  ave_atac = AverageExpression(obj_plot,slot = 'data')$activity
  ave_atac = ave_atac[rowSums(ave_atac)>1,]
  cor.exp = reshape2::melt(cor(ave_exp,method = 'spearman'))
  cor.atac = reshape2::melt(cor(ave_atac,method = 'spearman'))
  cell_fate = obj_plot@meta.data[,c('cell_fate','cell_type')]
  cell_fate = cell_fate[!is.na(cell_fate[,1]),]
  colnames(cell_fate) = c('barcodes','celltype')
  cor.barcode = reshape2::melt(as.matrix(cell_type_fate_similartiy(cell_fate,plot=F,out_similar_mt = T,method = 'spearman')))
  multi_rela = cbind(cor.barcode,cor.exp[,3],cor.atac[,3])
  multi_rela = multi_rela[multi_rela$Var1!=multi_rela$Var2,]
  colnames(multi_rela)[3:5] = c('lineage','rna','atac')
  min_num = min(multi_rela[,4],multi_rela[,5])
  if (min_num-0.05 < 0) {
    min_num = 0
  }else{
    min_num = min_num-0.05
  }
  rl_slope = round(coef(lm(rna ~ lineage, data = multi_rela))[2],3)
  al_slope = round(coef(lm(atac ~ lineage, data = multi_rela))[2],3)
  p1=ggplot(multi_rela, aes(x = atac, y = lineage)) +
    geom_point(size = 4) +
    labs(x = "ATAC pair correlation", y = "Barcode pair correlation") +
    theme_minimal()+xlim(c(min_num,1))+ylim(c(min_num,1))+geom_smooth(method = 'lm')+
    ggtitle(paste0('atac+lineage slope: ',al_slope))
  p2=ggplot(multi_rela, aes(x = rna, y = lineage)) +
    geom_point(size = 4) +
    labs(x = "RNA pair correlation", y = "Barcode pair correlation") +
    theme_minimal()+xlim(c(min_num,1))+ylim(c(min_num,1))+geom_smooth(method = 'lm')+
    ggtitle(paste0('rna+lineage slope: ',rl_slope))
  print(p1|p2)
  return(multi_rela)
}

rna_atac_relationship <- function(obj){
  library(ggplot2)
  metadata = data.frame(colnames(obj),obj$cell_fate)
  metadata = metadata[!is.na(metadata[,2]),]
  metadata[,2] = as.character(metadata[,2])
  colnames(metadata) = c('cellname','barcodes')
  obj_use = subset(obj,cells = metadata[,1])
  ave_exp = AverageExpression(obj_use,slot = 'data',group.by = 'cell_fate')$RNA
  ave_atac = AverageExpression(obj_use,slot = 'data',group.by = 'cell_fate')$activity
  share_gene = intersect(rownames(ave_exp),rownames(ave_atac))
  ave_exp = ave_exp[share_gene,]
  ave_atac = ave_atac[share_gene,]
  clone_cor = cor(x=ave_exp,y=ave_atac,method = 'pearson')
  clone_cor = reshape2::melt(clone_cor)
  clone_cor$cor_type = 'interclonal'
  clone_cor[clone_cor[,1]==clone_cor[,2],4] = 'intraclonal'
  p.val = wilcox.test(clone_cor[clone_cor[,1]==clone_cor[,2],3],
              clone_cor[clone_cor[,1]!=clone_cor[,2],3])
  ggplot(clone_cor,aes(x=cor_type,y=value,fill=cor_type))+geom_boxplot()+
    theme_minimal()+xlab('')+ylab('Expression-accessibilty corr')+
    theme(text = element_text(size=12),axis.text = element_text(size = 12))
}
