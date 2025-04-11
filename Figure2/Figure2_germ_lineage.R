library(pheatmap)
library(ggtree)
metadata1 = all_sclt_V3@meta.data
metadata1 = metadata1[!is.na(metadata1$barcodes),]
metadata1 = metadata1[metadata1$pgen<0.001,]
table(metadata1$germ)
which(colnames(metadata1)=='germ')
metadata1 = metadata1[,c('germ','barcodes','time')]
colnames(metadata1)[1] = 'celltype'
f1 = FateMapper::fate_mapping(metadata1)
write.csv(f1,'/data/jiangjunyao/polyATAC/germ_lineage/all_germ_lineage.csv')
heatmap_col = c('#F7FBFF','#DEEBF7','#C6DBEF','#9ECAE1','#6BAED6','#4292C6',
                '#2171B5','#08519C','#08306B')
germ_col = c('#98DF8AFF',
             '#DBDB8DFF',
             '#AEC7E8FF',
             '#C5B0D5FF',
             '#FFBB78FF')
germ = c('NeuroEctoderm','SurfaceEctoderm','Mesoderm',
         'Endoderm','Extraembryonic')
for (i in unique(metadata1$time)) {
  a1= metadata1[metadata1$time==i,]
  colnames(a1)[1] = 'celltype'
  f1 = FateMapper::fate_mapping(a1)
  write.csv(f1,
            paste0('/data/jiangjunyao/polyATAC/germ_lineage/',i,'.csv'))
}

dir1 = dir('/data/jiangjunyao/polyATAC/germ_lineage/')
dir1 = dir1[grep('coup',dir1)]
for (i in dir1) {
  
  coup = read.csv(paste0('/data/jiangjunyao/polyATAC/germ_lineage/',i),row.names = 1)
  coup <- coup[rownames(coup) != "Extraembryonic", ]
  coup <- coup[, colnames(coup) != "Extraembryonic"]
  lineage_simi = pheatmap(coup,clustering_method = 'mcquitty')$tree_row
  p1 = ggtree(lineage_simi,size=0.7)
  label_use <- data.frame(label =lineage_simi$labels )
  label_use$celltype = label_use[,1]
  p1=p1%<+%label_use+geom_tiplab(offset = 0.0005)+xlim(NA,0.001)+geom_tippoint(aes(col=celltype))+
    scale_color_manual(values = germ_col,breaks = germ)
  print(p1)
  ggsave(paste0('/data/jiangjunyao/polyATAC/polyatac_figure/single cell fig/germ_tree/',
                i,'.pdf'),width = 8,height = 6)

  
  # coup = read.csv(paste0('/data/jiangjunyao/polyATAC/germ_lineage/',i),row.names = 1)
  # pheatmap(coup,clustering_method = 'single',
  #          color = colorRampPalette(heatmap_col)(50),
  #          border_color = NA,filename =paste0('/data/jiangjunyao/polyATAC/polyatac_figure/single cell fig/germ_tree/tree_',
  #                                             i,'.pdf'),width = 18,height = 16,cellwidth = 20,cellheight = 20)
}

