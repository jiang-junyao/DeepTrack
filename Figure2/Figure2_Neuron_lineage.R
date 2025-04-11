library(FateMapper)
library(ggtree)
cns_col = c(rgb(29/255,75/255,42/255),
            rgb(74/255,105/255,75/255),
            rgb(45/255,133/255,96/255),
            rgb(79/255,146/255,69/255),
            rgb(149/255,162/255,85/255),
            rgb(51/255,168/255,117/255),
            rgb(189/255,206/255,123/255),
            rgb(172/255,219/255,47/255),
            rgb(216/255,186/255,106/255))
names(cns_col) = c('Retinal primordium','Di telencephalon','Mesencephalon MHB',
                   'Anterior floor plate','Posterior floor plate','Hindbrain',
                   "Spinal cord (anterior)","Spinal cord (posterior)",'Neuromesodermal progenitors')
brain_sub <- readRDS("~/polyATAC/E9_10 visium/brain_sub.rds")
metadata = brain_sub@meta.data
metadata = metadata[!is.na(metadata$barcodes),]
metadata = metadata[metadata$pgen<0.001,]
UMAPPlot(brain_sub,group.by='celltype2')+theme_void()+
  scale_color_manual(values = cns_col)+theme(text = element_text(size=16))
ggsave('/data/jiangjunyao/polyATAC/polyatac_figure/single cell fig/subcluster lineage analysis/brain_AP_annotation.pdf',
       width = 12,height = 8)

for (i in unique(brain_sub$time)) {
  meta_use = metadata[metadata$time==i,]
  clone_mt=fate_mapping(meta_use)
  write.csv(clone_mt,paste0('/data/jiangjunyao/polyATAC/subcelltype_coupling/CNS_time/',i,'.csv'))
}
dir1 = dir('/data/jiangjunyao/polyATAC/subcelltype_coupling/CNS_time')
dir1 = dir1[grep('coupling',dir1)]
for (i in dir1[2:4]) {

  coup = read.csv(paste0('/data/jiangjunyao/polyATAC/subcelltype_coupling/CNS_time/',i),row.names = 1)
  lineage_simi = pheatmap(coup,clustering_method = 'mcquitty',display_numbers = T)$tree_row
  p1 = ggtree(lineage_simi,size=0.7)
  label_use <- data.frame(label =lineage_simi$labels )
  label_use$celltype = label_use[,1]
  p1=p1%<+%label_use+geom_tiplab(offset = 0.009)+xlim(NA,0.001)+geom_tippoint(aes(col=celltype))+
    scale_color_manual(values = cns_col)
  print(p1)
  ggsave(paste0('/data/jiangjunyao/polyATAC/polyatac_figure/single cell fig/subcluster lineage analysis/CNS_tree_',
                i,'.pdf'),width = 8,height = 6)
}
