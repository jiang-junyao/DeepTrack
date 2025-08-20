library(Seurat)
library(tidyverse)
library(BuenColors)
obj <- readRDS("~/polyATAC/polyATAC_20250616.rds")
proj = readRDS('/data/jiangjunyao/polyATAC/multiomi/archr_subset_20250604/Save-ArchR-Project.rds')
atac_qua = data.frame(proj$nFrags,proj$Sample)
colnames(atac_qua) = c('nFrag','sample')
obj$mouse = 'mouse2'
obj@meta.data[obj$Batch=='mouse1','mouse']='mouse1'
obj@meta.data[obj$Batch%in%c('tail1','tail2'),'mouse']='mouse3 (tail)'
cell_inter = rownames(proj@cellColData)[rownames(proj@cellColData)%in% colnames(obj)]
obj$nFrag = proj@cellColData[cell_inter,'nFrags']

ggplot(atac_qua,aes(x=sample,y=nFrag,fill=sample))+geom_violin()+
  theme_classic()+theme(text = element_text(size=16))

VlnPlot(obj,group.by = 'mouse',features = c('nCount_RNA'),pt.size=0)+
  scale_fill_manual(values = jdb_palette('corona'))
ggsave('/data/jiangjunyao/polyATAC/figure_after_science/multiomics fig/cell type statistic and fraction plot/nCount_RNA_per_mouse.pdf',width = 10,height = 8)
VlnPlot(obj,group.by = 'mouse',features = c('nFrag'),pt.size=0)+
  scale_fill_manual(values = jdb_palette('corona'))
ggsave('/data/jiangjunyao/polyATAC/figure_after_science/multiomics fig/cell type statistic and fraction plot/nFrags_ATAC_per_mouse.pdf',width = 10,height = 8)


ggplot(meta1,aes(x=Batch,y=nCount_RNA,fill=Batch))+geom_violin()+
  theme_classic()+theme(text = element_text(size=16))
