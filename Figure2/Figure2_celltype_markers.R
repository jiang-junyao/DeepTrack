.libPaths('/data/jiangjunyao/miniconda3/envs/r421/lib/R/library')
library(tidyverse)
library(Seurat)
### load data
rna = readRDS('E:\\polyATAC\\all_sc/all_sclt_V3.rds')

### set cell type order
ct_neworder <- rev(c('Di/telencephalon','Mesencephalon/MHB','Retinal primordium',
                     'Motor neurons','Neuron progenitor cells','Hindbrain',
                     'Roof plate','Anterior floor plate','Posterior floor plate',
                     'Spinal cord','Neural crest','Neuromesodermal progenitors',
                     'Otic epithelium','Placodal area','Pre-epidermal keratinocytes',
                     'Notochord','Osteoblast progenitors',
                     'Skeletal muscle progenitors','Chondrocyte and osteoblast progenitors',
                     'Pharyngeal mesoderm','Paraxial mesoderm','Intermediate mesoderm',
                     'Renal epithelium','Limb mesenchyme progenitors','Mesenchyme',
                     'Cardiomyocytes','Splanchnic mesoderm',
                     'Amniochorionic mesoderm','Extraembryonic mesoderm','Allantois',
                     'Endothelium','Hematoendothelial progenitors','Blood progenitors',
                     'Primitive erythroid cells','Gut','Hepatocyte','Extraembryonic ectoderm',
                     'Parietal endoderm','Visceral endoderm','Extraembryonic visceral endoderm'))

### set marker for each cell type
### The order of markers must be consistent with the order of cell types
marker_use = c('Wnt8b','Pax6','Foxg1','En1','Pax5','Rax',
               'Six3','Nlgn1','Olig2','Isl2','Lhx4','Mnx1',
               'Lhx3','Chat','Mybl1','Prmt8','Ascl1','Sox2','Crabp1',
               'Egr2','Lmx1a','Msx1','Wnt1','Ntn1','Shh','Lmx1b',
               'Foxa2','Shh','Nkx2-2','Foxa2','Nkx6-1','Foxb1',
               'Pax6','Hes5','Sox10','Dlx2','Foxd3','Tfap2b','Hoxb9',
               'Epha5','Hes3','Cdx4','Fgf8','Six1','Dlx5','Gata3',
               'Eya1','Trp63','Krt5','Noto','Chrd','T','Prrx1','Myf5',
               'Fap','Lamc3','Pax1','Pax9','Tbx1','Tbx6','Dll1',
               'Osr1','Wt1','Emx2','Sim1','Mecom','Cpa2','Fgf10',
               'Tbx5','Col3a1','Tnc','Pdgfra','Nkx2-5','Smarcd3',
               'Gata4','Hcn4','Myh7','Tcf21','Gata4','Foxf1','Bmp2',
               'Postn','Col1a1','Hand1','Tbx4','Hoxa10','Twist2','Bmp4',
               'Cdh5','Pecam1','Kdr','Eng','Lmo2','Fgd5','Etv2','Runx1',
               'Gata1','Gata1','Hba-a2','Apela','Foxa1','Alb','A1cf',
               'Hhex','Hnf4a','Elf5','Tfap2c','Ascl2','Afp',
               'Sparc','Glis1','Krt8','Krt18','Apoe','Emb','Ttr')
marker_use=unique(marker_use)
rna$celltype[rna$celltype=='Caudal mesoderm'] = 'Pharyngeal mesoderm'
### dot plot
rna$ct_plot <- factor(rna$celltype,levels = ct_neworder)
DotPlot(rna,features = marker_use,group.by = 'ct_plot',dot.scale=5,dot.min=0,
        scale.by='size',col.max = 1,col.min=-1) + 
  scale_colour_gradient2(low='#0072BD', mid='#fffffb', high='#d93a49')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(text = element_text(size=13))

### save figure
ggsave('E:\\polyATAC\\single cell fig\\cell type statistic and fraction plot/celltype_markers_expression.pdf',width = 24,
       height = 8)
