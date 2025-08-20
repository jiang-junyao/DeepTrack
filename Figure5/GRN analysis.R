library(tidyverse)
library(venn)
library(ggrepel)
library(BuenColors)
library(org.Mm.eg.db)
library(reshape2)
grn_score = read.csv('/data/jiangjunyao/polyATAC/nmp_multi_20250717/grn/grn_score.csv')
nmp_deg_two_fate_30 <- read.csv("~/polyATAC/nmp_multi_20250717/differential analysis/nmp_deg_two_fate_30.csv", row.names=1)
pm_bias_grn <- read.csv("/data/jiangjunyao/polyATAC/nmp_multi_20250717/grn/pxmnmp_grn.csv", row.names=1)
pm_bias_grn = pm_bias_grn[!is.na(pm_bias_grn$p),]
sp_bias_grn <- read.csv("/data/jiangjunyao/polyATAC/nmp_multi_20250717/grn/scnmp_grn.csv", row.names=1)
sp_bias_grn = sp_bias_grn[!is.na(sp_bias_grn$p),]
pm_bias_grn = pm_bias_grn[pm_bias_grn$p<0.001,]
sp_bias_grn = sp_bias_grn[sp_bias_grn$p<0.001,]
pm_bias_grn = pm_bias_grn[order(pm_bias_grn$coef_abs,decreasing = T),]
sp_bias_grn = sp_bias_grn[order(sp_bias_grn$coef_abs,decreasing = T),]
pm_bias_grn_use = pm_bias_grn[1:2000,]
sc_bias_grn_use = sp_bias_grn[1:2000,]

### out pxmnmp specific grn
for_pxm = nmp_deg_two_fate_30
#pxmnmp_mean_exp = AverageExpression(nmp_multi,assays = 'RNA',slot = 'scale.data',group.by = 'celltype')$RNA
for_pxm[for_pxm$cluster!='PXM NMP',2] = -for_pxm[for_pxm$cluster!='PXM NMP',2]
rownames(for_pxm) = for_pxm$gene

grn_score_pxmnmp  = grn_score[grn_score$cluster=='PXM NMP',]
pm_bias_grn_use_purple = pm_bias_grn_use[pm_bias_grn_use$source%in% c('Ets1','Foxc1','Cdx2','Id2','Tfap2b','Ebf3','Rreb1','Lmo2','Sp7','Twist1','Foxf2'),]
pm_bias_grn_use_purple$regulation_type = ifelse(pm_bias_grn_use_purple$coef_mean>0,'activate','repress')
pm_bias_grn_use_purple$regulation_source = pm_bias_grn_use_purple$source
pm_bias_grn_use_purple$regulation_source_num = pm_bias_grn_use_purple$source
pm_bias_grn_use_purple$regulation_source_num = as.numeric(as.factor(pm_bias_grn_use_purple$regulation_source_num))
pm_bias_grn_use_purple$source_degree_centrality = grn_score_pxmnmp[match(pm_bias_grn_use_purple$source,grn_score_pxmnmp$X),3]
pm_bias_grn_use_purple = pm_bias_grn_use_purple[pm_bias_grn_use_purple$target %in% for_pxm$gene,]
pm_bias_grn_use_purple$target_log2fc = for_pxm[pm_bias_grn_use_purple$target,2]

write.csv(pm_bias_grn_use_purple,'/data/jiangjunyao/polyATAC/nmp_multi_20250717/grn/pm_bias_grn_use_purple.csv',
            quote = F,row.names = F)
###
for_sc = nmp_deg_two_fate_30
#pxmnmp_mean_exp = AverageExpression(nmp_multi,assays = 'RNA',slot = 'scale.data',group.by = 'celltype')$RNA
for_sc[for_sc$cluster=='PXM NMP',2] = -for_sc[for_sc$cluster=='PXM NMP',2]
rownames(for_sc) = for_sc$gene

scnmp_deg = nmp_deg_two_fate_30[nmp_deg_two_fate_30$cluster=='SC NMP','gene']
grn_score_scnmp  = grn_score[grn_score$cluster=='SC NMP',]
sc_bias_grn_use_purple = sc_bias_grn_use[sc_bias_grn_use$source%in% c('Rfx4','Meis1','Zic1','Neurog2','Nr2f2','Pax3','Pax6','Zic4','Ebf1','Neurod4','Glis3','Rarb','Pou3f2'),]
sc_bias_grn_use_purple$regulation_type = ifelse(sc_bias_grn_use_purple$coef_mean>0,'activate','repress')
sc_bias_grn_use_purple$regulation_source = sc_bias_grn_use_purple$source
sc_bias_grn_use_purple$regulation_source_num = sc_bias_grn_use_purple$source
sc_bias_grn_use_purple$regulation_source_num = as.numeric(as.factor(sc_bias_grn_use_purple$regulation_source_num))
sc_bias_grn_use_purple$source_degree_centrality = grn_score_scnmp[match(sc_bias_grn_use_purple$source,grn_score_scnmp$X),3]
sc_bias_grn_use_purple = sc_bias_grn_use_purple[sc_bias_grn_use_purple$target %in% for_sc$gene,]
sc_bias_grn_use_purple$target_log2fc = for_sc[sc_bias_grn_use_purple$target,2]
write.csv(sc_bias_grn_use_purple,'/data/jiangjunyao/polyATAC/nmp_multi_20250717/grn/sc_bias_grn_use_green.csv',
          quote = F,row.names = F)


### T
pm_t = pm_bias_grn_use[pm_bias_grn_use$source=='T',]
pm_t$regulation_type = ifelse(pm_t$coef_mean>0,'activate','repress')
write.table(pm_t,
            '/data/jiangjunyao/polyATAC/nmp_multi_20250717/grn/pm_T_grn.txt',
            quote = F,sep = '\t',row.names = F)
### meis1
sp_meis1 = sp_bias_grn_use[sp_bias_grn_use$source=='Meis1',]
sp_meis1$regulation_type = ifelse(sp_meis1$coef_mean>0,'activate','repress')
active_enrich = gene_enrich(sp_meis1[sp_meis1$coef_mean>0,'target'],org.Mm.eg.db,'GO')
active_enrich = setReadable(active_enrich,OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
repress_enrich = gene_enrich(sp_meis1[sp_meis1$coef_mean<0,'target'],org.Mm.eg.db,'GO')
repress_enrich = setReadable(repress_enrich,OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
write.table(sp_meis1,
            '/data/jiangjunyao/polyATAC/nmp_multi_20250616/grn_final/sp_meis1.txt',
            quote = F,sep = '\t',row.names = F)
selected_func = rbind(active_enrich@result[c(1,2,5,9),],
                      repress_enrich@result[c(4,35,36,40),])
selected_func$logp = -log10(selected_func$qvalue)
selected_func$rank = 8:1
selected_func$group = c(rep('activated genes',4),rep('repressed genes',4))
ggplot(selected_func,aes(x=logp,y=reorder(Description,rank),fill=group))+
  geom_col()+scale_x_continuous(expand=c(0,0))+theme_classic()+
  xlab('-log10 qvalue')+ylab('GO terms')+
  theme(text = element_text(size=16))+
  scale_fill_manual(values = c('#FEB24C','#BDBDBD'))
ggsave('/data/jiangjunyao/polyATAC/polyatac_figure/multiomics fig/NMP_fate_analysis/nmp_grn/GO_for_Meis1_downstream.pdf',
       width = 10,height = 8)
sp_meis1_ac = sp_meis1[sp_meis1$regulation_type=='activate',]
sp_meis1_re = sp_meis1[sp_meis1$regulation_type!='activate',]
write.csv(sp_meis1_ac,'/data/jiangjunyao/polyATAC/nmp_final/celloracle_grn/filtered_grn/Meis1_activate.csv',quote = F,row.names = F)
write.csv(sp_meis1_re,'/data/jiangjunyao/polyATAC/nmp_final/celloracle_grn/filtered_grn/Meis1_repress.csv',quote = F,row.names = F)
write.csv(selected_func,'/data/jiangjunyao/polyATAC/nmp_final/celloracle_grn/Meis1_target_GO_selected.csv')
### cdx2
pm_cdx2 = pm_bias_grn_use[pm_bias_grn_use$source=='Cdx2',]
sc_cdx2 = sc_bias_grn_use[sc_bias_grn_use$source=='Cdx2',]
# pm_cdx2$index = paste0(pm_cdx2$source,'-',pm_cdx2$target)
# sp_cdx2$index = paste0(sp_cdx2$source,'-',sp_cdx2$target)
pm_cdx2$regulation_type = ifelse(pm_cdx2$coef_mean>0,'activate','repress')
write.table(pm_cdx2,
            '/data/jiangjunyao/polyATAC/nmp_multi_20250717/grn/pm_Cdx2_grn.txt',
            quote = F,sep = '\t',row.names = F)
write.csv(pm_cdx2[pm_cdx2$coef_mean>0,],'/data/jiangjunyao/polyATAC/nmp_multi_20250717/grn/pm_Cdx2_activate.csv',quote = F)
write.csv(pm_cdx2[pm_cdx2$coef_mean<0,],'/data/jiangjunyao/polyATAC/nmp_multi_20250717/grn/pm_Cdx2_repress.csv',quote = F)
active_enrich = gene_enrich(pm_cdx2[pm_cdx2$coef_mean>0,'target'],org.Mm.eg.db,'GO')
active_enrich = setReadable(active_enrich,OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
active_enrich = active_enrich@result
active_enrich = active_enrich[active_enrich$qvalue<0.05 & active_enrich$Count>1,]
repress_enrich = gene_enrich(pm_cdx2[pm_cdx2$coef_mean<0,'target'],org.Mm.eg.db,'GO',use_internal_data = T)
repress_enrich = setReadable(repress_enrich,OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
repress_enrich = repress_enrich@result
repress_enrich = repress_enrich[repress_enrich$qvalue<=0.05 & repress_enrich$Count>1,]
repress_enrich$idx = 1:nrow(repress_enrich)

selected_func = rbind(active_enrich[c(4,62,63),],
                      repress_enrich[c(13,40,52),])
selected_func$logp = -log10(selected_func$qvalue)
selected_func$rank = 6:1
selected_func$group = c(rep('activated genes',3),rep('repressed genes',3))
ggplot(selected_func,aes(x=logp,y=reorder(Description,rank),fill=group))+
  geom_col()+scale_x_continuous(expand=c(0,0))+theme_classic()+
  xlab('-log10 qvalue')+ylab('GO terms')+
  theme(text = element_text(size=16))+
  scale_fill_manual(values = c('#FEB24C','#BDBDBD'))
ggsave('/data/jiangjunyao/polyATAC/figure_after_science/multiomics fig/NMP_fate_analysis/GO_for_Cdx2_downstream.pdf',
       width = 10,height = 8)
write.csv(selected_func,'/data/jiangjunyao/polyATAC/nmp_multi_20250717/grn/Cdx2_target_GO_selected.csv')


pm_bias_grn_use_tf = pm_bias_grn_use[pm_bias_grn_use$target %in% pm_bias_grn_use$source,]
sp_bias_grn_use_tf = sp_bias_grn_use[sp_bias_grn_use$target %in% sp_bias_grn_use$source,]
pm_bias_grn_use_tf$regulation_type = ifelse(pm_bias_grn_use_tf $coef_mean>0,'activate','repress')
sp_bias_grn_use_tf$regulation_type = ifelse(sp_bias_grn_use_tf $coef_mean>0,'activate','repress')
write.table(pm_bias_grn_use_tf,
          '/data/jiangjunyao/polyATAC/nmp_final/celloracle_grn/filtered_grn/pm_bias_tf_network.txt',
          quote = F,sep = '\t',row.names = F)
write.table(sp_bias_grn_use_tf,
          '/data/jiangjunyao/polyATAC/nmp_final/celloracle_grn/filtered_grn/sp_bias_tf_network.txt',
          quote = F,sep = '\t',row.names = F)


###GRN SCORE analysis
# grn_score_pxmnmp = read.csv('/data/jiangjunyao/polyATAC/nmp_multi_20250616/grn_pxm_tunning/30_10_1500_grn_score_pxmnmp.csv')
# grn_score_scnmp = read.csv('/data/jiangjunyao/polyATAC/nmp_multi_20250616/grn_sc_tunning/30_10_1500_grn_score_scnmp.csv')
grn_score = read.csv('/data/jiangjunyao/polyATAC/nmp_multi_20250717/grn/grn_score.csv')
#grn_score = rbind(grn_score_pxmnmp,grn_score_scnmp)
grn_group_score = list()
group_tf_list = list()
top_number = 20
all_tf = c()
for (i in c('PXM NMP','SC NMP')) {
  grn_score_use = grn_score[grn_score$cluster==i,]
  grn_score_use = grn_score_use[order(grn_score_use$degree_centrality_all,decreasing = T),]
  group_tf_list[[i]] = grn_score_use[1:top_number,1]
  all_tf = c(all_tf,grn_score_use[1:top_number,1])
  grn_score_use = grn_score_use[1:top_number,]
  grn_group_score[[i]] = grn_score_use
  p1=ggplot(grn_score_use,aes(y=reorder(X,degree_centrality_all),x=degree_centrality_all,color='red'))+
    geom_point(size=4)+scale_color_manual(values = '#145b7d')+theme_minimal()+
    xlab('degree centrality')+ylab('')+theme(legend.position = "none",
                                             text = element_text(size=14))+
    ggtitle(paste0('Top 20 in ',i))+theme(plot.title = element_text(hjust = 0.5))
  print(p1)
  ggsave(paste0('/data/jiangjunyao/polyATAC/figure_after_science/multiomics fig/NMP_fate_analysis/grn/Node_degree_',i,'.pdf'),width = 6,height = 5)
}
venn_result=venn(group_tf_list,plotsize = T,zcolor = jdb_palette('corona'),ilabels = "counts",box=F)

grn_group_score_df = do.call(bind_rows,grn_group_score)
ggplot(grn_group_score_df,aes(y=X,x=cluster,color=degree_centrality_all))+
  geom_point(aes(size=betweenness_centrality))


### plot pm GRN VS SP GRN
pm_sp = grn_score[grn_score$cluster %in% c('PXM NMP',
                                           'SC NMP'),]
pm_sp_plot = dcast(pm_sp,X~cluster,value.var = 'degree_centrality_all',
                    fill = 0)
colnames(pm_sp_plot)[2:3] = c('PM','SP')

pm_sp_plot$type = 'other genes'
thr_use = 0.1
pm_sp_plot[pm_sp_plot$PM>thr_use  & pm_sp_plot$SP<=thr_use ,'type'] = 'Paraxial mesoderm bias specific'
pm_sp_plot[pm_sp_plot$SP>thr_use  & pm_sp_plot$PM<=thr_use ,'type'] = 'Spinal cord bias specific'
pm_sp_plot[pm_sp_plot$SP>thr_use  & pm_sp_plot$PM>thr_use ,'type'] = 'both'

top_pm <- pm_sp_plot[pm_sp_plot$PM>thr_use ,]
top_sc <- pm_sp_plot[pm_sp_plot$SP>thr_use ,]
dual_tf = intersect(top_pm[,1],top_sc[,1])
top_combined <- bind_rows(top_pm, top_sc) %>%
  distinct()

ggplot(pm_sp_plot,aes(x=PM,y=SP,color=type))+geom_point()+theme_classic()+
  geom_text_repel(data = top_combined,
                  aes(label = X),
                  size = 4,
                  max.overlaps = 15,
                  box.padding = unit(1, "lines"),
                  point.padding = unit(0.5, "lines"),
                  segment.color = 'grey')+
  scale_color_manual(values = c('#f47a55','grey','#6950a1','#84bf96'))+
  xlab('degree centrality for gens in GRN of \n Paraxial mesoderm bias NMP')+
  ylab('degree centrality for gens in GRN of \n Spinal cord bias NMP')+
  geom_vline(xintercept = thr_use, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = thr_use, linetype = "dashed", color = "grey")
ggsave('/data/jiangjunyao/polyATAC/figure_after_science/multiomics fig/NMP_fate_analysis/grn/differential_network_node.pdf',
       width = 10,height = 8)


#         z = ~BA,
#         color = ~type,  # 根据type列设置颜色
#           # 选择颜色范围，也可以自定义其他颜色
#         text = ~X,
#         type = 'scatter3d',
#         mode = 'markers') %>%
#   layout(scene = list(xaxis = list(title = 'Degree in Paraxial mesoderm bias NMP'),
#                       yaxis = list(title = 'Degree in Spinal cord bias NMP'),
#                       zaxis = list(title = 'Degree in Balance NMP')),title = "")


## venn gene module
library(venn)
group_df_20241115 <- read.csv("~/polyATAC/nmp_final/nmp_pseudotime_module/only_rna/group_df_20250709.csv",row.names = 1)
module_tf_pm = rownames(group_df_20241115[group_df_20241115$module=="early mesoderm regulate",])
module_tf_sp = rownames(group_df_20241115[group_df_20241115$module=="early neuron regulate",])
grn_pm = pm_sp_plot[pm_sp_plot$type=='Paraxial mesoderm bias specific',1]
grn_sp = pm_sp_plot[pm_sp_plot$type=='Spinal cord bias specific',1]
list_venn_pm = list('PXM NMP specific'=grn_pm,
                 'early mesoderm regulate gene'=module_tf_pm)
list_venn_sp = list('SC NMP specific'=grn_sp,
                 'early neuron regulate gene'=module_tf_sp)
venn(list_venn_pm,ilabels = "counts",box=F,zcolor = jdb_palette('corona'),ggplot = T,
     ilcs=1,sncs=1.2)
ggsave('/data/jiangjunyao/polyATAC/polyatac_figure/multiomics fig/NMP_fate_analysis/nmp_grn/venn_early_meso_grn.pdf',
       width = 10,height = 8)
venn(list_venn_sp,ilabels = "counts",box=F,zcolor = jdb_palette('corona'),ggplot = T,
     ilcs=1,sncs=1.2)
ggsave('/data/jiangjunyao/polyATAC/polyatac_figure/multiomics fig/NMP_fate_analysis/nmp_grn/venn_early_neuo_grn.pdf',
       width = 10,height = 8)
intersect(list_venn_pm$`PXM NMP specific`,list_venn_pm$`early mesoderm regulate gene`)
intersect(list_venn_sp$`SC NMP specific`,list_venn_sp$`early neuron regulate gene`)

write.table(grn_pm,'/data/jiangjunyao/polyATAC/polyatac_figure/multiomics fig/NMP_fate_analysis/grn_venn_plot_gene_list/PM GRN.csv',
          quote = F,row.names = F,col.names = F)
write.table(grn_sp,'/data/jiangjunyao/polyATAC/polyatac_figure/multiomics fig/NMP_fate_analysis/grn_venn_plot_gene_list/SP GRN.csv',
          quote = F,row.names = F,col.names = F)
write.table(module_tf_pm,'/data/jiangjunyao/polyATAC/polyatac_figure/multiomics fig/NMP_fate_analysis/grn_venn_plot_gene_list/early mesoderm gene.csv',
          quote = F,row.names = F,col.names = F)
write.table(module_tf_sp,'/data/jiangjunyao/polyATAC/polyatac_figure/multiomics fig/NMP_fate_analysis/grn_venn_plot_gene_list/early neuron gene.csv',
          quote = F,row.names = F,col.names = F)

### overlap deg
deg_pm = nmp_deg_two_fate_30[nmp_deg_two_fate_30$cluster=='PXM NMP','gene']
deg_sp = nmp_deg_two_fate_30[nmp_deg_two_fate_30$cluster!='PXM NMP','gene']
list_venn_pm = list('PXM NMP specific'=grn_pm,
                    'early mesoderm regulate gene'=deg_pm)
list_venn_sp = list('SC NMP specific'=grn_sp,
                    'early neuron regulate gene'=deg_sp)
venn(list_venn_pm,ilabels = "counts",box=F,zcolor = jdb_palette('corona'),ggplot = T,
     ilcs=1,sncs=1.2)

venn(list_venn_sp,ilabels = "counts",box=F,zcolor = jdb_palette('corona'),ggplot = T,
     ilcs=1,sncs=1.2)
intersect(list_venn_pm$`PXM NMP specific`,list_venn_pm$`early mesoderm regulate gene`)
intersect(list_venn_sp$`SC NMP specific`,list_venn_sp$`early neuron regulate gene`)
VlnPlot(nmp_multi,intersect(list_venn_pm$`PXM NMP specific`,list_venn_pm$`early mesoderm regulate gene`),group.by='celltype')
VlnPlot(nmp_multi,intersect(list_venn_sp$`SC NMP specific`,list_venn_sp$`early neuron regulate gene`),group.by='celltype')
###
pm_gene_all = pm_gene_all[,c(1,3)]
sc_gene_all = sc_gene_all[,c(1,3)]
colnames(pm_gene_all) = c('gene','degree centrality')
colnames(sc_gene_all) = c('gene','degree centrality')
pm_gene_all = pm_gene_all[order(pm_gene_all$`degree centrality`,decreasing = T),]
sc_gene_all = sc_gene_all[order(sc_gene_all$`degree centrality`,decreasing = T),]
write.table(pm_gene_all,'/data/jiangjunyao/polyATAC/polyatac_figure/multiomics fig/NMP_fate_analysis/grn_venn_plot_gene_list/PM GRN gene degree.csv',
            quote = F,row.names = F,col.names = F)
write.table(sc_gene_all,'/data/jiangjunyao/polyATAC/polyatac_figure/multiomics fig/NMP_fate_analysis/grn_venn_plot_gene_list/SC GRN gene degree.csv',
            quote = F,row.names = F,col.names = F)
pm_tf_all=pm_gene_all[pm_gene_all$gene %in% pm_bias_grn_use$source,]
sc_tf_all=sc_gene_all[sc_gene_all$gene %in% sp_bias_grn_use$source,]
write.table(pm_tf_all,'/data/jiangjunyao/polyATAC/polyatac_figure/multiomics fig/NMP_fate_analysis/grn_venn_plot_gene_list/PM GRN tf degree.csv',
            quote = F,row.names = F,col.names = F)
write.table(sc_tf_all,'/data/jiangjunyao/polyATAC/polyatac_figure/multiomics fig/NMP_fate_analysis/grn_venn_plot_gene_list/SC GRN tf degree.csv',
            quote = F,row.names = F,col.names = F)

### overlap venn
load("~/public/transfac/Tranfac201803_Mm_MotifTFsF.rda")
pm_gene_all = pm_sp[pm_sp$cluster %in% 'NMP with Paraxial mesoderm bias fate',]
sc_gene_all = pm_sp[pm_sp$cluster %in% 'NMP with Spinal cord bias fate',]
pm_gene_all = pm_gene_all[,c(1,3)]
sc_gene_all = sc_gene_all[,c(1,3)]
colnames(pm_gene_all) = c('gene','degree centrality')
colnames(sc_gene_all) = c('gene','degree centrality')
pm_gene_all = pm_gene_all[order(pm_gene_all[,2],decreasing = T),][1:20,]
sc_gene_all = sc_gene_all[order(sc_gene_all[,2],decreasing = T),][1:20,]
pm_tf_all=pm_gene_all[pm_gene_all$gene %in% unlist(strsplit(Tranfac201803_Mm_MotifTFsF$TFs,';')),]
sc_tf_all=sc_gene_all[sc_gene_all$gene %in% unlist(strsplit(Tranfac201803_Mm_MotifTFsF$TFs,';')),]
list_venn_pm = list('early mesoderm regulate gene'=module_tf_pm,
                    'Top genes in GRN'=pm_tf_all$gene)
list_venn_sp = list('early neuron regulate gene'=module_tf_sp,
                    'Top genes in GRN'=sc_tf_all$gene)
venn(list_venn_sp,ilabels = "counts",box=F,zcolor = jdb_palette('corona'),ggplot = T,
     ilcs=1,sncs=1.2)
venn(list_venn_pm,ilabels = "counts",box=F,zcolor = jdb_palette('corona'),ggplot = T,
     ilcs=1,sncs=1.2)
intersect(list_venn_pm$`early mesoderm regulate gene`,list_venn_pm$`Top genes in GRN`)
intersect(list_venn_sp$`early neuron regulate gene`,list_venn_sp$`Top genes in GRN`)


pm_gene_all = pm_sp_plot[pm_sp_plot$type=='Paraxial mesoderm bias specific',1]
sc_gene_all = pm_sp_plot[pm_sp_plot$type=='Spinal cord bias specific',1]
pm_tf_all=pm_gene_all[pm_gene_all %in% unlist(strsplit(Tranfac201803_Mm_MotifTFsF$TFs,';'))]
sc_tf_all=sc_gene_all[sc_gene_all %in% unlist(strsplit(Tranfac201803_Mm_MotifTFsF$TFs,';'))]
list_venn_pm = list('early mesoderm regulate gene'=module_tf_pm,
                    'PXM NMP TFs'=pm_tf_all)
list_venn_sp = list('early neuron regulate gene'=module_tf_sp,
                    'SC NMP TFs'=sc_tf_all)
venn(list_venn_sp,ilabels = "counts",box=F,zcolor = jdb_palette('corona'),ggplot = T,
     ilcs=1,sncs=1.2)
venn(list_venn_pm,ilabels = "counts",box=F,zcolor = jdb_palette('corona'),ggplot = T,
     ilcs=1,sncs=1.2)
intersect(list_venn_pm$`early mesoderm regulate gene`,list_venn_pm$`PXM NMP TFs`)
intersect(list_venn_sp$`early neuron regulate gene`,list_venn_sp$`SC NMP TFs`)
