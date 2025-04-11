library(ggplot2)
library(tidyverse)
library(BuenColors)
cag_polyexpress_meta <- read.csv("E:/polyATAC/all_sc/cag_polyexpress_meta.csv", row.names=1)
cag_polyexpress_meta = cag_polyexpress_meta[!is.na(cag_polyexpress_meta$ct_leiden),]
plot_sclt_statistic = function(metadata,idx='ct_leiden'){
  all_diversity = c()
  all_detection = c()
  all_number = c()
  for (i in unique(metadata[,idx])) {
    meta_use = metadata[metadata[,idx]==i,]
    all_number = c(all_number,nrow(meta_use))
    all_detection = c(all_detection,
                      nrow(meta_use[meta_use$detection=='with barcode',])/nrow(meta_use))
    meta_use = meta_use[!is.na(meta_use$barcodes),]
    all_diversity = c(all_diversity,length(unique(meta_use$barcodes)))
  }
  barcode_statistic = data.frame(all_diversity,all_detection,all_number,
                                 unique(metadata[,idx]))
  colnames(barcode_statistic) = c('unique_barcode_number','detection_rate',
                                  'cell_number','celltype')
  barcode_statistic = barcode_statistic[order(barcode_statistic$unique_barcode_number),]
  barcode_statistic$order = 1:nrow(barcode_statistic)
  ggplot(barcode_statistic,aes(x=unique_barcode_number,y=reorder(celltype,order)))+
    geom_point(aes(size=cell_number,color=detection_rate))+
    theme_minimal()+scale_color_gradientn(colors = jdb_palette("solar_basic"))+
    xlab('Unique barcode number')+ylab('Cell types')+
    theme(text = element_text(size=15))
}
plot_sclt_statistic(cag_polyexpress_meta,'celltype')
ggsave('E:\\polyATAC\\single cell fig\\cell type statistic and fraction plot/celltype barcode statistic.pdf',
       width = 10,height = 8)
plot_sclt_statistic(polyATAC_v2_metadata,'celltype')
###
plot_clone_size = function(metadata,idx){
  size_list = list()
  for (i in unique(metadata[,idx])) {
    meta_use = metadata[metadata[,idx]==i,]
    meta_use = meta_use[!is.na(meta_use$barcodes),]
    meta_use = meta_use[meta_use$pgen<0.001,]
    clone_size = as.data.frame(table(meta_use$barcodes))
    clone_size[,1] = as.character(clone_size[,1])
    clone_size$fraction = clone_size[,2]/sum(clone_size[,2])
    if (max(clone_size[,3])<0.05 | max(clone_size[,2])<3) {
      clone_size[,1] = 'other clones'

    }else{
      clone_size[clone_size$fraction<0.05,1] = 'other clones'
      clone_size[clone_size$Freq<5,1] = 'other clones'

    }
    clone_size$celltype = i
    size_list[[i]] = clone_size
  }
  size_df = do.call(bind_rows,size_list)
  ggplot(size_df,aes(y=celltype,x=fraction,fill=Var1))+geom_col()+theme_minimal()+
    scale_fill_manual(values = c(jdb_palette('corona')[1:(length(unique(size_df$Var1))-1)],'grey'))+
    theme(text = element_text(size=15))
}
plot_clone_size(cag_polyexpress_meta,'celltype')
### top clone size
singlecell_clone_embedding = function(rna,barcode_use){
  germ=c(rgb(25/255,90/255,53/255),
         rgb(119/255,42/255,34/255),
         rgb(41/255,112/255,160/255),
         rgb(153/255,125/255,37/255))
  names(germ)=c('NeuroEctoderm','SurfaceEctoderm','Mesoderm',
                'Endoderm')
  library(tidyverse)
  rna@meta.data = rna@meta.data[!is.na(rna@meta.data$celltype),]
  barcode_anno = rep('other cells',nrow(rna@meta.data))
  idx <- which(rna$barcodes %in% barcode_use)
  barcode_anno[idx] <- rna$barcodes[idx]
  coor = as.data.frame(rna@reductions$umap@cell.embeddings)
  coor$barcode_type = barcode_anno
  colnames(coor)[1:2] = c('UMAP_1','UMAP_2')
  coor$border = NA
  coor1 = coor[coor$barcode_type ==c('other cells'),]
  coor1$barcode_type = NA
  coor2 = coor[!coor$barcode_type %in% c('other cells'),]
  coor2$border=barcode_use
  coor2$barcode_type = rna@meta.data[rownames(coor2),'germ']
  coor2$barcode_type[coor2$barcode_type=='mesoderm'] = 'Mesoderm'
  coor2$barcode_type[coor2$barcode_type=='Surface ectoderm'] = 'SurfaceEctoderm'
  coor2$barcode_type[coor2$barcode_type=='Neuroectoderm'] = 'NeuroEctoderm'
  ggplot(coor1,aes(x=UMAP_1,y=UMAP_2,fill=barcode_type,color=border))+geom_point(data=coor1,size=1.5,shape=21)+
    scale_color_manual(values =  rgb(215/255,215/255,215/255))+
    geom_point(data=coor2,aes(x=UMAP_1,y=UMAP_2),shape=21,size=4)+theme_void()+
    scale_color_manual(values = 'black', na.value = rgb(215/255,215/255,215/255))+
    theme(text = element_text(size=16))+
    scale_fill_manual(values = germ, na.value = rgb(215/255,215/255,215/255))+
    guides(color="none")+ggtitle(barcode_use)+theme(plot.title = element_text(hjust = 0.5))
}

cellfate = meta[!is.na(meta$barcodes),]
cellfate = cellfate[cellfate$pgen<0.001,]
clone_size = as.data.frame(table(cellfate$barcodes))
clone_size = clone_size[order(clone_size$Freq,decreasing = T),]
for (i in 1:20) {
  singlecell_clone_embedding(rna,as.character(clone_size[i,1]))
  ggsave(paste0('E:\\polyATAC\\single cell fig\\single cell embedding\\clone embedding/',
                as.character(clone_size[i,1]),'.pdf'),width = 10,height = 8)
}

##
sc_group_df = readRDS('E:\\polyATAC\\all_barcode_kmeans_V2.rds')
sc_group_df = sc_group_df[!is.na(sc_group_df$mouse),]
sc_group_df[sc_group_df$fate=='Endoderm-restricted','fate'] = 'Endoderm-bias'
sc_group_df[sc_group_df$fate=='Mesoderm-restricted','fate'] = 'Mesoderm-bias'
sc_group_df[sc_group_df$fate=='NeuroEctoderm-restricted','fate'] = 'NeuroEctoderm-bias'
sc_group_df[sc_group_df$fate=='SurfaceEctoderm-restricted','fate'] = 'SurfaceEctoderm-bias'
sc_group_df$type = '3.1< clone size <20'
sc_group_df$type[sc_group_df$size>20] = '2.20 < clone szie < 40'
sc_group_df$type[sc_group_df$size>40] = '1.clone szie > 40'
sc_group_df$type[sc_group_df$size==1] = '4.clone szie =1'
ggplot(sc_group_df,aes(y=fate,x=size,fill=type))+geom_col()+theme_minimal()+
  theme(text = element_text(size=15))+xlab('cell number')+scale_fill_manual(values = jdb_palette('corona'))
ggsave('E:\\polyATAC\\all_sc\\figures/clone_size.pdf',width = 12,height = 6)

## per mouse clone analysis
mouse_clone_num = as.data.frame(table(sc_group_df$mouse))
ggplot(mouse_clone_num,aes(y=Var1,x=Freq))+geom_col()+theme_minimal()+
  theme(text = element_text(size=18))+ylab('mouse')+xlab('clone number')

### clone size
size_df=as.data.frame(table(sc_group_df$size))
size_df[,1] = as.numeric(size_df[,1])
ggplot(size_df,aes(x=Var1,y=Freq))+geom_col()+theme_minimal()+
  theme(text = element_text(size=18))+xlab('clone size')+ylab('frequency')+scale_y_log10()

### clone size cell type number
fatebias=c('Endoderm-bias'='#C5B0D5FF',
           'Mesoderm-bias'='#AEC7E8FF',
           'NeuroEctoderm-bias'='#98DF8AFF',
           'SurfaceEctoderm-bias'='#DBDB8DFF',
           'Multilineage'=rgb(242/255,150/255,155/255))
all_ct_num = c()
for (i in 1:nrow(sc_group_df)) {
  barcode_use = sc_group_df[i,'barcodes']
  ct_num = length(unique(cag_polyexpress_meta[cag_polyexpress_meta$barcodes==barcode_use,'celltype']))
  all_ct_num = c(all_ct_num,ct_num)
}
sc_group_df$celltype_num = all_ct_num
ggplot(sc_group_df,aes(x=celltype_num,y=size,color=fate))+geom_point(size=2)+theme_minimal()+
  xlab('clone related cell type number')+ylab('clone size')+
  scale_color_manual(values = fatebias)+theme(text = element_text(size=14))+
  ggtitle('spearman correlation 0.91')
ggplot(sc_group_df,aes(x=))
