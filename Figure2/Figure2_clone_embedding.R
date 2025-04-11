library(Seurat)
library(Signac)
library(ggrepel)
clone_embedding <- function(barcode_use,meta,coor,mode='single',barcode_key='barcodes',
                            colors=jdb_palette('corona')[-8],bg_size=1.5,barcode_size=4,
                            ct_key='celltype'){
  coor = as.data.frame(coor)
  barcode_anno = rep('other',nrow(meta))
  idx <- which(meta[,barcode_key] %in% barcode_use)
  barcode_anno[idx] <- meta[idx,barcode_key]
  coor$barcode_type = barcode_anno
  colnames(coor)[1:2] = c('UMAP_1','UMAP_2')
  coor$barcode_type = barcode_anno
  coor1 = coor[coor$barcode_type %in% c('other'),]
  coor1$barcode_type = NA
  coor1$border=NA
  coor2 = coor[!coor$barcode_type %in% c('other'),]
  coor2$border='with barcode'
  if (mode=='single') {
    coor2$barcode_type = meta[rownames(coor2),ct_key]
  }
  
  p1=ggplot(coor1,aes(x=UMAP_1,y=UMAP_2,fill=barcode_type,color=border))+geom_point(data=coor1,size=bg_size,shape=21)+
    geom_point(data=coor2,aes(x=UMAP_1,y=UMAP_2),shape=21,size=barcode_size)+theme_void()+
    scale_color_manual(values = 'black', na.value = rgb(215/255,215/255,215/255))+
    theme(text = element_text(size=16))+
    scale_fill_manual(values = colors, na.value = rgb(215/255,215/255,215/255))+
    guides(color="none")+theme(legend.position = "none")
  return(p1)
}
all_sclt_V3 <- readRDS("~/polyATAC/all_sclt_V3.rds")
germ_col = c('#98DF8AFF',
             '#DBDB8DFF',
             '#AEC7E8FF',
             '#C5B0D5FF',
             '#FFBB78FF')
names(germ_col) = c('NeuroEctoderm','SurfaceEctoderm','Mesoderm',
         'Endoderm','Extraembryonic')
top_barcodes = group_df[group_df$size>=10 & group_df$size<=20,]
top_barcodes = top_barcodes %>%
  group_by(fate)  %>%
  slice_max(order_by = size, n = 10) %>%
    ungroup()
for (i in 1:nrow(top_barcodes)) {
  id1 = paste(top_barcodes$barcodes[i],top_barcodes[i,2])
  p1=clone_embedding(top_barcodes$barcodes[i],all_sclt_V3@meta.data,
                     all_sclt_V3@reductions$umap@cell.embeddings,mode = 'single',
                     barcode_key = 'barcodes',barcode_size = 5,
                     ct_key = 'germ',colors = germ_col)+
    ggtitle(id1)
  print(p1)
  ggsave(paste0('/data/jiangjunyao/polyATAC/polyatac_figure/single cell fig/single cell embedding/clone embedding mid/',
                id1,'.pdf'),width = 10,height = 8)
}
