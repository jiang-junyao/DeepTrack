group_df = readRDS('E:\\polyATAC\\bulk\\barcode_kmeans_20240925.rds')
group_df$mouse = str_split(rownames(group_df),'_',simplify = T)[,2]
sc_group_df = readRDS('E:\\polyATAC\\all_barcode_kmeans_V2.rds')
sc_group_df = sc_group_df[!is.na(sc_group_df$mouse),]

downsample_clone_fate=function(group_df,outdir){
  library(BuenColors)
  anno_cols = c('Endoderm-restricted'=jdb_palette('brewer_spectra')[6],
                              'Endoderm-bias'=jdb_palette('brewer_spectra')[7],
                              'Mesoderm-bias'=jdb_palette('brewer_spectra')[2],
                              'NeuroEctoderm-bias'=jdb_palette('brewer_spectra')[4],
                              'SurfaceEctoderm-bias'=jdb_palette('brewer_spectra')[8],
                              "Extraembryonic-bias"='#f2eada',
                              'Multilineage'="#855C59")
  anno_cols = anno_cols[names(anno_cols) %in% group_df$fate]
  all_fate = names(anno_cols)
  sample_list = list()
  for (i in c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) {
    set.seed(114)
    idx = sample(1:nrow(group_df),nrow(group_df)*i)
    group_df_use = group_df[idx,]
    freq = as.data.frame(table(group_df_use$fate))
    exclude_fate = all_fate[!all_fate %in% freq[,1]]
    exclude_df = data.frame(exclude_fate,rep(0,length(exclude_fate)))
    colnames(exclude_df) = colnames(freq)
    freq = rbind(freq,exclude_df)
    freq$ratio = freq[,2]/sum(freq[,2])
    freq$sample_percent = as.character(i)
    sample_list[[as.character(i)]] = freq
  }
  sample_df = do.call(bind_rows,sample_list)
  colnames(sample_df)[1] = 'fate'
  ggplot(sample_df,aes(x=sample_percent,y=ratio,fill=fate))+geom_col()+
    theme_classic()+scale_y_continuous(expand = c(0,0))+xlab('clone percentage')+
    ylab('fate percentage')+theme(text = element_text(size=16))+
    scale_fill_manual(values = anno_cols)
  ggsave(paste0(outdir,'_barplot.pdf'),width = 10,height = 8)
  sample_df$sample_percent = as.numeric(sample_df$sample_percent)
  ggplot(sample_df,aes(x=sample_percent,y=ratio,color=fate))+geom_line(size=1)+
    theme_classic()+scale_color_manual(values = anno_cols)+xlab('clone percentage')+
    ylab('fate percentage')+theme(text = element_text(size=16))
  ggsave(paste0(outdir,'_lineplot.pdf'),width = 10,height = 8)
  write.csv(sample_df,paste0(outdir,'_summary_table.csv'))
}
### bulk
downsample_clone_fate(group_df,
                      'E:\\polyATAC\\downsample_epiblast_fatebias/bulk_mouse_all')
downsample_clone_fate(group_df[group_df$mouse=='2',],
                      'E:\\polyATAC\\downsample_epiblast_fatebias/bulk_mouse_2')
downsample_clone_fate(group_df[group_df$mouse=='3',],
                      'E:\\polyATAC\\downsample_epiblast_fatebias/bulk_mouse_3')

### single cell
downsample_clone_fate(sc_group_df,
                      'E:\\polyATAC\\downsample_epiblast_fatebias/sc_mouse_all')
for (i in unique(sc_group_df$mouse)) {
  group_df_use = sc_group_df[sc_group_df$mouse==i,]
  downsample_clone_fate(group_df_use,
                        paste0('E:\\polyATAC\\downsample_epiblast_fatebias/sc_',i))
}
