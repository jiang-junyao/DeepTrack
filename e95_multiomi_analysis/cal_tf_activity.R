library(chromVAR)
library(chromVARmotifs)
library(BSgenome.Mmusculus.UCSC.mm10)
run_chromvar <- function(atac_count,metadata,bsdb,species = 'Mus musculus',motifs){
  library(chromVAR)
  library(motifmatchr)
  library(SingleCellExperiment)
  sce <- SummarizedExperiment(assays = list(counts = atac_count),
                              colData = metadata,
                              rowRanges = GRanges(rownames(atac_count)))
  sce <- addGCBias(sce,genome = bsdb)
  motif_ix <- matchMotifs(motifs, sce,
                          genome = bsdb)
  BiocParallel::register(BiocParallel::SerialParam())
  dev <- computeDeviations(object = sce,
                           annotations = motif_ix)
  return(dev)
}
atac_count = obj@assays$Peak_ct@counts
metadata = obj@meta.data
peakname=str_split(rownames(atac_count),'-',simplify = T)
peakname2 = paste0(peakname[,1],':',peakname[,2],'-',peakname[,3])
rownames(atac_count) = peakname2
atac_count = atac_count[peakname[,1]%in%paste0('chr',1:22),]
data("mouse_pwms_v2")
dev_mt_cisbp = run_chromvar(atac_count,metadata,BSgenome.Mmusculus.UCSC.mm10,
                      motifs=mouse_pwms_v2)
obj[['TfActCisbp']] = CreateAssayObject(dev_mt_cisbp@assays@data$deviations)
saveRDS(dev_mt,'/data/jiangjunyao/polyATAC/multiomi/dev_mt_cisbp.rds')
