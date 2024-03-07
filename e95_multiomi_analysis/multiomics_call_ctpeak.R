library(Signac)
library(Seurat)
library(GenomicRanges)
obj <- readRDS("~/polyATAC/multiomi/multiomics_seurat.rds")
annotation <- readRDS('/data/jiangjunyao/public/Mmv79_annotations.rds')
source('/data/jiangjunyao/polyATAC/tool/normalize_peak_width.R')
DefaultAssay(obj) = 'ATAC'
peaks <- CallPeaks(
  object = obj,
  group.by = "celltype",macs2.path = '/data/jiangjunyao/miniconda3/envs/normal/bin/macs2'
)
obj@active.ident = as.factor(obj$celltype)
peaks2 = subdivideGRanges(peaks,core = 8)
peaks2_counts <- FeatureMatrix(
  fragments = Fragments(obj),
  features = peaks2,
  cells = colnames(obj)
)
obj[["Peak_ct"]] <- CreateChromatinAssay(
  counts = peaks2_counts,
  fragments = "/data/jiangjunyao/polyATAC/multiomi/atac_fragments.tsv.gz",
  annotation = annotation
)
DefaultAssay(obj) = "Peak_ct"
obj = RunTFIDF(obj)
saveRDS(obj,"~/polyATAC/multiomi/multiomics_seurat.rds")
