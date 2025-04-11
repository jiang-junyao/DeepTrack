library(dplyr)
library(tidyr)
library(readr)
library(GenomicRanges)
library(reshape2)
library(Signac)
library(openxlsx)
library(ggplot2)
library(Matrix)
library(EpiTrace)
library(Seurat)
library(SeuratObject)
library(ggtree)
library(patchwork)
library(ArchR)
library(parallel)
library(ChIPseeker)
library(ArchR)
### set file
rna_metadata = read.csv('/data/jiangjunyao/polyATAC/multiomi/polyatac_rna_meta.csv',
                        row.names = 1)
setwd('/data/jiangjunyao/polyATAC/multiomi/archr')
addArchRThreads(threads = 20)
addArchRGenome("mm10")
frag_list = c('/data/jiangjunyao/polyATAC/multiomi/atac_fragments.tsv.gz',
              '/data/jiangjunyao/polyATAC/multiomi/mouse2_raw/mouse2_rep1/atac_fragments.tsv.gz',
              '/data/jiangjunyao/polyATAC/multiomi/mouse2_raw/mouse2_rep2/atac_fragments.tsv.gz',
              '/data/jiangjunyao/polyATAC/multiomi/mouse2_raw/mouse2_rep3/atac_fragments.tsv.gz',
              '/data/jiangjunyao/polyATAC/multiomi/mouse2_raw/mouse2_rep4/atac_fragments.tsv.gz')
names(frag_list) = c('mouse1','mouse2_rep1','mouse2_rep2','mouse2_rep3',
                     'mouse2_rep4')

### create arrow
ArrowFiles <- createArrowFiles(
  inputFiles = frag_list,
  sampleNames = names(frag_list),
  minTSS = 4,
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force=T
)
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "/data/jiangjunyao/polyATAC/multiomi/archr",
  copyArrows = TRUE
)

### subset obj & call peaks
cell_used = intersect(rna_metadata$cell,rownames(proj@cellColData))
proj_used = subsetArchRProject(proj,cells = cell_used,
                   outputDirectory = '/data/jiangjunyao/polyATAC/multiomi/archr_subset'
                   ,force = T)
rownames(rna_metadata) = rna_metadata$cell
proj_used@cellColData$celltype = rna_metadata[rownames(proj_used@cellColData),
                                               'ct_leiden']
proj_used <- addGroupCoverages(ArchRProj = proj_used, groupBy = 'celltype',
                               force = T)
proj_used <- addReproduciblePeakSet(
  ArchRProj = proj_used,
  groupBy = "celltype",
  pathToMacs2 = '/data/jiangjunyao/miniconda3/envs/normal/bin/macs2',
  maxPeaks = 100000,cutOff = 0.0001,force = T
)
proj_used <- addPeakMatrix(proj_used)
peak_mt = getMatrixFromProject(proj_used,"PeakMatrix",threads = 20)
peak_mt = peak_mt@assays@data$PeakMatrix
rownames(peak_mt) = paste0(proj_used@peakSet@seqnames,':',
                           proj_used@peakSet@ranges)
### add to object

obj = readRDS("~/polyATAC/polyatac_V2.rds")
obj = subset(obj,cell %in% rownames(proj_used@cellColData))
obj[['Peak_ct']] <- CreateChromatinAssay(
  counts = peak_mt,
  sep = c(":", "-"),
  annotation = readRDS('/data/jiangjunyao/public/Mmv79_annotations.rds')
)
### output
library(reticulate)
use_python("/data/jiangjunyao/miniconda3/envs/normal/bin/python")
saveArchRProject(ArchRProj = proj_used, outputDirectory = '/data/jiangjunyao/polyATAC/multiomi/archr_subset', load = FALSE)
saveRDS(proj_used,'/data/jiangjunyao/polyATAC/multiomi/merge_wolf_archr_subset/Save-ArchR-Project.rds')
saveRDS(obj@assays$Peak_ct,'/data/jiangjunyao/polyATAC/multiomi/merge_wolf/peak_assay.rds')
saveRDS(obj,"~/polyATAC/polyatac_V2.rds")
DefaultAssay(obj) = 'Peak_ct'
sceasy::convertFormat(obj, from="seurat", to="anndata",assay='Peak_ct',
                      outFile='/data/jiangjunyao/polyATAC/polyatac_peak.h5ad')
