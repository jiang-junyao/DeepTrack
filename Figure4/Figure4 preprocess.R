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
rna_metadata = read.csv('/data/jiangjunyao/polyATAC/multiomi/polyatac_rna_meta_20250604.csv',
                        row.names = 1)
setwd('/data/jiangjunyao/polyATAC/multiomi/archr_20250604')
addArchRThreads(threads = 40)
addArchRGenome("mm10")
frag_list = c('/data/jiangjunyao/polyATAC/multiomi/atac_fragments.tsv.gz',
              '/data/jiangjunyao/polyATAC/multiomi/mouse2_raw/mouse2_rep1/atac_fragments.tsv.gz',
              '/data/jiangjunyao/polyATAC/multiomi/mouse2_raw/mouse2_rep2/atac_fragments.tsv.gz',
              '/data/jiangjunyao/polyATAC/multiomi/mouse2_raw/mouse2_rep3/atac_fragments.tsv.gz',
              '/data/jiangjunyao/polyATAC/multiomi/mouse2_raw/mouse2_rep4/atac_fragments.tsv.gz',
              '/data/jiangjunyao/polyATAC/20250603_cellranger/tail1/atac_fragments.tsv.gz',
              '/data/jiangjunyao/polyATAC/20250603_cellranger/tail2/atac_fragments.tsv.gz')
names(frag_list) = c('mouse1','mouse2_rep1','mouse2_rep2','mouse2_rep3',
                     'mouse2_rep4','tail1','tail2')

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
  outputDirectory = "/data/jiangjunyao/polyATAC/multiomi/archr_20250604",
  copyArrows = TRUE
)

### subset obj & call peaks
cell_used = intersect(rna_metadata$cell,rownames(proj@cellColData))
proj_used = subsetArchRProject(proj,cells = cell_used,
                   outputDirectory = '/data/jiangjunyao/polyATAC/multiomi/archr_subset_20250604'
                   ,force = T)
rownames(rna_metadata) = rna_metadata$cell
proj_used@cellColData$celltype = rna_metadata[rownames(proj_used@cellColData),
                                               'celltype']
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

obj = readRDS("/data/jiangjunyao/polyATAC/polyATAC_20250604.rds")
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
saveRDS(proj_used,'/data/jiangjunyao/polyATAC/multiomi/archr_subset_20250604/Save-ArchR-Project.rds')
saveRDS(obj@assays$Peak_ct,'/data/jiangjunyao/polyATAC/multiomi/merge_wolf/peak_assay.rds')
saveRDS(obj,"/data/jiangjunyao/polyATAC/polyATAC_20250604.rds")
DefaultAssay(obj) = 'Peak_ct'
sceasy::convertFormat(obj, from="seurat", to="anndata",assay='Peak_ct',
                      outFile='/data/jiangjunyao/polyATAC/polyatac_peak_20250604.h5ad')
getGroupBW(
  ArchRProj = proj_used,
  groupBy = "celltype",
  normMethod = "ReadsInTSS",
  tileSize = 40,
  maxCells = 2000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)
