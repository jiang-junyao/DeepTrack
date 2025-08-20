library(tidyverse)
library(ggpubr)
library(BuenColors)
data1 <- read.csv("/data/jiangjunyao/polyATAC/seq Benchmark/Supplementary Table1. Summary of datasets in scLTdb 0915(1).csv", header=FALSE)
data2 = read.csv('/data/jiangjunyao/polyATAC/seq Benchmark/deeptrack_loxcode.csv',row.names = 1)
colnames(data1)[4] = 'tissue'
colnames(data1)[10] = 'barcode_number'
colnames(data1)[11] = 'barcoded_cell_rate'
colnames(data1)[6] = 'tech'
data1 = data1[,colnames(data2)]
data1 = rbind(data1,data2)
data1 = data1[!is.na(data1$barcoded_cell_rate),]



ggplot(data1, aes(y = reorder(tech,barcoded_cell_rate), x = barcoded_cell_rate)) +

  geom_jitter(alpha = 0.8, height = 0,color = "blue") +
  stat_summary(
    fun = "mean", geom = "bar", aes(fill = tissue),
    alpha = 0.7
  ) +
  stat_summary(
    fun.data = "mean_se", geom = "errorbar", width = 0.2,
    position = position_dodge(width = 0.9)
  ) +

  theme_minimal() +xlab('barcoded cell rate')+ylab('tech')+
  theme(text = element_text(size = 16))+
  scale_fill_manual(values = jdb_palette('brewer_spectra')[c(1,2,3,6,7,8)])
ggsave('/data/jiangjunyao/polyATAC/figure_after_science/bulk fig/seq tech benchmark.pdf',width = 10,height = 8)


