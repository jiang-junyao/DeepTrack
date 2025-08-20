library(BuenColors)
library(ggplot2)
library(stats)
prob_final = list()
for (i in seq(100, 10000, by = 100)) {
  N <- i
  k <- 1
  threshold <- 0.995
  p_list <- seq(1e-6, 0.01, length.out = 10000)

  prob_list <- sapply(p_list, function(p) {
    pbinom(k, size = N, prob = p)

  })
  df1 <- data.frame(
    Probability = p_list,

    Cumulative_Prob = prob_list

  )
  target_values <- c(0.9, 0.95, 0.975)
  closest_indices <- sapply(target_values, function(target) {
    which.min(abs(df1$Cumulative_Prob - target))
  })
  df1$Cumulative_Prob[closest_indices] <- target_values
  df1$cell_number = i
  prob_final[[as.character(i)]] = df1
}
prob_df = do.call(bind_rows,prob_final)


prob_df_use = prob_df[prob_df$Cumulative_Prob %in% c(0.9,0.95,0.975),]
prob_df_use$Cumulative_Prob = as.character(prob_df_use$Cumulative_Prob)
ggplot(prob_df_use, aes(x = cell_number, y = Probability, color = Cumulative_Prob)) +
  geom_smooth(method = "lm", se = FALSE,size=1.5) +
  scale_x_log10() +
  scale_y_log10() +
  geom_vline(xintercept = 500, linetype = "dashed", color = "black",size=1) +
  geom_hline(yintercept = 0.001, linetype = "dashed", color = "black",size=1) +
  theme_classic()+theme(text = element_text(size=16))+
  scale_color_manual(values = jdb_palette('corona'))+
  ylab('Barcode generation probability')+
  xlab('Number of precursor cells barcoded')
ggsave('/data/jiangjunyao/polyATAC/figure_after_science/bulk fig/pgen prob line plot 500.pdf',
       width = 10,height = 8)
