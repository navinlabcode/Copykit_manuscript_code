library(devtools)
library(copykit)
library(BiocParallel)
library(glue)
library(dplyr)
library(ggplot2)

register(MulticoreParam(workers = 50), default = T)
list.files("/mnt/USR1/junke/Projects/copykit/rerun_manu/singu_0424/")
mat_gt <- matrix(nrow = 50, ncol = 50)
mat_result <- matrix(nrow = 50, ncol = 50)

for (i in 1:50) {

  obj <- mock_bincounts(ncells = 50,
                     ncells_diploid = i,
                     position_gain = 4900:5493,
                     position_del = 6523:7056,
                     resolution = '220kb')

  obj <- runVst(obj)
  obj <- runSegmentation(obj)
  obj <- findAneuploidCells(obj)

  mat_gt[,i] <- colData(obj)$ground_truth
  mat_result[,i] <- colData(obj)$is_aneuploid

}



precision <- vector()
recall <- vector()
sensitivity <- vector()
specificity <- vector()

for (i in 1:ncol(mat_gt)) {
  conf <- NULL

  mat_gt_res <- factor(mat_gt[,i])
  levels(mat_gt_res) <- c(TRUE, FALSE)

  mat_result_res <- factor(mat_result[,i])
  levels(mat_result_res) <- c(TRUE, FALSE)

  conf <- table(mat_gt_res, mat_result_res)

  precision[i] <- conf[1]/(conf[1]+conf[2])
  recall[i] <- conf[1]/(conf[1]+conf[3])
  sensitivity[i] <- conf[1]/(conf[1]+conf[3])
  specificity[i] <- conf[4]/(conf[4]+conf[2])


}

precision_whole_chr <- precision
recall_whole_chr <- recall

# one chr gain

mat_gt <- matrix(nrow = 50, ncol = 50)
mat_result <- matrix(nrow = 50, ncol = 50)

for (i in 1:50) {

  obj <- mock_bincounts(ncells = 50,
                        ncells_diploid = i,
                        position_gain = 4900:5132,
                        position_del = NULL,
                        resolution = '220kb')

  obj <- runVst(obj)
  obj <- runSegmentation(obj)
  obj <- findAneuploidCells(obj)

  mat_gt[,i] <- colData(obj)$ground_truth
  mat_result[,i] <- colData(obj)$is_aneuploid

}

precision <- vector()
recall <- vector()
sensitivity <- vector()
specificity <- vector()

for (i in 1:ncol(mat_gt)) {
  conf <- NULL

  mat_gt_res <- factor(mat_gt[,i])
  levels(mat_gt_res) <- c(TRUE, FALSE)

  mat_result_res <- factor(mat_result[,i])
  levels(mat_result_res) <- c(TRUE, FALSE)

  conf <- table(mat_gt_res, mat_result_res)

  precision[i] <- conf[1]/(conf[1]+conf[2])
  recall[i] <- conf[1]/(conf[1]+conf[3])
  sensitivity[i] <- conf[1]/(conf[1]+conf[3])
  specificity[i] <- conf[4]/(conf[4]+conf[2])


}

precision_parm <- precision
recall_parm <- recall


# plot

df_whole_chr  <- data.frame(precision = precision_whole_chr+0.015,
                            recall = recall_whole_chr+0.015,
                            cell = 1:50,
                            data = "gain7del10")

df_parm <- data.frame(precision = precision_parm-0.015,
                      recall = recall_parm-0.015,
                      cell = 1:50,
                      data = "gain7p")

df <- bind_rows(df_whole_chr,
                df_parm)

save(list = ls(all.names = T), file = "/mnt/USR1/junke/Projects/copykit/rerun_manu/singu_0424/simu.Rdata")

## Visualizations
library(cowplot)
load("/mnt/USR1/junke/Projects/copykit/rerun_manu/singu_0424/simu.Rdata")
p_precision <- ggplot(df, aes(cell, precision, color = data)) +
  geom_point() +
  cowplot::theme_cowplot() +
  scale_color_viridis_d(option = 'turbo') +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  labs(x = 'number of euploid cells',
       color = 'copy number') +
  ylim(c(0, 1.02))

p_recall <- ggplot(df, aes(cell, recall, color = data)) +
  geom_point() +
  cowplot::theme_cowplot() +
  scale_color_viridis_d(option = 'turbo') +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  labs(x = 'number of euploid cells',
       color = 'copy number') +
  ylim(c(0, 1.02))

p_patch <- p_precision + p_recall + patchwork::plot_layout(guides = 'collect')


cowplot::ggsave2("/mnt/USR1/junke/Projects/copykit/rerun_manu/singu_0424/p_precision_recall_findAneu.pdf", width = 10, height = 3)

col_fun=circlize::colorRamp2(breaks = c(-1,0,1),
                             c("blue", "white", "red"))
pdf("/mnt/USR1/junke/Projects/copykit/rerun_manu/singu_0424/p_findAneuSimul_parm.pdf", width = 8, height = 8)
plotHeatmap(obj, label = 'is_aneuploid', col = col_fun)
dev.off()
