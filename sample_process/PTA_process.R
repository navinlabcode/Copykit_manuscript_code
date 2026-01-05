library(copykit)
library(BiocParallel)
set.seed(17)
register(MulticoreParam(workers = 100), default = T)


ck_md231 <- runVarbin("~/Projects/copykit/MC_revision/PTA_data/MDA231/marked/", remove_Y = T, is_paired_end = T)

pdf("~/Projects/copykit/MC_revision/PTA_data/md231pta_ht.pdf", width = 7, height = 5)
plotHeatmap(ck_md231, order_cells = 'hclust', n_threads = 50)
dev.off()

ck_md231 <- runMetrics(ck_md231)

colData(ck_md231)

ck_pta <- runVarbin("~/Projects/copykit/MC_revision/PTA_data/pnas_paper/marked/", remove_Y = T, is_paired_end = T)
pdf("~/Projects/copykit/MC_revision/PTA_data/gmpta_ht.pdf", width = 7, height = 5)
plotHeatmap(ck_pta, order_cells = 'hclust', n_threads = 50)
dev.off()

saveRDS(ck_pta, "~/Projects/copykit/MC_revision/PTA_data/pnas_paper/ckobj_raw.rds")

ck_pta <- findOutliers(ck_pta, resolution = 0.25)
# Marked 6 cells as outliers.
ck_pta <- runMetrics(ck_pta)
pdf("~/Projects/copykit/MC_revision/PTA_data/gmpta_isoutlier_ht.pdf", width = 6, height = 3)
plotHeatmap(ck_pta, order_cells = 'hclust', label = "outlier", row_split = "outlier", n_threads = 50)
dev.off()



