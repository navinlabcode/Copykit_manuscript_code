load("~/Projects/copykit/benchmark/copykit/20230217_knn/20230217_knn_copykitobj.rds")

mean(brcan_filt_knn$reads_total)
mean(brcan_filt_knn$percentage_duplicates)
# 9.35%
mean(co5_filt_knn$reads_total)
mean(co5_filt_knn$percentage_duplicates)
mean(tn7_filt_knn$reads_total)
mean(tn7_filt_knn$percentage_duplicates)
mean(co8_filt_knn$reads_total)
mean(co8_filt_knn$percentage_duplicates)
mean(pmtc1_filt_knn$reads_total)
mean(pmtc1_filt_knn$percentage_duplicates)
mean(pmtc6_filt_knn$reads_total)
mean(pmtc6_filt_knn$percentage_duplicates)
mean(pmtc7_filt_knn$reads_total)
mean(pmtc7_filt_knn$percentage_duplicates)


slected_name <- brcan_filt_knn$sample[brcan_filt_knn$reads_total - brcan_filt_knn$reads_duplicates>1000000]

pmtc1_breast <- runVarbin("/volumes/seq/projects/CNA_projects/DT_CNA/FFPE/Breast/PMTC1_PM2040_PM2054/output/sort/", genome = "hg19", is_paired_end = F)
pmtc1_breast <- findAneuploidCells(pmtc1_breast)
pmtc1_breast_filt <- pmtc1_breast[, colData(pmtc1_breast)$is_aneuploid == TRUE]
# pmtc1_breast_filt_knn <- knnSmooth(pmtc1_breast_filt)
# pmtc1_breast_filt_knn <- calcInteger(pmtc1_breast_filt_knn, method = 'scquantum', assay = 'smoothed_bincounts')

mean(pmtc1_breast_filt$reads_total)
mean(pmtc1_breast_filt$percentage_duplicates)
slected_name_pmtc1 <- pmtc1_breast_filt$sample[pmtc1_breast_filt$reads_total - pmtc1_breast_filt$reads_duplicates>750000]

### making files for downsampling
sapply(slected_name, function(x){
  cmd <- glue("ln -s ~/Projects/copykit/ginkgo/brcan/bam/{x}.bam ~/Projects/copykit/downsampling/data/brcan/fulldata/")
  system(cmd)
})

### making files for downsampling
sapply(slected_name_pmtc1, function(x){
  cmd <- glue("ln -s /volumes/seq/projects/CNA_projects/DT_CNA/FFPE/Breast/PMTC1_PM2040_PM2054/output/sort/{x}.bam ~/Projects/copykit/downsampling/data/pmtc1breast/fulldata/")
  system(cmd)
})
