load("~/Projects/copykit/benchmark/ginkgo/mda231/20230828_knn_gingko_compare.Rdata")
library(tidyverse)
library(BiocParallel)
library(glue)
library("distances")
set.seed(123)
register(MulticoreParam(workers = 100, progressbar = F), default = T)
outputdir <- "/volumes/USR1/junke/Projects/copykit/benchmark/sc_hmmcopy/mda231/"
chrlevels <- sapply(c(1:22,"X","Y"), function(x){paste0("chr",x)})

## mda231
########
mda231_seg_schmm <-read_csv("~/Projects/copykit/benchmark/sc_hmmcopy/mda231_c50/results/reads.csv.gz")
mda231_hmm_cnmat <- mda231_seg_schmm %>% dplyr::select(start, end, chr, state, cell_id) %>%
  spread(key = cell_id, value = state)
mda231_hmm_cnmat$chr <- factor(mda231_hmm_cnmat$chr, levels = chrlevels)
mda231_hmm_cnmat <- mda231_hmm_cnmat %>% arrange(chr, start)

mda231_hmm_ploidy <- colMeans(mda231_hmm_cnmat[,-c(1:3)])

segment_range <- mda231_hmm_cnmat %>% dplyr::select(chr, start, end)
colnames(segment_range) <- c("seqnames", "start","end")
segment_range$seqnames <- as.character(segment_range$seqnames)

mda231_hmm_cnmat_filt <- mda231_hmm_cnmat[, colnames(mda231_hmm_cnmat)%in%sapply(colData(mda231_filt_knn)$sample, function(x){unlist(strsplit(x,".", fixed=T))[1]})]
hc <- fastcluster::hclust(distance_matrix(distances(t(mda231_hmm_cnmat_filt))), method = "ward.D2")
pdf(glue("{outputdir}/mda231_hmm_cn_filt.pdf"), height = 4, width = 8)
print(plotHeatmap_gk(mda231, df=mda231_hmm_cnmat_filt[,hc$order], chr_ranges = segment_range))
dev.off()

#### Compare CV
cv <- function(seg){vapply(
  seg, function(z) {
    sd(z) / mean(z)
  },
  numeric(1)
)}

mda231_cv <- data.frame(
  ginkgo = cv(mda231_cnmat[,colnames(mda231_hmm_cnmat_filt)]), 
  hmmcopy = cv(mda231_hmm_cnmat_filt),
  copykit = cv(assay(mda231_filt_knn,"integer"))
) %>%
  rownames_to_column("cellname") %>%
  mutate(sample = "mda231")

cv_plot_df <- mda231_cv %>% dplyr::select(cellname, sample, ginkgo,copykit,hmmcopy) %>%
  gather(value="cv", key="tool", -sample, -cellname)

my_comparisons <- list( c("ginkgo", "copykit"), c("ginkgo", "hmmcopy"), c("copykit", "hmmcopy") )

ggpubr::ggviolin(cv_plot_df,
                 x = "tool", y = "cv",
                 color = "sample", palette = "jco", add="boxplot") + 
  ggpubr::stat_compare_means(comparisons = my_comparisons,label = "p.signif")+ # Add pairwise comparisons p-value
  # ggpubr::stat_compare_means(label.y = 50)  +
  theme(axis.text.x = element_text(angle = 45, hjust=0.8),
        legend.position = "none") +
  facet_grid(~sample) 
cowplot::ggsave2(glue("{outputdir}/cv_compare_violin_mda231.pdf"),width = 3, height = 4)

bm_buk <- round(assay(mda231_merged, "segment_ratios") * 2.5, digits = 0)
brkpt <- rle(round(assay(mda231_merged, "segment_ratios") * 2.5, digits = 0) %>% unlist() %>% as.numeric())

index <- sapply(1:length(brkpt[["lengths"]]), function(x){
  sum(brkpt[["lengths"]][1:x])
})[-length(brkpt[["lengths"]])]

brkpt_rg <- hg19_rg[index,]
brkpt_rg$start_up <- ifelse(brkpt_rg$start-220000>1, brkpt_rg$start-220000, 1)
brkpt_rg$end_down <- brkpt_rg$end+220000

library(GenomicRanges)
brkpt_gr <- makeGRangesFromDataFrame(brkpt_rg%>%dplyr::select(chr, start_up, end_down, strand, width), seqnames.field = "chr", start.field = "start_up", end.field = "end_down", keep.extra.columns = T)

brkpt_ck <- apply(assay(mda231_filt_knn, "integer"), 2, rle)
index_ck <- lapply(brkpt_ck, function(brkpt){
  sapply(1:length(brkpt[["lengths"]]), function(x){
    sum(brkpt[["lengths"]][1:x])
  })[-length(brkpt[["lengths"]])]
})

acc_ck <- lapply(index_ck, function(x){
  length(intersect(c(index+1,index, index-1), x))/length(index)
})

brkpt_hmm <- apply(mda231_hmm_cnmat_filt, 2, rle)
index_hmm <- lapply(brkpt_hmm, function(brkpt){
  sapply(1:length(brkpt[["lengths"]]), function(x){
    sum(brkpt[["lengths"]][1:x])
  })[-length(brkpt[["lengths"]])]
})

acc_hmm <- lapply(index_hmm, function(ind){
  hmm_gr <- makeGRangesFromDataFrame(segment_range[ind,])
  ct <- as_tibble(findOverlaps(hmm_gr, brkpt_gr)) %>% group_by(subjectHits) %>% summarise(n=n())
  nrow(ct)/length(index)
})


brkpt_gk <- apply(mda231_cnmat[,colnames(mda231_hmm_cnmat_filt)], 2, rle)
index_gk <- lapply(brkpt_gk, function(brkpt){
  sapply(1:length(brkpt[["lengths"]]), function(x){
    sum(brkpt[["lengths"]][1:x])
  })[-length(brkpt[["lengths"]])]
})

acc_gk <- lapply(index_gk, function(ind){
  gk_gr <- makeGRangesFromDataFrame(ranges[ind,])
  ct <- as_tibble(findOverlaps(gk_gr, brkpt_gr)) %>% group_by(subjectHits) %>% summarise(n=n())
  nrow(ct)/length(index)
})

acc_mat <- data.frame(copykit=unlist(acc_ck), hmmcopy=unlist(acc_hmm), ginkgo=unlist(acc_gk))
ggpubr::ggviolin(acc_mat %>% gather(value="acc_220k", key="tool"),
                 x = "tool", y = "acc_220k",
                 palette = "jco", add="boxplot") + 
  ggpubr::stat_compare_means(comparisons = my_comparisons,label = "p.signif")+ # Add pairwise comparisons p-value
  # ggpubr::stat_compare_means(label.y = 50)  +
  theme(axis.text.x = element_text(angle = 45, hjust=0.8),
        legend.position = "none") 
cowplot::ggsave2(glue("{outputdir}/acc_220k_compare_violin_mda231.pdf"),width = 3, height = 4)

save(list = ls(all.names = TRUE), file = glue("{outputdir}/20230828_mda231_compare.RData"))



