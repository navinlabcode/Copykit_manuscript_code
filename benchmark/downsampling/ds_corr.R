library(tidyverse)
library(copykit)
library(BiocParallel)
library(glue)
register(MulticoreParam(workers = 100, progressbar = F), default = T)
source("~/Projects/copykit/downsampling/ds_corr_func.R")
outdir <- "~/Projects/copykit/downsampling/ck_obj/"


## copykit
# brcan_ds <- run_ds_ck(path = "/volumes/USR1/junke/Projects/copykit/downsampling/data/brcan/")
# # names(brcan_ds) <- c("50k", "75k", "125k", "250k", "500k", "750k", "1M", "fulldata")
# wafer231p_ds <- run_ds_ck(path = "/volumes/USR1/junke/Projects/copykit/downsampling/data/wafer231p/")
# # names(wafer231p_ds) <- c("50k", "75k", "125k", "250k", "500k", "750k", "1M", "fulldata")
# P4P_ds <- run_ds_ck(path = "/volumes/USR1/junke/Projects/copykit/downsampling/data/P4P/", is_pairend = T)
# # names(P4P_ds) <- c("50k", "75k", "125k", "250k", "500k", "750k", "1M", "fulldata")
# 
# save(list = ls(all.names = T), file = glue("{outdir}/20230904_ds_ck_objs.Rdata"))
# 
# ### add one ffpe sample - 20230925
# load(glue("{outdir}/20230904_ds_ck_objs.Rdata"))
# source("~/Projects/copykit/downsampling/ds_corr_func.R")
# 
# pmtc1breast_ds <-  run_ds_ck(path = "/volumes/USR1/junke/Projects/copykit/downsampling/data/pmtc1breast/", res = c("downsampled_sort_50k", "downsampled_sort_75k", "downsampled_sort_125k", "downsampled_sort_250k", "downsampled_sort_500k", "downsampled_sort_750k", "fulldata"))
# save(list = ls(all.names = T), file = glue("{outdir}/20230925_ds_ck_objs.Rdata"))

load(glue("{outdir}/20230925_ds_ck_objs.Rdata"))

brcan_matlist <- lapply(brcan_ds, segment_ratios)
wafer231p_matlist <- lapply(wafer231p_ds, segment_ratios)
P4P_matlist <- lapply(P4P_ds, segment_ratios)
pmtc1breast_matlist <- lapply(pmtc1breast_ds, segment_ratios)


brcan_cordf_ck <- create_corr_df(brcan_matlist)
wafer231p_cordf_ck <- create_corr_df(wafer231p_matlist)
P4P_cordf_ck <- create_corr_df(P4P_matlist)
pmtc1breast_cordf_ck <- create_corr_df(pmtc1breast_matlist)


## scHMM
brcan_hmm_list <- read_ds_hmm("~/Projects/copykit/benchmark/sc_hmmcopy/brcan")
# names(brcan_hmm_list) <- c("50k", "75k", "125k", "250k", "500k", "750k", "1M", "fulldata")
wafer231p_hmm_list <- read_ds_hmm("~/Projects/copykit/benchmark/sc_hmmcopy/wafer231p")
P4P_hmm_list <- read_ds_hmm("~/Projects/copykit/benchmark/sc_hmmcopy/P4P")
pmtc1breast_hmm_list <- read_ds_hmm("~/Projects/copykit/benchmark/sc_hmmcopy/pmtc1breast", res = c("downsampled_sort_50k", "downsampled_sort_75k", "downsampled_sort_125k", "downsampled_sort_250k", "downsampled_sort_500k", "downsampled_sort_750k", "fulldata"))

brcan_cordf_hmm <- create_corr_df(brcan_hmm_list)
wafer231p_cordf_hmm <- create_corr_df(wafer231p_hmm_list)
P4P_cordf_hmm <- create_corr_df(P4P_hmm_list)
pmtc1breast_cordf_hmm <- create_corr_df(pmtc1breast_hmm_list)


## Ginkgo
brcan_gk_list <- read_ds_gk("/volumes/seq/code/PIPELINES/copy_number_pipeline_hg38/dev_lib/code/copy_number_pipeline_hg38/ginkgo/uploads/brcan")
wafer231p_gk_list <- read_ds_gk("/volumes/seq/code/PIPELINES/copy_number_pipeline_hg38/dev_lib/code/copy_number_pipeline_hg38/ginkgo/uploads/wafer231p")
P4P_gk_list <- read_ds_gk("/volumes/seq/code/PIPELINES/copy_number_pipeline_hg38/dev_lib/code/copy_number_pipeline_hg38/ginkgo/uploads/P4P")
pmtc1breast_gk_list <- read_ds_gk("/volumes/seq/code/PIPELINES/copy_number_pipeline_hg38/dev_lib/code/copy_number_pipeline_hg38/ginkgo/uploads/pmtc1breast", res = c("downsampled_sort_50k", "downsampled_sort_75k", "downsampled_sort_125k", "downsampled_sort_250k", "downsampled_sort_500k", "downsampled_sort_750k", "fulldata"))

brcan_cordf_gk <- create_corr_df(brcan_gk_list)
wafer231p_cordf_gk <- create_corr_df(wafer231p_gk_list)
P4P_cordf_gk <- create_corr_df(P4P_gk_list)
pmtc1breast_cordf_gk <- create_corr_df(pmtc1breast_gk_list)

p1 <- ggplot(do.call(rbind, list(
  brcan_cordf_ck %>% mutate(tool="copykit"),
  brcan_cordf_gk %>% mutate(tool="ginkgo"),
  brcan_cordf_hmm %>% mutate(tool="hmmcopy")
)) %>% mutate(n_reads = fct_relevel(n_reads, n_reads_levels)) , 
             aes(x = n_reads, y = correlation, fill = tool)) + 
  ylim(c(0, 1)) + 
  # geom_point(aes(color = tool),position = position_jitterdodge(seed = 1), size = 0.05) +
  # geom_boxplot(aes(color = tool),alpha=0, size = 0.5) + 
  stat_summary(aes(group = tool, color = tool),fun = median, fun.min = function(z) quantile(z, 0.25),
               fun.max = function(z) quantile(z, 0.75), position = position_dodge(width = 0.3)) +
  stat_summary(
    fun = median,
    geom = 'line',
    aes(group = tool, color = tool),
    position = position_dodge(width = 0.3) 
  ) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  cowplot::theme_cowplot() +
  xlab("number of reads (downsampled)") + ylab("correlation with original profile")

p1

p1_j <- ggplot(do.call(rbind, list(
  brcan_cordf_ck %>% mutate(tool="copykit"),
  brcan_cordf_gk %>% mutate(tool="ginkgo"),
  brcan_cordf_hmm %>% mutate(tool="hmmcopy")
)) %>% mutate(n_reads = fct_relevel(n_reads, n_reads_levels)) , 
aes(x = n_reads, y = correlation, fill = tool)) + 
  ylim(c(0, 1)) + 
  geom_point(aes(color = tool),position = position_jitterdodge(seed = 1), size = 0.05, alpha=0.5) +
  # geom_boxplot(aes(color = tool),alpha=0, size = 0.5) + 
  stat_summary(fun = median, fun.min = function(z) quantile(z, 0.25),
               fun.max = function(z) quantile(z, 0.75), position = position_dodge(width = 0.75)) +
  stat_summary(
    fun = median,
    geom = 'line',
    aes(group = tool, color = tool),
    position = position_dodge(width = 0.75) 
  ) +
  cowplot::theme_cowplot() +
  xlab("number of reads (downsampled)") + ylab("correlation with original profile")

p1_j
cowplot::ggsave2(glue("{outdir}/downsample_brcan_cor_line.pdf"), p1, width = 5.5, height = 4)  
cowplot::ggsave2(glue("{outdir}/downsample_brcan_cor_jitline.pdf"), p1_j, width = 5.5, height = 4)  
saveRDS(do.call(rbind, list(
  brcan_cordf_ck %>% mutate(tool="copykit"),
  brcan_cordf_gk %>% mutate(tool="ginkgo"),
  brcan_cordf_hmm %>% mutate(tool="hmmcopy")
)) %>% mutate(n_reads = fct_relevel(n_reads, n_reads_levels)), glue("{outdir}/downsample_brcan_cor.rds"))

  
p2 <- ggplot(do.call(rbind, list(
  wafer231p_cordf_ck %>% mutate(tool="copykit"),
  wafer231p_cordf_gk %>% mutate(tool="ginkgo"),
  wafer231p_cordf_hmm %>% mutate(tool="hmmcopy")
)) %>% mutate(n_reads = fct_relevel(n_reads, n_reads_levels)) , 
aes(x = n_reads, y = correlation, fill = tool)) + 
  ylim(c(0, 1)) + 
  # geom_point(aes(color = tool),position = position_jitterdodge(seed = 1), size = 0.05) +
  # geom_boxplot(aes(color = tool),alpha=0, size = 0.5) + 
  stat_summary(aes(group = tool, color = tool),fun = median, fun.min = function(z) quantile(z, 0.25),
               fun.max = function(z) quantile(z, 0.75), position = position_dodge(width = 0.3)) +
  stat_summary(
    fun = median,
    geom = 'line',
    aes(group = tool, color = tool),
    position = position_dodge(width = 0.3) 
  ) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  cowplot::theme_cowplot() +
  xlab("number of reads (downsampled)") + ylab("correlation with original profile")

p2

p2_j <- ggplot(do.call(rbind, list(
  wafer231p_cordf_ck %>% mutate(tool="copykit"),
  wafer231p_cordf_gk %>% mutate(tool="ginkgo"),
  wafer231p_cordf_hmm %>% mutate(tool="hmmcopy")
)) %>% mutate(n_reads = fct_relevel(n_reads, n_reads_levels)) , 
aes(x = n_reads, y = correlation, fill = tool)) + 
  ylim(c(0, 1)) + 
  geom_point(aes(color = tool),position = position_jitterdodge(seed = 1), size = 0.05, alpha=0.5) +
  # geom_boxplot(aes(color = tool),alpha=0, size = 0.5) + 
  stat_summary(fun = median, fun.min = function(z) quantile(z, 0.25),
               fun.max = function(z) quantile(z, 0.75), position = position_dodge(width = 0.75)) +
  stat_summary(
    fun = median,
    geom = 'line',
    aes(group = tool, color = tool),
    position = position_dodge(width = 0.75) 
  ) +
  cowplot::theme_cowplot() +
  xlab("number of reads (downsampled)") + ylab("correlation with original profile")

p2_j
cowplot::ggsave2(glue("{outdir}/downsample_wafer231p_cor_line.pdf"), p2, width = 5.5, height = 4)  
cowplot::ggsave2(glue("{outdir}/downsample_wafer231p_cor_jitline.pdf"), p2_j, width = 5.5, height = 4)  
saveRDS(do.call(rbind, list(
  wafer231p_cordf_ck %>% mutate(tool="copykit"),
  wafer231p_cordf_gk %>% mutate(tool="ginkgo"),
  wafer231p_cordf_hmm %>% mutate(tool="hmmcopy")
)) %>% mutate(n_reads = fct_relevel(n_reads, n_reads_levels)), glue("{outdir}/downsample_wafer231p_cor.rds"))


p3 <- ggplot(do.call(rbind, list(
  P4P_cordf_ck %>% mutate(tool="copykit"),
  P4P_cordf_gk %>% mutate(tool="ginkgo"),
  P4P_cordf_hmm %>% mutate(tool="hmmcopy")
)) %>% mutate(n_reads = fct_relevel(n_reads, n_reads_levels)) , 
aes(x = n_reads, y = correlation, fill = tool)) + 
  ylim(c(0, 1)) + 
  # geom_point(aes(color = tool),position = position_jitterdodge(seed = 1), size = 0.05) +
  # geom_boxplot(aes(color = tool),alpha=0, size = 0.5) + 
  stat_summary(aes(group = tool, color = tool),fun = median, fun.min = function(z) quantile(z, 0.25),
               fun.max = function(z) quantile(z, 0.75), position = position_dodge(width = 0.3)) +
  stat_summary(
    fun = median,
    geom = 'line',
    aes(group = tool, color = tool),
    position = position_dodge(width = 0.3) 
  ) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  cowplot::theme_cowplot() +
  xlab("number of reads (downsampled)") + ylab("correlation with original profile")

p3
cowplot::ggsave2(glue("{outdir}/downsample_P4P_cor_line.pdf"), p3, width = 5.5, height = 4)  
saveRDS(do.call(rbind, list(
  P4P_cordf_ck %>% mutate(tool="copykit"),
  P4P_cordf_gk %>% mutate(tool="ginkgo"),
  P4P_cordf_hmm %>% mutate(tool="hmmcopy")
)) %>% mutate(n_reads = fct_relevel(n_reads, n_reads_levels)), glue("{outdir}/downsample_P4P_cor.rds"))



p4 <- ggplot(do.call(rbind, list(
  pmtc1breast_cordf_ck %>% mutate(tool="copykit"),
  pmtc1breast_cordf_gk %>% mutate(tool="ginkgo"),
  pmtc1breast_cordf_hmm %>% mutate(tool="hmmcopy")
)) %>% mutate(n_reads = fct_relevel(n_reads, n_reads_levels)) , 
aes(x = n_reads, y = correlation, fill = tool)) + 
  ylim(c(0, 1)) + 
  # geom_point(aes(color = tool),position = position_jitterdodge(seed = 1), size = 0.05) +
  # geom_boxplot(aes(color = tool),alpha=0, size = 0.5) + 
  stat_summary(aes(group = tool, color = tool),fun = median, fun.min = function(z) quantile(z, 0.25),
               fun.max = function(z) quantile(z, 0.75), position = position_dodge(width = 0.3)) +
  stat_summary(
    fun = median,
    geom = 'line',
    aes(group = tool, color = tool),
    position = position_dodge(width = 0.3) 
  ) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  cowplot::theme_cowplot() +
  xlab("number of reads (downsampled)") + ylab("correlation with original profile")

p4
cowplot::ggsave2(glue("{outdir}/downsample_pmtc1breast_cor_line.pdf"), p4, width = 5, height = 4)  
saveRDS(do.call(rbind, list(
  pmtc1breast_cordf_ck %>% mutate(tool="copykit"),
  pmtc1breast_cordf_gk %>% mutate(tool="ginkgo"),
  pmtc1breast_cordf_hmm %>% mutate(tool="hmmcopy")
)) %>% mutate(n_reads = fct_relevel(n_reads, n_reads_levels)), glue("{outdir}/downsample_pmtc1breast_cor.rds"))

# save(list = ls(all.names = T), file = glue("{outdir}/20230904_cor.rds"))

## copykit w smoothing

# brcan_ds <- run_ds_ck_smooth(path = "/volumes/USR1/junke/Projects/copykit/downsampling/data/brcan/")
# # names(brcan_ds) <- c("50k", "75k", "125k", "250k", "500k", "750k", "1M", "fulldata")
# wafer231p_ds <- run_ds_ck_smooth(path = "/volumes/USR1/junke/Projects/copykit/downsampling/data/wafer231p/")
# # names(wafer231p_ds) <- c("50k", "75k", "125k", "250k", "500k", "750k", "1M", "fulldata")
# P4P_ds <- run_ds_ck_smooth(path = "/volumes/USR1/junke/Projects/copykit/downsampling/data/P4P/", is_pairend = T)
# # names(P4P_ds) <- c("50k", "75k", "125k", "250k", "500k", "750k", "1M", "fulldata")
# 
# pmtc1breast_ds <-  run_ds_ck_smooth(path = "/volumes/USR1/junke/Projects/copykit/downsampling/data/pmtc1breast/", res = c("downsampled_sort_50k", "downsampled_sort_75k", "downsampled_sort_125k", "downsampled_sort_250k", "downsampled_sort_500k", "downsampled_sort_750k", "fulldata"))
# 
# save(list = ls(all.names = T), file = glue("{outdir}/20240304_ds_ck_smooth_objs.Rdata"))

load(glue("{outdir}/20240304_ds_ck_smooth_objs.Rdata"))

brcan_matlist <- lapply(brcan_ds, segment_ratios)
wafer231p_matlist <- lapply(wafer231p_ds, segment_ratios)
P4P_matlist <- lapply(P4P_ds, segment_ratios)
pmtc1breast_matlist <- lapply(pmtc1breast_ds, segment_ratios)


brcan_cordf_ck_sm <- create_corr_df(brcan_matlist)
wafer231p_cordf_ck_sm <- create_corr_df(wafer231p_matlist)
P4P_cordf_ck_sm <- create_corr_df(P4P_matlist)
pmtc1breast_cordf_ck_sm <- create_corr_df(pmtc1breast_matlist)

brcan_cordf <- readRDS(glue("{outdir}/downsample_brcan_cor.rds"))
brcan_cordf <- rbind(brcan_cordf, brcan_cordf_ck_sm %>% mutate(tool="ck_w_sm"))

p1 <- ggplot(brcan_cordf %>% mutate(n_reads = fct_relevel(n_reads, n_reads_levels)) , 
aes(x = n_reads, y = correlation, fill = tool)) + 
  ylim(c(0, 1)) + 
  # geom_point(aes(color = tool),position = position_jitterdodge(seed = 1), size = 0.05) +
  # geom_boxplot(aes(color = tool),alpha=0, size = 0.5) + 
  stat_summary(aes(group = tool, color = tool),fun = median, fun.min = function(z) quantile(z, 0.25),
               fun.max = function(z) quantile(z, 0.75), position = position_dodge(width = 0.3)) +
  stat_summary(
    fun = median,
    geom = 'line',
    aes(group = tool, color = tool),
    position = position_dodge(width = 0.3) 
  ) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  cowplot::theme_cowplot() +
  xlab("number of reads (downsampled)") + ylab("correlation with original profile")

p1
cowplot::ggsave2(glue("{outdir}/downsample_brcan_wsm_cor_line.pdf"), p1, width = 5.5, height = 4)  


wafer231p_cordf <- readRDS(glue("{outdir}/downsample_wafer231p_cor.rds"))
wafer231p_cordf <- rbind(wafer231p_cordf, wafer231p_cordf_ck_sm %>% mutate(tool="ck_w_sm"))

p2 <- ggplot(wafer231p_cordf %>% mutate(n_reads = fct_relevel(n_reads, n_reads_levels)) , 
             aes(x = n_reads, y = correlation, fill = tool)) + 
  ylim(c(0, 1)) + 
  # geom_point(aes(color = tool),position = position_jitterdodge(seed = 1), size = 0.05) +
  # geom_boxplot(aes(color = tool),alpha=0, size = 0.5) + 
  stat_summary(aes(group = tool, color = tool),fun = median, fun.min = function(z) quantile(z, 0.25),
               fun.max = function(z) quantile(z, 0.75), position = position_dodge(width = 0.3)) +
  stat_summary(
    fun = median,
    geom = 'line',
    aes(group = tool, color = tool),
    position = position_dodge(width = 0.3) 
  ) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  cowplot::theme_cowplot() +
  xlab("number of reads (downsampled)") + ylab("correlation with original profile")

p2
cowplot::ggsave2(glue("{outdir}/downsample_wafer231p_wsm_cor_line.pdf"), p2, width = 5.5, height = 4)  

pmtc1breast_cordf <- readRDS(glue("{outdir}/downsample_pmtc1breast_cor.rds"))
pmtc1breast_cordf <- rbind(pmtc1breast_cordf, pmtc1breast_cordf_ck_sm %>% mutate(tool="ck_w_sm"))

p3 <- ggplot(pmtc1breast_cordf %>% mutate(n_reads = fct_relevel(n_reads, n_reads_levels)) , 
             aes(x = n_reads, y = correlation, fill = tool)) + 
  ylim(c(0, 1)) + 
  # geom_point(aes(color = tool),position = position_jitterdodge(seed = 1), size = 0.05) +
  # geom_boxplot(aes(color = tool),alpha=0, size = 0.5) + 
  stat_summary(aes(group = tool, color = tool),fun = median, fun.min = function(z) quantile(z, 0.25),
               fun.max = function(z) quantile(z, 0.75), position = position_dodge(width = 0.3)) +
  stat_summary(
    fun = median,
    geom = 'line',
    aes(group = tool, color = tool),
    position = position_dodge(width = 0.3) 
  ) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  cowplot::theme_cowplot() +
  xlab("number of reads (downsampled)") + ylab("correlation with original profile")

p3
cowplot::ggsave2(glue("{outdir}/downsample_pmtc1breast_wsm_cor_line.pdf"), p3, width = 5, height = 4)  

save(list = ls(all.names = T), file = glue("{outdir}/20240305_cor.rds"))
