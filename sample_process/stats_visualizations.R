library(devtools)
library(copykit, lib.loc = "/mnt/USR1/junke/Projects/copykit/rerun_manu/lib/")
library(ggdist)
library(tidyverse)
library(jcolors)
library(BiocParallel)
library(glue)
library(janitor)
library(boot)
library(cowplot)
set.seed(17)
register(MulticoreParam(workers = 50), default = T)

# 
load("/mnt/USR1/junke/Projects/copykit/rerun_manu/singu_0424/20240424_ck_singul.Rdata")
outdir="/mnt/USR1/junke/Projects/copykit/rerun_manu/singu_0424"

calcSmoothOver <- function(scCNA,
                           BPPARAM = bpparam()) {
  
  dat_bin <- assay(scCNA, 'smoothed_bincounts')
  overdisp <- BiocParallel::bplapply(dat_bin,
                                     overdispersion,
                                     BPPARAM = BPPARAM
  )
  
  overdisp <- unlist(overdisp)
  
  SummarizedExperiment::colData(scCNA)$overdispersion_smooth <- overdisp
  
  return(scCNA)
  
}

shan <- function(data, indices) {
  data_ind <- data[indices]
  prop <- janitor::tabyl(data_ind) %>% pull(percent)
  div <- -sum(prop*log(prop))
  return(div)
}

shan_sample <- function(meta, group) {
  meta_group <- meta %>%
    filter(L_foci == !!group)
  group_perc <- janitor::tabyl(as.character(meta_group$subclones)) %>% pull(percent)
  group_div <- -sum(group_perc*log(group_perc))
  
  if (group_div == 0) {
    df <- data.frame(foci = group,
                     div = 0,
                     lci = 0,
                     uci = 0)
    return(df)
  }
  
  shan_group_boot <-
    boot::boot(
      as.character(meta_group$subclones),
      statistic = shan,
      R = 2000
    )
  shan_group_boot_ci <- suppressWarnings(boot::boot.ci(shan_group_boot))
  df <- data.frame(foci = group,
                   div = group_div,
                   lci = shan_group_boot_ci$normal[2],
                   uci = shan_group_boot_ci$normal[3])
  
}

plotOverdispersionComparison <- function(scCNA) {
  
  coldata_df <- as.data.frame(colData(scCNA)) %>%
    dplyr::select(sample,
                  overdispersion,
                  overdispersion_smooth) %>%
    dplyr::rename(original = 'overdispersion',
                  knnSmooth = 'overdispersion_smooth') %>%
    tidyr::gather(key = 'metric',
                  value = 'overdispersion',
                  -sample) %>%
    dplyr::mutate(metric = forcats::fct_relevel(metric, c('original','knnSmooth')))
  
  ggplot() +
    geom_line(data = coldata_df, aes(x = metric, y = overdispersion, group = sample), alpha = .2, color = 'gray85') +
    ggdist::stat_slab(
      data = coldata_df %>% filter(metric == 'original'),
      aes(x = metric, y = overdispersion),
      adjust = .3,
      width = .4,
      fill = '#286886',
      side = 'left') +
    ggdist::stat_slab(
      data = coldata_df %>% filter(metric == 'knnSmooth'),
      aes(x = metric, y = overdispersion),
      adjust = .3,
      width = .4,
      fill = '#7B321C',
      side = 'right') +
    geom_point(data = coldata_df, aes(x = metric, y = overdispersion), color = 'black',
               size = .2, shape = 21) +
    cowplot::theme_cowplot() +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 9)) +
    xlab('')
  
}




## pmtc6 plots

pmtc6_metrics %>%
  summarise(reads_total_mean = mean(reads_total),
            reads_total_sd = sd(reads_total),
            percent_dup_mean = mean(percentage_duplicates),
            percent_dup_sd = sd(percentage_duplicates)
  )
# reads_total_mean reads_total_sd percent_dup_mean percent_dup_sd
#          1036403       245796.5       0.08145516    0.009694156

pmtc6_metrics %>% as.data.frame() %>%
  rstatix::kruskal_test(overdispersion ~ outlier)
# .y.                n statistic    df      p method        
# overdispersion   892      6.13     1 0.0133 Kruskal-Wallis

pmtc6_metrics %>% as.data.frame() %>%
  rstatix::kruskal_test(breakpoint_count ~ outlier)
# .y.                  n statistic    df        p method        
# breakpoint_count   892      54.8     1 1.35e-13 Kruskal-Wallis

##### WHYYYYYY
cowplot::ggsave2(glue("{outdir}/pmtc6_p_metrics.pdf"), pmtc6_p_metrics, width = 6, height = 4)


wilcox.test(as.data.frame(colData(pmtc6_liver_knn))$overdispersion, as.data.frame(colData(pmtc6_liver_knn))$overdispersion_smooth, paired = TRUE )
# V = 314820, p-value < 2.2e-16


pmtc6_seg_ratios_dist <- amap::Dist(t(segment_ratios(pmtc6_liver_knn)), method = 'manhattan', nbproc = 30)

pmtc6_liver_knn_hclust <- fastcluster::hclust(pmtc6_seg_ratios_dist, method = 'ward.D2')

pdf(glue("{outdir}/pmtc6_liver_filt_knnhclust.pdf"), width = 7, height = 8)
plotHeatmap(pmtc6_liver_filt[,pmtc6_liver_knn_hclust$order], raster_quality = 5)
dev.off()

pdf(glue("{outdir}/pmtc6_liver_knn_knnhclust.pdf"), width = 7, height = 8)
plotHeatmap(pmtc6_liver_knn[,pmtc6_liver_knn_hclust$order], raster_quality = 5)
dev.off()

pmtc6_over_comparison_p <- plotOverdispersionComparison(pmtc6_liver_knn)
cowplot::ggsave2(glue("{outdir}/pmtc6_overdisp_comparison.pdf"), pmtc6_over_comparison_p, width =3.5, height = 5)

pmtc6_metrics_aftfilt %>%
  summarise(reads_total_mean = mean(reads_total),
            reads_total_sd = sd(reads_total),
            percent_dup_mean = mean(percentage_duplicates),
            percent_dup_sd = sd(percentage_duplicates)
  )
# reads_total_mean reads_total_sd percent_dup_mean percent_dup_sd
#          1064625       220056.2       0.08174202    0.009733249

cowplot::ggsave2(glue("{outdir}/pmtc6_liver_knn_pl_umap_p_ploidysubmedian.pdf"), pmtc6_liver_knn_pl_umap_p_ploidysubmedian, width = 4, height = 3)


cowplot::ggsave2(glue("{outdir}/pmtc6_primary_suggestedK_p.pdf"), pmtc6_primary_suggestedK_p, width = 5, height = 3)

janitor::tabyl(colData(pmtc6_liver_knn_pl)$subclones)
 # colData(pmtc6_liver_knn_pl)$subclones   n    percent
 #                                    c1  64 0.08876560
 #                                    c2  25 0.03467406
 #                                    c3  62 0.08599168
 #                                    c4 570 0.79056865

cowplot::ggsave2(glue("{outdir}/pmtc6_primary_umap.pdf"), pmtc6_liver_knn_pl_umap_p, width = 4, height = 3)

pmtc6_liver_knn_ploidy_umap_p <- plotUmap(pmtc6_liver_knn_pl, label = 'ploidy') + scale_fill_jcolors_contin('rainbow')
cowplot::ggsave2(glue("{outdir}/pmtc6_ploidy_umap.pdf"), pmtc6_liver_knn_ploidy_umap_p, width = 4, height = 3)

cowplot::ggsave2(glue("{outdir}/pmtc6_primary_umap_l_foci.pdf"), pmtc6_liver_knn_pl_umap_l_foci_p, width = 4.1, height = 3)

cowplot::ggsave2(glue("{outdir}/pmtc6_fociS1_pie.pdf"), pmtc6_foci_pies[[1]], width = 3, height = 3)
cowplot::ggsave2(glue("{outdir}/pmtc6_fociS4S5_pie.pdf"), pmtc6_foci_pies[[2]], width = 3, height = 3)
cowplot::ggsave2(glue("{outdir}/pmtc6_fociS6S7_pie.pdf"), pmtc6_foci_pies[[3]], width = 3, height = 3)
cowplot::ggsave2(glue("{outdir}/pmtc6_fociS3_pie.pdf"), pmtc6_foci_pies[[4]], width = 3, height = 3)
cowplot::ggsave2(glue("{outdir}/pmtc6_fociS2_pie.pdf"), pmtc6_foci_pies[[5]], width = 3, height = 3)
cowplot::ggsave2(glue("{outdir}/pmtc6_fociS8S9_pie.pdf"), pmtc6_foci_pies[[6]], width = 3, height = 3)


library(ggrepel)
library(ggnewscale)
pmtc6_lfoci_umap_p <- ggplot(pmtc6_lfoci_umap, aes(V1, V2, label = subclones)) +
  geom_point(aes(fill = L_foci), shape = 21,
             stroke = 0.1,
             size = 2.5) +
  ggnewscale::new_scale_fill() +
  geom_text_repel(aes(color = subclones),
                  min.segment.length = 0,
                  box.padding = 0.01,
                  size = 4,
                  max.overlaps = Inf,
                  na.rm = T) +
  scale_color_manual(values = subclones_pal(),
                     limits = force) +
  theme_classic() +
  my_theme

cowplot::ggsave2(glue("{outdir}/pmtc6_lfoci_umap.pdf"), pmtc6_lfoci_umap_p, width = 4, height = 3)

# * foci diversity ----
pmtc6_meta <- colData(pmtc6_liver_knn_pl)

pmtc6_meta <- as.data.frame(pmtc6_meta)
pmtc6_meta <- pmtc6_meta[c('subclones', 'L_foci')]

pmtc6_L1_div <- shan_sample(pmtc6_meta, "S1")
pmtc6_L2_div <- shan_sample(pmtc6_meta, "S2")
pmtc6_L3_div <- shan_sample(pmtc6_meta, "S3")
pmtc6_L4L5_div <- shan_sample(pmtc6_meta, "S4S5")
pmtc6_L6L7_div <- shan_sample(pmtc6_meta, "S6S7")
pmtc6_L8L9_div <- shan_sample(pmtc6_meta, "S8S9")

pmtc6_div <- bind_rows(
  pmtc6_L1_div,
  pmtc6_L2_div,
  pmtc6_L3_div,
  pmtc6_L4L5_div,
  pmtc6_L6L7_div,
  pmtc6_L8L9_div
)

pmtc6_div
# foci        div         lci       uci
# 1   S1 0.36766228  0.30044742 0.4384327
# 2   S2 0.67749440  0.62095736 0.7644487
# 3   S3 0.00000000  0.00000000 0.0000000
# 4 S4S5 0.05465882 -0.02931279 0.1489054
# 5 S6S7 0.75680931  0.67865928 0.8513240
# 6 S8S9 0.27118937  0.03129424 0.5562124
table(colData(pmtc6_liver_knn_pl)$subclones, colData(pmtc6_liver_knn_pl)$L_foci)
#     S1  S2  S3 S4S5 S6S7 S8S9
# c1   0   0   0    0   64    0
# c2   0   0   0    0    1   24
# c3  39  20   0    1    2    0
# c4 285  14  66  102  101    2
pmtc6_div_p <- pmtc6_div %>%
  ggplot() +
  geom_errorbar(aes(
    x = fct_reorder(foci, div),
    ymin = lci,
    ymax = uci
  ),
  width = .1,
  size = 1) +
  geom_point(aes(
    x = fct_reorder(foci, div),
    y = div,
    color = foci
  ), size = 5) +
  cowplot::theme_cowplot() +
  theme(legend.position = "none",
        axis.text.x = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1
        )) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 12)) +
  ylab("shannon diversity index") +
  xlab("") +
  coord_flip()

cowplot::ggsave2(glue("{outdir}/pmtc6_foci_shan_diversity.pdf"), pmtc6_div_p, height = 3, width = 4.5)

pmtc7_metrics %>%
  summarise(reads_total_mean = mean(reads_total),
            reads_total_sd = sd(reads_total),
            percent_dup_mean = mean(percentage_duplicates),
            percent_dup_sd = sd(percentage_duplicates)
  )
# reads_total_mean reads_total_sd percent_dup_mean percent_dup_sd
#        975057.5       303911.9       0.08932727     0.01065848

pmtc7_metrics %>% as.data.frame() %>%
  rstatix::kruskal_test(breakpoint_count ~ outlier)
# .y.                  n statistic    df      p method        
# breakpoint_count  1375      3.39     1 0.0655 Kruskal-Wallis

pmtc7_metrics %>% as.data.frame() %>%
  rstatix::kruskal_test(overdispersion ~ outlier)
# .y.                n statistic    df             p method        
# overdispersion  1375      32.9     1 0.00000000958 Kruskal-Wallis

pmtc7_metrics %>% as.data.frame() %>%
  rstatix::kruskal_test(reads_total ~ outlier)
# .y.             n statistic    df        p method        
# reads_total  1375      55.7     1 8.36e-14 Kruskal-Wallis

cowplot::ggsave2(glue("{outdir}/pmtc7_p_metrics.pdf"), pmtc7_p_metrics, width = 6, height = 4)


pmtc7_over_comparison_p <- plotOverdispersionComparison(pmtc7_liver_knn)
cowplot::ggsave2(glue("{outdir}/pmtc7_overdisp_comparison.pdf"), pmtc7_over_comparison_p, width = 2, height = 3)

pmtc7_metrics_aftfilt %>%
  summarise(reads_total_mean = mean(reads_total),
            reads_total_sd = sd(reads_total),
            percent_dup_mean = mean(percentage_duplicates),
            percent_dup_sd = sd(percentage_duplicates)
  )
# reads_total_mean reads_total_sd percent_dup_mean percent_dup_sd
#          997439.9       296605.5       0.08932845     0.01070951

cowplot::ggsave2(glue("{outdir}/pmtc7_liver_knn_pl_umap_p_ploidysubmedian.pdf"), pmtc7_liver_knn_pl_umap_p_ploidysubmedian, width = 4, height = 3)

cowplot::ggsave2(glue("{outdir}/pmtc7_primary_suggestedK_p.pdf"), pmtc7_primary_suggestedK_p, width = 5, height = 3)
cowplot::ggsave2(glue("{outdir}/pmtc7_primary_umap.pdf"), pmtc7_liver_knn_pl_umap_p, width = 4, height = 3)

cowplot::ggsave2(glue("{outdir}/pmtc7_primary_umap_l_foci.pdf"), pmtc7_liver_knn_pl_umap_l_foci_p, width = 4, height = 3)

cowplot::ggsave2(glue("{outdir}/pmtc7_fociL1_pie.pdf"), pmtc7_foci_pies[[1]], width = 3, height = 3)
cowplot::ggsave2(glue("{outdir}/pmtc7_fociL5_pie.pdf"), pmtc7_foci_pies[[2]], width = 3, height = 3)
cowplot::ggsave2(glue("{outdir}/pmtc7_fociL6_pie.pdf"), pmtc7_foci_pies[[3]], width = 3, height = 3)
cowplot::ggsave2(glue("{outdir}/pmtc7_fociL7_pie.pdf"), pmtc7_foci_pies[[4]], width = 3, height = 3)
cowplot::ggsave2(glue("{outdir}/pmtc7_fociL4_pie.pdf"), pmtc7_foci_pies[[5]], width = 3, height = 3)

pmtc7_lfoci_umap_p <- ggplot(pmtc7_lfoci_umap, aes(V1, V2, label = subclones)) +
  geom_point(aes(fill = L_foci), shape = 21,
             stroke = 0.1,
             size = 2.5) +
  ggnewscale::new_scale_fill() +
  geom_text_repel(aes(color = subclones),
                  min.segment.length = 0,
                  box.padding = 0.01,
                  label.padding = 0.15,
                  size = 4,
                  max.overlaps = Inf,
                  na.rm = T) +
  scale_color_manual(values = subclones_pal(),
                     limits = force) +
  theme_classic() +
  my_theme

cowplot::ggsave2(glue("{outdir}/pmtc7_lfoci_umap.pdf"), pmtc7_lfoci_umap_p, width = 4, height = 3)


pmtc7_L1_div <- shan_sample(pmtc7_meta, "S1")
pmtc7_L4_div <- shan_sample(pmtc7_meta, "S4")
pmtc7_L5_div <- shan_sample(pmtc7_meta, "S5")
pmtc7_L6_div <- shan_sample(pmtc7_meta, "S6")
pmtc7_L7_div <- shan_sample(pmtc7_meta, "S7")

pmtc7_div <- bind_rows(
  pmtc7_L1_div,
  pmtc7_L4_div,
  pmtc7_L5_div,
  pmtc7_L6_div,
  pmtc7_L7_div
)

pmtc7_div
# foci       div       lci       uci
# 1   S1 0.6529304 0.6097658 0.7040394
# 2   S4 0.8683769 0.7521261 1.0072650
# 3   S5 1.9577063 1.8912295 2.0517739
# 4   S6 1.6562524 1.5394296 1.8164097
# 5   S7 1.0340151 0.9247150 1.1687441


table(colData(pmtc7_liver_knn_pl)$subclones, colData(pmtc7_liver_knn_pl)$L_foci)
  #      S1  S4  S5  S6  S7
  # c0    0   0   0   0   0
  # c1    0 234  25   2   0
  # c2    0   0 112   8   0
  # c3    0  43   0   0   0
  # c4    0   0  11  14   0
  # c5    0   0  54   0   0
  # c6    0   0  42   0   0
  # c7    0   2  11  15  10
  # c8    0   0  25   5   0
  # c9    0   3  64   5   3
  # c10   0  20  10   0   0
  # c11   0   2   1  14  63
  # c12   0   6   4  41   3
  # c13   0   1   0  84  81
  # c14  51   0   0   0   0
  # c15  91   0   0   0   0

pmtc7_div_p <- pmtc7_div %>%
  ggplot() +
  geom_errorbar(aes(
    x = fct_reorder(foci, div),
    ymin = lci,
    ymax = uci
  ),
  width = .1,
  size = 1) +
  geom_point(aes(
    x = fct_reorder(foci, div),
    y = div,
    color = foci
  ), size = 5) +
  cowplot::theme_cowplot() +
  theme(legend.position = "none",
        axis.text.x = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1
        )) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 12)) +
  ylab("shannon diversity index") +
  xlab("") +
  expand_limits(y = 0) +
  coord_flip()


cowplot::ggsave2(glue("{outdir}/pmtc7_foci_shan_diversity.pdf"), pmtc7_div_p, height = 3, width = 4.5)


pmtc1_breast_overd_plot <- plotOverdispersionComparison(pmtc1_breast_knn_reserv)

cowplot::ggsave2(glue("{outdir}/pmtc1_breast_overd_plot.pdf"),
                 pmtc1_breast_overd_plot,
                 width = 4.5,
                 height = 5)


pmtc1_seg_ratios_dist <- amap::Dist(t(segment_ratios(pmtc1_breast_knn)), method = 'manhattan', nbproc = 30)

pmtc1_breast_knn_hclust <- fastcluster::hclust(pmtc1_seg_ratios_dist, method = 'ward.D2')

pdf(glue("{outdir}/pmtc1_breast_prior_smomothing_knnhclust.pdf"), width = 7, height = 8)
plotHeatmap(pmtc1_breast_filt[,pmtc1_breast_knn_hclust$order], raster_quality = 5)
dev.off()


pdf(glue("{outdir}/pmtc1_breast_after_smomothing_knnhclust.pdf"), width = 7, height = 8)
plotHeatmap(pmtc1_breast_knn[,pmtc1_breast_knn_hclust$order], raster_quality = 5)
dev.off()


pmtc1_metrics_aftfilt %>%
  summarise(reads_total_mean = mean(reads_total),
            reads_total_sd = sd(reads_total),
            percent_dup_mean = mean(percentage_duplicates),
            percent_dup_sd = sd(percentage_duplicates)
  )
# reads_total_mean reads_total_sd percent_dup_mean percent_dup_sd
# 1          1270697       857154.4        0.1709941      0.1656996

cowplot::ggsave2(glue("{outdir}/pmtc1_merged_umap_timepoint.pdf"), pmtc1_merged_tp_umap_p, width = 4.2, height = 3)
cowplot::ggsave2(glue("{outdir}/pmtc1_merged_umap.pdf"), pmtc1_merged_umap_p, width = 4, height = 3)

cowplot::ggsave2(glue("{outdir}/pmtc1_merged_nod_suggestedk.pdf"), pmtc1_merged_nod_suggestedk, width = 5, height = 3)

cowplot::ggsave2(glue("{outdir}/pmtc1_merged_suggestedk.pdf"), pmtc1_merged_suggestedk, width = 5, height = 3)
cowplot::ggsave2(glue("{outdir}/pmtc1_merged_umap_withdoublets.pdf"), pmtc1_merged_umap_p, width = 4, height = 3)
(ape::cophenetic.phylo(consensusPhylo(pmtc1_merged_knn_pl_nod_noc0)))
#          c1       c2       c3       c4       c5
# c1    0.000  511.930 1057.374 1246.656 1690.491
# c2  511.930    0.000 1185.714 1374.996 1818.831
# c3 1057.374 1185.714    0.000  526.400 1490.542
# c4 1246.656 1374.996  526.400    0.000 1679.825
# c5 1690.491 1818.831 1490.542 1679.825    0.000
cowplot::ggsave2(glue("{outdir}/pmtc1_merged_consensus_phylo.pdf"), pmtc1_merged_consensus_phylo, width = 5, height = 5)

cowplot::ggsave2(glue("{outdir}/pmtc1_merged_knn_pl_nod_noc0_umap_p_ploidysubmedian.pdf"), pmtc1_merged_knn_pl_nod_noc0_umap_p_ploidysubmedian, width = 4, height = 3)
cowplot::ggsave2(glue("{outdir}/pmtc1_merged_hvg_gc.pdf"), pmtc1_merged_hvg_gc, width = 5, height = 8)

co8_metrics_aftfilt %>%
  summarise(reads_total_mean = mean(reads_total),
            reads_total_sd = sd(reads_total),
            percent_dup_mean = mean(percentage_duplicates),
            percent_dup_sd = sd(percentage_duplicates)
  )
# reads_total_mean reads_total_sd percent_dup_mean percent_dup_sd
#         1050636       411230.5       0.08980761     0.01427085
cowplot::ggsave2(glue("{outdir}/co8_merged_knn_pl_suggestedk.pdf"), co8_merged_knn_pl_suggestedk, width = 5, height = 3)
cowplot::ggsave2(glue("{outdir}/co8_merged_knn_pl_noc0_umap_timepoint.pdf"), co8_merged_knn_pl_noc0_tp_umap_p, width = 4.62, height = 3)
cowplot::ggsave2(glue("{outdir}/co8_ploidy_umap.pdf"), co8_merged_knn_pl_p, width = 4, height = 3)
cowplot::ggsave2(glue("{outdir}/co8_merged_knn_pl_noc0_umap.pdf"), co8_merged_knn_pl_noc0_umap_p, width = 4.7, height = 3)
cowplot::ggsave2(glue("{outdir}/co8_merged_knn_pl_noc0_selected_hvg_gc.pdf"), co8_merged_knn_pl_noc0_selected_hvg_gc, width = 5, height = 7)


cowplot::ggsave2(glue("{outdir}/co8_primary_umap.pdf"), co8_tumor_primary_umap_p, width = 5, height = 5)

try(co8_merged_knn_pl_noc0_consensus_phylo <- plotPhylo(co8_merged_knn_pl_noc0, label = 'subclones', consensus = TRUE, group = 'timepoint'))
cowplot::ggsave2(glue("{outdir}/co8_merged_knn_pl_noc0_consensus_phylo.pdf"), co8_merged_knn_pl_noc0_consensus_phylo, width = 5, height = 5)
cowplot::ggsave2(glue("{outdir}/co8_merged_knn_pl_noc0_umap_p_ploidysubmedian.pdf"), co8_merged_knn_pl_noc0_umap_p_ploidysubmedian, width = 4, height = 3)


cowplot::ggsave2(glue("{outdir}/co5_merged_knn_pl_suggestedk.pdf"), co5_merged_knn_pl_suggestedk, width = 5, height = 3)
cowplot::ggsave2(glue("{outdir}/co5_merged_knn_pl_noc0_umap_timepoint.pdf"), co5_merged_knn_pl_noc0_tp_umap_p, width = 4.62, height = 3)
cowplot::ggsave2(glue("{outdir}/co5_merged_knn_pl_noc0_umap.pdf"), co5_merged_knn_pl_noc0_umap_p, width = 4.2, height = 3)
cowplot::ggsave2(glue("{outdir}/co5_merged_knn_pl_noc0_selected_hvg_gc.pdf"), co5_merged_knn_pl_noc0_selected_hvg_gc, width = 5, height = 7)
cowplot::ggsave2(glue("{outdir}/co5_primary_umap.pdf"), co5_tumor_primary_umap_p, width = 5, height = 5)



consensusPhylo(co5_merged_knn_pl_noc0) <- phytools::rotateNodes(consensusPhylo(co5_merged_knn_pl_noc0), c(19,21:27, 29:33))

try(co5_merged_knn_pl_noc0_consensus_phylo <- plotPhylo(co5_merged_knn_pl_noc0, label = 'subclones', consensus = TRUE, group = 'timepoint'))
cowplot::ggsave2(glue("{outdir}/co5_merged_knn_pl_noc0_consensus_phylo.pdf"), co5_merged_knn_pl_noc0_consensus_phylo, width = 5, height = 5)

pdf(glue("{outdir}/co5_merged_knn_pl_noc0_subclones.pdf"), width = 7, height = 8)
plotHeatmap(co5_merged_knn_pl_noc0, label = c('subclones', 'timepoint'), order_cells = 'consensus_tree', raster_quality = 5)
dev.off()

pdf(glue("{outdir}/co5_merged_knn_pl_noc0_subclones_integer.pdf"), width = 7, height = 8)
plotHeatmap(co5_merged_knn_pl_noc0, label = c('subclones', 'timepoint'), order_cells = 'consensus_tree', raster_quality = 5, assay = 'integer')
dev.off()



co5_merged_knn_pl_noc0_metrics_aftfilt %>%
  summarise(reads_total_mean = mean(reads_total),
            reads_total_sd = sd(reads_total),
            percent_dup_mean = mean(percentage_duplicates),
            percent_dup_sd = sd(percentage_duplicates)
  )
# reads_total_mean reads_total_sd percent_dup_mean percent_dup_sd
# 1          1043909       265948.3       0.08673491     0.01112969
cowplot::ggsave2(glue("{outdir}/co5_merged_knn_pl_noc0_umap_p_ploidysubmedian.pdf"), co5_merged_knn_pl_noc0_umap_p_ploidysubmedian, width = 4, height = 3)

pdf(glue("{outdir}/co5_merged_subclones_integer_integer_submedian.pdf"), width = 7, height = 8)
plotHeatmap(co5_merged_knn_pl_noc0, label = c('subclones', 'timepoint'), order_cells = 'consensus_tree', raster_quality = 4, assay = 'integer')
dev.off()


pdf(glue("{outdir}/co5_merged_knn_pl_noc0_consensus_ht.pdf"), width = 10, height = 5)
plotHeatmap(co5_merged_knn_pl_noc0, label = c('subclones'), consensus = TRUE, genes = co5_selected_genes)
dev.off()

pdf(glue("{outdir}/co5_merged_knn_pl_noc0_consensus_ht_integer.pdf"), width = 10, height = 5)
plotHeatmap(co5_merged_knn_pl_noc0, label = c('subclones'), consensus = TRUE, genes = co5_selected_genes, assay = 'integer', order_cells = 'consensus_tree')
dev.off()




brcan_metrics %>%
  summarise(reads_total_mean = mean(reads_total),
            reads_total_sd = sd(reads_total),
            percent_dup_mean = mean(percentage_duplicates),
            percent_dup_sd = sd(percentage_duplicates)
  )
# reads_total_mean reads_total_sd percent_dup_mean percent_dup_sd
#          1047467       304357.9       0.08537978      0.0168745
brcan_metrics %>% as.data.frame() %>%
  rstatix::kruskal_test(overdispersion ~ outlier)
#   .y.                n statistic    df        p method        
# * <chr>          <int>     <dbl> <int>    <dbl> <chr>         
# 1 overdispersion   890      62.4     1 2.77e-15 Kruskal-Wallis
brcan_metrics %>% as.data.frame() %>%
  rstatix::kruskal_test(breakpoint_count ~ outlier)
# 1 breakpoint_count   890      30.6     1 0.0000000312 Kruskal-Wallis

cowplot::ggsave2(glue("{outdir}/brcan_p_metrics.pdf"), brcan_p_metrics, width = 6, height = 4)

brcan_over_comparison_p <- plotOverdispersionComparison(brcan_knn)

cowplot::ggsave2(glue("{outdir}/brcan_overdisp_comparison.pdf"), brcan_over_comparison_p, width = 7, height = 5)
brcan_metrics_aftfilt %>%
  summarise(reads_total_mean = mean(reads_total),
            reads_total_sd = sd(reads_total),
            percent_dup_mean = mean(percentage_duplicates),
            percent_dup_sd = sd(percentage_duplicates)
  )
# reads_total_mean reads_total_sd percent_dup_mean percent_dup_sd
# 1          1071694       282527.2       0.08531321     0.01588034

cowplot::ggsave2(glue("{outdir}/brcan_primary_suggestedK_p.pdf"), brcan_primary_suggestedK_p, width = 5, height = 3)
janitor::tabyl(colData(brcan_knn)$subclones)
#  colData(brcan_knn)$subclones   n    percent
#                            c0   0 0.00000000
#                            c1  81 0.14972274
#                            c2  55 0.10166359
#                            c3  65 0.12014787
#                            c4  56 0.10351201
#                            c5  42 0.07763401
#                            c6  24 0.04436229
#                            c7  74 0.13678373
#                            c8 144 0.26617375
cowplot::ggsave2(glue("{outdir}/brcan_primary_umap.pdf"), brcan_knn_umap_p, width = 4, height = 3)
cowplot::ggsave2(glue("{outdir}/brcan_ploidy_umap.pdf"), brcan_knn_umap_ploidy_p, width = 4, height = 3)


cowplot::ggsave2(glue("{outdir}/tn17_10x_p_metrics.pdf"), tn17_10x_p_metrics, width = 6, height = 4)
tn17_10x_metrics %>%
  summarise(reads_total_mean = mean(reads_total),
            reads_total_sd = sd(reads_total),
            percent_dup_mean = mean(percentage_duplicates),
            percent_dup_sd = sd(percentage_duplicates)
  )
#   reads_total_mean reads_total_sd percent_dup_mean percent_dup_sd
# 1           628310       251573.4        0.1417568    0.003473852
tn17_10x_metrics %>% as.data.frame() %>%
  rstatix::kruskal_test(overdispersion ~ outlier)
# .y.                n statistic    df            p method        
# 1 overdispersion  1168      30.7     1 0.0000000297 Kruskal-Wallis
tn17_10x_metrics %>% as.data.frame() %>%
  rstatix::kruskal_test(breakpoint_count ~ outlier)
# 1 breakpoint_count  1168      40.9     1 1.62e-10 Kruskal-Wallis
tn17_10x_over_comparison_p <- plotOverdispersionComparison(tn17_10x_knn)
cowplot::ggsave2(glue("{outdir}/tn17_10x_overdisp_comparison.pdf"), tn17_10x_over_comparison_p, width = 7, height = 5)


tn17_10x_knn_umap_ploidy_p <- plotUmap(tn17_10x_knn_pl_reserv, label = 'ploidy') + scale_fill_jcolors_contin('rainbow')


tn17_10x_metrics_aftfilt %>%
  summarise(reads_total_mean = mean(reads_total),
            reads_total_sd = sd(reads_total),
            percent_dup_mean = mean(percentage_duplicates),
            percent_dup_sd = sd(percentage_duplicates)
  )
#   reads_total_mean reads_total_sd percent_dup_mean percent_dup_sd
# 1           629628       249971.6        0.1416812    0.003556581
cowplot::ggsave2(glue("{outdir}/tn17_10x_primary_suggestedK_p.pdf"), tn17_10x_primary_suggestedK_p, width = 5, height = 3)
janitor::tabyl(colData(tn17_10x_knn)$subclones)
# colData(tn17_10x_knn)$subclones   n    percent
# c0   0 0.00000000
# c1  67 0.06291080
# c2  30 0.02816901
# c3  57 0.05352113
# c4  90 0.08450704
# c5 378 0.35492958
# c6  99 0.09295775
# c7 344 0.32300469
cowplot::ggsave2(glue("{outdir}/tn17_10x_primary_umap.pdf"), tn17_10x_knn_umap_p, width = 4, height = 3)
cowplot::ggsave2(glue("{outdir}/tn17_10x_ploidy_umap.pdf"), tn17_10x_knn_umap_ploidy_p, width = 4, height = 3)
cowplot::ggsave2(glue("{outdir}/tn17_10x_knn_pl_umap_p_ploidysubmedian.pdf"), tn17_10x_knn_pl_umap_p_ploidysubmedian, width = 4, height = 3)


dlp_sa928_a90553_metrics %>%
  summarise(reads_total_mean = mean(reads_total),
            reads_total_sd = sd(reads_total),
            percent_dup_mean = mean(percentage_duplicates),
            percent_dup_sd = sd(percentage_duplicates)
  )
# reads_total_mean reads_total_sd percent_dup_mean percent_dup_sd
# 1         318533.1        81400.7       0.08774115    0.008828223
dlp_sa928_a90553_metrics %>% as.data.frame() %>%
  rstatix::kruskal_test(overdispersion ~ outlier)
# .y.                n statistic    df     p method        
# * <chr>          <int>     <dbl> <int> <dbl> <chr>         
#   1 overdispersion  6189  0.000429     1 0.983 Kruskal-Wallis
dlp_sa928_a90553_metrics %>% as.data.frame() %>%
  rstatix::kruskal_test(breakpoint_count ~ outlier)
# breakpoint_count  6189      249.     1 5.48e-56 Kruskal-Wallis

cowplot::ggsave2(glue("{outdir}/dlp_sa928_a90553_p_metrics.pdf"), dlp_sa928_a90553_p_metrics, width = 7, height = 4)

dlp_sa928_a90553_metrics_aftfilt %>%
  summarise(reads_total_mean = mean(reads_total),
            reads_total_sd = sd(reads_total),
            percent_dup_mean = mean(percentage_duplicates),
            percent_dup_sd = sd(percentage_duplicates)
  )
# reads_total_mean reads_total_sd percent_dup_mean percent_dup_sd
#         317681.8       80525.91       0.08772711    0.008838496

cowplot::ggsave2(glue("{outdir}/co8_ploidy_umap_rev.pdf"), co8_merged_knn_pl_p, width = 4, height = 3)
cowplot::ggsave2(glue("{outdir}/brcan_ploidy_umap_rev.pdf"), brcan_knn_umap_ploidy_p, width = 4, height = 3)

