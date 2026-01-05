library(tidyverse)
library(BiocParallel)
library(glue)
library("distances")
set.seed(123)
register(MulticoreParam(workers = 100, progressbar = F), default = T)
load("~/Projects/copykit/benchmark/cp_3tools/20230222_3tools_ploidy_comparison.RData")
cr_cutoff=0.3

ploidy_comp <- rbind(brcan_ploidy_comp %>% dplyr::filter(confidence_ratio>cr_cutoff) %>% dplyr::select(cellname, confidence_ratio,uber_id, facs_ploidy, SoS_ploidy, ploidy, hmm_ploidy),
                     pmtc6_ploidy_comp %>% dplyr::filter(confidence_ratio>cr_cutoff) %>% dplyr::select(cellname, confidence_ratio,uber_id, facs_ploidy, SoS_ploidy, ploidy, hmm_ploidy),
                     pmtc7_ploidy_comp %>% dplyr::filter(confidence_ratio>cr_cutoff) %>% dplyr::select(cellname, confidence_ratio,uber_id, facs_ploidy, SoS_ploidy, ploidy, hmm_ploidy),
                     co5_ploidy_comp %>% dplyr::filter(confidence_ratio>cr_cutoff) %>% dplyr::select(cellname, confidence_ratio,uber_id, facs_ploidy, SoS_ploidy, ploidy, hmm_ploidy),
                     co8_ploidy_comp %>% dplyr::filter(confidence_ratio>cr_cutoff) %>% dplyr::select(cellname, confidence_ratio,uber_id, facs_ploidy, SoS_ploidy, ploidy, hmm_ploidy)) 

ploidy_comp$ginkgo <- (ploidy_comp$SoS_ploidy-ploidy_comp$facs_ploidy)^2
ploidy_comp$scquantum <- (ploidy_comp$ploidy-ploidy_comp$facs_ploidy)^2
ploidy_comp$hmmcopy <- (ploidy_comp$hmm_ploidy-ploidy_comp$facs_ploidy)^2

ploidy_plot_df <- ploidy_comp %>% 
  gather(value="squared_error", key="tool", -cellname, -confidence_ratio,-uber_id, -facs_ploidy, -SoS_ploidy, -ploidy, -hmm_ploidy)


my_comparisons <- list( c("ginkgo", "scquantum"), c("ginkgo", "hmmcopy"), c("scquantum", "hmmcopy") )

ggpubr::ggviolin(ploidy_plot_df %>%
                   mutate(log_se=log(squared_error, base=0.5)),
                 x = "tool", y = "log_se",
                 color = "uber_id", palette = "jco",
                 add = "boxplot") + 
  ggpubr::stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  #ggpubr::stat_compare_means(label.y = 50)  +
  theme(axis.text.x = element_text(angle = 45, hjust=0.8),
        legend.position = "none") +
  facet_grid(~uber_id) 
cowplot::ggsave2(glue("{outputdir}/ploidy_compare_bysample_cr{cr_cutoff}_violin.pdf"),width = 14, height = 4)

ggpubr::ggviolin(ploidy_plot_df%>%dplyr::select(-uber_id)%>%
                   mutate(log_se=log(squared_error, base=0.5)), x = "tool", y = "log_se",
                 add = "boxplot") + 
  ggpubr::stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  ggpubr::stat_compare_means(label.y = 40) 
cowplot::ggsave2(glue("{outputdir}/ploidy_compare_bytool_cr{cr_cutoff}_violin.pdf"),width = 5, height = 4)

ploidy_comp$ginkgo <- abs(ploidy_comp$SoS_ploidy-ploidy_comp$facs_ploidy)
ploidy_comp$scquantum <- abs(ploidy_comp$ploidy-ploidy_comp$facs_ploidy)
ploidy_comp$hmmcopy <- abs(ploidy_comp$hmm_ploidy-ploidy_comp$facs_ploidy)

ploidy_plot_df2 <- ploidy_comp %>% 
  gather(value="absolute_err", key="tool", -cellname, -confidence_ratio,-uber_id, -facs_ploidy, -SoS_ploidy, -ploidy, -hmm_ploidy)


my_comparisons <- list( c("ginkgo", "scquantum"), c("ginkgo", "hmmcopy"), c("scquantum", "hmmcopy") )

ggpubr::ggviolin(ploidy_plot_df2,
                 x = "tool", y = "absolute_err",
                 color = "uber_id", palette = "jco",
                 add = "boxplot") + 
  ggpubr::stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  theme(axis.text.x = element_text(angle = 45, hjust=0.8),
        legend.position = "none") +
  facet_grid(~uber_id) 
cowplot::ggsave2(glue("{outputdir}/ploidy_compare_bysample_cr{cr_cutoff}_abs_violin.pdf"),width = 10, height = 4)
ggpubr::ggviolin(ploidy_plot_df2%>%dplyr::select(-uber_id), x = "tool", y = "absolute_err",
                 add = "boxplot") + 
  ggpubr::stat_compare_means(comparisons = my_comparisons)
cowplot::ggsave2(glue("{outputdir}/ploidy_compare_bytool_cr{cr_cutoff}_abs_violin.pdf"),width = 5, height = 4)

write_tsv(ploidy_plot_df %>%
  mutate(log_se=log(squared_error, base=0.5)) %>%
  inner_join(ploidy_plot_df2) %>%
  dplyr::rename(scquantum_ploidy=ploidy), glue("{outputdir}/0222_ploidy_compare_bytool.tsv"))

#############
############# try abs error normalized
#############
ploidy_comp$ginkgo <- abs((ploidy_comp$SoS_ploidy-ploidy_comp$facs_ploidy))
ploidy_comp$scquantum <- abs((ploidy_comp$ploidy-ploidy_comp$facs_ploidy))
ploidy_comp$hmmcopy <- abs((ploidy_comp$hmm_ploidy-ploidy_comp$facs_ploidy))

m <- as.matrix(ploidy_comp %>% dplyr::select(ginkgo, scquantum, hmmcopy))
m <- m/colSums(t(m))
m <- cbind(m, ploidy_comp$uber_id)
colnames(m)[4] <- "uber_id"

ploidy_plot_df <- m %>% 
  as_tibble() %>%
  gather(value="abs_error_norm", key="tool", -uber_id)
ploidy_plot_df$abs_error_norm <- as.numeric(ploidy_plot_df$abs_error_norm)

my_comparisons <- list( c("ginkgo", "scquantum"), c("ginkgo", "hmmcopy"), c("scquantum", "hmmcopy") )

ggpubr::ggviolin(ploidy_plot_df,
                 x = "tool", y = "abs_error_norm",
                 color = "uber_id", palette = "jco",
                 add = "boxplot") + 
  ggpubr::stat_compare_means(comparisons = my_comparisons,label = "p.signif")+ # Add pairwise comparisons p-value
  #ggpubr::stat_compare_means(label.y = 50)  +
  theme(axis.text.x = element_text(angle = 45, hjust=0.8),
        legend.position = "none") +
  facet_grid(~uber_id) 
cowplot::ggsave2(glue("{outputdir}/ploidy_compare_abs_error_norm_bysample_cr{cr_cutoff}_violin.pdf"),width = 8, height = 4)


#### Compare CV
cv <- function(seg){vapply(
  seg, function(z) {
    sd(z) / mean(z)
  },
  numeric(1)
)}

brcan_cv <- data.frame(
                       ginkgo = cv(brcan_cnmat[,colnames(brcan_hmm_cnmat_filt)]), 
                       hmmcopy = cv(brcan_hmm_cnmat_filt),
                       copykit = cv(assay(brcan_filt_knn,"integer"))
                       ) %>%
            rownames_to_column("cellname") %>%
            mutate(sample = "BRCAN")

co5_cv <- data.frame(
  hmmcopy = cv(co5_hmm_cnmat_filt),
  ginkgo = cv(co5_cnmat[,sapply(colnames(co5_hmm_cnmat_filt), function(c){gsub(x=c,pattern = "-",replacement = ".", fixed = T)})]), 
  copykit = cv(assay(co5_filt_knn,"integer"))
) %>%
  rownames_to_column("cellname") %>%
  mutate(sample = "CO5")

co8_cv <- data.frame(
  hmmcopy = cv(co8_hmm_cnmat_filt),
  ginkgo = cv(co8_cnmat[,sapply(colnames(co8_hmm_cnmat_filt), function(c){gsub(x=c,pattern = "-",replacement = ".", fixed = T)})]),
  copykit = cv(assay(co8_filt_knn,"integer"))
) %>%
  rownames_to_column("cellname") %>%
  mutate(sample = "CO8")

pmtc6_cv <- data.frame(
  hmmcopy = cv(pmtc6_hmm_cnmat_filt),
  ginkgo = cv(pmtc6_cnmat[,colnames(pmtc6_hmm_cnmat_filt)]), 
  copykit = cv(assay(pmtc6_filt_knn,"integer"))
) %>%
  rownames_to_column("cellname") %>%
  mutate(sample = "PMTC6")

pmtc7_cv <- data.frame(
  hmmcopy = cv(pmtc7_hmm_cnmat_filt),
  ginkgo = cv(pmtc7_cnmat[,colnames(pmtc7_hmm_cnmat_filt)]), 
  copykit = cv(assay(pmtc7_filt_knn,"integer"))
) %>%
  rownames_to_column("cellname") %>%
  mutate(sample = "PMTC7")

cv_plot_df <- rbind(brcan_cv, 
                    co5_cv, 
                    co8_cv, 
                    pmtc6_cv, 
                    pmtc7_cv) %>% dplyr::select(cellname, sample, ginkgo,copykit,hmmcopy) %>%
  gather(value="cv", key="tool", -sample, -cellname)

my_comparisons <- list( c("ginkgo", "copykit"), c("ginkgo", "hmmcopy"), c("copykit", "hmmcopy") )

ggpubr::ggboxplot(cv_plot_df,
                  x = "tool", y = "cv",
                  color = "sample", palette = "jco",
                  add = "jitter") + 
  ggpubr::stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  # ggpubr::stat_compare_means(label.y = 50)  +
  theme(axis.text.x = element_text(angle = 45, hjust=0.8),
        legend.position = "none") +
  facet_grid(~sample) 
cowplot::ggsave2(glue("{outputdir}/cv_compare_bysample.pdf"),width = 14, height = 4)

ggpubr::ggboxplot(cv_plot_df%>%dplyr::select(-sample), x = "tool", y = "cv",
                  add = "jitter") + 
  ggpubr::stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
  ggpubr::stat_compare_means(label.y = 1.1) 
cowplot::ggsave2(glue("{outputdir}/cv_compare_bytool.pdf"),width = 5, height = 4)

ggpubr::ggviolin(cv_plot_df,
                 x = "tool", y = "cv",
                 color = "sample", palette = "jco", add="boxplot") + 
  ggpubr::stat_compare_means(comparisons = my_comparisons,label = "p.signif")+ # Add pairwise comparisons p-value
  # ggpubr::stat_compare_means(label.y = 50)  +
  theme(axis.text.x = element_text(angle = 45, hjust=0.8),
        legend.position = "none") +
  facet_grid(~sample) 
cowplot::ggsave2(glue("{outputdir}/cv_compare_bysample_violin.pdf"),width = 8, height = 4)


ggpubr::ggviolin(cv_plot_df%>%dplyr::select(-sample), x = "tool", y = "cv",
                 add = "boxplot") + 
  ggpubr::stat_compare_means(comparisons = my_comparisons,label = "p.signif") + # Add pairwise comparisons p-value
  ggpubr::stat_compare_means(label.y = 1.5) 
cowplot::ggsave2(glue("{outputdir}/cv_compare_bytool_violin.pdf"),width = 5, height = 4)
write_tsv(cv_plot_df, glue("{outputdir}/0222_cv_compare_bytool.tsv"))

#### Compare overdispersion
overdispersion <- function(dat_bin) {
  unlist(BiocParallel::bplapply(dat_bin, function(v){
    # 3 elements, 2 differences, can find a standard deviation
    stopifnot(length(v) >= 3)
    # Differences between pairs of values
    y <- v[-1]
    x <- v[-length(v)]
    # Normalize the differences using the sum. The result should be around zero,
    # plus or minus square root of the index of dispersion
    vals.unfiltered <- (y - x) / sqrt(y + x)
    # Remove divide by zero cases, and--considering this is supposed to be count
    # data--divide by almost-zero cases
    vals <- vals.unfiltered[y + x >= 1]
    # Check that there's anything left
    stopifnot(length(vals) >= 2)
    # Assuming most of the normalized differences follow a normal distribution,
    # estimate the standard deviation
    val.sd <- l2e.normal.sd(vals)
    # Square this standard deviation to obtain an estimate of the index of
    # dispersion
    iod <- val.sd^2
    # subtract one to get the overdispersion criteria
    iod.over <- iod - 1
    # normalizing by mean bincounts
    iod.norm <- iod.over / mean(v)
    return(iod.norm)
  }))
}

brcan_ovdp <- data.frame(
  ginkgo = overdispersion(brcan_cnmat[,colnames(brcan_hmm_cnmat_filt)]), 
  hmmcopy = overdispersion(brcan_hmm_cnmat_filt),
  copykit = overdispersion(assay(brcan_filt_knn,"integer"))
) %>%
  rownames_to_column("cellname") %>%
  mutate(sample = "BRCAN")

co5_ovdp <- data.frame(
  hmmcopy = overdispersion(co5_hmm_cnmat_filt),
  ginkgo = overdispersion(co5_cnmat[,sapply(colnames(co5_hmm_cnmat_filt), function(c){gsub(x=c,pattern = "-",replacement = ".", fixed = T)})]), 
  copykit = overdispersion(assay(co5_filt_knn,"integer"))
) %>%
  rownames_to_column("cellname") %>%
  mutate(sample = "CO5")

co8_ovdp <- data.frame(
  hmmcopy = overdispersion(co8_hmm_cnmat_filt),
  ginkgo = overdispersion(co8_cnmat[,sapply(colnames(co8_hmm_cnmat_filt), function(c){gsub(x=c,pattern = "-",replacement = ".", fixed = T)})]),
  copykit = overdispersion(assay(co8_filt_knn,"integer"))
) %>%
  rownames_to_column("cellname") %>%
  mutate(sample = "CO8")

pmtc6_ovdp <- data.frame(
  hmmcopy = overdispersion(pmtc6_hmm_cnmat_filt),
  ginkgo = overdispersion(pmtc6_cnmat[,colnames(pmtc6_hmm_cnmat_filt)]), 
  copykit = overdispersion(assay(pmtc6_filt_knn,"integer"))
) %>%
  rownames_to_column("cellname") %>%
  mutate(sample = "PMTC6")

pmtc7_ovdp <- data.frame(
  hmmcopy = overdispersion(pmtc7_hmm_cnmat_filt),
  ginkgo = overdispersion(pmtc7_cnmat[,colnames(pmtc7_hmm_cnmat_filt)]), 
  copykit = overdispersion(assay(pmtc7_filt_knn,"integer"))
) %>%
  rownames_to_column("cellname") %>%
  mutate(sample = "PMTC7")

ovdp_plot_df <- rbind(brcan_ovdp, 
                    co5_ovdp, 
                    co8_ovdp, 
                    pmtc6_ovdp, 
                    pmtc7_ovdp) %>% dplyr::select(cellname, sample, ginkgo,copykit,hmmcopy) %>%
  gather(value="overdispersion", key="tool", -sample, -cellname)

my_comparisons <- list( c("ginkgo", "copykit"), c("ginkgo", "hmmcopy"), c("copykit", "hmmcopy") )

ggpubr::ggboxplot(ovdp_plot_df,
                  x = "tool", y = "overdispersion",
                  color = "sample", palette = "jco",
                  add = "jitter") + 
  ggpubr::stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  # ggpubr::stat_compare_means(label.y = 50)  +
  theme(axis.text.x = element_text(angle = 45, hjust=0.8),
        legend.position = "none") +
  facet_grid(~sample) 
cowplot::ggsave2(glue("{outputdir}/ovdp_compare_bysample.pdf"),width = 14, height = 4)

ggpubr::ggboxplot(ovdp_plot_df%>%dplyr::select(-sample), x = "tool", y = "overdispersion",
                  add = "jitter") + 
  ggpubr::stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
  ggpubr::stat_compare_means(label.y = 0.3) 
cowplot::ggsave2(glue("{outputdir}/ovdp_compare_bytool.pdf"),width = 5, height = 4)

ggpubr::ggviolin(ovdp_plot_df,
                 x = "tool", y = "overdispersion",
                 color = "sample", palette = "jco", add="boxplot") + 
  ggpubr::stat_compare_means(comparisons = my_comparisons,label = "p.signif")+ # Add pairwise comparisons p-value
  # ggpubr::stat_compare_means(label.y = 50)  +
  theme(axis.text.x = element_text(angle = 45, hjust=0.8),
        legend.position = "none") +
  facet_grid(~sample) 
cowplot::ggsave2(glue("{outputdir}/ovdp_compare_bysample_violin.pdf"),width = 8, height = 4)

write_tsv(ovdp_plot_df, glue("{outputdir}/0222_ovdp_compare_bytool.tsv"))

