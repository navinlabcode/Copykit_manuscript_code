# copykit manuscript samples v0.1.2
#SETUP

library(devtools)
library(copykit)
library(BiocParallel)
library(glue)
library(dplyr)
library(ggplot2)

set.seed(17)
register(MulticoreParam(workers = 50), default = T)
 
varbin=T
outdir="/mnt/USR1/junke/Projects/copykit/rerun_manu/singu_0424"
sink(glue("{outdir}/output.log"))

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

#BL1 -----


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PMTC6 -----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Running samples with runVarbin, commented out to avoid accidental running
if(varbin){
  pmtc6_liver <- runVarbin("/mnt/lab/users/dminussi/projects/copy_number_pipeline_hg38/test_tumor_data/PMTC6/marked/", remove_Y = T)
}
# 
# aneuploidy cells and plot
pmtc6_liver_filt <- findAneuploidCells(pmtc6_liver)

pdf(glue("{outdir}/pmtc6_liver_knn_is_aneuploid.pdf"), width = 7, height = 8)
plotHeatmap(pmtc6_liver_filt, label = c('is_aneuploid'),
            row_split = 'is_aneuploid',
            order_cells = 'hclust',
            n_threads = 50)
dev.off()
pmtc6_liver_filt <- pmtc6_liver_filt[,colData(pmtc6_liver_filt)$is_aneuploid == TRUE]

# filtering outliers and plot
pmtc6_liver_filt <- findOutliers(pmtc6_liver_filt)

pdf(glue("{outdir}/pmtc6_liver_knn_outlier.pdf"), width = 7, height = 8)
plotHeatmap(pmtc6_liver_filt, label = c('outlier'),
            row_split = 'outlier',
            order_cells = 'hclust',
            n_threads = 50)
dev.off()


pmtc6_liver_filt <- runMetrics(pmtc6_liver_filt)

pmtc6_p_metrics <- plotMetrics(
  pmtc6_liver_filt,
  metric = c(
    "reads_total",
    "percentage_duplicates",
    "overdispersion",
    "breakpoint_count"
  ), label = 'outlier',
  dodge.width = 0.8,
  ncol = 2
) + scale_fill_manual(values = c("TRUE" = "#DA614D",
                                 "FALSE" = "#5F917A"))

pmtc6_metrics <- as.data.frame(colData(pmtc6_liver_filt))
#########
# plot metrics here
#########
pmtc6_liver_filt <- pmtc6_liver_filt[,colData(pmtc6_liver_filt)$outlier == FALSE]

pmtc6_liver_knn <- knnSmooth(pmtc6_liver_filt)
pmtc6_liver_knn <- calcSmoothOver(pmtc6_liver_knn)

pmtc6_liver_knn <- calcInteger(pmtc6_liver_knn, method = 'scquantum', assay = 'smoothed_bincounts')

table(colData(pmtc6_liver_knn)$ploidy_score < 0.2)
# FALSE  TRUE 
# 70   723 
pmtc6_liver_knn_pl <- pmtc6_liver_knn[,colData(pmtc6_liver_knn)$ploidy_score < 0.2]

plotMetrics(pmtc6_liver_knn_pl, metric = 'ploidy', label = 'ploidy_score')

pmtc6_metrics_aftfilt <- as.data.frame(colData(pmtc6_liver_knn_pl))


pmtc6_liver_knn_pl <- runUmap(pmtc6_liver_knn_pl)
pmtc6_liver_knn_pl <- findSuggestedK(pmtc6_liver_knn_pl, k_range = 10:20)
pmtc6_primary_suggestedK_p <- plotSuggestedK(pmtc6_liver_knn_pl)

pmtc6_liver_knn_pl <- findClusters(pmtc6_liver_knn_pl)

# * Umap ----
pmtc6_liver_knn_pl_umap_p <- plotUmap(pmtc6_liver_knn_pl, label = 'subclones')


pmtc6_liver_knn_pl <- calcConsensus(pmtc6_liver_knn_pl)
pmtc6_liver_knn_pl <- runConsensusPhylo(pmtc6_liver_knn_pl, root = 'mrca')

pdf(glue("{outdir}/pmtc6_primary_subclones.pdf"), width = 7, height = 8)
plotHeatmap(pmtc6_liver_knn_pl, label = c('subclones'), order_cells = 'consensus_tree', raster_quality = 4)
dev.off()

pdf(glue("{outdir}/pmtc6_primary_subclones_integer.pdf"), width = 7, height = 8)
plotHeatmap(pmtc6_liver_knn_pl, label = c('subclones'), order_cells = 'consensus_tree', raster_quality = 4, assay = 'integer')
dev.off()

pmtc6_metrics_aftfilt <- as.data.frame(colData(pmtc6_liver_knn_pl))

pmtc6_median_sub_ploidy <- pmtc6_metrics_aftfilt %>%
  group_by(subclones) %>%
  dplyr::summarise(median_pl_subclones = median(ploidy))

colData(pmtc6_liver_knn_pl)$ploidy <-
  pmtc6_median_sub_ploidy$median_pl_subclones[match(colData(pmtc6_liver_knn_pl)$subclones, pmtc6_median_sub_ploidy$subclones)]

pmtc6_liver_knn_pl <- calcInteger(pmtc6_liver_knn_pl, method = 'metadata')

pdf(glue("{outdir}/pmtc6_primary_subclones_integer_submedian.pdf"), width = 7, height = 8)
plotHeatmap(pmtc6_liver_knn_pl, label = c('subclones'), order_cells = 'consensus_tree', raster_quality = 4, assay = 'integer')
dev.off()

pmtc6_liver_knn_pl_umap_p_ploidysubmedian <-plotUmap(pmtc6_liver_knn_pl, label = 'ploidy')

# * L-sectors ----
colData(pmtc6_liver_knn_pl)$L_foci <- stringr::str_extract(colData(pmtc6_liver_knn_pl)$sample, "(L[0-9]+L[0-9]+|L[0-9]+)")

colData(pmtc6_liver_knn_pl)$L_foci <- stringr::str_replace_all(colData(pmtc6_liver_knn_pl)$L_foci,
                                                               "L", "S")

pmtc6_liver_knn_pl_umap_l_foci_p <- plotUmap(pmtc6_liver_knn_pl, label = 'L_foci')

# pie charts of locations
pmtc6_focis <- unique(colData(pmtc6_liver_knn_pl)$L_foci)

pmtc6_foci_pies <- list()

for (i in seq_along(pmtc6_focis)) {
  
  df <- as.data.frame(colData(pmtc6_liver_knn_pl)) %>%
    dplyr::filter(L_foci == pmtc6_focis[i])
  
  pmtc6_foci_pies[[i]] <- ggplot(df) +
    geom_bar(aes(x = "", y = L_foci, fill = subclones),
             stat = 'identity') +
    theme_void() +
    scale_fill_manual(values = subclones_pal(),
                      limits = force) +
    coord_polar(theta = "y") +
    ggtitle(pmtc6_focis[i])
  
}



my_theme <- list(
  ggplot2::theme(
    axis.title.x = element_text(colour = "gray28", size = 20),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(colour = "gray28", size = 20),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    legend.text = element_text(size = 14),
    panel.border = element_rect(fill = NA, color = 'black')
  ),
  xlab("umap1"),
  ylab("umap2")
)

pmtc6_lfoci_umap <- as.data.frame(reducedDim(pmtc6_liver_knn_pl, 'umap'))

pmtc6_lfoci_umap$L_foci <- colData(pmtc6_liver_knn_pl)$L_foci
pmtc6_lfoci_umap$subclones <- as.character(colData(pmtc6_liver_knn_pl)$subclones)

pmtc6_lfoci_umap$subclones <- ifelse(!duplicated(pmtc6_lfoci_umap$subclones), pmtc6_lfoci_umap$subclones, "")


pdf(glue("{outdir}/pmtc6_liver_knn_pl_cs_ht.pdf"), width = 7,height = 4)
plotHeatmap(pmtc6_liver_knn_pl, label = 'subclones', consensus = T, group = 'L_foci',
            order_cells = 'consensus_tree', raster_quality = 4,
            genes = c(
              "FHIT",
              "CUX1",
              "WWOX",
              "FOXO1",
              "APC",
              "BCL2",
              "KRAS",
              "PGR",
              "PDGFRA",
              "ROS1",
              "PTPN11",
              "CCND1",
              "BTG1",
              "FGFR3",
              "PTEN",
              "FGFR2"
            ))
dev.off()

pmtc6_liver_knn_pl <- calcConsensus(pmtc6_liver_knn_pl, assay = 'integer')

pdf(glue("{outdir}/pmtc6_liver_knn_pl_cs_ht_integer.pdf"), width = 7,height = 4)
plotHeatmap(pmtc6_liver_knn_pl, label = 'subclones', consensus = T, group = 'L_foci',
            order_cells = 'consensus_tree', raster_quality = 4, assay = 'integer',
            genes = c(
              "FHIT",
              "CUX1",
              "WWOX",
              "FOXO1",
              "APC",
              "BCL2",
              "KRAS",
              "PGR",
              "PDGFRA",
              "ROS1",
              "PTPN11",
              "CCND1",
              "BTG1",
              "FGFR3",
              "PTEN",
              "FGFR2"
            ))
dev.off()





# BL2 ----

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PMTC7 -----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(varbin){
 pmtc7_liver <- runVarbin("/mnt/lab/users/dminussi/projects/copy_number_pipeline_hg38/test_tumor_data/PMTC7/marked/", remove_Y = T)
}

pmtc7_liver_filt <- findAneuploidCells(pmtc7_liver)

pdf(glue("{outdir}/pmtc7_liver_is_aneuploid.pdf"), width = 7, height = 8)
plotHeatmap(pmtc7_liver_filt, label = c('is_aneuploid'),
            row_split = 'is_aneuploid',
            n_threads = 50,
            order_cells = 'hclust',
            raster_quality = 4)
dev.off()
pmtc7_liver_filt <- pmtc7_liver_filt[,colData(pmtc7_liver_filt)$is_aneuploid == TRUE]

pmtc7_liver_filt <- findOutliers(pmtc7_liver_filt)

pdf(glue("{outdir}/pmtc7_liver_filt_outlier.pdf"), width = 7, height = 8)
plotHeatmap(pmtc7_liver_filt, label = c('outlier'),
            row_split = 'outlier',
            n_threads = 50,
            order_cells = 'hclust',
            raster_quality = 4)
dev.off()

pmtc7_liver_filt <- runMetrics(pmtc7_liver_filt)

pmtc7_p_metrics <- plotMetrics(
  pmtc7_liver_filt,
  metric = c(
    "reads_total",
    "percentage_duplicates",
    "overdispersion",
    "breakpoint_count"
  ), label = 'outlier',
  ncol = 2,
  dodge.width = .8
) + scale_fill_manual(values = c("TRUE" = "#DA614D",
                                 "FALSE" = "#5F917A"))


pmtc7_metrics <- as.data.frame(colData(pmtc7_liver_filt))


pmtc7_liver_filt <- pmtc7_liver_filt[,colData(pmtc7_liver_filt)$outlier == FALSE]

pmtc7_liver_knn <- knnSmooth(pmtc7_liver_filt)

pmtc7_liver_knn <- calcSmoothOver(pmtc7_liver_knn)

pmtc7_liver_knn <- calcInteger(pmtc7_liver_knn, method = 'scquantum', assay = 'smoothed_bincounts')

table(colData(pmtc7_liver_knn)$ploidy_score < 0.2)
# FALSE  TRUE 
# 78  1189 
pmtc7_liver_knn_pl <- pmtc7_liver_knn[,colData(pmtc7_liver_knn)$ploidy_score < 0.2]

plotMetrics(pmtc7_liver_knn_pl, metric = 'ploidy', label = 'ploidy_score')

pmtc7_metrics_aftfilt <- as.data.frame(colData(pmtc7_liver_knn_pl))


pmtc7_liver_knn_pl <- runUmap(pmtc7_liver_knn_pl)
pmtc7_liver_knn_pl <- findSuggestedK(pmtc7_liver_knn_pl)
pmtc7_primary_suggestedK_p <- plotSuggestedK(pmtc7_liver_knn_pl)


pmtc7_liver_knn_pl <- findClusters(pmtc7_liver_knn_pl)

pmtc7_liver_knn_pl <- pmtc7_liver_knn_pl[,colData(pmtc7_liver_knn_pl)$subclones != 'c0']

pmtc7_liver_knn_pl_umap_p <- plotUmap(pmtc7_liver_knn_pl, label = 'subclones')

pmtc7_liver_knn_pl <- calcConsensus(pmtc7_liver_knn_pl)
pmtc7_liver_knn_pl <- runConsensusPhylo(pmtc7_liver_knn_pl, root = 'mrca')

pdf(glue("{outdir}/pmtc7_primary_subclones.pdf"), width = 7, height = 8)
plotHeatmap(pmtc7_liver_knn_pl, label = c('subclones'), order_cells = 'consensus_tree', raster_quality = 4)
dev.off()

pdf(glue("{outdir}/pmtc7_primary_subclones_integer.pdf"), width = 7, height = 8)
plotHeatmap(pmtc7_liver_knn_pl, label = c('subclones'), assay = 'integer',order_cells = 'consensus_tree', raster_quality = 4)
dev.off()


pmtc7_metrics_aftfilt <- as.data.frame(colData(pmtc7_liver_knn_pl))

pmtc7_median_sub_ploidy <- pmtc7_metrics_aftfilt %>%
  group_by(subclones) %>%
  dplyr::summarise(median_pl_subclones = median(ploidy))

colData(pmtc7_liver_knn_pl)$ploidy <-
  pmtc7_median_sub_ploidy$median_pl_subclones[match(colData(pmtc7_liver_knn_pl)$subclones, pmtc7_median_sub_ploidy$subclones)]

pmtc7_liver_knn_pl <- calcInteger(pmtc7_liver_knn_pl, method = 'metadata')

pmtc7_liver_knn_pl_umap_p_ploidysubmedian <-plotUmap(pmtc7_liver_knn_pl, label = 'ploidy')

pdf(glue("{outdir}/pmtc7_primary_subclones_integer_submedian.pdf"), width = 7, height = 8)
plotHeatmap(pmtc7_liver_knn_pl, label = c('subclones'), order_cells = 'consensus_tree', raster_quality = 4, assay = 'integer')
dev.off()

# * L-sectors ----
colData(pmtc7_liver_knn_pl)$L_foci <- stringr::str_extract(colData(pmtc7_liver_knn_pl)$sample, "L[0-9]+")
colData(pmtc7_liver_knn_pl)$L_foci <- stringr::str_replace(colData(pmtc7_liver_knn_pl)$L_foci, "L", "S")

pmtc7_liver_knn_pl_umap_l_foci_p <- plotUmap(pmtc7_liver_knn_pl, label = 'L_foci')

# pie charts of locations
pmtc7_focis <- unique(colData(pmtc7_liver_knn_pl)$L_foci)

pmtc7_foci_pies <- list()

for (i in seq_along(pmtc7_focis)) {
  
  df <- as.data.frame(colData(pmtc7_liver_knn_pl)) %>%
    dplyr::filter(L_foci == pmtc7_focis[i])
  
  pmtc7_foci_pies[[i]] <- ggplot(df) +
    geom_bar(aes(x = "", y = L_foci, fill = subclones),
             stat = 'identity') +
    theme_void() +
    scale_fill_manual(values = subclones_pal(),
                      limits = force) +
    coord_polar(theta = "y") +
    ggtitle(pmtc7_focis[i])
  
}


pmtc7_lfoci_umap <- as.data.frame(reducedDim(pmtc7_liver_knn_pl, 'umap'))

pmtc7_lfoci_umap$L_foci <- colData(pmtc7_liver_knn_pl)$L_foci
pmtc7_lfoci_umap$subclones <- as.character(colData(pmtc7_liver_knn_pl)$subclones)

pmtc7_lfoci_umap$subclones <- ifelse(!duplicated(pmtc7_lfoci_umap$subclones), pmtc7_lfoci_umap$subclones, "")


pdf(glue("{outdir}/pmtc7_liver_knn_pl_cs_ht.pdf"), width = 8,height = 5)
plotHeatmap(pmtc7_liver_knn_pl, label = 'subclones', consensus = TRUE, group = 'L_foci',
            order_cells = 'consensus_tree', raster_quality = 4,
            genes = c(
              "CCND1",
              "TP53",
              "GATA1",
              "SOX2",
              "ERBB2",
              "MYC",
              "RB1",
              "PIK3CA",
              "MMP3",
              "BRCA1",
              "BRCA2",
              "FHIT"
            ))
dev.off()

pmtc7_liver_knn_pl <- calcConsensus(pmtc7_liver_knn_pl, assay = 'integer')

pdf(glue("{outdir}/pmtc7_liver_knn_pl_cs_ht_integer.pdf"), width = 8,height = 5)
plotHeatmap(pmtc7_liver_knn_pl, label = 'subclones', consensus = TRUE, group = 'L_foci',
            order_cells = 'consensus_tree', raster_quality = 4, assay = 'integer',
            genes = c(
              "CCND1",
              "TP53",
              "GATA1",
              "SOX2",
              "ERBB2",
              "MYC",
              "RB1",
              "PIK3CA",
              "MMP3",
              "BRCA1",
              "BRCA2",
              "FHIT"
            ))
dev.off()


# * foci diversity ----
pmtc7_meta <- colData(pmtc7_liver_knn_pl)

pmtc7_meta <- as.data.frame(pmtc7_meta)
pmtc7_meta <- pmtc7_meta[c('subclones', 'L_foci')]




# BM1 ----

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PMTC1 ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# * Breast FFPE ----

if(varbin){
 pmtc1_breast <- runVarbin("/mnt/lab/users/dminussi/projects/copy_number_pipeline_hg38/test_tumor_data/PMTC1/breast/marked/", remove_Y = TRUE)
}

pmtc1_breast <- findOutliers(pmtc1_breast, resolution = 0.8)

pmtc1_breast <- findAneuploidCells(pmtc1_breast)

pdf(glue("{outdir}/pmtc1_breast_outlier_is_aneuploid.pdf"), width = 7, height = 8)
plotHeatmap(pmtc1_breast, label = c("outlier", 'is_aneuploid'),
            row_split = 'outlier', n_threads = 50, order_cells = 'hclust', raster_quality = 4)
dev.off()

pmtc1_breast_filt <- pmtc1_breast[,colData(pmtc1_breast)$outlier == FALSE]
pmtc1_breast_filt <- pmtc1_breast_filt[,colData(pmtc1_breast_filt)$is_aneuploid == TRUE]

pmtc1_breast_filt <- runMetrics(pmtc1_breast_filt)

pmtc1_breast_knn <- knnSmooth(pmtc1_breast_filt)
colData(pmtc1_breast_knn)$overdispersion <- colData(pmtc1_breast_filt)$overdispersion
pmtc1_breast_knn <- calcSmoothOver(pmtc1_breast_knn)


colData(pmtc1_breast_knn)$timepoint <- 'breast'
pmtc1_breast_knn_reserv <- pmtc1_breast_knn
# ~~~~~~~~~~~~~~~~~~~~~~~


# * liver ----
if(varbin){
 pmtc1_liver <- runVarbin("/mnt/lab/users/dminussi/projects/copy_number_pipeline_hg38/test_tumor_data/PMTC1/liver/marked/", remove_Y = TRUE)
}

pmtc1_liver <- findOutliers(pmtc1_liver)

pmtc1_liver <- findAneuploidCells(pmtc1_liver)

pdf(glue("{outdir}/pmtc1_liver_outlier_is_aneuploid.pdf"), width = 7, height = 8)
plotHeatmap(pmtc1_liver, label = c("outlier", 'is_aneuploid'),
            row_split = 'outlier', n_threads = 50, order_cells = 'hclust', raster_quality = 4)
dev.off()

pmtc1_liver_filt <- pmtc1_liver[,colData(pmtc1_liver)$outlier == FALSE]
pmtc1_liver_filt <- pmtc1_liver_filt[,colData(pmtc1_liver_filt)$is_aneuploid == TRUE]

pmtc1_liver_knn <- knnSmooth(pmtc1_liver_filt)

colData(pmtc1_liver_knn)$timepoint <- 'liver'

# * pleural ----
if(varbin){
 pmtc1_pleural <- runVarbin("/mnt/lab/users/dminussi/projects/copy_number_pipeline_hg38/test_tumor_data/PMTC1/pleural_effusion/marked/", remove_Y = TRUE)
}
 
pmtc1_pleural <- findOutliers(pmtc1_pleural)

pmtc1_pleural <- findAneuploidCells(pmtc1_pleural)

pdf(glue("{outdir}/pmtc1_pleural_outlier_is_aneuploid.pdf"), width = 7, height = 8)
plotHeatmap(pmtc1_pleural, label = c("outlier", 'is_aneuploid'),
            row_split = 'outlier', n_threads = 50, order_cells = 'hclust', raster_quality = 4)
dev.off()

pmtc1_pleural_filt <- pmtc1_pleural[,colData(pmtc1_pleural)$outlier == FALSE]
pmtc1_pleural_filt <- pmtc1_pleural_filt[,colData(pmtc1_pleural_filt)$is_aneuploid == TRUE]

pmtc1_pleural_knn <- knnSmooth(pmtc1_pleural_filt)

colData(pmtc1_pleural_knn)$timepoint <- 'pleural'


# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# pmtc1 merged ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~

# setting NULL to the extra data frame info
colData(pmtc1_breast_knn)$overdispersion_smooth <- NULL
colData(pmtc1_breast_knn)$overdispersion <- NULL
colData(pmtc1_breast_knn)$breakpoint_count <- NULL

pmtc1_merged_knn <- SingleCellExperiment::cbind(pmtc1_breast_knn,
                                                pmtc1_liver_knn,
                                                pmtc1_pleural_knn)


pmtc1_merged_knn <- calcInteger(pmtc1_merged_knn, method = 'scquantum', assay = 'smoothed_bincounts')
plotMetrics(pmtc1_merged_knn, 'ploidy', label = 'timepoint')
plotMetrics(pmtc1_merged_knn, 'ploidy', label = 'ploidy_score')
table(colData(pmtc1_merged_knn)$ploidy_score < 0.4)
pmtc1_merged_knn_pl <- pmtc1_merged_knn[,colData(pmtc1_merged_knn)$ploidy_score < 0.4]

pmtc1_metrics_aftfilt <- as.data.frame(colData(pmtc1_merged_knn_pl))


pmtc1_merged_knn_pl <- runUmap(pmtc1_merged_knn_pl)

pmtc1_merged_knn_pl <- findSuggestedK(pmtc1_merged_knn_pl)
pmtc1_merged_suggestedk <- plotSuggestedK(pmtc1_merged_knn_pl)


pmtc1_merged_knn_pl <- findClusters(pmtc1_merged_knn_pl)

pmtc1_merged_knn_pl <- calcConsensus(pmtc1_merged_knn_pl)
pmtc1_merged_knn_pl <- runConsensusPhylo(pmtc1_merged_knn_pl)

pdf(glue("{outdir}/pmtc1_merged_subclones_withdoublets.pdf"), width = 7, height = 8)
plotHeatmap(pmtc1_merged_knn_pl, label = c('subclones', 'timepoint'), order_cells = 'consensus_tree')
dev.off()

pmtc1_merged_umap_p <- plotUmap(pmtc1_merged_knn_pl, label = 'subclones')

# removing cluster of doublets
pmtc1_merged_knn_pl_nod <- pmtc1_merged_knn_pl[,colData(pmtc1_merged_knn_pl)$subclones != 'c1']

# reclustering after doublet removal
pmtc1_merged_knn_pl_nod <- runUmap(pmtc1_merged_knn_pl_nod, n_neighbors = 15)

pmtc1_merged_knn_pl_nod <- findSuggestedK(pmtc1_merged_knn_pl_nod, k_range = 5:15)
pmtc1_merged_nod_suggestedk <- plotSuggestedK(pmtc1_merged_knn_pl)

pmtc1_merged_knn_pl_nod <- findClusters(pmtc1_merged_knn_pl_nod)

pmtc1_merged_knn_pl_nod_noc0 <- pmtc1_merged_knn_pl_nod[,colData(pmtc1_merged_knn_pl_nod)$subclones != 'c0']

pmtc1_merged_tp_umap_p <- plotUmap(pmtc1_merged_knn_pl_nod_noc0, label = 'timepoint')

pmtc1_merged_umap_p <- plotUmap(pmtc1_merged_knn_pl_nod_noc0, label = 'subclones')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PMTC1 tree ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# pmtc1_breast <- runUmap(pmtc1_breast)
#
# pmtc1_breast <- findSuggestedK(pmtc1_breast)
# pmctc1_breast <- findClusters(pmtc1_breast)

colData(pmtc1_breast)$subclones <- 'c1'

pmtc1_breast <- calcConsensus(pmtc1_breast)

pmtc1_breast <- inferMrca(pmtc1_breast)

pmtc1_merged_knn_pl_nod_noc0 <- calcConsensus(pmtc1_merged_knn_pl_nod_noc0)

pmtc1_merged_knn_pl_nod_noc0 <- runConsensusPhylo(pmtc1_merged_knn_pl_nod_noc0, root = 'user',
                                                  root_user = metadata(pmtc1_breast)$inferred_mrca)


pdf(glue("{outdir}/pmtc1_merged_subclones.pdf"), width = 7, height = 8)
plotHeatmap(pmtc1_merged_knn_pl_nod_noc0, label = c('subclones', 'timepoint'), order_cells = 'consensus_tree', raster_quality = 4)
dev.off()

pdf(glue("{outdir}/pmtc1_merged_subclones_integer.pdf"), width = 7, height = 8)
plotHeatmap(pmtc1_merged_knn_pl_nod_noc0, label = c('subclones', 'timepoint'), order_cells = 'consensus_tree', raster_quality = 4, assay = 'integer')
dev.off()

table(colData(pmtc1_merged_knn_pl_nod_noc0)$subclones, colData(pmtc1_merged_knn_pl_nod_noc0)$timepoint)


pmtc1_merged_consensus_phylo <- plotPhylo(pmtc1_merged_knn_pl_nod_noc0, label = 'subclones', consensus = TRUE)


pdf(glue("{outdir}/pmtc1_merged_consensus_ht.pdf"), width = 8, height = 4)
plotHeatmap(pmtc1_merged_knn_pl_nod_noc0, label = 'subclones', consensus = TRUE, order_cells = 'consensus_tree', raster_quality = 4,
            genes = c("MYC", "MYB", "BRCA1", "ERBB2", "CDH1", "FGFR1", "AKT2", "CDK4", "CCNE1", "CCND1", "MTOR", "FGF10", "BRAF", "AURKA"))
dev.off()

pmtc1_merged_knn_pl_nod_noc0 <- calcConsensus(pmtc1_merged_knn_pl_nod_noc0,
                                              assay = 'integer')

pdf(glue("{outdir}/pmtc1_merged_consensus_ht_integer.pdf"), width = 8, height = 4)
plotHeatmap(pmtc1_merged_knn_pl_nod_noc0, label = 'subclones', consensus = TRUE, order_cells = 'consensus_tree', raster_quality = 4, assay = 'integer',
            genes = c("MYC", "MYB", "BRCA1", "ERBB2", "CDH1", "FGFR1", "AKT2", "CDK4", "CCNE1", "CCND1", "MTOR", "FGF10", "BRAF", "AURKA"))
dev.off()



pmtc1_metrics_aftfilt <- as.data.frame(colData(pmtc1_merged_knn_pl_nod_noc0))

pmtc1_median_sub_ploidy <- pmtc1_metrics_aftfilt %>%
  group_by(subclones) %>%
  dplyr::summarise(median_pl_subclones = median(ploidy))

colData(pmtc1_merged_knn_pl_nod_noc0)$ploidy <-
  pmtc1_median_sub_ploidy$median_pl_subclones[match(colData(pmtc1_merged_knn_pl_nod_noc0)$subclones, pmtc1_median_sub_ploidy$subclones)]

pmtc1_merged_knn_pl_nod_noc0 <- calcInteger(pmtc1_merged_knn_pl_nod_noc0, method = 'metadata')

pmtc1_merged_knn_pl_nod_noc0_umap_p_ploidysubmedian <-plotUmap(pmtc1_merged_knn_pl_nod_noc0, label = 'ploidy')

pdf(glue("{outdir}/pmtc1_merged_subclones_integer_integer_submedian.pdf"), width = 7, height = 8)
plotHeatmap(pmtc1_merged_knn_pl_nod_noc0, label = c('subclones', 'timepoint'), order_cells = 'consensus_tree', raster_quality = 4, assay = 'integer')
dev.off()


pmtc1_merged_hvg_gc <- plotGeneCopy(pmtc1_merged_knn_pl_nod_noc0, genes = c("MYC", "MYB", "BRCA1", "ERBB2", "CDH1", "FGFR1", "AKT2", "CDK4", "CCNE1", "CCND1", "MTOR", "FGF10", "BRAF", "AURKA"),
                                    label = 'timepoint', dodge.width = .8) +
  scale_fill_hue() + coord_flip()






# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CO8 primary -----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(varbin){
 co8_tumor_primary <- runVarbin("/mnt/lab/users/dminussi/projects/copy_number_pipeline_hg38/test_tumor_data/CO8/primary/marked", remove_Y = TRUE)
}

co8_tumor_primary <- findOutliers(co8_tumor_primary)

co8_tumor_primary <- findAneuploidCells(co8_tumor_primary)

pdf(glue("{outdir}/co8_primary_filtered_is_normal.pdf"), width = 7, height = 8)
plotHeatmap(co8_tumor_primary, label = c("outlier", 'is_aneuploid'), order_cells = 'hclust',
            row_split = 'outlier', n_threads = 30, raster_quality = 4)
dev.off()

co8_tumor_primary_filt <- co8_tumor_primary[,colData(co8_tumor_primary)$outlier == FALSE]
co8_tumor_primary_filt <- co8_tumor_primary_filt[,colData(co8_tumor_primary_filt)$is_aneuploid == TRUE]

colData(co8_tumor_primary_filt)$timepoint <- 'primary'

co8_tumor_primary_knn <- knnSmooth(co8_tumor_primary_filt)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CO8 met -----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(varbin){
 co8_tumor_met <- runVarbin("/mnt/lab/users/dminussi/projects/copy_number_pipeline_hg38/test_tumor_data/CO8/met/marked", remove_Y = TRUE)
}

co8_tumor_met <- findOutliers(co8_tumor_met)

co8_tumor_met <- findAneuploidCells(co8_tumor_met)

pdf(glue("{outdir}/co8_met_filtered_is_normal.pdf"), width = 7, height = 8)
plotHeatmap(co8_tumor_met, label = c("outlier", 'is_aneuploid'),
            row_split = 'outlier', n_threads = 30, raster_quality = 4)
dev.off()

co8_tumor_met_filt <- co8_tumor_met[,colData(co8_tumor_met)$outlier == FALSE]
co8_tumor_met_filt <- co8_tumor_met_filt[,colData(co8_tumor_met_filt)$is_aneuploid == TRUE]

colData(co8_tumor_met_filt)$timepoint <- 'metastasis'

co8_tumor_met_knn <- knnSmooth(co8_tumor_met_filt)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# CO8 merged ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~


co8_merged_knn <- cbind(co8_tumor_primary_knn,
                        co8_tumor_met_knn)

co8_merged_knn <- calcInteger(co8_merged_knn, method = 'scquantum', assay = 'smoothed_bincounts')
plotMetrics(co8_merged_knn, 'ploidy', label = 'timepoint')
plotMetrics(co8_merged_knn, 'ploidy', label = 'ploidy_score')
table(colData(co8_merged_knn)$ploidy_score < 0.2)
co8_merged_knn_pl <- co8_merged_knn[,colData(co8_merged_knn)$ploidy_score < 0.2]

co8_metrics_aftfilt <- as.data.frame(colData(co8_merged_knn_pl))


co8_merged_knn_pl <- runUmap(co8_merged_knn_pl)

co8_merged_knn_pl <- findSuggestedK(co8_merged_knn_pl)
co8_merged_knn_pl_suggestedk <- plotSuggestedK(co8_merged_knn_pl)

co8_merged_knn_pl <- findClusters(co8_merged_knn_pl)

co8_merged_knn_pl_noc0 <- co8_merged_knn_pl[,colData(co8_merged_knn_pl)$subclones != 'c0']

co8_merged_knn_pl_noc0_tp_umap_p <- plotUmap(co8_merged_knn_pl_noc0, label = 'timepoint')

co8_merged_knn_pl_noc0 <- calcInteger(co8_merged_knn_pl_noc0, method = 'scquantum', assay = 'smoothed_bincounts')
co8_merged_knn_pl_p <- plotUmap(co8_merged_knn_pl_noc0, label = 'ploidy') + scale_fill_viridis_c(breaks = scales::pretty_breaks(n = 6))
# * scquantum ----
co8_merged_knn_pl_p <- plotUmap(co8_merged_knn_pl_noc0, label = 'ploidy') + scale_fill_viridis_c(option = 'turbo', breaks = scales::pretty_breaks(n = 6))

co8_merged_knn_pl_noc0 <- calcConsensus(co8_merged_knn_pl_noc0)
co8_merged_knn_pl_noc0 <- runConsensusPhylo(co8_merged_knn_pl_noc0)


co8_lb_colors = list(subclones = subclones_pal(),
                     ploidy = circlize::colorRamp2(
                       seq(min(colData(co8_merged_knn_pl_noc0)$ploidy),
                           max(colData(co8_merged_knn_pl_noc0)$ploidy),
                           length = 300),
                       viridis::viridis(300, option = 'turbo')
                     ))

pdf(glue("{outdir}/co8_merged_knn_pl_noc0_subclones_integer_ploidyanno_raster.pdf"), width = 7, height = 8)
plotHeatmap(co8_merged_knn_pl_noc0, label = c('subclones', 'ploidy'), order_cells = 'consensus_tree', use_raster = TRUE, raster_quality = 10, assay = 'integer', label_colors = co8_lb_colors)
dev.off()

co8_merged_knn_pl_noc0_umap_p <- plotUmap(co8_merged_knn_pl_noc0, label = 'subclones')

co8_merged_knn_pl_noc0 <- calcConsensus(co8_merged_knn_pl_noc0)
co8_merged_knn_pl_noc0 <- runConsensusPhylo(co8_merged_knn_pl_noc0)

colData(co8_merged_knn_pl_noc0)$timepoint <-
  forcats::fct_relevel(colData(co8_merged_knn_pl_noc0)$timepoint,
                       c("primary", "metastasis"))

co8_selected_genes = c(
  "SMAD3",
  "FHIT",
  "APC",
  "SOX4",
  "IGFBP7",
  "CDK8",
  "PIK3CA",
  "MYC",
  "TP53",
  "GATA4",
  "CHEK1",
  "TGFB1",
  "TIAM1"
)

co8_merged_knn_pl_noc0_selected_hvg_gc <- plotGeneCopy(co8_merged_knn_pl_noc0, genes = co8_selected_genes, label = 'timepoint', dodge.width = .8) + scale_fill_hue(direction = -1) + coord_flip()

table(colData(co8_merged_knn_pl_noc0)$subclones, colData(co8_merged_knn_pl_noc0)$timepoint)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CO8 tree ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

co8_tumor_primary_knn <- runUmap(co8_tumor_primary_knn)
co8_tumor_primary_knn <- findSuggestedK(co8_tumor_primary_knn)
co8_tumor_primary_knn_p <- plotSuggestedK(co8_tumor_primary_knn)

co8_tumor_primary_knn <- findClusters(co8_tumor_primary_knn)

co8_tumor_primary_knn_noc0 <- co8_tumor_primary_knn[,colData(co8_tumor_primary_knn)$subclones != 'c0']


co8_tumor_primary_umap_p <- plotUmap(co8_tumor_primary_knn_noc0, label = 'subclones')

co8_tumor_primary_knn_noc0 <- calcConsensus(co8_tumor_primary_knn_noc0)

co8_tumor_primary_knn_noc0 <- inferMrca(co8_tumor_primary_knn_noc0)

co8_merged_knn_pl_noc0 <- runConsensusPhylo(co8_merged_knn_pl_noc0, root = 'user',
                                            root_user = metadata(co8_tumor_primary_knn_noc0)$inferred_mrca)

colData(co8_merged_knn_pl_noc0)$timepoint <- forcats::fct_relevel(colData(co8_merged_knn_pl_noc0)$timepoint, c("metastasis", "primary"))

try(co8_merged_knn_pl_noc0_consensus_phylo <- plotPhylo(co8_merged_knn_pl_noc0, label = 'subclones', consensus = TRUE, group = 'timepoint'))

pdf(glue("{outdir}/co8_merged_knn_pl_noc0_subclones.pdf"), width = 7, height = 8)
plotHeatmap(co8_merged_knn_pl_noc0, label = c('subclones', 'timepoint'), order_cells = 'consensus_tree', raster_quality = 5)
dev.off()

pdf(glue("{outdir}/co8_merged_knn_pl_noc0_subclones_integer.pdf"), width = 7, height = 8)
plotHeatmap(co8_merged_knn_pl_noc0, label = c('subclones', 'timepoint'), order_cells = 'consensus_tree', raster_quality = 5, assay = 'integer')
dev.off()

co8_merged_knn_pl_noc0_metrics_aftfilt <- as.data.frame(colData(co8_merged_knn_pl_noc0))

co8_median_sub_ploidy <- co8_merged_knn_pl_noc0_metrics_aftfilt %>%
  group_by(subclones) %>%
  dplyr::summarise(median_pl_subclones = median(ploidy))

colData(co8_merged_knn_pl_noc0)$ploidy <-
  co8_median_sub_ploidy$median_pl_subclones[match(colData(co8_merged_knn_pl_noc0)$subclones, co8_median_sub_ploidy$subclones)]

co8_merged_knn_pl_noc0 <- calcInteger(co8_merged_knn_pl_noc0, method = 'metadata')

co8_merged_knn_pl_noc0_umap_p_ploidysubmedian <-plotUmap(co8_merged_knn_pl_noc0, label = 'ploidy')


pdf(glue("{outdir}/co8_merged_subclones_integer_integer_submedian.pdf"), width = 7, height = 8)
plotHeatmap(co8_merged_knn_pl_noc0, label = c('subclones', 'timepoint'), order_cells = 'consensus_tree', raster_quality = 4, assay = 'integer')
dev.off()


pdf(glue("{outdir}/co8_merged_knn_pl_noc0_consensus_ht.pdf"), width = 10, height = 5)
plotHeatmap(co8_merged_knn_pl_noc0, label = c('subclones'), consensus = TRUE, genes = co8_selected_genes)
dev.off()

co8_merged_knn_pl_noc0 <- calcConsensus(co8_merged_knn_pl_noc0, assay = 'integer')

pdf(glue("{outdir}/co8_merged_knn_pl_noc0_consensus_ht_integer.pdf"), width = 10, height = 5)
plotHeatmap(co8_merged_knn_pl_noc0, label = c('subclones'), consensus = TRUE, genes = co8_selected_genes, assay = 'integer', order_cells = 'consensus_tree')
dev.off()



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CO5 primary -----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(varbin){
 co5_tumor_primary <- runVarbin("/mnt/lab/users/dminussi/projects/copy_number_pipeline_hg38/test_tumor_data/CO5/primary/marked/", remove_Y = TRUE)
}
 
co5_tumor_primary <- findOutliers(co5_tumor_primary)

co5_tumor_primary <- findAneuploidCells(co5_tumor_primary)

pdf(glue("{outdir}/co5_primary_filtered_is_normal.pdf"), width = 7, height = 8)
plotHeatmap(co5_tumor_primary, label = c("outlier", 'is_aneuploid'), order_cells = 'hclust',
            row_split = 'outlier', n_threads = 30, raster_quality = 4)
dev.off()

co5_tumor_primary_filt <- co5_tumor_primary[,colData(co5_tumor_primary)$outlier == FALSE]
co5_tumor_primary_filt <- co5_tumor_primary_filt[,colData(co5_tumor_primary_filt)$is_aneuploid == TRUE]

colData(co5_tumor_primary_filt)$timepoint <- 'primary'

co5_tumor_primary_knn <- knnSmooth(co5_tumor_primary_filt)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# co5 met -----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(varbin){
 co5_tumor_met <- runVarbin("/mnt/lab/users/dminussi/projects/copy_number_pipeline_hg38/test_tumor_data/CO5/met/marked", remove_Y = TRUE)
}
 
co5_tumor_met <- findOutliers(co5_tumor_met)

co5_tumor_met <- findAneuploidCells(co5_tumor_met)

pdf(glue("{outdir}/co5_met_filtered_is_normal.pdf"), width = 7, height = 8)
plotHeatmap(co5_tumor_met, label = c("outlier", 'is_aneuploid'),
            row_split = 'outlier', n_threads = 30, raster_quality = 4)
dev.off()

co5_tumor_met_filt <- co5_tumor_met[,colData(co5_tumor_met)$outlier == FALSE]
co5_tumor_met_filt <- co5_tumor_met_filt[,colData(co5_tumor_met_filt)$is_aneuploid == TRUE]

colData(co5_tumor_met_filt)$timepoint <- 'metastasis'

co5_tumor_met_knn <- knnSmooth(co5_tumor_met_filt)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# co5 merged ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~


co5_merged_knn <- cbind(co5_tumor_primary_knn,
                        co5_tumor_met_knn)

co5_merged_knn <- calcInteger(co5_merged_knn, method = 'scquantum', assay = 'smoothed_bincounts')
plotMetrics(co5_merged_knn, 'ploidy', label = 'timepoint')
plotMetrics(co5_merged_knn, 'ploidy', label = 'ploidy_score')
table(colData(co5_merged_knn)$ploidy_score < 0.2)
co5_merged_knn_pl <- co5_merged_knn[,colData(co5_merged_knn)$ploidy_score < 0.2]

co5_merged_knn_pl <- runUmap(co5_merged_knn_pl)

co5_merged_knn_pl <- findSuggestedK(co5_merged_knn_pl)
co5_merged_knn_pl_suggestedk <- plotSuggestedK(co5_merged_knn_pl)

co5_merged_knn_pl <- findClusters(co5_merged_knn_pl)

co5_merged_knn_pl_noc0 <- co5_merged_knn_pl[,colData(co5_merged_knn_pl)$subclones != 'c0']

co5_merged_knn_pl_noc0_tp_umap_p <- plotUmap(co5_merged_knn_pl_noc0, label = 'timepoint')

co5_merged_knn_pl_noc0_umap_p <- plotUmap(co5_merged_knn_pl_noc0, label = 'subclones')

co5_merged_knn_pl_noc0 <- calcConsensus(co5_merged_knn_pl_noc0)
co5_merged_knn_pl_noc0 <- runConsensusPhylo(co5_merged_knn_pl_noc0)


colData(co5_merged_knn_pl_noc0)$timepoint <-
  forcats::fct_relevel(colData(co5_merged_knn_pl_noc0)$timepoint,
                       c("primary", "metastasis"))

co5_selected_genes = c(
  "SMAD3",
  "FHIT",
  "APC",
  "SOX4",
  "IGFBP7",
  "CDK8",
  "PIK3CA",
  "MYC",
  "TP53",
  "GATA4",
  "CHEK1",
  "TGFB1",
  "TIAM1"
)

co5_merged_knn_pl_noc0_selected_hvg_gc <- plotGeneCopy(co5_merged_knn_pl_noc0, genes = co5_selected_genes, label = 'timepoint', dodge.width = .8) + scale_fill_hue(direction = -1) + coord_flip()

table(colData(co5_merged_knn_pl_noc0)$subclones, colData(co5_merged_knn_pl_noc0)$timepoint)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# co5 tree ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

co5_tumor_primary_knn <- runUmap(co5_tumor_primary_knn)
co5_tumor_primary_knn <- findSuggestedK(co5_tumor_primary_knn)
co5_tumor_primary_knn_p <- plotSuggestedK(co5_tumor_primary_knn)

co5_tumor_primary_knn <- findClusters(co5_tumor_primary_knn)

co5_tumor_primary_knn_noc0 <- co5_tumor_primary_knn[,colData(co5_tumor_primary_knn)$subclones != 'c0']


co5_tumor_primary_umap_p <- plotUmap(co5_tumor_primary_knn_noc0, label = 'subclones')

co5_tumor_primary_knn_noc0 <- calcConsensus(co5_tumor_primary_knn_noc0)

co5_tumor_primary_knn_noc0 <- inferMrca(co5_tumor_primary_knn_noc0)

co5_merged_knn_pl_noc0 <- runConsensusPhylo(co5_merged_knn_pl_noc0, root = 'user',
                                            root_user = metadata(co5_tumor_primary_knn_noc0)$inferred_mrca)
# consensusPhylo(co5_merged_knn_pl_noc0) <- phytools::rotateNodes(consensusPhylo(co5_merged_knn_pl_noc0), c(19,21:27, 29:33))
colData(co5_merged_knn_pl_noc0)$timepoint <- forcats::fct_relevel(colData(co5_merged_knn_pl_noc0)$timepoint, c("metastasis", "primary"))

# co5_merged_knn_pl_noc0_flip <- co5_merged_knn_pl_noc0
# consensusPhylo(co5_merged_knn_pl_noc0_flip) <- phytools::rotateNodes(consensusPhylo(co5_merged_knn_pl_noc0_flip), c(19,21:27, 29:33))
# co5_merged_knn_pl_noc0_consensus_phylo <- plotPhylo(co5_merged_knn_pl_noc0_flip, label = 'subclones', consensus = TRUE, group = 'timepoint')
# cowplot::ggsave2(glue("{outdir}/co5_merged_knn_pl_noc0_consensus_phylo_test.pdf"), co5_merged_knn_pl_noc0_consensus_phylo, width = 5, height = 5)

# try(co5_merged_knn_pl_noc0_consensus_phylo <- plotPhylo(co5_merged_knn_pl_noc0, label = 'subclones', consensus = TRUE, group = 'timepoint'))
# 
# pdf(glue("{outdir}/co5_merged_knn_pl_noc0_subclones.pdf"), width = 7, height = 8)
# plotHeatmap(co5_merged_knn_pl_noc0, label = c('subclones', 'timepoint'), order_cells = 'consensus_tree', raster_quality = 5)
# dev.off()
# 
# pdf(glue("{outdir}/co5_merged_knn_pl_noc0_subclones_integer.pdf"), width = 7, height = 8)
# plotHeatmap(co5_merged_knn_pl_noc0, label = c('subclones', 'timepoint'), order_cells = 'consensus_tree', raster_quality = 5, assay = 'integer')
# dev.off()

co5_merged_knn_pl_noc0_metrics_aftfilt <- as.data.frame(colData(co5_merged_knn_pl_noc0))


co5_median_sub_ploidy <- co5_merged_knn_pl_noc0_metrics_aftfilt %>%
  group_by(subclones) %>%
  dplyr::summarise(median_pl_subclones = median(ploidy))

colData(co5_merged_knn_pl_noc0)$ploidy <-
  co5_median_sub_ploidy$median_pl_subclones[match(colData(co5_merged_knn_pl_noc0)$subclones, co5_median_sub_ploidy$subclones)]

co5_merged_knn_pl_noc0 <- calcInteger(co5_merged_knn_pl_noc0, method = 'metadata')
co5_merged_knn_pl_noc0_umap_p_ploidysubmedian <-plotUmap(co5_merged_knn_pl_noc0, label = 'ploidy')


# pdf(glue("{outdir}/co5_merged_subclones_integer_integer_submedian.pdf"), width = 7, height = 8)
# plotHeatmap(co5_merged_knn_pl_noc0, label = c('subclones', 'timepoint'), order_cells = 'consensus_tree', raster_quality = 4, assay = 'integer')
# dev.off()
# 
# 
# pdf(glue("{outdir}/co5_merged_knn_pl_noc0_consensus_ht.pdf"), width = 10, height = 5)
# plotHeatmap(co5_merged_knn_pl_noc0, label = c('subclones'), consensus = TRUE, genes = co5_selected_genes)
# dev.off()

co5_merged_knn_pl_noc0 <- calcConsensus(co5_merged_knn_pl_noc0, assay = 'integer')

# pdf(glue("{outdir}/co5_merged_knn_pl_noc0_consensus_ht_integer.pdf"), width = 10, height = 5)
# plotHeatmap(co5_merged_knn_pl_noc0, label = c('subclones'), consensus = TRUE, genes = co5_selected_genes, assay = 'integer', order_cells = 'consensus_tree')
# dev.off()



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# brcan -----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(varbin){
 brcan <- runVarbin("/mnt/lab/users/dminussi/projects/copy_number_pipeline_hg38/test_tumor_data/BRCAN/marked/", remove_Y = T)
}


brcan_filt <- findAneuploidCells(brcan)

pdf(glue("{outdir}/brcan_knn_is_aneuploid.pdf"), width = 7, height = 8)
plotHeatmap(brcan_filt, label = c('is_aneuploid'),
            row_split = 'is_aneuploid',
            n_threads = 50)
dev.off()
brcan_filt <- brcan_filt[,colData(brcan_filt)$is_aneuploid == TRUE]


brcan_filt <- findOutliers(brcan_filt)

pdf(glue("{outdir}/brcan_knn_outlier.pdf"), width = 7, height = 8)
plotHeatmap(brcan_filt, label = c('outlier'),
            row_split = 'outlier',
            n_threads = 50)
dev.off()

brcan_filt <- runMetrics(brcan_filt)
brcan_metrics <- as.data.frame(colData(brcan_filt))
brcan_p_metrics <- plotMetrics(
  brcan_filt,
  metric = c(
    "reads_total",
    "percentage_duplicates",
    "overdispersion",
    "breakpoint_count"
  ), label = 'outlier',
  dodge.width = 0.8,
  ncol = 2
) + scale_fill_manual(values = c("TRUE" = "#DA614D",
                                 "FALSE" = "#5F917A"))


brcan_filt <- brcan_filt[,colData(brcan_filt)$outlier == FALSE]

brcan_knn <- knnSmooth(brcan_filt)

brcan_knn <- runMetrics(brcan_knn)
brcan_knn <- calcSmoothOver(brcan_knn)


brcan_metrics_aftfilt <- as.data.frame(colData(brcan_knn))


brcan_knn <- runUmap(brcan_knn)
brcan_knn <- findSuggestedK(brcan_knn)
brcan_primary_suggestedK_p <- plotSuggestedK(brcan_knn)


brcan_knn <- findClusters(brcan_knn)
brcan_knn <- brcan_knn[,colData(brcan_knn)$subclones != 'c0']


# * Umap ----
brcan_knn_umap_p <- plotUmap(brcan_knn, label = 'subclones')

brcan_knn <- calcConsensus(brcan_knn)
brcan_knn <- runConsensusPhylo(brcan_knn, root = 'mrca')

pdf(glue("{outdir}/brcan_primary_subclones.pdf"), width = 7, height = 8)
plotHeatmap(brcan_knn, label = c('subclones'), order_cells = 'consensus_tree')
dev.off()


pdf(glue("{outdir}/brcan_knn_cs_ht.pdf"), width = 7,height = 4)
plotHeatmap(brcan_knn, label = 'subclones', consensus = T,
            genes = c(
              "FHIT",
              "CUX1",
              "WWOX",
              "FOXO1",
              "APC",
              "BCL2",
              "KRAS",
              "PGR",
              "PDGFRA",
              "ROS1",
              "PTPN11",
              "CCND1",
              "BTG1",
              "FGFR3",
              "PTEN",
              "FGFR2"
            ))
dev.off()

# *scquantum ----
brcan_knn <- calcInteger(brcan_knn, method = 'scquantum', assay = 'smoothed_bincounts')

table(colData(brcan_knn)$ploidy_score < 0.2)
brcan_knn_pl <- brcan_knn[,colData(brcan_knn)$ploidy_score < 0.2]

plotMetrics(brcan_knn_pl, metric = 'ploidy', label = 'ploidy_score')

pdf(glue("{outdir}/brcan_subclones_integer.pdf"), width = 7, height = 8)
plotHeatmap(brcan_knn_pl, label = c('subclones'), order_cells = 'consensus_tree', raster_quality = 5, assay = 'integer')
dev.off()
brcan_lb_colors = list(subclones = subclones_pal(),
                       ploidy = circlize::colorRamp2(
                         seq(min(colData(brcan_knn_pl)$ploidy),
                             max(colData(brcan_knn_pl)$ploidy),
                             length = 300),
                         viridis::viridis(300, option = 'turbo')
                       ))
pdf(glue("{outdir}/brcan_primary_subclones_ploidyanno.pdf"), width = 7, height = 8)
plotHeatmap(brcan_knn_pl, label = c('subclones', 'ploidy'), order_cells = 'consensus_tree', use_raster = TRUE, raster_quality = 10, assay = 'integer', label_colors = brcan_lb_colors)
dev.off()


# brcan_knn_umap_ploidy_p <- plotUmap(brcan_knn_pl, label = 'ploidy') + scale_fill_jcolors_contin('rainbow')
brcan_knn_umap_ploidy_p <- plotUmap(brcan_knn_pl, label = 'ploidy') + scale_fill_viridis_c(breaks = scales::pretty_breaks(n = 6))
brcan_knn_umap_ploidy_p <- plotUmap(brcan_knn_pl, label = 'ploidy') + scale_fill_viridis_c(breaks = scales::pretty_breaks(n = 6), option = 'turbo')

brcan_metrics_aftfilt <- as.data.frame(colData(brcan_knn_pl))

brcan_median_sub_ploidy <- brcan_metrics_aftfilt %>%
  group_by(subclones) %>%
  dplyr::summarise(median_pl_subclones = median(ploidy))

colData(brcan_knn_pl)$ploidy <-
  brcan_median_sub_ploidy$median_pl_subclones[match(colData(brcan_knn_pl)$subclones, brcan_median_sub_ploidy$subclones)]

brcan_knn_pl <- calcInteger(brcan_knn_pl, method = 'metadata')

pdf(glue("{outdir}/brcan_subclones_integer_submedian.pdf"), width = 7, height = 8)
plotHeatmap(brcan_knn_pl, label = c('subclones'), order_cells = 'consensus_tree', raster_quality = 4, assay = 'integer')
dev.off()



# 10X TN7 -----

if(varbin){
 tn17_10x <- runVarbin("/mnt/lab/users/dminussi/projects/copy_number_pipeline_hg38/test_tumor_data/TN17_10XCNA/marked/", remove_Y = T, is_paired_end = TRUE)
}

tn17_10x_filt <- findAneuploidCells(tn17_10x)

pdf(glue("{outdir}/tn17_10x_knn_is_aneuploid.pdf"), width = 7, height = 8)
plotHeatmap(tn17_10x_filt, label = c('is_aneuploid'),
            row_split = 'is_aneuploid',
            n_threads = 50)
dev.off()
tn17_10x_filt <- tn17_10x_filt[,colData(tn17_10x_filt)$is_aneuploid == TRUE]


tn17_10x_filt <- findOutliers(tn17_10x_filt)

pdf(glue("{outdir}/tn17_10x_knn_outlier.pdf"), width = 7, height = 8)
plotHeatmap(tn17_10x_filt, label = c('outlier'),
            row_split = 'outlier',
            n_threads = 50)
dev.off()


tn17_10x_filt <- runMetrics(tn17_10x_filt)
tn17_10x_filt <- calcSmoothOver(tn17_10x_filt)

tn17_10x_p_metrics <- plotMetrics(
  tn17_10x_filt,
  metric = c(
    "reads_total",
    "percentage_duplicates",
    "overdispersion",
    "breakpoint_count"
  ), label = 'outlier',
  dodge.width = 0.8,
  ncol = 2
) + scale_fill_manual(values = c("TRUE" = "#DA614D",
                                 "FALSE" = "#5F917A"))


tn17_10x_metrics <- as.data.frame(colData(tn17_10x_filt))





tn17_10x_filt <- tn17_10x_filt[,colData(tn17_10x_filt)$outlier == FALSE]

tn17_10x_knn <- knnSmooth(tn17_10x_filt)

colData(tn17_10x_knn)$overdispersion <- colData(tn17_10x_filt)$overdispersion
tn17_10x_knn <- calcSmoothOver(tn17_10x_knn)






tn17_10x_metrics_aftfilt <- as.data.frame(colData(tn17_10x_knn))


tn17_10x_knn <- runUmap(tn17_10x_knn, n_neighbors = 40)
tn17_10x_knn <- findSuggestedK(tn17_10x_knn, k_range = 15:30)
tn17_10x_primary_suggestedK_p <- plotSuggestedK(tn17_10x_knn)


tn17_10x_knn <- findClusters(tn17_10x_knn)
tn17_10x_knn <- tn17_10x_knn[,colData(tn17_10x_knn)$subclones != 'c0']


# * Umap ----
tn17_10x_knn_umap_p <- plotUmap(tn17_10x_knn, label = 'subclones')

tn17_10x_knn <- calcConsensus(tn17_10x_knn)
tn17_10x_knn <- runConsensusPhylo(tn17_10x_knn, root = 'mrca')

pdf(glue("{outdir}/tn17_10x_primary_subclones.pdf"), width = 7, height = 8)
plotHeatmap(tn17_10x_knn, label = c('subclones'), order_cells = 'consensus_tree')
dev.off()


pdf(glue("{outdir}/tn17_10x_knn_cs_ht.pdf"), width = 7,height = 4)
plotHeatmap(tn17_10x_knn, label = 'subclones', consensus = T,
            genes = c(
              "FHIT",
              "CUX1",
              "WWOX",
              "FOXO1",
              "APC",
              "BCL2",
              "KRAS",
              "PGR",
              "PDGFRA",
              "ROS1",
              "PTPN11",
              "CCND1",
              "BTG1",
              "FGFR3",
              "PTEN",
              "FGFR2"
            ))
dev.off()

# *scquantum ----
tn17_10x_knn <- calcInteger(tn17_10x_knn, method = 'scquantum', assay = 'smoothed_bincounts')

table(colData(tn17_10x_knn)$ploidy_score < 0.2)
tn17_10x_knn_pl <- tn17_10x_knn[,colData(tn17_10x_knn)$ploidy_score < 0.2]

plotMetrics(tn17_10x_knn_pl, metric = 'ploidy', label = 'ploidy_score')

pdf(glue("{outdir}/tn17_10x_subclones_integer.pdf"), width = 7, height = 8)
plotHeatmap(tn17_10x_knn_pl, label = c('subclones'), order_cells = 'consensus_tree', raster_quality = 5, assay = 'integer')
dev.off()

pdf(glue("{outdir}/tn17_10x_primary_subclones_ploidyanno.pdf"), width = 7, height = 8)
plotHeatmap(tn17_10x_knn_pl, label = c('subclones', 'ploidy'), order_cells = 'consensus_tree', use_raster = TRUE, raster_quality = 10, assay = 'integer')
dev.off()

tn17_10x_knn_pl_reserv <- tn17_10x_knn_pl

tn17_10x_metrics_aftfilt <- as.data.frame(colData(tn17_10x_knn_pl))

tn17_10x_median_sub_ploidy <- tn17_10x_metrics_aftfilt %>%
  group_by(subclones) %>%
  dplyr::summarise(median_pl_subclones = median(ploidy))

colData(tn17_10x_knn_pl)$ploidy <-
  tn17_10x_median_sub_ploidy$median_pl_subclones[match(colData(tn17_10x_knn_pl)$subclones, tn17_10x_median_sub_ploidy$subclones)]

tn17_10x_knn_pl <- calcInteger(tn17_10x_knn_pl, method = 'metadata')
tn17_10x_knn_pl_umap_p_ploidysubmedian <-plotUmap(tn17_10x_knn_pl, label = 'ploidy')



pdf(glue("{outdir}/tn17_10x_subclones_integer_submedian.pdf"), width = 7, height = 8)
plotHeatmap(tn17_10x_knn_pl, label = c('subclones'), order_cells = 'consensus_tree', raster_quality = 4, assay = 'integer')
dev.off()

# ~~~~~~~~~~~
# DLP -----
# ~~~~~~~~~~~

if(varbin){
 dlp_sa928_a90553c <- runVarbin("/mnt/lab/users/dminussi/projects/copykit_manuscript/dlp_plus/marked/", remove_Y = T, is_paired_end = TRUE)
}
 
dlp_sa928_a90553c <- findAneuploidCells(dlp_sa928_a90553c)

# pdf(glue("{outdir}/dlp_sa928_a90553c_is_aneuploid.pdf"), width = 7, height = 8)
# plotHeatmap(dlp_sa928_a90553c, label = c('is_aneuploid'),
#             row_split = 'is_aneuploid',
#             n_threads = 50)
# dev.off()
# dlp_sa928_a90553c <- dlp_sa928_a90553c[,colData(dlp_sa928_a90553c)$is_aneuploid == TRUE]

dlp_sa928_a90553c <- findOutliers(dlp_sa928_a90553c, resolution = 0.8)

pdf(glue("{outdir}/dlp_sa928_a90553c_filtered.pdf"), width = 7, height = 8)
plotHeatmap(dlp_sa928_a90553c, label = c('outlier'),
            row_split = 'outlier',
            col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
            n_threads = 50)
dev.off()

dlp_sa928_a90553c <- runMetrics(dlp_sa928_a90553c)

dlp_sa928_a90553_p_metrics <- plotMetrics(
  dlp_sa928_a90553c,
  metric = c(
    "reads_total",
    "percentage_duplicates",
    "overdispersion",
    "breakpoint_count"
  ), label = 'outlier',
  dodge.width = 0.8,
  ncol = 2
) + scale_fill_manual(values = c("TRUE" = "#DA614D",
                                 "FALSE" = "#5F917A"))

dlp_sa928_a90553_metrics <- as.data.frame(colData(dlp_sa928_a90553c))


dlp_sa928_a90553c <- dlp_sa928_a90553c[,colData(dlp_sa928_a90553c)$outlier == FALSE]

dlp_sa928_a90553_metrics_aftfilt <- as.data.frame(colData(dlp_sa928_a90553c))


# dlp_sa928_a90553c <- runUmap(dlp_sa928_a90553c)
# dlp_sa928_a90553c <- findSuggestedK(dlp_sa928_a90553c)
# dlp_sa928_a90553_primary_suggestedK_p <- plotSuggestedK(dlp_sa928_a90553c)

# cowplot::ggsave2(glue("{outdir}/dlp_sa928_a90553_primary_suggestedK_p.pdf"), dlp_sa928_a90553_primary_suggestedK_p, width = 5, height = 3)

# dlp_sa928_a90553c <- findClusters(dlp_sa928_a90553c)

## for revision
co8_merged_knn_pl_p <- plotUmap(co8_merged_knn_pl_noc0, label = 'ploidy') + scale_fill_viridis_c(limits = c(2.2,4.2),option = "turbo" )

co8_lb_colors = list(subclones = subclones_pal(),
                     ploidy = circlize::colorRamp2(
                       seq(2.2,
                           4.2,
                           length = 300),
                       viridis::viridis(300,option = "turbo")
                     ))

pdf(glue("{outdir}/co8_merged_knn_pl_noc0_subclones_integer_ploidyanno_raster_rev.pdf"), width = 7, height = 8)
plotHeatmap(co8_merged_knn_pl_noc0, label = c('subclones', 'ploidy'), order_cells = 'consensus_tree', use_raster = TRUE, raster_quality = 10, assay = 'integer', label_colors = co8_lb_colors, n_threads = 100)
dev.off()

brcan_lb_colors = list(subclones = subclones_pal(),
                       ploidy = circlize::colorRamp2(
                         seq(1.3,
                             3.3,
                             length = 300),
                         viridis::viridis(300, option = 'turbo')
                       ))
pdf(glue("{outdir}/brcan_primary_subclones_ploidyanno_rev.pdf"), width = 7, height = 8)
plotHeatmap(brcan_knn_pl, label = c('subclones', 'ploidy'), order_cells = 'consensus_tree', use_raster = TRUE, raster_quality = 10, assay = 'integer', label_colors = brcan_lb_colors)
dev.off()


# brcan_knn_umap_ploidy_p <- plotUmap(brcan_knn_pl, label = 'ploidy') + scale_fill_jcolors_contin('rainbow')
brcan_knn_umap_ploidy_p <- plotUmap(brcan_knn_pl, label = 'ploidy') + scale_fill_viridis_c(limits = c(1.3,3.3), option = 'turbo')

save(list = ls(all.names = T), file = glue("{outdir}/20240424_ck_singul.Rdata"))


sink()