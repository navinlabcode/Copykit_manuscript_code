# BiocManager::install("SCOPE")
# BiocManager::install("WGSmapp")


library(SCOPE)
library(WGSmapp)
library(BSgenome.Hsapiens.UCSC.hg19)
workdir <- "/mnt/USR1/junke/Projects/copykit/benchmark/SCOPE/"
bamfolder <- "/mnt/lab/users/dminussi/projects/copy_number_pipeline_hg38/test_tumor_data/benchmark_cells_hg19"
bamFile <- list.files(bamfolder, pattern = '*.sort.markdup.bam$')[200:300]
ncores=1

## read in bam files and QC
setwd(workdir)
bamdir <- file.path(bamfolder, bamFile)
sampname_raw <- sapply(strsplit(bamFile, ".", fixed = TRUE), "[", 1)
bambedObj <- get_bam_bed(bamdir = bamdir, sampname = sampname_raw, 
                         hgref = "hg19", resolution = 220)
ref_raw <- bambedObj$ref

mapp <- get_mapp(ref_raw, hgref = "hg19")

gc <- get_gc(ref_raw, hgref = "hg19")
values(ref_raw) <- cbind(values(ref_raw), DataFrame(gc, mapp))
ref_raw

coverageObj <- get_coverage_scDNA(bambedObj, mapqthres = 30, 
                                  seq = 'single-end', hgref = "hg19")

Y_raw <- coverageObj$Y

QCmetric_raw <- get_samp_QC(bambedObj)
qcObj <- perform_qc(Y_raw = Y_raw, 
                    sampname_raw = sampname_raw, ref_raw = ref_raw, 
                    QCmetric_raw = QCmetric_raw)
Y <- qcObj$Y
sampname <- qcObj$sampname
ref <- qcObj$ref
QCmetric <- qcObj$QCmetric

## Infers normal cells by Gini coef and use them as negative control to do normalization
# get gini coefficient for each cell
Gini <- get_gini(Y)
normObj <- normalize_codex2_ns_noK(Y_qc = Y,
                                       gc_qc = ref$gc,
                                       norm_index = which(Gini<=0.12))
ploidy <- initialize_ploidy(Y = Y, Yhat = normObj$Yhat, ref = ref)

# If using high performance clusters, parallel computing is 
# easy and improves computational efficiency. Simply use 
# normalize_scope_foreach() instead of normalize_scope(). 
# All parameters are identical. 
normObj.scope <- normalize_scope_foreach(Y_qc = Y, gc_qc = ref$gc,
                                         K = 1, ploidyInt = ploidy,
                                         norm_index = which(Gini<=0.12), T = 1:10,
                                         beta0 = normObj$beta.hat, nCores = ncores)

Yhat <- normObj.scope$Yhat[[which.max(normObj.scope$BIC)]]
fGC.hat <- normObj.scope$fGC.hat[[which.max(normObj.scope$BIC)]]

plot_EM_fit(Y_qc = Y, gc_qc = ref$gc, norm_index = which(Gini<=0.12), 
            T = 1:10,
            ploidyInt = ploidy, beta0 = normObj$beta.hat,
            filename = "plot_EM_fit_demo.pdf")


## cross sample segmentation
chrs <- unique(as.character(seqnames(ref)))
segment_cs <- vector('list',length = length(chrs))
names(segment_cs) <- chrs
# this part can be parallelized by chromosomes but it is not part of the feature from the package
for (chri in chrs) {
  message('\n', chri, '\n')
  segment_cs[[chri]] <- segment_CBScs(Y = Y,
                                      Yhat = Yhat,
                                      sampname = colnames(Y),
                                      ref = ref,
                                      chr = chri,
                                      mode = "integer", max.ns = 1)
}

iCN <- do.call(rbind, lapply(segment_cs, function(z){z[["iCN"]]}))


## plot Heatmap
plot_iCN(iCNmat = iCN, ref = ref, Gini = Gini, 
         filename = "plot_iCN_demo")


## save segment ranges, ploidy, copy number matrices and original object
# saveRDS(list(ref, ploidy, iCN, normObj.scope), "brcan_100cells_example.rds")

