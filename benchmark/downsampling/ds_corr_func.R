library(tidyverse)
library(copykit)
library(BiocParallel)
library(glue)
register(MulticoreParam(workers = 100, progressbar = F), default = T)
n_reads_levels <-  c("10k", "25k", "50k", "75k", "125k", "250k", "500k", "750k", "1M")
n_reads_levels <- rev(n_reads_levels)

## read in function for copykit
run_ds_ck <- function(path, is_pairend=F, res=c("downsampled_sort_50k", "downsampled_sort_75k", "downsampled_sort_125k", "downsampled_sort_250k", "downsampled_sort_500k", "downsampled_sort_750k", "downsampled_sort_1M", "fulldata")){
  
  a <- lapply(res, function(x){
    obj <- runVarbin(dir = glue("{path}/{x}"), genome = "hg19", is_paired_end = is_pairend, min_bincount = 0, remove_Y = T)
    obj <- calcInteger(obj, method = "scquantum")
    obj
  })
  names(a) <- str_remove(string = res, pattern = "downsampled_sort_")
  a
  
}

####### !!!!!! ADD IN smoothing module 20240304
## read in function for copykit
run_ds_ck_smooth <- function(path, is_pairend=F, res=c("downsampled_sort_50k", "downsampled_sort_75k", "downsampled_sort_125k", "downsampled_sort_250k", "downsampled_sort_500k", "downsampled_sort_750k", "downsampled_sort_1M", "fulldata")){
  
  a <- lapply(res, function(x){
    obj <- runVarbin(dir = glue("{path}/{x}"), genome = "hg19", is_paired_end = is_pairend, min_bincount = 0, remove_Y = T)
    obj <- knnSmooth(obj)
    obj <- calcInteger(obj, method = "scquantum", assay = 'smoothed_bincounts')
    obj
  })
  names(a) <- str_remove(string = res, pattern = "downsampled_sort_")
  a
  
}
## read in function for ginkgo
read_ds_gk <- function(path, res = c("downsampled_sort_50k", "downsampled_sort_75k", "downsampled_sort_125k", "downsampled_sort_250k", "downsampled_sort_500k", "downsampled_sort_750k", "downsampled_sort_1M", "fulldata")){
  a <- lapply(res, function(x){
    seg <- read_tsv(glue("{path}.{x}/SegNorm"))
    seg[,-c(1:3)]
  })
  names(a) <- str_remove(string = res, pattern = "downsampled_sort_")
  a
}

  
## read in function for sc_hmmcopy
read_ds_hmm <- function(path, res=c("downsampled_sort_50k", "downsampled_sort_75k", "downsampled_sort_125k", "downsampled_sort_250k", "downsampled_sort_500k", "downsampled_sort_750k", "downsampled_sort_1M", "fulldata") ){
  a <- lapply(res, function(x){
    seg_schmm <- read_csv(glue("{path}.{x}/results/reads.csv.gz"))
    seg <- seg_schmm %>% dplyr::select(start, end, chr, state, cell_id) %>%
      spread(key = cell_id, value = state)
    seg[,-c(1:3)]
  })
  names(a) <- str_remove(string = res, pattern = "downsampled_sort_")
  a
}

## calc correlations
corr_cells <- function(original,
                       downsampled) {
  
  ###
  # Loops across all cells matching by their names and calculates spearman
  # correlation
  
  
  # original: uber.sample.seg.txt matrix of segment ratios from original files
  # downsampled: uber.sample.seg.txt matrix of segment ratios from downsampled files
  ###
  
  # keeping only cells that were downsampled
  original <- original[,names(downsampled)]
  
  # sanity check
  stopifnot(identical(names(downsampled), names(original)))
  
  # cell names vector
  cell_names <- names(downsampled)
  
  # running correlation
  cor_list <- BiocParallel::bplapply(cell_names, function(x)
    cor(downsampled[,x], original[,x], method = 'spearman')
  )
  
  # binding to a vector
  cor_vector <- do.call(c, cor_list)
  names(cor_vector) <- cell_names
  cor_vector
  
}
create_corr_df <- function(matlist) {
  
  original = matlist$fulldata
  corlist <- lapply(matlist[-which(names(matlist)=="fulldata")], function(x){corr_cells(original = original, downsampled = x)})
  corlist_mod <- do.call(rbind, lapply(names(corlist), function(x){data.frame(cell=names(corlist[[x]]),
                                                                              correlation=corlist[[x]],
                                                                              n_reads=x)}))
  
  corlist_mod %>% mutate(n_reads = fct_relevel(n_reads, n_reads_levels)) 
  
}
