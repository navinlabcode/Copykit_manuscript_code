library(microbenchmark)
library(copykit)
library(BiocParallel)
register(MulticoreParam(progressbar = F, workers = 1), default = T)

run_copykit_benchmark <- function() {
  ck_benchmark_obj <- runVarbin("/mnt/lab/users/dminussi/projects/copy_number_pipeline_hg38/test_tumor_data/benchmark_cells_hg19", genome = 'hg19', remove_Y = T)
  ck_benchmark_obj <- calcInteger(ck_benchmark_obj, method = 'scquantum')
}

ck <- microbenchmark(
  run_copykit_benchmark(),
  times = 6
)

run_hmm_benchmark <- function() {
  system("singularity exec --bind /rsrch4/ --cleanenv -H /rsrch4/home/genetics/jwang48/projects/copykit/benchmark/ /rsrch4/home/genetics/jwang48/pipeline/single_cell_pipeline_hmmcopy_v0.8.26.sif single_cell hmmcopy --input_yaml input/bams.yaml  --tmpdir temp/  --pipelinedir pipeline/  --library_id benchmark --submit local --maxjobs 1 --nocleanup  --sentinel_only --loglevel DEBUG  --config_file input/config.yaml --output_prefix result")
  system("rm -rf result temp/* pipeline/*")

}
hmm <- microbenchmark(
  run_hmm_benchmark(),
  times = 6
)

run_scope_benchmark <- function() {
  system("Rscript scope.r")

}
sco <- microbenchmark(
  run_scope_benchmark(),
  times = 6
)

run_ginkgo_benchmark <- function() {
  system("/mnt/seq/code/PIPELINES/copy_number_pipeline_hg38/dev_lib/code/copy_number_pipeline_hg38/ginkgo/scripts/analyze.sh benchmark_cells")

}
gk <- microbenchmark(
  run_ginkgo_benchmark(),
  times = 6
)


bench <- data.frame(copykit = ck$time/1e9/3600,
                    gingko = gk$time/1e9/3600+2.166667,
                    hmmcopy = hmm$time/1e9/3600,
                    scope = sco$time/1e9/3600)
bench$gingko[1]<-bench$gingko[1]-13

ggplot(bench %>% gather(key = "tool", value = "Time (hrs)") , aes(x=tool, y=`Time (hrs)`))+
  geom_boxplot(fill = "cornflowerblue", 
               alpha = .7) +
  ylim(c(0,10)) +
  theme_bw()
ggsave("~/Projects/copykit/benchmark/sc_hmmcopy/benchmark/benchmark.pdf", width = 3, height = 3)
plotdata <- bench %>% gather(key = "tool", value = "Time (hrs)") %>%
  group_by(tool) %>%
  summarize(n = n(),
            mean = mean(`Time (hrs)`),
            sd = sd(`Time (hrs)`),
            se = sd / sqrt(n),
            ci = qt(0.95, df = n - 1) * sd / sqrt(n))
ggplot(plotdata, 
       aes(x = tool, 
           y = mean, 
           group = 1)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = mean - se, 
                    ymax = mean + se), 
                width = .1)
