library(tidyverse)
library(glue)

# SAMPLEID = "co8"
SAMPLEID=commandArgs(trailingOnly=TRUE)[1]
# cat(SAMPLEID)
SAMPLEDIR <- str_replace(SAMPLEID, "\\.", "/") 
  workdir <- glue("/volumes/USR1/junke/Projects/copykit/benchmark/sc_hmmcopy/{SAMPLEID}")
  system(glue("mkdir -p {workdir}"))
  setwd(workdir)
  system(glue("mkdir -p input pipeline temp"))
  system(glue("cp /volumes/USR1/junke/Projects/copykit/benchmark/sc_hmmcopy/test/input/config.yaml input/"))

  outfile <- "input/bams.yaml"
  bamfolder <- glue("/volumes/USR1/junke/Projects/copykit/downsampling/data/{SAMPLEDIR}")
  bamFile <- list.files(bamfolder, pattern = '*.bam$')
  for(i in 1:length(bamFile)){
    bami <- bamFile[i]
    system(glue("samtools index {bamfolder}/{bami}"))
    cellname <- strsplit(bami, ".", fixed=T)[[1]][1]
    write_lines(glue("{cellname}: "), outfile,append = T)
    write_lines(glue("  bam: {bamfolder}/{bami}"), outfile,append = T)
    write_lines(glue("  column: 01"), outfile,append = T)
    write_lines(glue("  condition: A"), outfile,append = T)
    write_lines(glue("  img_col: {i}"), outfile,append = T)
    write_lines(glue("  index_i5: i5-{i}"), outfile,append = T)
    write_lines(glue("  index_i7: i7-{i}"), outfile,append = T)
    write_lines(glue("  pick_met: CELL"), outfile,append = T)
    write_lines(glue("  primer_i5: ACTACTATT"), outfile,append = T)
    write_lines(glue("  primer_i7: AGTAGTACT"), outfile,append = T)
    write_lines(glue("  row: {i}"), outfile,append = T)
    write_lines(glue("  sample_type: C"), outfile,append = T)
    write_lines(glue("  sample_id: {SAMPLEID}"), outfile,append = T)
    write_lines(glue("  library_id: {SAMPLEID}"), outfile,append = T)
    write_lines(glue("  is_control: FALSE"), outfile,append = T)
  }
  cmd <- glue("nice -n 10 singularity exec --bind /volumes/ --cleanenv -H /volumes/USR1/junke/Projects/copykit/benchmark/sc_hmmcopy/{SAMPLEID}/ /volumes/USR1/junke/software/sc_hmmcopy/single_cell_pipeline_hmmcopy_v0.8.26.sif single_cell hmmcopy --input_yaml input/bams.yaml  --tmpdir temp/  --pipelinedir pipeline/  --library_id {SAMPLEID} --submit local --maxjobs 1 --nocleanup  --sentinel_only --loglevel DEBUG  --config_file input/config.yaml --output_prefix results")
  cat(cmd)

  