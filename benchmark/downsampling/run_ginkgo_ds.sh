## run ginkgo for ds
samp="pmtc1breast"
export samp
foo(){
  mkdir -p /volumes/USR1/junke/Projects/copykit/ginkgo/ginkgo/uploads/${samp}.${1}
  cd /volumes/USR1/junke/Projects/copykit/ginkgo/ginkgo/uploads/${samp}.${1}
  cp /volumes/USR1/junke/Projects/copykit/downsampling/data/${samp}/${1}_bed/* .
  ls *.bed>list
  cp ../simul_clone1/config .
  cd ..
  /volumes/seq/code/PIPELINES/copy_number_pipeline_hg38/dev_lib/code/copy_number_pipeline_hg38/ginkgo/scripts/analyze.sh ${samp}.${1}
}
export -f foo

cd /volumes/USR1/junke/Projects/copykit/downsampling/data/${samp}/
  ls --ignore="*_bed" | parallel -j 8 foo 

## run schmmcopy for ds
for samp in `echo "fulldata" "downsampled_sort_50k" "downsampled_sort_75k" "downsampled_sort_125k" "downsampled_sort_250k" "downsampled_sort_500k" "downsampled_sort_750k" "downsampled_sort_10k" "downsampled_sort_25k"`;
do
SAMPID="pmtc1breast.${samp}";
cmd=`Rscript ~/Projects/copykit/downsampling/run_schmmcopy_ds.R ${SAMPID}`;
cd /volumes/USR1/junke/Projects/copykit/benchmark/sc_hmmcopy/${SAMPID};
$cmd;
done


