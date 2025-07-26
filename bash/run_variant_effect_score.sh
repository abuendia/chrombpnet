export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

CUDA_VISIBLE_DEVICES=0 python /oak/stanford/groups/akundaje/abuen/personal_genome/abuen_fork/personal_genome/chrombpnet/evaluation/variant_effect_prediction/snp_scoring.py \
  -snps /oak/stanford/groups/akundaje/abuen/personal_genome/abuen_fork/personal_genome/data/variant_benchmark/caqtls.african.lcls.asb.benchmarking.all.tsv.gz \
  -g /oak/stanford/groups/akundaje/abuen/personal_genome/abuen_fork/personal_genome/data/hg38.fa \
  -m /oak/stanford/groups/akundaje/abuen/personal_genome/chrombpnet/pipeline_train_output/models/chrombpnet.h5 \
  -op /oak/stanford/groups/akundaje/abuen/personal_genome/abuen_fork/personal_genome/output/ref_perf \
  -bs 64
