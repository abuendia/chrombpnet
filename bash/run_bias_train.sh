sample_id=GM18498

CUDA_VISIBLE_DEVICES=0 chrombpnet bias pipeline \
  --bigwig /oak/stanford/groups/akundaje/ziwei75/african-omics/data/processed/GM18498/GM18498_unstranded.bw \
  -d "ATAC" \
  -g /oak/stanford/groups/akundaje/abuen/personal_genome/chrombpnet/data/hg38.fa \
  -c /oak/stanford/groups/akundaje/abuen/personal_genome/chrombpnet/data/hg38.chrom.sizes \
  -p /oak/stanford/groups/akundaje/ziwei75/african-omics/data/encode_output/${sample_id}/peak/overlap_reproducibility/overlap.conservative_peak.narrowPeak.gz \
  -n /oak/stanford/groups/akundaje/ziwei75/african-omics/data/processed/GM18498/GM18498_negatives.bed \
  -fl /oak/stanford/groups/akundaje/abuen/personal_genome/chrombpnet/data/fold_0.json \
  -b 0.5 \
  -o /oak/stanford/groups/akundaje/abuen/personal_genome/chrombpnet/output
