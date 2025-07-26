sample_id=GM18498
sample_id_in_vcf=NA18498

rm -rf /oak/stanford/groups/akundaje/abuen/personal_genome/personal_genome/personalized_output_${sample_id}

CUDA_VISIBLE_DEVICES=0 chrombpnet pipeline \
  --bigwig /oak/stanford/groups/akundaje/ziwei75/african-omics/data/processed/${sample_id}/${sample_id}_unstranded.bw \
  -d "ATAC" \
  -g /oak/stanford/groups/akundaje/abuen/personal_genome/personal_genome/data/hg38.fa \
  -c /oak/stanford/groups/akundaje/abuen/personal_genome/personal_genome/data/hg38.chrom.sizes \
  -p /oak/stanford/groups/akundaje/ziwei75/african-omics/data/encode_output/${sample_id}/peak/overlap_reproducibility/overlap.conservative_peak.narrowPeak.gz \
  -n /oak/stanford/groups/akundaje/ziwei75/african-omics/data/processed/${sample_id}/${sample_id}_negatives.bed \
  -fl /oak/stanford/groups/akundaje/abuen/personal_genome/personal_genome/data/fold_0.json \
  -b /oak/stanford/groups/akundaje/abuen/personal_genome/personal_genome/bias_train_output/models/bias.h5 \
  -o /oak/stanford/groups/akundaje/abuen/personal_genome/personal_genome/personalized_output_${sample_id} \
  --vcf-file /oak/stanford/groups/akundaje/ziwei75/african-omics/data/genotype_new/hwepass01.withgnomadfreq.var.merged.processedMKK.split1kG30x.bcf \
  --sample-id ${sample_id_in_vcf}
