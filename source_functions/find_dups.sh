#!/bin/bash

srun -c 12 --mem=24G plink --vcf data/derived_data/joint_genotyping/select_variants/select_variants.28.vcf.gz --make-bed --double-id --cow -set-missing-var-ids @:#\$1\$2 --threads 12 --out select_variants.28 &> select_variants.28.plink.log
srun -c 12 --mem=24G source_functions/king -b select_variants.28.bed --duplicate --prefix select_variants.28 --sexchr 30 --cpus 12

srun -c 12 --mem=24G plink --vcf data/derived_data/joint_genotyping/select_variants/select_variants.29.vcf.gz --make-bed --double-id --cow -set-missing-var-ids @:#\$1\$2 --threads 12 --out select_variants.29 &> select_variants.29.plink.log
srun -c 12 --mem=24G source_functions/king -b select_variants.29.bed --duplicate --prefix select_variants.29 --sexchr 30 --cpus 12
