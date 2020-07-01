# Demographic and selection history of cattle and related species

## Genotype calling

Evaluation of computing resources used for each step of genotype calling can be found in `notebooks/psrecord.Rmd` with results in `html/psrecord.html`

1. Ancient sample BAM files pre-processed and haplotypes called in `source_functions/ancient_preprocess.snakefile`
2. All samples combined and genotypes called in `source_functions/joint_genotyping.snakefile`
    + Sample cohorts combined using `GATK CombineGVCFs` (cohorts determined in `source_functions/genotyping_cohorts.R`)
    + Joint genotypes called using `GATK GenotypeGVCFs`
3. After joint genotype calling, INFO field filter values and depth of coverage at each variant on chromosome 28 extracted in `source_functions/filter_eval.snakefile` using `GATK VariantsToTable` and `vcftools --site-mean-depth`. Descriptive statistics & distribution of these values explored in `notebooks/qc_eval.Rmd`. Results can be found in `html/qc_eval.html` and were used to inform filtering cutoffs in the next step
4. Variant callset filtered in `source_functions/joint_genotyping.snakefile`
    + Variants restricted to biallelic SNPs using `GATK SelectVariants`
    + Site-level and genotype-level filtering annotated using `GATK VariantFiltration`. Then failing sites removed + failing genotypes set to missing using `GATK SelectVariants`
5. Chromosome-by-chromosome files concatenated in order to one whole-genome file using `Picard GatherVcfs`
6. Create summary of the resulting callset using `Picard CollectVariantCallingMetrics` then evaluate in `source_functions/joint_genotyping.Rmd`
7. Duplicate samples identified using `king` then removed in `source_functions/find_dups.snakefile`

## Meta-data processing

1. SRA metadata scraped from downloaded XML files then initially tidied in `notebooks/xml_scraping.Rmd`
2. Population labels futher categorized by region, continent, and species in `notebooks/categorizing.Rmd`
3. Coverage of samples in tidied dataset assessed `notebooks/coverage.Rmd`
4. Duplicate samples removed from metadata in `source_functions/duplicate_samples.R`

**Resulting sample metadata file: `data/derived_data/sample_metadata.csv`**