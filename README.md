# Demographic and selection history of cattle and related species

## Meta-data processing

1. SRA metadata scraped from downloaded XML files then initially tidied in `notebooks/xml_scraping.Rmd`
2. Population labels futher categorized by region, continent, and species in `notebooks/categorizing.Rmd`
3. Coverage of samples in tidied dataset assessed `notebooks/coverage.Rmd`
4. Duplicate and low quality samples removed from metadata in `source_functions/remove_samples.R`

**Resulting sample metadata file: `data/derived_data/metadata/bovine_demo.sample_metadata.csv`**

### Notes about population label assignment

* I know little to nothing about yak breeds, but where I could I tried to use the same population designations for yak samples that were used in the papers they came from. These include:
    * Wild yaks (*Bos mutus*)
    * Datong yaks which were recently developed as a cross between wild and domesticated yaks [to be hornless](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0158642)
    * [Tianzhu white yaks](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0158642), bred in the Qilian mountains of Gansu province
    * Jinchuan yaks, which typically have an [additional thoracic vertebra](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-6598-9) and are found in Sichuan province
    * [Qinghai-Tibet Plateau (QTP)](https://www.nature.com/articles/ncomms10283) yaks
* Cattle samples with their species listed as "composite" are known taurus/indicus hybrids. These include:
    * Breeds intentionally developed within the last 100 years in America and Australia (Beefmaster, Droughtmaster, Santa Gertrudis)
    * African sanga & zenga breeds (Ankole, Boran, Fulani)
    * Asian advanced generation composites (all others)

## Genotype calling

Evaluation of computing resources used for each step of genotype calling can be found in `notebooks/psrecord.Rmd` with results in `html/psrecord.html`

1. Ancient samples pre-processed and haplotypes called in `source_functions/ancient_preprocess.snakefile`
    * Trim 5 bp from ends of reads using `bamUtil trimBam`
    * Realign indels using `GATK IndelRealigner`
    * Call haplotypes using `GATK HaplotypeCaller`
2. All samples combined and genotypes called in `source_functions/joint_genotyping.snakefile`
    * Sample cohorts combined using `GATK CombineGVCFs` (cohorts determined in `source_functions/genotyping_cohorts.R`)
    * Joint genotypes called using `GATK GenotypeGVCFs`
3. After joint genotype calling, INFO field filter values and depth of coverage at each variant on chromosome 28 extracted in `source_functions/filter_eval.snakefile` using `GATK VariantsToTable` and `vcftools --site-mean-depth`. Descriptive statistics & distribution of these values explored in `notebooks/filter_eval.Rmd`. Results can be found in `html/filter_eval.html` and were used to inform filtering cutoffs in the next step
4. Callset filtered in `source_functions/joint_genotyping.snakefile`
    * Variants restricted to biallelic SNPs using `GATK SelectVariants`
    * Site-level and genotype-level filtering annotated using `GATK VariantFiltration`. Then failing sites removed and failing genotypes set to missing using `GATK SelectVariants`
5. Summary stats for each chromosome generated using `Picard CollectVariantCallingMetrics` then evaluated in `source_functions/joint_genotyping.Rmd`, VCF format checked using `GATK ValidateVariants`. `CollectVariantCallingMetrics` results:
    * **data/derived_data/joint_genotyping/bovine_demo.variant_metrics.summary_chr.csv** contains a summary by chromosome
    * **data/derived_data/joint_genotyping/bovine_demo.variant_metrics.detail_wg.csv** contains a summary by sample averaged/summed across all chromosomes
    * **data/derived_data/joint_genotyping/bovine_demo.variant_metrics.detail_chr.csv** contains a summary by sample and by chromosome
6. Duplicate samples identified based on chromosome 28 and chromosome 29 variants using [`king`](http://people.virginia.edu/~wc9c/KING/manual.html) in `source_functions/find_dups.sh`. Duplicates and other low quality samples + sites that are no longer variant after sample removal removed in `source_functions/reheader.snakemake`

## Phasing

See `source_functions/phasing.snakefile` and `notebooks/phasing.Rmd`

1. In order to phase X chromosome, missing sexes imputed and incorrectly assigned sexes fixed 
    * Ended up using the ratio of average coverage on the X chromosome/average coverage on all autosomes to determine cutoffs. Of all other tested metrics, I think this the only one that should be similar across all species in the dataset + agnostic to effective population size
2. Genetic map inferred using several published cattle recombination maps *TODO*
3. Pre-phasing QC
    * For all chromosomes, sites with > 10% missingness removed
    * For all chromosomes, listed sex updated to imputed sex
    * Pseudo-autosomal region removed from X chromosome
    * Heterozygous genotypes set to missing on Y chromosome
4. Phase autosomes and sex chromosomes separately using `SHAPEIT` *TODO*

## Exploratory analyses

* `fastStructure`
    * Using output of pre-phasing QC in `source_functions/phasing.snakefile`, variants removed from each chromosome with an X% probability of being retained, downsample individuals
    * See `source_functions/faststructure.bovine_demo.snakefile` for analysis and `notebooks/faststructure.Rmd` for thinning/downsampling dataset designations & results
* `EIGENSOFT smartpca`
* [`SMC++`](https://github.com/popgenmethods/smcpp)