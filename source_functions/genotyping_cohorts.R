#split list of lab ids into chunks of 200 for genotypeGVCFs

library(dplyr)
library(glue)
library(purrr)
library(stringr)

genotyping_cohorts <- function(df, chromosome, n = 50) {
  
  chr_dir <- glue::glue("data/derived_data/joint_genotyping/genotyping_cohorts/{chromosome}")
  
  # make a directory for the chromosome if it doesn't already exist
  if (!dir.exists(here::here(chr_dir))) {
    dir.create(here::here(chr_dir))
  }
  
  chunk <- n
  n <- nrow(df)
  r  <- rep(1:ceiling(n / chunk), each = chunk)[1:n]
  
  df %>%
    arrange(region, species, population, lab_id) %>% 
    mutate(cohort = r) %>% 
    mutate(
      path = as.character(glue::glue(
        "/storage/htc/schnabellab/results/9913/wgs/1kbulls_ars1.2/gvcf/{chromosome}/{lab_id}.{chromosome}.g.vcf.gz"
      ))
    ) %>% 
    select(path, cohort) %>%
    # Add ancient samples manually
    dplyr::add_row(
      path = as.character(glue::glue(
        "/storage/hpc/group/UMAG/WORKING/hjdzpd/bovine_diversity/data/derived_data/ancient_preprocess/derbyshire/haplotype_caller.derbyshire.{chromosome}.g.vcf.gz")),
      cohort = 1) %>% 
    dplyr::add_row(
      path = as.character(glue::glue(
        "/storage/hpc/group/UMAG/WORKING/hjdzpd/bovine_diversity/data/derived_data/ancient_preprocess/siberian/haplotype_caller.siberian.{chromosome}.g.vcf.gz")),
      cohort = 1) %>% 
    group_by(cohort) %>%
    group_walk( ~ write_tsv(.x,
                            here::here(
                              glue::glue("{chr_dir}/cohort_{.y$cohort}.list")
                            ),
                            col_names = FALSE))
  
}

sorted <- read_rds(here::here("data/derived_data/categorizing/sorted.sample_metadata.rds"))

cohorts <- 
  # List of samples that have gVCFs at /storage/htc/schnabellab/results/9913/wgs/1kbulls_ars1.2/gvcf/
  c("data/raw_data/cattle_gvcf.txt") %>% 
  purrr::set_names("cattle") %>% 
  purrr::map_dfr(~ read_table2(here::here(.x), col_names = "path"), .id = "group") %>% 
  mutate(lab_id = as.numeric(str_extract(path, "^[[:digit:]]+(?=\\.)")))%>% 
  distinct(group, lab_id) %>% 
  right_join(sorted %>% 
               select(lab_id, species, region, population)) %>% 
  filter(!is.na(group) & population != "ancient") %>% 
  select(lab_id, species, region, population)

walk(.x = c(1:29, "X", "Y"), ~ genotyping_cohorts(df = cohorts, chromosome = .x, n = 85))