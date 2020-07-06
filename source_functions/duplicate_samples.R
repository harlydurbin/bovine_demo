
library(readr) 
library(dplyr)
library(tidyr)

# Notes & questions

# IBS0 -- the number of sites where one sample is hom-ref and another is hom-alt
# IBS2 -- the number of sites where the samples have the same genotype

sample_metadata <- read_csv(here::here("data/derived_data/metadata/coverage/coverage.sample_metadata.csv"))

dups <-
  c(28, 29) %>% 
  purrr::set_names() %>% 
  purrr::map_dfr(~ read_table2(here::here(glue::glue("data/derived_data/joint_genotyping/find_dups/select_variants.{.x}.con"))), .id = "chr") %>% 
  janitor::clean_names() %>% 
  select(-contains("fid")) %>% 
  arrange(id1, id2) %>% 
  # Remove only if identified as a duplicate on both chr28 and chr29
  group_by(id1, id2) %>% 
  filter(n_distinct(chr) == 2) %>% 
  ungroup() %>% 
  left_join(sample_metadata %>% 
              select(id1 = international_id, cov1 = avg_coverage)) %>% 
  left_join(sample_metadata %>% 
              select(id2 = international_id, cov2 = avg_coverage)) %>% 
  select(id1, cov1, id2, cov2) %>% 
  distinct() %>% 
  mutate_at(vars(contains("cov")), ~ replace_na(., 0)) %>% 
  mutate(drop = if_else(cov1 > cov2 , id2, id1)) %>% 
  select(drop) %>% 
  distinct() 

# This is lazy but couldn't figure out how to programatically remove these samples that appear in both id1 and id2
multiples <- tibble::tibble(drop = c("UMCUSAM000000196764", "UMCUSAU000000194604", "HOLGBRM000000598172"))

dups %>% 
  bind_rows(multiples) %>% 
  distinct() %>% 
  write_tsv(here::here("data/derived_data/joint_genotyping/find_dups/drop_dups.txt"), col_names = FALSE)
