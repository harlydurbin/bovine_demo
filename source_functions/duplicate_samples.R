
library(readr) 
library(dplyr)
library(tidyr)

# Notes & questions

# IBS0 -- the number of sites where one sample is hom-ref and another is hom-alt
# IBS2 -- the number of sites where the samples have the same genotype

# Setup

sample_metadata <-
  read_csv(here::here("data/derived_data/coverage/coverage.sample_metadata.csv"))

dups <-
  read_table2(here::here("data/derived_data/joint_genotyping/filter_eval/filter_eval.with_indels.28.con")) %>% 
  janitor::clean_names() %>% 
  select(-contains("fid")) 

# Choose which of each duplicate pair to keep 

drop <-
  dups %>% 
  left_join(sample_metadata %>% 
              select(id1 = international_id, cov1 = avg_coverage)) %>% 
  left_join(sample_metadata %>% 
              select(id2 = international_id, cov2 = avg_coverage)) %>% 
  mutate_at(vars(contains("cov")), ~ replace_na(., 0)) %>% 
  mutate(drop = if_else(cov1 > cov2 , id2, id1)) %>% 
  select(id1, cov1, id2, cov2, drop, everything()) %>% 
  pull(drop)

sample_metadata %>% 
  filter(!international_id %in% drop) %>% 
  write_csv(here::here("data/derived_data/sample_metadata.csv"), na = "")
