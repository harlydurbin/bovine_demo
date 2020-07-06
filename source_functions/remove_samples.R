library(readr)
library(dplyr)
library(tidyr)
library(tidylog)

sample_metadata <- read_csv(here::here("data/derived_data/metadata/coverage/coverage.sample_metadata.csv"))

# List of duplicates, unknown pops, low coverage to be removed

remove <-
  read_tsv(
    here::here(
      "data/derived_data/joint_genotyping/find_dups/drop_dups.txt"
    ),
    col_names = "international_id"
  ) %>%
  bind_rows(read_csv(
    here::here("data/derived_data/joint_genotyping/remove_unknown.csv")
  ) %>%
    select(international_id)) %>%
  bind_rows(
    sample_metadata %>%
      filter(5 >= avg_coverage) %>%
      filter(species %in% c("bos taurus", "bos indicus", "composite")) %>%
      select(international_id)
  ) %>% 
  distinct()


remove %>% 
  write_tsv(here::here("data/derived_data/joint_genotyping/remove_samples/remove.txt"), col_names = FALSE)

sample_metadata %>% 
  filter(!international_id %in% remove$international_id) %>% 
  write_csv(here::here("data/derived_data/metadata/bovine_demo.sample_metadata.csv"))
