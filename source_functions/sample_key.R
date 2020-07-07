library(readr)
library(dplyr)

sample_metadata <- read_csv(here::here("data/derived_data/metadata/bovine_demo.sample_metadata.csv"))

read_table2(here::here("data/derived_data/joint_genotyping/reheader/sample_list.txt"), col_names = "international_id") %>%
  left_join(sample_metadata %>% 
              select(international_id, lab_id)) %>% 
  write_tsv(here::here("data/derived_data/joint_genotyping/reheader/sample_key.txt"), col_names = FALSE)
