library(dplyr)
library(readr)
library(stringr)

fam <- as.character(commandArgs(trailingOnly = TRUE)[1])

out <- as.character(commandArgs(trailingOnly = TRUE)[2])

sample_metadata <-
  read_csv(here::here("data/derived_data/metadata/bovine_demo.sample_metadata.csv")) %>%
  mutate(international_id = if_else(international_id == "CPC98_Bprimigenius_EnglandDerbyshire_5936", "ancient_derbyshire", international_id))


read_table2(here::here(fam), col_names = FALSE) %>%
  select(international_id = X1, X2) %>%
  left_join(sample_metadata %>%
              select(international_id, population, region, species)) %>%
  mutate(population = case_when(
    international_id == "ancient_derbyshire" ~ "AncientDerbyshire",
    international_id == "A2494_Yenisei_River" ~ "AncientSiberian",
    species == "bos grunniens" ~ "YakDomestic",
    species == "bos mutus" ~ "YakWild",
    species == "bos frontalis" ~ "Gayal",
    species == "bos gaurus" ~ "Gaur",
    species == "bos javanicus" ~ "Banteng",
    population == "hereford miniature" ~ "Hereford",
    population == "angus lowline" ~ "Angus",
    region == "podolian steppe" ~ "PodolianSteppe",
    TRUE ~ stringr::str_to_title(population)
  ),
  population = str_remove_all(population, " |[[:punct:]]")) %>%
  #select(-species) %>%
  group_by(population, species) %>% 
  mutate(
    drop = 
      case_when(
        population %in% c("Tibetan", "Mongolian", "ChaidamuYellow", "PodolianSteppe") ~ FALSE,
        population %in% c("Beefmaster", "EastAfricanZebu") ~ TRUE,
        n() < 5 ~ TRUE,
        TRUE ~ FALSE
      )
  ) %>% 
  ungroup() %>% 
  filter(drop == FALSE) %>%
  write_tsv(here::here(out), col_names = FALSE)
