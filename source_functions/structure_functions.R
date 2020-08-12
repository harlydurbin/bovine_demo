source(here::here("source_functions/coalesce_join.R"))

read_structure <-
  function(dataset_thinning, k) {
    # population ID key
    pop_key <-
      readr::read_table2(
        here::here(
          glue::glue(
            "data/derived_data/plink_qc/thin_variants/{dataset_thinning}/merge_thinned.{dataset_thinning}.fam"
          )
        ),
        col_names = FALSE,
        col_types = readr::cols(.default = "c")
      ) %>%
      dplyr::select(international_id = X2) %>%
      dplyr::left_join(
        readr::read_csv(
          here::here(
            "data/derived_data/metadata/bovine_demo.sample_metadata.csv"
          )
        ) %>%
          dplyr::mutate(
            international_id = dplyr::if_else(
              international_id == "CPC98_Bprimigenius_EnglandDerbyshire_5936",
              "ancient_derbyshire",
              international_id
            )
          ) %>%
          dplyr::select(international_id, population, region, continent, species)
      )
    
    # Read in mean Q file, append population label
    dat <-
      readr::read_table2(here::here(
        glue::glue(
          "data/derived_data/faststructure/structure/{dataset_thinning}/structure.{dataset_thinning}.{k}.meanQ"
        )
      ), col_names = FALSE) %>%
      dplyr::bind_cols(pop_key) %>%
      dplyr::filter(!is.na(population))
    
    
    # Establish column order to sort and arrange individuals on on
    col_order <-
      # Take the column names of `dat` (one column for each ancestry component)
      colnames(dat) %>%
      dplyr::as_tibble() %>%
      # exclude pop label
      dplyr::filter(str_detect(value, "X")) %>%
      dplyr::arrange(value) %>%
      dplyr::pull(value)
    
    # Reordering mean Q file
    dat %<>%
      # Arrange by X1 then X2 then X3 and so on so on
      dplyr::arrange(!!!rlang::syms(col_order)) %>%
      # Based on rearrangement, set sample ID factor levels, capitalize population labels
      dplyr::mutate(international_id = forcats::fct_inorder(international_id))
    
    dat %<>%
      mutate(
        species = 
          case_when(
            species %in% c("bos grunniens", "bos mutus") ~ "Wild & domestic yak",
            species %in% c("bos frontalis", "bos gaurus") ~ "Gaur & gayal",
            species == "bos javanicus" ~ "Banteng",
            species == "bos primigenius" ~ "Ancient",
            species == "bison bison" ~ "American bison",
            TRUE ~ stringr::str_to_sentence(species)),
        species = forcats::fct_relevel(
          as.factor(species),
          "Ancient",
          "Bos taurus", 
          "Composite", 
          "Bos indicus",
          "Wild & domestic yak",
          "Gaur & gayal",
          "Banteng",
          "American bison"),
        continent = case_when(
          continent %in% c("americas", "australia") ~ "Americas & Australia",
          population == "ancient" ~ "Ancient",
          TRUE ~ str_to_title(continent)
        ),
        continent = forcats::fct_relevel(
          as.factor(continent),
          "Ancient",
          "Americas & Australia",
          "Europe",
          "Africa",
          "Asia"
        ),
        region = case_when(
          region %in% c("americas", "australia") ~ "Americas & Australia",
          region == "india and pakistan" ~ "India & Pakistan",
          region == "japan and korea" ~ "Japan & Korea",
          population == "ancient" ~ "Ancient",
          TRUE ~ str_to_title(region)
        ),
        region = forcats::fct_relevel(
          as.factor(region),
          "Ancient",
          "Americas & Australia",
          "British Isles",
          "France",
          "Continental Europe",
          "Scandinavia",
          "Iberia",
          "Italy",
          "Podolian Steppe",
          "West Africa",
          "East Africa",
          "Middle East",
          "India & Pakistan",
          "Japan & Korea",
          "Southeast Asia",
          "Asian Steppe",
          "Northeast China",
          "East Central China",
          "South Central China",
          "Southeast China"
        ),
        population = case_when(
          population == "yak qtp" ~ "Yak QTP",
          TRUE ~ stringr::str_to_title(population)
        )
      )
    
    dat %<>%
      # Go from wide df to long df: results in data frame with columns below plus a column with X1...Xn and a column with the corresponding assignment value
      reshape2::melt(id = c(
        "international_id",
        "population",
        "region",
        "continent",
        "species"
      )) %>%
      # Reorder population factor levels
      dplyr::arrange(continent, region, population) %>%
      dplyr::mutate(population = forcats::fct_inorder(population))
    
    return(dat)
    
  }

# Function to calculate the mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
  
}


# Janky function foor keeping colors consistent between structure plots
# color_key <-
#   function(df) {
#     extra_colors <- c("#8e4866", "#84d6d3")
#     
#     color_key <-
#       df %>%
#       group_by(international_id) %>%
#       # This part's the important part: giving each id a "likely assignment"
#       # (i.e, this sample is most likely mostly African etc)
#       mutate(likely_assignment = variable[which.max(value)]) %>%
#       ungroup() %>%
#       group_by(region) %>%
#       # Find which ancestry component likely corresponds to each con
#       summarise(variable = Mode(likely_assignment)) %>%
#       ungroup() %>%
#       mutate(
#         color_var =
#           case_when(
#             region %in% c("British Isles") ~ "#126180",
#             region %in% c("India & Pakistan") ~ "#ff7436",
#             region %in% c("Southeast China") ~ "#629d62"
#           )
#       ) %>%
#       filter(!is.na(color_var)) %>%
#       distinct(variable, color_var)
#     
#     
#     
#     color_key <-
#       df %>%
#       distinct(variable) %>%
#       mutate(color_var = NA_character_) %>%
#       coalesce_join(color_key,
#                     by = "variable",
#                     join = dplyr::left_join) %>%
#       # Need to fix this
#       mutate(color_var = if_else(is.na(color_var),
#                                  sample(extra_colors, 1),
#                                  color_var))
#     
#     df <-
#       df %>%
#       left_join(color_key)
#     
#     return(df)
#   }

ggstructure <- function(df, wrap_var, custom_color = FALSE) {
  plot <-
    if (custom_color == FALSE) {
      df  %>%
        ggplot2::ggplot(ggplot2::aes(x = international_id,
                                     y = value,
                                     fill = variable))
    } else {
      df  %>%
        ggplot2::ggplot(ggplot2::aes(x = international_id,
                                     y = value,
                                     fill = color_var)) +
        ggplot2::scale_fill_identity()
      
    }
  
  
  # Make a stacked bar plot that's scaled by ancestry percentage using the scales package
  plot <-
    plot +
    ggplot2::geom_bar(position = "fill",
                    stat = "identity",
                    width = 3.5) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::guides(fill = FALSE) +
    ggplot2::theme_minimal() +
    # Remove all x axis title info (i.e., international IDs)
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.spacing.y = ggplot2::unit(0.6, "lines"),
      panel.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(angle = 90,
                                           size = 11),
      strip.text.y = ggplot2::element_text(size = 28, face = "bold")
    ) +
    ggforce::facet_row(
      facets = wrap_var,
      scales = "free_x",
      space = "free",
      labeller = ggplot2::label_wrap_gen()
    )
  
  return(plot)
  
}