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
      # Based on rearrangement, set sample ID factor levels
      dplyr::mutate(international_id = forcats::fct_inorder(international_id)) %>%
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

ggstructure <- function(df, wrap_var) {
  df  %>%
    ggplot2::ggplot(ggplot2::aes(x = international_id,
                                 y = value,
                                 fill = variable)) +
    # Make a stacked bar plot that's scaled by ancestry percentage using the scales package
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
  
  
}