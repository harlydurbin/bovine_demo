plot_pca <-
  
  function(dataset_thinning,
           pc1,
           pc2,
           color_var,
           shape_var = NULL,
           title_string = NULL) {
    col1 <- rlang::sym(str_c("PC", as.character(pc1)))
    col2 <- rlang::sym(str_c("PC", as.character(pc2)))
    
    shape_var <- rlang::enquo(shape_var)
    
    color_var <- rlang::enquo(color_var)
    
    
    path <-
      glue::glue(
        "data/derived_data/smartpca/{dataset_thinning}/smartpca.{dataset_thinning}.evec"
      )
    
    val <-
      read_lines(here::here(path), n_max = 1) %>%
      str_remove("#eigvals:") %>%
      str_squish() %>%
      str_split(pattern = "[[:space:]]") %>%
      flatten_chr() %>%
      map_dbl( ~ as.numeric(.x))
    
    vec <-
      read_table2(here::here(path), skip = 1, col_names = FALSE) %>%
      select(-X12) %>%
      set_names(c("international_id", map_chr(1:10, ~ str_c("PC", .x)))) %>%
      mutate(international_id = str_remove(international_id, ":[[:graph:]]+")) %>%
      left_join(
        read_csv(here::here("data/derived_data/metadata/bovine_demo.sample_metadata.csv")) %>% 
          mutate(international_id = if_else(international_id == "CPC98_Bprimigenius_EnglandDerbyshire_5936", "ancient_derbyshire", international_id)) %>%
          select(international_id, population, region, continent, species)
      ) %>%
      mutate(
        species = stringr::str_to_sentence(species),
        region = stringr::str_to_title(region),
        continent = stringr::str_to_title(continent),
        population = stringr::str_to_title(population)
      )
    
    plot <-
      if (is.null(shape_var)) {
        vec %>%
          ggplot(aes(
            x = !!col1,
            y = !!col2
          )) +
          geom_point(aes(color = !!color_var),
                     alpha = 0.6)
      } else{
        vec %>%
          ggplot(aes(
            x = !!col1,
            y = !!col2
          )) +
          geom_point(aes(color = !!color_var, 
                         shape = !!shape_var),
                     alpha = 0.6)
      }
    
    plot <-
      plot +
      bakeoff::scale_color_bakeoff() +
      guides(colour = guide_legend(override.aes = list(alpha = 1))) +
      theme_classic() +
      labs(
        x = str_c("PC ", pc1, ": ", scales::percent(val[pc1] / sum(val))),
        y = str_c("PC ", pc2, ": ", scales::percent(val[pc2] / sum(val))),
        color = stringr::str_to_sentence(rlang::quo_name(rlang::enquo(color_var)))
      )
    
    plot <-
      if (is.null(shape_var)) {
        plot
      } else{
        plot +
          labs(shape = stringr::str_to_sentence(rlang::quo_name(shape_var)))
      }
    
    plot <-
      if(is.null(title_string)){
        plot +
          labs(title = glue::glue("{dataset_thinning}: PC {pc1} vs. PC {pc2}"))
      } else{
        plot +
          labs(title = glue::glue("{title_string}: PC {pc1} vs. PC {pc2}"))
      }
    
    return(plot)
  }

ggscree <- function(dataset_thinning, title_string = NULL) {
  
  path <-
    glue::glue(
      "data/derived_data/smartpca/{dataset_thinning}/smartpca.{dataset_thinning}.evec"
    )
  
  val <-
    read_lines(here::here(path), n_max = 1) %>%
    str_remove("#eigvals:") %>%
    str_squish() %>%
    str_split(pattern = "[[:space:]]") %>%
    flatten_chr() %>%
    as_tibble() %>%
    mutate(
      PC = as.character(glue::glue("PC{row_number()}")),
      PC = forcats::fct_inorder(PC),
      value = as.double(value),
      varexp = round((value / sum(value))*100, digits = 2)
    )
  
  plot <-
    val %>%
    ggplot2::ggplot() +
    ggplot2::geom_col(mapping = ggplot2::aes(x = PC, y = varexp), stat = "identity") +
    #scale_y_continuous(limits = c(0, 15)) +
    geom_text(
      aes(
        x = PC,
        y = varexp,
        label = str_c(varexp, "%")
      ),
      position = position_dodge(width = 0.9),
      vjust = -0.25
    ) +
    theme_classic() +
    labs(x = NULL,
         y = "% Variation explained") 
  
  plot <-
    if(is.null(title_string)){
      plot +
        labs(title = dataset_thinning)
    } else{
      plot +
        labs(title = title_string)
    }
  
  return(plot)
  
}

read_vec <- 
  
  function(dataset_thinning) {
    
    path <-
      glue::glue(
        "data/derived_data/smartpca/{dataset_thinning}/smartpca.{dataset_thinning}.evec"
      )
    
    vec <-
      read_table2(here::here(path), skip = 1, col_names = FALSE) %>%
      select(-X12) %>%
      set_names(c("international_id", map_chr(1:10, ~ str_c("PC", .x)))) %>%
      mutate(international_id = str_remove(international_id, ":[[:graph:]]+")) %>%
      left_join(
        sample_metadata %>%
          select(international_id, population, region, continent, species)
      ) 
    
    return(vec)
    
  }