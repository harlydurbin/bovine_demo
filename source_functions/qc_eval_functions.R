


plot_filter_density <-
  
  function(df, var) {
    fvar <- rlang::enquo(var)
    
    df %>%
      ggplot(aes(x = !!fvar)) +
      geom_density(alpha = 0.3) +
      labs(
        x = rlang::quo_name(fvar),
        y = "Kernel density",
        title = glue::glue("Density of {rlang::quo_name(fvar)} values, chromosome 28")
      )
    
  }


var_sum <-
  function(var) {
    var <- rlang::enquo(var)
    
    df %>%
      filter(!is.infinite(!!var)) %>%
      summarise(
        !!glue::glue("{rlang::quo_name(var)}_min") := min(!!var, na.rm = TRUE),
        !!glue::glue("{rlang::quo_name(var)}_max") := max(!!var, na.rm = TRUE),
        !!glue::glue("{rlang::quo_name(var)}_mean") := mean(!!var, na.rm = TRUE)
      )
    
  }
