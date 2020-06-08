flatten_biosample <-
  
  function(biosample_list, index_num, approved_ids) {
    # print(glue::glue("processing index {index_num}"))
    
    ids <-
      purrr::set_names(biosample_list[[index_num]][["Ids"]],
                       # For the nth BioSample,
                       # Set the name of each element in "Ids" to its attribute name
                       purrr::map(
                         .x = seq_along(biosample_list[[index_num]][["Ids"]]),
                         ~ attributes(biosample_list[[index_num]][["Ids"]][[.x]])[1] %>%
                           unlist() %>%
                           unname()
                       )) %>%
      purrr::flatten() %>%
      # Janitor takes care of duplicate column names
      dplyr::as_tibble(.name_repair = janitor::make_clean_names)
    
    
    long_ids <-
      ids %>%
      tidyr::pivot_longer(everything())
    
    if (any(long_ids$value %in% approved_ids)) {
      attribute <-
        purrr::set_names(
          biosample_list[[index_num]][["Attributes"]],
          purrr::map(
            .x = seq_along(biosample_list[[index_num]][["Attributes"]]),
            ~ attributes(biosample_list[[index_num]][["Attributes"]][[.x]])[1] %>%
              unlist() %>%
              unname()
          )
        ) %>%
        purrr::flatten() %>%
        dplyr::as_tibble(.name_repair = janitor::make_clean_names)
      
      # Pluck out owner if one is listed
      contact <-
        biosample_list[[index_num]][["Owner"]] %>%
        purrr::pluck("Name", .default = NA_character_) %>%
        unlist()
      
      # Pluck out species name from Description list if one is listed
      species <-
        biosample_list[[index_num]][["Description"]] %>%
        purrr::pluck("Organism", "OrganismName", .default = NA_character_) %>%
        unlist()
      
      # Pluck out project title from Description if one is listed
      desc <-
        biosample_list[[index_num]][["Description"]] %>%
        purrr::pluck("Title", .default = NA_character_) %>%
        unlist()
      
      flattened <-
        dplyr::bind_cols(ids) %>%
        dplyr::mutate(organism = species,
                      title = desc,
                      owner = contact)
      # Pivot attributes wide ---> long to reduce the number of columns
      
      flattened <-
        if (length(attribute) > 0) {
          flattened %>%
            bind_cols(attribute) %>% 
            tidyr::pivot_longer(cols = colnames(attribute), names_to = "attribute")
        } else {flattened}
      
      return(flattened)
    }
    
  }