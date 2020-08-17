# f4 test of admixture ratios using TreeMix fourpop output

admix_prop <-
  function(df, o, x, a, b, c) {
    
    numerator <-
      df %>% 
      filter_at(vars(contains("combo")), any_vars(str_detect(., glue("{o},{a}|{a},{o}")))) %>% 
      filter_at(vars(contains("combo")), any_vars(str_detect(., glue("{x},{c}|{c},{x}")))) %>% 
      mutate(position = "numerator")
    # pull(f4)
    
    denominator <-
      df %>% 
      filter_at(vars(contains("combo")), any_vars(str_detect(., glue("{o},{a}|{a},{o}")))) %>% 
      filter_at(vars(contains("combo")), any_vars(str_detect(., glue("{b},{c}|{c},{b}")))) %>% 
      mutate(position = "denominator")
    #  pull(f4)
    
    
    alpha <- (numerator$f4_abs[[1]])/(denominator$f4_abs[[1]])
    
    tribble(
      ~pop, ~position, ~estimate, ~test,
      b, "B", alpha, glue("({numerator$combo_1[[1]]};{numerator$combo_2[[1]]})/({denominator$combo_1[[1]]};{denominator$combo_2[[1]]})"),
      c, "C", (1-alpha), glue("({numerator$combo_1[[1]]};{numerator$combo_2[[1]]})/({denominator$combo_1[[1]]};{denominator$combo_2[[1]]})")
    )
  }