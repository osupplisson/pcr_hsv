# Function for Posterior predictive check----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pp_check_function <- function(fit_input, Nsample = 100) {
  dfbinomial <- dfbinomial %>%
    mutate(
      ageclass = cut(age, seq(0, 100, 10), include.lowest = TRUE),
      ageclassbis = cut(
        age,
        breaks = c(0, 24, 100),
        include.lowest = TRUE
      )
    )
  
  if (!exists("state")) {
    draw_pp <- generate(
      object = fit_input,
      newdata = dfbinomial,
      n.samples = Nsample,
      seed = 123,
      formula = fml
    )
  } else{
    print("Using state")
    draw_pp <- inlabru::evaluate_model(
      model = fit_input$bru_info$model,
      state = state,
      data = dfbinomial,
      predictor = fml
    )
  }
  
  draw_pp <- draw_pp %>%
    imap(~ cbind(.x, dfbinomial)) %>%
    imap(~ .x %>% mutate(sampleid = .y)) %>%
    do.call(rbind, .)
  
  
  draw_pp <-  draw_pp %>%
    mutate(
      ageclass = cut(age, seq(0, 100, 10), include.lowest = TRUE),
      ageclassbis = cut(
        age,
        breaks = c(0, 24, 100),
        include.lowest = TRUE
      )
    )
  
  
  fun_sum <- function(vect = c("virus", "ageclass", "sexe", "year", "nom"),
                      type) {
    obs_exp <- dfbinomial %>%
      group_by_at(vect) %>%
      summarise(obs = sum(result), Ntot = sum(N)) %>%
      ungroup() %>%
      mutate(prop = obs / Ntot)
    
    obs_prop <- dfbinomial %>%
      group_by_at(vect) %>%
      summarise(obs = sum(result), Ntot = sum(N)) %>%
      ungroup() %>%
      group_by_at(c(vect[-which(vect == "virus")])) %>%
      mutate(x = ifelse(virus == "hsv2", 1, 0)) %>%
      arrange(x) %>%
      mutate(Ntot = obs + obs[1],
             obs = obs,
             prop = obs / (obs + obs[1])) %>%
      filter(virus == "hsv2") %>%
      select(-x)
    
    if (type == "pp") {
      tmp <- draw_pp   %>%
        rowwise() %>%
        mutate(value = rbinom(1, size = N, prob = p)) %>%
        ungroup() %>%
        group_by_at(c(vect, "sampleid")) %>%
        summarise(value = sum(value), Ntot = sum(N))
      
      proportion <- tmp  %>%
        mutate(value = value / Ntot) %>%
        group_by_at(vect) %>%
        summary_draw() %>%
        left_join(obs_exp, by = vect)
      
      proportion_relative <- tmp %>%
        group_by_at(c(vect[-which(vect == "virus")], "sampleid")) %>%
        mutate(x = ifelse(virus == "hsv2", 1, 0)) %>%
        arrange(x) %>%
        mutate(value = value / (value + value[1])) %>%
        filter(virus == "hsv2") %>%
        group_by_at(vect) %>%
        summary_draw() %>%
        left_join(obs_prop, by = vect)
      
    } else{
      tmp <- draw_pp    %>%
        mutate(expected_n = p * N) %>%
        group_by_at(c(vect, "sampleid")) %>%
        summarise(value = sum(expected_n), Ntot = sum(N))
      
      proportion <- tmp %>%
        mutate(value = value / Ntot) %>%
        group_by_at(vect) %>%
        summary_draw() %>%
        left_join(obs_exp, by = vect)
      
      proportion_relative <- tmp %>%
        group_by_at(c(vect[-which(vect == "virus")], "sampleid")) %>%
        mutate(x = ifelse(virus == "hsv2", 1, 0)) %>%
        arrange(x) %>%
        mutate(value = value / (value + value[1])) %>%
        group_by_at(vect) %>%
        summary_draw() %>%
        filter(virus == "hsv2") %>%
        left_join(obs_prop, by = vect)
    }
    return(list(
      "proportion" = proportion,
      "proportion_relative" = proportion_relative
    ))
  }
  out <- list()
  out[["pp"]][["v"]] <- fun_sum(vect = c("virus"), type = "pp")
  out[["pp"]][["va"]] <- fun_sum(vect = c("virus", "ageclass"), type = "pp")
  out[["pp"]][["vabis"]] <- fun_sum(vect = c("virus", "ageclassbis"), type = "pp")
  out[["pp"]][["vs"]] <- fun_sum(vect = c("virus", "sexe"), type = "pp")
  out[["pp"]][["vy"]] <- fun_sum(vect = c("virus", "year"), type = "pp")
  out[["pp"]][["vn"]] <- fun_sum(vect = c("virus", "nom"), type = "pp")
  out[["pp"]][["vsa"]] <- fun_sum(vect = c("virus", "sexe", "ageclass"),
                                  type = "pp")
  
  out[["exp"]][["v"]] <- fun_sum(vect = c("virus"), type = "exp")
  out[["exp"]][["va"]] <- fun_sum(vect = c("virus", "ageclass"), type = "exp")
  out[["exp"]][["vabis"]] <- fun_sum(vect = c("virus", "ageclassbis"), type = "exp")
  out[["exp"]][["vs"]] <- fun_sum(vect = c("virus", "sexe"), type = "exp")
  out[["exp"]][["vy"]] <- fun_sum(vect = c("virus", "year"), type = "exp")
  out[["exp"]][["vn"]] <- fun_sum(vect = c("virus", "nom"), type = "exp")
  out[["exp"]][["vsa"]] <- fun_sum(vect = c("virus", "sexe", "ageclass"),
                                   type = "exp")
  
  
  return(out)
}

# Function for formating --------------------------------------------------
numeric_to_text <- function(df) {
  space <- dfbinomial %>% distinct(insee_dep, .keep_all = T) %>% select(idspace, insee_dep, nom, insee_reg, nom_region)
  df %>%
    mutate(
      virus = ifelse(hsv2 == 1, 'HSV-2', 'HSV-1'),
      groupsex = ifelse(females == 1, "Females", "Males")
    ) %>%
    select(-hsv2) %>%
    left_join(space, by = "idspace") %>%
    left_join(dfdate, by = "idtime")
}

# Function to compute AME -------------------------------------------------
fun_sum_ame <- function(d, vect, input) {
  ##Proportion of GUD due to HSV-1 and HSV-2
  tmp <- d %>%
    group_by_at(c(vect, "sampleid", "group")) %>%
    summarise(expected_n = sum(expected_n), N = sum(N)) %>%
    group_by_at(c(vect, "sampleid")) %>%
    arrange(group) %>%
    mutate(value = (expected_n - expected_n[1]) / N)
  
  ame <- tmp %>%
    filter(group == 2) %>%
    group_by_at(vect) %>%
    summary_draw()
  
  prev <- tmp  %>%
    group_by_at(c(vect, "group")) %>%
    mutate(value = expected_n / N) %>%
    summary_draw()
  
  #Relative proportion of anogenital herpes due to HSV-2
  prop <- d %>%
    group_by_at(c(vect, "sampleid", "group")) %>%
    summarise(expected_n = sum(expected_n)) %>%
    group_by_at(c(vect[-which(vect == "virus")], "sampleid", "group")) %>%
    mutate(x = ifelse(virus == "hsv2", 1, 0)) %>%
    arrange(x) %>%
    mutate(value = expected_n / (expected_n + expected_n[1])) %>%
    filter(virus == "hsv2") %>%
    group_by_at(c(vect, "group")) %>%
    summary_draw()
  
  
  prop_ame <- d %>%
    group_by_at(c(vect, "sampleid", "group")) %>%
    summarise(expected_n = sum(expected_n)) %>%
    group_by_at(c(vect[-which(vect == "virus")], "sampleid", "group")) %>%
    mutate(x = ifelse(virus == "hsv2", 1, 0)) %>%
    arrange(x) %>%
    mutate(value = expected_n / (expected_n + expected_n[1])) %>%
    filter(virus == "hsv2") %>%
    arrange(group) %>%
    group_by_at(c(vect, "sampleid")) %>%
    mutate(value = (value - value[1])) %>%
    group_by_at(c(vect, "group")) %>%
    summary_draw()
  
  ame <- ame %>%
    mutate(
      type_ame = input$type_ame,
      ref_ame = input$ref_ame,
      ref = input$ref,
      new = input$new
    )
  
  prev <- prev %>%
    mutate(
      type_ame = input$type_ame,
      ref_ame = input$ref_ame,
      ref = input$ref,
      new = input$new
    ) %>%
    mutate(group = ifelse(group == 1, "Reference", "New"))
  
  prop <- prop %>%
    mutate(
      type_ame = input$type_ame,
      ref_ame = input$ref_ame,
      ref = input$ref,
      new = input$new
    ) %>%
    mutate(group = ifelse(group == 1, "Reference", "New"))
  
  prop_ame <- prop_ame %>%
    mutate(
      type_ame = input$type_ame,
      ref_ame = input$ref_ame,
      ref = input$ref,
      new = input$new
    )
  
  
  list(
    "ame" = ame,
    "expected_prev" = prev,
    "prop" = prop,
    "prop_ame" = prop_ame
  )
}

ame_function <- function(fit_input,
                         Nsample = 100,
                         type_ame = "district",
                         ref_ame = "observed",
                         ref = NA,
                         new = 1) {
  if (ref_ame == "observed") {
    ref <- NA
  }
  input_list <- list(
    "type_ame" = type_ame,
    "ref_ame" = ref_ame,
    "ref" = ref,
    "new" = new
  )
  
  print(input_list)
  
  if (type_ame == "district") {
    if (ref_ame == "observed") {
      data_set <- rbind(
        dfbinomial %>% mutate(group = 1),
        dfbinomial %>% mutate(group = 2) %>% mutate(idspace = new)
      )
      dep_ref <- NA
    } else{
      data_set <- rbind(
        dfbinomial %>% mutate(group = 1) %>% mutate(idspace = ref),
        dfbinomial %>% mutate(group = 2) %>% mutate(idspace = new)
      )
    }
  } else if (type_ame == "age") {
    if (ref_ame == "observed") {
      data_set <- rbind(
        dfbinomial %>% mutate(group = 1),
        dfbinomial %>% mutate(group = 2) %>% mutate(age = new)
      )
      age_ref <- NA
    } else{
      data_set <- rbind(
        dfbinomial %>% mutate(group = 1) %>% mutate(age = ref),
        dfbinomial %>% mutate(group = 2) %>% mutate(age = new)
      )
    }
  } else if (type_ame == "sexe") {
    if (ref == "observed") {
      data_set <- rbind(
        dfbinomial %>% mutate(group = 1),
        dfbinomial %>% mutate(group = 2) %>% mutate(females = new)
      )
    } else{
      data_set <- rbind(
        dfbinomial %>% mutate(group = 1) %>% mutate(females = ref),
        dfbinomial %>% mutate(group = 2) %>% mutate(females = new)
      )
    }
  } else if (type_ame == "time")  {
    if (ref == "observed") {
      data_set <- rbind(
        dfbinomial %>% mutate(group = 1),
        dfbinomial %>% mutate(group = 2) %>% mutate(idtime = new)
      )
    } else{
      data_set <- rbind(
        dfbinomial %>% mutate(group = 1) %>% mutate(idtime = ref),
        dfbinomial %>% mutate(group = 2) %>% mutate(idtime = new)
      )
    }
  }
  
  print(unique(data_set$idspace))
  
  
  if (!exists("state")) {
    draw <- generate(
      object = fit_input,
      newdata = data_set,
      n.samples = Nsample,
      seed = 123,
      formula = fml
    )
  } else{
    print("Using state")
    draw <- inlabru::evaluate_model(
      model = fit_input$bru_info$model,
      state = state,
      data = data_set,
      predictor = fml
    )
  }
  
  
  draw <- draw %>%
    imap( ~ cbind(.x, data_set)) %>%
    imap( ~ .x %>% mutate(sampleid = .y)) %>%
    do.call(rbind, .) %>%
    mutate(expected_n = p * N)
  
  out <- list()
  out[["v"]] <- fun_sum_ame(d = draw,
                            vect = c("virus"),
                            input = input_list)
  if (type_ame == "age") {
    out[["vs"]] <- fun_sum_ame(d = draw,
                               vect = c("virus", "sexe"),
                               input = input_list)
  }
  out[["input"]] <- input_list
  rm("draw")
  gc()
  return(out)
}
