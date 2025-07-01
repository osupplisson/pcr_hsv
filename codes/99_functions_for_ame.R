# Function for Posterior predictive check----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pp_check_function <- function(draws = draw_pp) {
  library(tidytable)
  setDTthreads(threads = cpudt)
  draw_pp <- draws %>%
    imap(function(.x, .y) {
      m <- .y
      .x %>%
        imap(
          ~ .x %>%
            mutate(
              sample = paste(m, .y, sep = "_")
            ) %>%
            cbind(., dfbinomial %>% st_drop_geometry())
        ) %>%
        do.call(rbind, .)
    }) %>%
    do.call(rbind, .)
  detach("package:tidytable", unload = TRUE)
  
  
  draw_pp <-  draw_pp %>%
    mutate(
      ageclass = cut(age, seq(0, 100, 10), include.lowest = TRUE),
      ageclassbis = cut(
        age,
        breaks = c(0, 24, 100),
        include.lowest = TRUE
      )
    )
  
  dfbinomial <- dfbinomial %>%
    mutate(
      ageclass = cut(age, seq(0, 100, 10), include.lowest = TRUE),
      ageclassbis = cut(
        age,
        breaks = c(0, 24, 100),
        include.lowest = TRUE
      )
    )
  
  print(colnames(dfbinomial))
  print(colnames(draw_pp))
  
  fun_sum <- function(vect = c("virus", "ageclass", "sexe", "year", "nom"),
                      type) {
    print(vect)
    print(type)
    library(tidytable)
    setDTthreads(threads = cpudt)
    obs_exp <- dfbinomial %>%
      group_by(all_of(vect)) %>%
      summarise(obs = sum(result), Ntot = sum(N)) %>%
      ungroup() %>%
      mutate(prop = obs / Ntot)
    
    obs_prop <- dfbinomial %>%
      group_by(all_of(vect)) %>%
      summarise(obs = sum(result), Ntot = sum(N)) %>%
      ungroup() %>%
      group_by(all_of(c(vect[-which(vect == "virus")]))) %>%
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
        group_by(all_of(c(vect, "sample"))) %>%
        summarise(y = sum(value), Ntot = sum(N))
      
      proportion <- tmp  %>%
        mutate(value = y / Ntot) %>%
        group_by(all_of(vect)) %>%
        summary_draw() %>%
        left_join(obs_exp, by = vect)
      
      proportion_relative <- tmp %>%
        group_by(all_of(c(vect[-which(vect == "virus")], "sample"))) %>%
        mutate(x = ifelse(virus == "hsv2", 1, 0)) %>%
        arrange(x) %>%
        mutate(value = y / (y + y[1])) %>%
        filter(virus == "hsv2") %>%
        group_by(vect) %>%
        summary_draw() %>%
        left_join(obs_prop, by = vect)
      
    } else{
      tmp <- draw_pp    %>%
        mutate(expected_n = p * N) %>%
        group_by(all_of(c(vect, "sample"))) %>%
        summarise(y = sum(expected_n), Ntot = sum(N))
      
      proportion <- tmp %>%
        mutate(value = y / Ntot) %>%
        group_by(all_of(vect)) %>%
        summary_draw() %>%
        left_join(obs_exp, by = vect)
      
      proportion_relative <- tmp %>%
        group_by(all_of(c(vect[-which(vect == "virus")], "sample"))) %>%
        mutate(x = ifelse(virus == "hsv2", 1, 0)) %>%
        arrange(x) %>%
        mutate(value = y / (y + y[1])) %>%
        group_by(vect) %>%
        summary_draw() %>%
        filter(virus == "hsv2") %>%
        left_join(obs_prop, by = vect)
    }
    detach("package:tidytable", unload = TRUE)
    return(
      list(
        "proportion" = proportion %>% as_tibble(),
        "proportion_relative" = proportion_relative %>% as_tibble()
      )
    )
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
  print(length(unique(d$sample)))
  library(tidytable)
  setDTthreads(threads = cpudt)
  ##Proportion of GUD due to HSV-1 and HSV-2
  tmp <- d %>%
    group_by(all_of(c(vect, "sample", "group"))) %>%
    summarise(expected_n = sum(expected_n), N = sum(N)) %>%
    group_by(all_of(c(vect, "sample"))) %>%
    arrange(group) %>%
    mutate(value = (expected_n - expected_n[1]) / N)
  
  print(length(unique(tmp$sample)))
  
  ame <- tmp %>%
    filter(group == 2) %>%
    group_by(all_of(vect)) %>%
    summary_draw()
  
  #print(ame$Ndraws)
  
  prev <- tmp  %>%
    group_by(all_of(c(vect, "group"))) %>%
    mutate(value = expected_n / N) %>%
    summary_draw()
  
  #print(prev$Ndraws)
  
  #Relative proportion of anogenital herpes due to HSV-2
  prop <- d %>%
    group_by(all_of(c(vect, "sample", "group"))) %>%
    summarise(expected_n = sum(expected_n)) %>%
    group_by(all_of(c(vect[-which(vect == "virus")], "sample", "group"))) %>%
    mutate(x = ifelse(virus == "hsv2", 1, 0)) %>%
    arrange(x) %>%
    mutate(value = expected_n / (expected_n + expected_n[1])) %>%
    filter(virus == "hsv2") %>%
    group_by(all_of(c(vect, "group"))) %>%
    summary_draw()
  
  
  prop_ame <- d %>%
    group_by(all_of(c(vect, "sample", "group"))) %>%
    summarise(expected_n = sum(expected_n)) %>%
    group_by(all_of(c(vect[-which(vect == "virus")], "sample", "group"))) %>%
    mutate(x = ifelse(virus == "hsv2", 1, 0)) %>%
    arrange(x) %>%
    mutate(value = expected_n / (expected_n + expected_n[1])) %>%
    filter(virus == "hsv2") %>%
    arrange(group) %>%
    group_by(all_of(c(vect, "sample"))) %>%
    mutate(value = (value - value[1])) %>%
    group_by(all_of(c(vect, "group"))) %>%
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
  
  detach("package:tidytable", unload = TRUE)
  list(
    "ame" = ame %>% as_tibble(),
    "expected_prev" = prev %>% as_tibble(),
    "prop" = prop %>% as_tibble(),
    "prop_ame" = prop_ame %>% as_tibble()
  )
}

ame_function <- function(model_set,
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
  } else if (type_ame == "districttime")  {
    data_set <- rbind(
      dfbinomial %>% mutate(group = 1) %>% mutate(idspace = new, idtime = 1),
      dfbinomial %>% mutate(group = 2) %>% mutate(idspace = new, idtime = max(dfdate$idtime))
    )
  }
  
  draws <- list()
  for (i in 1:nrow(model_set)) {
    print(i)
    setting_paths(i)
    print("Ndraws:")
    print(ndrawsame)
    print("Using state")
    state_to_use <- state[seq(1, ndrawsame, by = 1L)]
    print(length(state_to_use))
    
    draws[[i]] <- inlabru::evaluate_model(
      model = info,
      state = state_to_use,
      data = data_set,
      predictor = fml
    )
    #draws[[i]] %>% imap(~is.null(.x))  %>% discard(., isFALSE) %>% print
    #print(draws[[i]])
    rm("info")
    rm("state_to_use")
    rm("state")
    gc()
  }
  gc()
  ndraw_sim <- draws %>%
    imap( ~ length(.x))
  print(model_set$draws_ame)
  print(Reduce(c, ndraw_sim))
  print(sum(Reduce(c, ndraw_sim)))
  print(model_set$draws_ame == Reduce(c, ndraw_sim))

  library(tidytable)
  setDTthreads(threads = cpudt)
  draws <-  draws %>%
    imap(function(.x, .y) {
      m <- .y
      .x %>%
        imap(~ .x %>%
               mutate(sample = paste(m, .y, sep = "_")) %>%
               cbind(., data_set %>% st_drop_geometry())) %>%
        do.call(rbind, .)
    }) %>%
    do.call(rbind, .) %>%
    mutate(expected_n = p * N)
  
  print(length(unique(draws$sample)))
  X <- tibble(s = unique(draws$sample)) %>%
    separate(col = s, into = c("model", "draw"), sep = "_")
  print(table(as.numeric(X$model)))
  print(model_set$draws_ame)
  gc()
  detach("package:tidytable", unload = TRUE)
  
  out <- list()
  out[["v"]] <- fun_sum_ame(d = draws,
                            vect = c("virus"),
                            input = input_list)
  gc()
  # if (type_ame == "age") {
  #   out[["vs"]] <- fun_sum_ame(d = draws,
  #                              vect = c("virus", "sexe"),
  #                              input = input_list)
  #   gc()
  # }
  out[["input"]] <- input_list
  rm("draws")
  gc()
  
  
  test_hpc(hpc)
  print(out$v$ame$Ndraws)
  print(out$v$ame$checkNA)
  return(out)
}
