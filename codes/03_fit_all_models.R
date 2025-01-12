print('Cleaning ent')
rm(list = ls())
# Packages ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
required_packages <- c(
  'tidybayes',
  'ggdist',
  'ggExtra',
  'ggpubr',
  'MatrixModels',
  'Matrix',
  'readr',
  'spdep',
  'sf',
  'janitor',
  'plyr',
  'INLA',
  'loo',
  'inlabru',
  'tidylog',
  'tidyverse'
)

# LOAD ALL PACKAGES
lapply(required_packages, library, character.only = TRUE)
set.seed(123)


sd_to_prec <- function(sigma) {
  tibble('sd' = sigma,
         'var' = sigma^2,
         'prec' = 1 / sigma^2)
}

# Function for summary  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
summary_draw <- function(df) {
  df %>%
    summarise(
      mean = mean(value, na.rm = T),
      median = median(value, na.rm = T),
      qi_lb = quantile(value, probs = 1 - 0.975, na.rm = T),
      qi_ub = quantile(value, probs = 0.975, na.rm = T),
      qi_lb_099 = quantile(value, probs = 1 - 0.99, na.rm = T),
      qi_ub_099 = quantile(value, probs = 0.99, na.rm = T),
      qi_lb_095 = quantile(value, probs = 1 - 0.95, na.rm = T),
      qi_ub_095 = quantile(value, probs = 0.95, na.rm = T),
      qi_lb_09 = quantile(value, probs = 1 - 0.9, na.rm = T),
      qi_ub_09 = quantile(value, probs = 0.9, na.rm = T),
      qi_lb_08 = quantile(value, probs = 1 - 0.8, na.rm = T),
      qi_ub_08 = quantile(value, probs = 0.8, na.rm = T),
      qi_lb_07 = quantile(value, probs = 1 - 0.7, na.rm = T),
      qi_ub_07 = quantile(value, probs = 0.7, na.rm = T),
      qi_lb_06 = quantile(value, probs = 1 - 0.6, na.rm = T),
      qi_ub_06 = quantile(value, probs = 0.6, na.rm = T),
      qi_lb_05 = quantile(value, probs = 1 - 0.5, na.rm = T),
      qi_ub_05 = quantile(value, probs = 0.5, na.rm = T),
      min = min(value, na.rm = T),
      max = max(value, na.rm = T),
      N = max(row_number()),
      N_higher_0 = sum(value > 0, na.rm = T),
      N_higher_05 = sum(value > 0.5, na.rm = T),
      N_higher_1 = sum(value > 1, na.rm = T),
      N_lower_0 = sum(value < 0, na.rm = T),
      N_lower_05 = sum(value < 0.5, na.rm = T),
      N_lower_1 = sum(value < 1, na.rm = T)
    ) %>%
    ungroup() %>%
    mutate(
      higher_0 = N_higher_0 / N,
      higher_05 = N_higher_05 / N,
      higher_1 = N_higher_1 / N,
      lower_0 = N_lower_0 / N,
      lower_05 = N_lower_05 / N,
      lower_1 = N_lower_1 / N
    )
}

source("pcr_hsv/codes/99_functions_for_ame.R", echo = T)


# Data --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('pcr_hsv/clean_data/input_models/data_for_fit.rda')


# Path --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
path_data <- 'pcr_hsv/clean_data/'
path_fit <- paste0(path_data, 'inla_fit/')



# Modifying variables -----------------------------------------------------
aggregated <- aggregated %>%
  mutate(
    intercept = 1,
    sexe = factor(sexe),
    year = factor(year),
    year2023 = ifelse(year == '2023', 1, 0)
  )


# Pivoting ----------------------------------------------------------------
dfbinomial <- aggregated   %>%
  pivot_longer(c(hsv1, hsv2), names_to = 'virus', values_to = 'result') %>%
  mutate(hsv1 = ifelse(virus == 'hsv1', 1, 0),
         hsv2 = ifelse(virus == 'hsv2', 1, 0)) %>%
  mutate(
    levelintercept = case_when(
      hsv2 == 0 & females == 0 ~ 1,
      hsv2 == 0 & females == 1 ~ 2,
      hsv2 == 1 & females == 0 ~ 3,
      hsv2 == 1 & females == 1 ~ 4
    )
  )
rm('aggregated')

# Neighborhood matrix ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
list_w <-
  c(names(spatial_objects$nb),
    names(spatial_objects$nb$nb_list))



sd_to_prec <- function(sigma) {
  tibble('sd' = sigma,
         'var' = sigma^2,
         'prec' = 1 / sigma^2)
}

inla.pardiso.check()
matrices <- c('soi', 'queen', 'delaunay')



w_set <- expand.grid(whsv1 = matrices, whsv2 = matrices) %>%
  mutate_at(vars(whsv1, whsv2), ~ factor(., levels = matrices))

print(w_set)

w_set <- w_set %>%
  mutate(link = "logit", type = "type1")

w_set <- w_set %>%
  rbind(w_set %>%
          mutate(link = "logit", type = "type2")) %>%
  rbind(w_set %>%
          mutate(link = "logit", type = "type3")) %>%
  rbind(w_set %>%
          mutate(link = "logit", type = "type4")) %>%
  arrange(row_number())


# ggplot(data.frame(x = seq(-30,30)), aes(x)) +
#   geom_function(fun = plogis, colour = "red")+
#   geom_function(fun = pnorm, colour = "blue")+
#   geom_function(fun = pchisq, colour = "green", args = list(df = 7))

return_matrix <- function(w) {
  test <- ifelse(str_detect(w, "nb"), 1, 0)
  test <- ifelse(str_detect(w, "gab"), 0, test)
  test <- ifelse(str_detect(w, "relative"), 0, test)
  print("Should be 1 if NN list:")
  print(w)
  print(test)
  if (test == 0) {
    W <- as(nb2mat(
      spatial_objects$nb[[w]],
      style = 'B',
      zero.policy = TRUE
    ),
    'Matrix')
  } else{
    W <- as(nb2mat(
      spatial_objects$nb$nb_list[[w]],
      style = 'B',
      zero.policy = TRUE
    ),
    'Matrix')
  }
  W
}

# Fit ---------------------------------------------------------------------
for (x in 1:nrow(w_set)) {
  print(x)
  w_type <- w_set %>%
    filter(row_number() == x) %>%
    mutate(whsv1 = as.character(whsv1), whsv2 = as.character(whsv2))
  print(w_type)
  name_fit <-
    paste0('fit_',
           paste(w_type$type, w_type$whsv1, w_type$whsv2, sep = '_'))
  
  path_fit <- 'pcr_hsv/clean_data/inla_fit/'
  list_fit <- base::list.files(
    path_fit,
    pattern = '.RDS',
    all.files = T,
    full.names = FALSE
  )
  
  
  if (!file.exists(paste0(path_fit, name_fit, '.RDS'))) {
    print('Fitting model')
    print(name_fit)
    
    strategy_input <- 'auto'
    int.strategy_input <- 'auto'
    init_bru <- NULL
    
    
    # https://stats.stackexchange.com/questions/454647/relation-between-gaussian-processes-and-gaussian-markov-random-fields
    prec_to_var <- function(prec) {
      tibble(
        'prec' = prec,
        'var' = 1 / prec,
        'std' = sqrt(1 / prec)
      )
    }
    
    
    Whsv1 <- return_matrix(w_type$whsv1)
    Whsv2 <- return_matrix(w_type$whsv2)
    type <- w_type$type
    
    #Baseline
    cmp <- ~ 0 + intercept(
      intercept,
      model = "linear",
      mean.linear = 0,
      prec.linear = 0.01
    ) +
      interceptadditional(levelintercept,
                          model = "factor_contrast",
                          hyper = list(prec = list(
                            prior = 'pc.prec', param = c(u = 1, a = 0.5)
                          )))
    
    if (type %in% c("type1", "type2", "type3")) {
      cmp <- update(
        cmp,
        ~ . +
          agehsv1(
            age + 1,
            model = "rw2",
            scale.model = TRUE,
            constr = TRUE,
            group = females + 1,
            ngroup = 2,
            control.group = list(model = "exchangeable", hyper = list(rho = list(
              prior = "normal", param = c(0, 0.2)
            ))),
            hyper = list(prec = list(
              prior = 'pc.prec', param = c(u = 1, a = 0.01)
            ))
          ) +
          agehsv2(
            age + 1,
            hsv2,
            model = "rw2",
            scale.model = TRUE,
            constr = TRUE,
            group = females + 1,
            ngroup = 2,
            control.group = list(model = "exchangeable", hyper = list(rho = list(
              prior = "normal", param = c(0, 0.2)
            ))),
            hyper = list(prec = list(
              prior = 'pc.prec', param = c(u = 1, a = 0.01)
            ))
          )
      )
    }
    
    
    if (type == "type1") {
      cmp <- update(
        cmp,
        ~ . +
          spacehsv1(
            idspace,
            group = females + 1,
            ngroup = 2,
            model = "bym2",
            graph = Whsv1,
            scale.model = TRUE,
            adjust.for.con.comp = TRUE,
            constr = TRUE,
            hyper = list(
              prec = list(prior = 'pc.prec', param = c(u = 1, a = 0.01)),
              phi = list(prior = 'pc', param =  c(u = 0.5, a = 0.5))
            ),
            control.group = list(model = "exchangeable", hyper = list(rho = list(
              prior = "normal", param = c(0, 0.2)
            )))
          ) +
          timehsv1(
            idtime,
            group = females + 1,
            ngroup = 2,
            model = "rw2",
            scale.model = TRUE,
            constr = TRUE,
            control.group = list(model = "exchangeable", hyper = list(rho = list(
              prior = "normal", param = c(0, 0.2)
            ))),
            hyper = list(prec = list(
              prior = 'pc.prec', param = c(u = 1, a = 0.01)
            ))
          ) +
          spacehsv2(
            idspace,
            hsv2,
            group = females + 1,
            ngroup = 2,
            model = "bym2",
            graph = Whsv2,
            scale.model = TRUE,
            adjust.for.con.comp = TRUE,
            constr = TRUE,
            hyper = list(
              prec = list(prior = 'pc.prec', param = c(u = 1, a = 0.01)),
              phi = list(prior = 'pc', param =  c(u = 0.5, a = 0.5))
            ),
            control.group = list(model = "exchangeable", hyper = list(rho = list(
              prior = "normal", param = c(0, 0.2)
            )))
          ) +
          timehsv2(
            idtime,
            hsv2,
            group = females + 1,
            ngroup = 2,
            model = "rw2",
            scale.model = TRUE,
            constr = TRUE,
            control.group = list(model = "exchangeable", hyper = list(rho = list(
              prior = "normal", param = c(0, 0.2)
            ))),
            hyper = list(prec = list(
              prior = 'pc.prec', param = c(u = 1, a = 0.01)
            ))
          )
      )
    }
    else if (type == "type2") {
      cmp <- update(
        cmp,
        ~ . +
          spacehsv1(
            idspace,
            group = idtime,
            ngroup = max(dfbinomial$idtime),
            model = "bym2",
            graph = Whsv1,
            scale.model = TRUE,
            adjust.for.con.comp = TRUE,
            constr = TRUE,
            hyper = list(
              prec = list(prior = 'pc.prec', param = c(u = 1, a = 0.01)),
              phi = list(prior = 'pc', param =  c(u = 0.5, a = 0.5))
            ),
            control.group = list(model = "rw2", scale.model = TRUE)
          ) +
          spacehsv2(
            idspace,
            hsv2,
            group = idtime,
            ngroup = max(dfbinomial$idtime),
            model = "bym2",
            graph = Whsv2,
            scale.model = TRUE,
            adjust.for.con.comp = TRUE,
            constr = TRUE,
            hyper = list(
              prec = list(prior = 'pc.prec', param = c(u = 1, a = 0.01)),
              phi = list(prior = 'pc', param =  c(u = 0.5, a = 0.5))
            ),
            control.group = list(model = "rw2", scale.model = TRUE)
          )
      )
    }
    else if (type == "type3") {
      cmp <- update(
        cmp,
        ~ . +
          spacehsv1(
            idspace,
            group = idtime,
            ngroup = max(dfbinomial$idtime),
            model = "bym2",
            graph = Whsv1,
            scale.model = TRUE,
            adjust.for.con.comp = TRUE,
            constr = TRUE,
            hyper = list(
              prec = list(prior = 'pc.prec', param = c(u = 1, a = 0.01)),
              phi = list(prior = 'pc', param =  c(u = 0.5, a = 0.5))
            ),
            control.group = list(model = "ar1", hyper = list(rho = list(
              prior = "pc.cor0", param = c(0.5, 0.5)
            )))
          ) +
          spacehsv2(
            idspace,
            hsv2,
            group = idtime,
            ngroup = max(dfbinomial$idtime),
            model = "bym2",
            graph = Whsv2,
            scale.model = TRUE,
            adjust.for.con.comp = TRUE,
            constr = TRUE,
            hyper = list(
              prec = list(prior = 'pc.prec', param = c(u = 1, a = 0.01)),
              phi = list(prior = 'pc', param =  c(u = 0.5, a = 0.5))
            ),
            control.group = list(model = "ar1", hyper = list(rho = list(
              prior = "pc.cor0", param = c(0.5, 0.5)
            )))
          ) +
          spacefemaleshsv1(
            idspace,
            females,
            group = idtime,
            ngroup = max(dfbinomial$idtime),
            model = "bym2",
            graph = Whsv1,
            scale.model = TRUE,
            adjust.for.con.comp = TRUE,
            constr = TRUE,
            hyper = list(
              prec = list(prior = 'pc.prec', param = c(u = 1, a = 0.01)),
              phi = list(prior = 'pc', param =  c(u = 0.5, a = 0.5))
            ),
            control.group = list(model = "ar1", hyper = list(rho = list(
              prior = "pc.cor0", param = c(0.5, 0.5)
            )))
          ) +
          spacefemaleshsv2(
            idspace,
            hsv2 * females,
            group = idtime,
            ngroup = max(dfbinomial$idtime),
            model = "bym2",
            graph = Whsv2,
            scale.model = TRUE,
            adjust.for.con.comp = TRUE,
            constr = TRUE,
            hyper = list(
              prec = list(prior = 'pc.prec', param = c(u = 1, a = 0.01)),
              phi = list(prior = 'pc', param =  c(u = 0.5, a = 0.5))
            ),
            control.group = list(model = "ar1", hyper = list(rho = list(
              prior = "pc.cor0", param = c(0.5, 0.5)
            )))
          )
      )
    }
    else if (type == "type4") {
      cmp <- update(
        cmp,
        ~ . +
          spacehsv1(
            idspace,
            group = idage,
            replicate = females + 1,
            model = "bym2",
            graph = Whsv1,
            scale.model = TRUE,
            adjust.for.con.comp = TRUE,
            constr = TRUE,
            hyper = list(
              prec = list(prior = 'pc.prec', param = c(u = 1, a = 0.01)),
              phi = list(prior = 'pc', param =  c(u = 0.5, a = 0.5))
            ),
            control.group = list(model = "rw2")
          ) +
          spacehsv2(
            idspace,
            hsv2,
            group = idage,
            replicate = females + 1,
            model = "bym2",
            graph = Whsv2,
            scale.model = TRUE,
            adjust.for.con.comp = TRUE,
            constr = TRUE,
            hyper = list(
              prec = list(prior = 'pc.prec', param = c(u = 1, a = 0.01)),
              phi = list(prior = 'pc', param =  c(u = 0.5, a = 0.5))
            ),
            control.group = list(model = "rw2")
          ) +
          timehsv1(
            idtime,
            group = females + 1,
            ngroup = 2,
            model = "rw2",
            scale.model = TRUE,
            constr = TRUE,
            control.group = list(model = "exchangeable", hyper = list(rho = list(
              prior = "normal", param = c(0, 0.2)
            ))),
            hyper = list(prec = list(
              prior = 'pc.prec', param = c(u = 1, a = 0.01)
            ))
          ) +
          timehsv2(
            idtime,
            group = females + 1,
            ngroup = 2,
            model = "rw2",
            scale.model = TRUE,
            constr = TRUE,
            control.group = list(model = "exchangeable", hyper = list(rho = list(
              prior = "normal", param = c(0, 0.2)
            ))),
            hyper = list(prec = list(
              prior = 'pc.prec', param = c(u = 1, a = 0.01)
            ))
          )
      )
      
    }
    
    
    
    fam <- list(control.link = list(model = "logit"))
    
    #Likelihood
    lik <- bru_obs(
      data = dfbinomial,
      Ntrials = N,
      control.family = fam,
      family = 'binomial',
      formula = result ~ .
    )
    
    #Bru call
    fit <- bru(
      components = cmp,
      lik,
      options = bru_options(
        bru_verbose = 4,
        verbose = T,
        bru_initial = NULL,
        control.compute = list(
          waic = TRUE,
          dic = TRUE,
          config = TRUE
        ),
        control.inla = list(
          strategy = "auto",
          int.strategy = "auto",
          optimise.strategy = "smart",
          fast = TRUE,
          compute.initial.values = T,
          stencil = 9,
          dz = 0.1,
          diff.logdens = 0.1,
          numint.maxfeval = 800000000,
          stupid.search = TRUE,
          restart = 4
        ),
        num.threads = "84:1",
        safe = T
      )
    )
    
    print("Adding graphs")
    fit[["Whsv1"]] <- w_type$whsv1
    fit[["Whsv2"]] <- w_type$whsv2
    fit[["graphWhsv1"]] <- Whsv1
    fit[["graphWhsv2"]] <- Whsv2
    fit[["type"]] <- w_type$type
    
    print('Computing CV')
    fit$loocv <-
      inla.group.cv(fit, num.level.sets = -1, strategy = 'posterior')
    for (x in c(-1, 5, 10, 15, 20, 25, 30)) {
      print(x)
      if (x > 0) {
        namex <- paste0('lgocv.m', x)
      } else{
        namex <- 'loo'
      }
      fit[['cv']][[namex]] <-
        inla.group.cv(
          fit,
          num.level.sets = x,
          strategy = 'posterior',
          size.max = 32
        )
    }
    print('SAVE')
    saveRDS(fit, paste0(path_fit, name_fit, '.RDS'))
  }
}


if (!file.exists('pcr_hsv/clean_data/assessment_metric.RDS')) {
  # Computing assessment metrics ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  list_fit <- base::list.files(
    path_fit,
    pattern = '.RDS',
    all.files = T,
    full.names = T
  ) %>%
    discard(str_detect(., 'init'))
  
  list_fit <- list_fit %>%
    purrr::set_names()
  names(list_fit) <-
    str_remove(str_remove(names(list_fit), 'pcr_hsv/clean_data/inla_fit/fit_'),
               '.RDS')
  
  
  summary <- list_fit %>%
    imap(function(.x, .y) {
      model <- readRDS(.x)
      cv <- model$cv %>%
        imap(~ tibble(type = .y, score = -sum(log(.x$cv)))) %>%
        do.call(rbind, .) %>%
        pivot_wider(names_from = type, values_from = score)
      
      
      metric <- tibble(
        fit = .y,
        "Whsv1" = model$Whsv1,
        "Whsv2" = model$Whsv2,
        "type" = model$type,
        waic = model$waic$waic,
        dic = model$dic$dic
      ) %>%
        cbind(cv)
      
      rm("model")
      gc()
      
      return(metric)
    })
  
  gc(full = T)
  assessment_metric <- summary %>%
    imap(function(.x, .y) {
      .x
    }) %>%
    do.call(rbind, .)
  
  #https://github.com/inlabru-org/inlabru/discussions/144
  
  
  print(assessment_metric)
  
  #Saving assessment
  saveRDS(assessment_metric,
          'pcr_hsv/clean_data/assessment_metric.RDS')
} else{
  assessment_metric <- readRDS('pcr_hsv/clean_data/assessment_metric.RDS')
}

# Importing best model ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
best_model <- assessment_metric  %>%
  arrange(lgocv.m30) %>%
  filter(row_number() == 1)
print(best_model)
# Structure of the model --------------------------------------------------
fit <- readRDS(paste0(best_model$fit, ".RDS"))
typeopt <- fit[["type"]]
Whsv1opt <- fit[["Whsv1"]]
Whsv2opt <- fit[["Whsv2"]]
Whsv1 <- fit[["graphWhsv1"]]
Whsv2 <- fit[["graphWhsv2"]]

saveRDS(Whsv1opt, 'pcr_hsv/clean_data/whsv1opt.RDS')
saveRDS(Whsv2opt, 'pcr_hsv/clean_data/whsv2opt.RDS')


if (typeopt == "type1") {
  fml <- ~ tibble(
    p = plogis(
      intercept + interceptadditional +  agehsv1 + spacehsv1 + timehsv1 + agehsv2 + spacehsv2 + timehsv2
    )
  )
} else if (typeopt %in% c("type2", "type3")) {
  fml <- ~ tibble(p = plogis(
    intercept + interceptadditional +  agehsv1 + spacehsv1 + agehsv2 + spacehsv2
  ))
}else if (typeopt %in% c("type4")) {
  fml <- ~ tibble(p = plogis(
    intercept + interceptadditional + spacehsv1 + spacehsv2 + timehsv1 + timehsv2
  ))
}


# Saving summary ----------------------------------------------------------
summary <- list(
  "fixed" = fit$summary.fixed,
  "random" = fit$summary.random,
  "hyper" = fit$summary.hyperpar,
  "marginal" = fit$marginals.random,
  "marginalhyper" = fit$marginals.hyperpar
)
saveRDS(summary, 'pcr_hsv/clean_data/summary_fixed_and_hyper.RDS')

rm("W")



# Generate draws ----------------------------------------------------------
if (!file.exists('pcr_hsv/clean_data/state.RDS')) {
  Nsample <- 3000
  data_grid <- expand.grid(
    age = seq(0, 100, 1),
    females = c(0, 1),
    hsv2 = c(0, 1),
    idtime = seq(1, max(dfbinomial$idtime)),
    idspace = sort(unique(dfbinomial$idspace))
  ) %>%
    mutate(
      levelintercept = case_when(
        hsv2 == 0 & females == 0 ~ 1,
        hsv2 == 0 & females == 1 ~ 2,
        hsv2 == 1 & females == 0 ~ 3,
        hsv2 == 1 & females == 1 ~ 4
      ),
      intercept = 1
    )
  
  
  state <- evaluate_state(
    model = fit$bru_info$model,
    result = fit,
    data = data_grid,
    property = 'sample',
    n = Nsample,
    seed = 123
  )
  
  saveRDS(state, 'pcr_hsv/clean_data/state.RDS')
} else{
  state <- readRDS('pcr_hsv/clean_data/state.RDS')
}

pass <- T
if (pass == F) {
  # PP ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  if (!file.exists('pcr_hsv/clean_data/pp_check.RDS')) {
    pp_check <- pp_check_function(fit_input = fit, Nsample = 2000)
    saveRDS(pp_check, 'pcr_hsv/clean_data/pp_check.RDS')
  }
  
  
  #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2661894/
  # Age ------------------------------------------------------------------
  if (!file.exists('pcr_hsv/clean_data/age_result.RDS')) {
    Nsample <- 2000
    # Table for draws ---------------------------------------------------------
    tible_pred <- expand.grid(
      age = seq(0, 100, 1),
      females = c(0, 1),
      hsv2 = c(0, 1),
      idtime = c(1, max(dfbinomial$idtime)),
      idspace = sort(unique(dfbinomial$idspace))
    ) %>%
      mutate(
        levelintercept = case_when(
          hsv2 == 0 & females == 0 ~ 1,
          hsv2 == 0 & females == 1 ~ 2,
          hsv2 == 1 & females == 0 ~ 3,
          hsv2 == 1 & females == 1 ~ 4
        ),
        intercept = 1
      )
    
    
    # Generate ----------------------------------------------------------------
    if (!exists("state")) {
      gen_age_sex <- generate(
        object = fit,
        n.samples = Nsample,
        newdata  = tible_pred,
        seed = 123L,
        formula  =  fml
      )
    } else{
      print("Using state")
      gen_age_sex <- inlabru::evaluate_model(
        model = fit$bru_info$model,
        state = state,
        data = tible_pred,
        predictor = fml
      )
    }
    
    
    # Prevalence --------------------------------------------------------------
    prev_age_sex <- gen_age_sex %>%
      imap(function(.x, .y) {
        .x %>% cbind(tible_pred) %>%
          mutate(sample = .y)
      }) %>%
      do.call(rbind, .) %>%
      numeric_to_text  %>%
      group_by(year_week,
               idtime,
               age,
               groupsex,
               virus,
               insee_dep,
               nom,
               insee_reg,
               nom_region) %>%
      dplyr::rename(value = p) %>%
      summary_draw()
    
    
    # Age 20---------------------------------------------------------------------
    tmp <- gen_age_sex %>%
      imap(function(.x, .y) {
        .x %>% cbind(tible_pred) %>%
          mutate(sample = .y)
      }) %>%
      do.call(rbind, .) %>%
      numeric_to_text
    
    fun_age_ref <- function(age_ref) {
      tmp %>%
        group_by(
          year_week,
          idtime,
          groupsex,
          virus,
          insee_dep,
          nom,
          insee_reg,
          nom_region,
          sample
        ) %>%
        mutate(x = ifelse(age == age_ref, 1, 0)) %>%
        arrange(desc(x)) %>%
        mutate(or = (p / (1 - p)) / (p[1] / (1 - p[1])),
               rr = p / p[1],
               rd = p - p[1]) %>%
        ungroup() %>%
        group_by(year_week,
                 idtime,
                 age,
                 groupsex,
                 virus,
                 insee_dep,
                 nom,
                 insee_reg,
                 nom_region)
    }
    diff_age_sex <- fun_age_ref(age_ref = 20)
    
    or_age <- diff_age_sex %>%
      mutate(value = or) %>%
      summary_draw
    
    rr_age <- diff_age_sex %>%
      mutate(value = rr) %>%
      summary_draw()
    
    rd_age <- diff_age_sex %>%
      mutate(value = rd) %>%
      summary_draw()
    
    # Age 55--------------------------------------------------------------------
    
    or_age_55 <- fun_age_ref(age_ref = 55) %>%
      mutate(value = or) %>%
      summary_draw
    
    # Age 70---------------------------------------------------------------------
    
    or_age_70 <- fun_age_ref(age_ref = 70) %>%
      mutate(value = or) %>%
      summary_draw
    
    
    # Sex ---------------------------------------------------------------------
    diff_age_sex <- tmp %>%
      group_by(year_week,
               idtime,
               age,
               virus,
               insee_dep,
               nom,
               insee_reg,
               nom_region,
               sample) %>%
      mutate(x = ifelse(groupsex == "Females", 1, 0)) %>%
      arrange(x) %>%
      mutate(or = (p / (1 - p)) / (p[1] / (1 - p[1])),
             rr = p / p[1],
             rd = p - p[1]) %>%
      ungroup() %>%
      group_by(year_week,
               groupsex,
               idtime,
               age,
               virus,
               insee_dep,
               nom,
               insee_reg,
               nom_region)
    
    
    or_sex <- diff_age_sex %>%
      mutate(value = or) %>%
      summary_draw
    
    rr_sex <- diff_age_sex %>%
      mutate(value = rr) %>%
      summary_draw()
    
    rd_sex <- diff_age_sex %>%
      mutate(value = rd) %>%
      summary_draw()
    
    
    # Virus -------------------------------------------------------------------
    diff_virus <- tmp %>%
      group_by(year_week,
               idtime,
               age,
               groupsex,
               insee_dep,
               nom,
               insee_reg,
               nom_region,
               sample) %>%
      mutate(x = ifelse(virus == "HSV-2", 1, 0)) %>%
      arrange(x) %>%
      mutate(
        or = (p / (1 - p)) / (p[1] / (1 - p[1])),
        rr = p / p[1],
        rd = p - p[1],
        prop = p / (p + p[1])
      ) %>%
      ungroup() %>%
      group_by(year_week,
               idtime,
               age,
               groupsex,
               virus,
               insee_dep,
               nom,
               insee_reg,
               nom_region)
    
    or_virus <- diff_virus %>%
      mutate(value = or) %>%
      summary_draw
    
    rr_virus <- diff_virus %>%
      mutate(value = rr) %>%
      summary_draw()
    
    rd_virus <- diff_virus %>%
      mutate(value = rd) %>%
      summary_draw()
    
    prop_virus <- diff_virus %>%
      mutate(value = prop) %>%
      summary_draw
    
    
    # District ---------------------------------------------------------------------
    diff_district <- tmp %>%
      group_by(year_week, idtime, groupsex, age, virus, sample) %>%
      mutate(x = ifelse(nom == "Paris", 1, 0)) %>%
      arrange(desc(x)) %>%
      mutate(or = (p / (1 - p)) / (p[1] / (1 - p[1])),
             rr = p / p[1],
             rd = p - p[1]) %>%
      ungroup() %>%
      group_by(year_week,
               idtime,
               age,
               groupsex,
               virus,
               insee_dep,
               nom,
               insee_reg,
               nom_region)
    
    or_district <- diff_district %>%
      mutate(value = or) %>%
      summary_draw
    
    rr_district <- diff_district %>%
      mutate(value = rr) %>%
      summary_draw()
    
    rd_district <- diff_district %>%
      mutate(value = rd) %>%
      summary_draw()
    
    
    
    # Time ---------------------------------------------------------------------
    diff_time <- tmp %>%
      group_by(groupsex,
               virus,
               age,
               insee_dep,
               nom,
               insee_reg,
               nom_region,
               sample) %>%
      mutate(x = ifelse(idtime == 1, 0, 1)) %>%
      arrange(x) %>%
      mutate(or = (p / (1 - p)) / (p[1] / (1 - p[1])),
             rr = p / p[1],
             rd = p - p[1]) %>%
      ungroup() %>%
      group_by(year_week,
               idtime,
               age,
               groupsex,
               virus,
               insee_dep,
               nom,
               insee_reg,
               nom_region)
    
    or_time <- diff_time %>%
      mutate(value = or) %>%
      summary_draw
    
    rr_time <- diff_time %>%
      mutate(value = rr) %>%
      summary_draw()
    
    rd_time <- diff_time %>%
      mutate(value = rd) %>%
      summary_draw()
    
    
    # Count sex ---------------------------------------------------------------
    diff_age_sex <- tmp %>%
      group_by(year_week,
               idtime,
               age,
               virus,
               insee_dep,
               nom,
               insee_reg,
               nom_region,
               sample) %>%
      mutate(x = ifelse(groupsex == "Females", 1, 0)) %>%
      arrange(x) %>%
      mutate(dummy = ifelse(p - p[1] > 0, 1, 0)) %>%
      filter(groupsex == "Females") %>%
      ungroup() %>%
      group_by(year_week, idtime, age, virus, sample) %>%
      summarise(value = sum(dummy))
    
    
    count_age_sex <- diff_age_sex %>%
      group_by(year_week, idtime, age, virus) %>%
      summary_draw
    
    # Count virus ---------------------------------------------------------------
    diff_age_sex <- tmp %>%
      group_by(year_week,
               idtime,
               groupsex,
               age,
               insee_dep,
               nom,
               insee_reg,
               nom_region,
               sample) %>%
      mutate(x = ifelse(virus == "HSV-2", 1, 0)) %>%
      arrange(x) %>%
      mutate(dummy = ifelse(p / (p + p[1]) > 0.5, 1, 0)) %>%
      filter(virus == "HSV-2") %>%
      ungroup() %>%
      group_by(year_week, idtime, age, groupsex, sample) %>%
      summarise(value = sum(dummy))
    
    
    count_age_virus <- diff_age_sex %>%
      group_by(year_week, idtime, age, groupsex) %>%
      summary_draw
    
    # OUT ---------------------------------------------------------------------
    age_result <- list(
      "prev_age_sex" = prev_age_sex,
      "rd_age" = rd_age,
      "rr_age" = rr_age,
      "or_age" = or_age,
      "or_age_55" = or_age_55,
      "or_age_70" = or_age_70,
      "rd_sex" = rd_sex,
      "rr_sex" = rr_sex,
      "or_sex" = or_sex,
      "rd_virus" = rd_virus,
      "rr_virus" = rr_virus,
      "or_virus" = or_virus,
      "prop_virus" = prop_virus,
      "rd_district" = rd_district,
      "rr_district" = rr_district,
      "or_district" = or_district,
      "rd_time" = rd_time,
      "rr_time" = rr_time,
      "or_time" = or_time,
      "count_age_sex" = count_age_sex,
      "count_age_virus" = count_age_virus
    )
    
    
    saveRDS(age_result, 'pcr_hsv/clean_data/age_result.RDS')
  }
}

# AME sex------------------------------------------------------------------
passsex <- T
passtime <- T
passdistrict <- F
passage <- T
Nchoice <- 2000
if (passsex == F) {
  ame_result <- list()
  ame_result <- ame_function(
    fit,
    Nsample = Nchoice,
    type_ame = "sexe",
    ref_ame = "Males",
    ref = 0,
    new = 1
  )
  gc()
  saveRDS(ame_result, 'pcr_hsv/clean_data/ame_result_ref_males.RDS')
}


# AME time------------------------------------------------------------------
if (passtime == F) {
  ame_result <- list()
  ame_result <- ame_function(
    fit,
    Nsample = Nchoice,
    type_ame = "time",
    ref_ame = "1st week of 2022",
    ref = 1,
    new = max(dfbinomial$idtime)
  )
  gc()
  saveRDS(ame_result,
          'pcr_hsv/clean_data/ame_result_ref_1stweek_studyperiod.RDS')
}


# AME district------------------------------------------------------------------
if (passdistrict == F) {
  ###May require several run if lack of RAM
  if (!file.exists('pcr_hsv/clean_data/ame_result_ref_paris.RDS')) {
    ame_result <- list()
  } else{
    ame_result <- readRDS('pcr_hsv/clean_data/ame_result_ref_paris.RDS')
  }
  
  for (i in spatial_objects$geometry_dep$idspace) {
    d <- spatial_objects$geometry_dep %>%
      filter(idspace == i) %>%
      .$nom
    if (i > length(ame_result)) {
      ame_result[[i]] <- ame_function(
        fit,
        Nsample = Nchoice,
        type_ame = "district",
        ref_ame = "Paris",
        ref = spatial_objects$geometry_dep %>% filter(nom == "Paris") %>% .$idspace,
        new = i
      )
      gc()
      saveRDS(ame_result,
              'pcr_hsv/clean_data/ame_result_ref_paris.RDS')
    }
  }
}

# AME age------------------------------------------------------------------
if (passage == F) {
  ###May require several run if lack of RAM
  if (!file.exists('pcr_hsv/clean_data/ame_result_ref_0.RDS')) {
    ame_result <- list()
  } else{
    ame_result <- readRDS('pcr_hsv/clean_data/ame_result_ref_0.RDS')
  }
  
  for (i in seq(1, 101)) {
    if (i > length(ame_result)) {
      ame_result[[i]] <- ame_function(
        fit,
        Nsample = Nchoice,
        type_ame = "age",
        ref_ame = "age 0",
        ref = 0,
        new = i - 1
      )
      gc()
      saveRDS(ame_result, 'pcr_hsv/clean_data/ame_result_ref_0.RDS')
    }
  }
}
