print("Cleaning ent")
rm(list = ls())
hpc <- FALSE
print(getwd())

# Function for detaching packages -----------------------------------------
detachAllPackages <- function() {
  basic.packages <- c(
    "package:stats",
    "package:graphics",
    "package:grDevices",
    "package:utils",
    "package:datasets",
    "package:methods",
    "package:base"
  )
  package.list <- search()[ifelse(unlist(gregexpr("package:", search())) == 1, TRUE, FALSE)]
  
  package.list <- setdiff(package.list, basic.packages)
  
  if (length(package.list) > 0)
    for (package in package.list)
      detach(package, character.only = TRUE)
}


# Detaching packages ------------------------------------------------------
detachAllPackages()

# Packages ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
required_packages <- c(
  "tidybayes",
  "MatrixModels",
  "Matrix",
  "readr",
  "spdep",
  "sf",
  "plyr",
  "INLA",
  "inlabru",
  "rSPDE",
  "tidylog",
  "tidyverse",
  "INLAspacetime"
)


extracting_model_part <- function(f) {
  list(
    "summary.fixed" = f$summary.fixed,
    "summary.random" = f$summary.random,
    "summary.hyperpar" = f$summary.hyperpar,
    'mode' = f$mode$theta
  )
}


# Function ----------------------------------------------------------------
function_group <- function(input) {
  #NB: Type 4 and 5 not consider b/c assumptions are stronger than type 1/2/3
  if (str_detect(input, "Type")) {
    if (str_detect(input, "Type1")) {
      out <- list(model = "rw1")
    } else if (str_detect(input, "Type2")) {
      out <- list(model = "rw2")
    } else if (str_detect(input, "Type3")) {
      out <- list(model = "ar1", hyper = list(rho = list(
        prior = "pc.cor1", param = c(0.75, 0.99)
      )))
    }
  } else{
    input <- case_when(
      str_detect(input, "220") ~ "220",
      str_detect(input, "202") ~ "202",
      str_detect(input, "121") ~ "121",
      str_detect(input, "102") ~ "102"
    )
    print(input)
    
    out <- stModel.define(
      smesh = mesh_space,
      tmesh = tmesh1,
      model = input,
      control.priors = list(
        prs = c(round(0.25 * map.size), 0.99),
        prt = c(round(0.25 * diff(
          range(vector_date)
        )), 0.99),
        psigma = c(1, 0.01)
      ),
      constr = TRUE
    )
  }
  return(out)
}

smoothedpmax <- function(x, y, type = 1) {
  if (type == 1) {
    #LogSumExp (Softmax)
    x <- 100^(-1) * log(exp(100 * x) + exp(100 * y))
  } else if (type == 2) {
    x <- 0.5 * (x + y + sqrt((x - y)^2 + 0.0000001))
  }
}


# LOAD ALL PACKAGES
lapply(required_packages, library, character.only = TRUE)
set.seed(123)
cpudt <- 15
CPUX <- "15:1"
inla.pardiso.check()

# Run on HPC or not ? -----------------------------------------------------
# This stuff was needed because one need to split the save between several partition on the HPC
# This is not needed if you can save ~600GB on the same place
if (str_detect(getwd(), "/projects/medium/phdthesis/phd_thesis")) {
  hpc <- TRUE
  path1hpc <- "/projects/medium/phdthesis/phd_thesis/"
  path2hpc <- "/users/supplisson/data/phd_thesis/"
}
print(hpc)
# expand.grid(x = seq(0, 1, 0.001), y = seq(0, 1, 0.001)) %>%
#   rowwise() %>%
#   mutate(pmax = pmax(x, y), smoothedpmax = smoothedpmax(x, y)) %>%
#   ggplot() +
#   geom_point(aes(x = pmax, y = smoothedpmax),
#             color = "red",
#             linetype = "dashed") +
#   geom_line(aes(x = x, y = x), color = "black", linetype = "solid") +
#   theme_minimal()

sd_to_prec <- function(sigma) {
  tibble(
    "sd" = sigma,
    "var" = round(sigma^2, 2),
    "prec" = round(1 / sigma^2, 2)
  )
}

# Function for summary  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
summary_draw <- function(df) {
  df %>%
    summarise(
      mean = mean(value, na.rm = T),
      median = as.double(median(value, na.rm = T)),
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
      Nrow = n(),
      Nrowbis = max(row_number()),
      Ndraws = sum(!is.na(value), na.rm = T),
      checkNA = sum(is.na(value)),
      N_higher_0 = sum(value > 0, na.rm = T),
      N_higher_05 = sum(value > 0.5, na.rm = T),
      N_higher_1 = sum(value > 1, na.rm = T),
      N_lower_0 = sum(value < 0, na.rm = T),
      N_lower_05 = sum(value < 0.5, na.rm = T),
      N_lower_1 = sum(value < 1, na.rm = T)
    ) %>%
    ungroup() %>%
    mutate(
      higher_0 = N_higher_0 / Ndraws,
      higher_05 = N_higher_05 / Ndraws,
      higher_1 = N_higher_1 / Ndraws,
      lower_0 = N_lower_0 / Ndraws,
      lower_05 = N_lower_05 / Ndraws,
      lower_1 = N_lower_1 / Ndraws
    )
}

source("pcr_hsv/codes/99_functions_for_ame.R", echo = T)


# Data --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load("pcr_hsv/clean_data/input_models/data_for_fit.rda")


# Path --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
path_data <- "pcr_hsv/clean_data/"
path_fit <- paste0(path_data, "inla_fit/")



# Modifying variables -----------------------------------------------------
aggregated <- aggregated %>%
  mutate(intercept = 1,
         sexe = factor(sexe),
         year = factor(year))

table(aggregated$year)


# Pivoting ----------------------------------------------------------------
dfbinomial <- aggregated   %>%
  pivot_longer(c(hsv1, hsv2), names_to = "virus", values_to = "result") %>%
  mutate(hsv1 = ifelse(virus == "hsv1", 1, 0),
         hsv2 = ifelse(virus == "hsv2", 1, 0)) %>%
  mutate(
    levelintercept = case_when(
      hsv2 == 0 & females == 0 ~ 1,
      hsv2 == 0 & females == 1 ~ 2,
      hsv2 == 1 & females == 0 ~ 3,
      hsv2 == 1 & females == 1 ~ 4
    )
  )
rm("aggregated")

# Neighborhood matrix ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
list_w <-
  c(names(spatial_objects$nb),
    names(spatial_objects$nb$nb_list))



sd_to_prec <- function(sigma) {
  tibble("sd" = sigma,
         "var" = sigma^2,
         "prec" = 1 / sigma^2)
}

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
      style = "B",
      zero.policy = TRUE
    ),
    "Matrix")
  } else{
    W <- as(nb2mat(
      spatial_objects$nb$nb_list[[w]],
      style = "B",
      zero.policy = TRUE
    ),
    "Matrix")
  }
  W
}

vector_age <- seq(0, 100)
vector_date <- dfdate %>% .$idtime
vector_date <- vector_date - (min(vector_date) - 1)
vector_space <- sort(unique(dfbinomial$idspace))
print(vector_date)

# Mesh space --------------------------------------------------------------
kmproj <-  st_crs(
  "+proj=lcc +lat_0=46.5 +lon_0=3 +lat_1=49 +lat_2=44 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs"
)
geometry_dep <- readRDS("pcr_hsv/clean_data/input_models/map_fr.RDS") %>% st_transform(kmproj)

bb <- matrix(st_bbox(geometry_dep), 2)

dbb <- apply(bb, 1, diff)
dbb

map.size <- mean(dbb)
map.size

map <- st_buffer(x = st_union(
  x = st_simplify(
    x = geometry_dep %>% summarise,
    preserveTopology = T,
    dTolerance = 1
  )
), dist = 0.1) %>% st_as_sf

if (!file.exists("pcr_hsv/clean_data/spatial_mesh.RDS")) {
  mesh_space <- fm_mesh_2d(
    offset = c(100, 100),
    max.edge = c(100, 200),
    cutoff = 10,
    min.angle = c(30, 21),
    loc = geometry_dep %>% st_centroid(),
    loc.domain	= map,
    boundary = map
  )
  
  #ggplot()+gg(mesh_space)+gg(map, alpha = 0.1)+gg(geometry_dep %>% st_centroid(), color = "red")
  mesh_space$n
  saveRDS(mesh_space, "pcr_hsv/clean_data/spatial_mesh.RDS")
} else{
  mesh_space <- readRDS("pcr_hsv/clean_data/spatial_mesh.RDS")
}
#NB: it is important the the last period has a knot because otherwise the variance is way too high
#as it will be computed as the barycentric average between something observed and something not observed but very far from the last knots with observations
min_date <- vector_date %>% min()
print(min_date)
max_date <- vector_date %>% max()
print(max_date)
interval <- c(min_date, max_date)
dt <- 5
knots <- seq(interval[1], interval[2], by = dt)
print(knots)
tmesh2 <- fm_mesh_1d(knots, degree = 2, boundary = "free")
tmesh1 <- fm_mesh_1d(seq(interval[1] - dt / 2, interval[2] + dt / 2, by = dt), degree = 1)
tmesh2$n * mesh_space$n
dfbinomial <- dfbinomial %>%
  left_join(geometry_dep %>% st_centroid() %>% select(insee_dep, geometry),
            by = "insee_dep")

dfpoisson <-  dfbinomial %>%
  distinct(age, females, idspace, idtime, .keep_all = T) %>%
  st_drop_geometry()


# Model set ---------------------------------------------------------------
matrices <- c("queen")

w_set <- expand.grid(whsv1 = matrices, whsv2 = matrices) %>%
  mutate_at(vars(whsv1, whsv2), ~ factor(., levels = matrices)) %>%
  mutate(link = "logit",
         type = "Type1",
         pathwd = "pathwd1")

w_set <- w_set %>%
  rbind(w_set %>% mutate(type = "Type2")) %>%
  rbind(w_set %>% mutate(type = "Type3"))
print(w_set)

w_set <- tibble(
  whsv1 = rep("queen", 4),
  whsv2 = rep("queen", 4),
  link = rep("logit", 4),
  type = c("202", "102", "220", "121"),
  pathwd = rep("pathwd1", 4)
) %>%
  add_row(w_set)


w_set <- w_set %>%
  rbind(w_set %>% filter(type == "Type3") %>%  mutate(link = 'cauchit', pathwd = "pathwd1")) %>%
  rbind(w_set %>% filter(type == "Type3") %>%  mutate(link = 'cloglog', pathwd = "pathwd1")) %>%
  rbind(w_set %>% filter(type == "Type3") %>%  mutate(link = 'probit', pathwd = "pathwd2")) %>%
  rbind(w_set %>% filter(type == "Type3") %>%  mutate(link = 'robit_3', pathwd = "pathwd2")) %>%
  rbind(w_set %>% filter(type == "Type3") %>%  mutate(link = 'robit_4', pathwd = "pathwd2")) %>%
  rbind(w_set %>% filter(type == "Type3") %>%  mutate(link = 'robit_5', pathwd = "pathwd2")) %>%
  rbind(w_set %>% filter(type == "Type3") %>%  mutate(link = 'robit_6', pathwd = "pathwd2"))

gc()

print(w_set)

tmp <- w_set %>%  filter(link == "robit_3")

w_set <- w_set %>%
  rbind(tmp %>%  mutate(type = 'Type3Specialbetaweightsmean', pathwd = "pathwd1")) 

#not there yet
  # rbind(tmp %>%  mutate(type = 'Type3Specialcorrelatedscalingmean', pathwd = "pathwd1")) %>% 
  # rbind(tmp %>%  mutate(type = 'Type3Specialbetaweightsmax', pathwd = "pathwd1")) %>%
  # rbind(tmp %>%  mutate(type = 'Type3Specialcorrelatedscalingmax', pathwd = "pathwd1"))

print(w_set)
rm("tmp")
print(w_set)


test_hpc <- function(hpc_test) {
  if (hpc_test == T) {
    setwd(path1hpc)
    print(getwd())
  } else{
    print(getwd())
  }
}



# Fit ---------------------------------------------------------------------
passfit <- F
if (passfit == F) {
  for (x in 1:(nrow(w_set))) {
    print(x)
    
    model <- w_set %>%
      filter(row_number() == x) %>%
      mutate(whsv1 = as.character(whsv1), whsv2 = as.character(whsv2))
    
    model <- w_set %>%
      filter(row_number() == x) %>%
      mutate(whsv1 = as.character(whsv1), whsv2 = as.character(whsv2))
    
    print(model)
    
    name_fit <-
      paste0("fit_",
             paste(
               ifelse(str_detect(model$type, "Type"), model$type, ""),
               model$link,
               model$whsv1,
               model$whsv2,
               sep = "_"
             )) %>%
      str_replace(., "__", "_")
    
    
    if (!str_detect(model$type, "Type")) {
      name_fit <- paste0(name_fit, "_spatiotemporal_", model$type)
    }
    
    print(name_fit)
    
    
    if (hpc == T & model$pathwd == "pathwd2") {
      setwd(path2hpc)
      print(getwd())
    } else if (hpc == T & model$pathwd == "pathwd1") {
      setwd(path1hpc)
      print(getwd())
    } else{
      
    }
    
    path_fit <- "pcr_hsv/clean_data/inla_fit/"
    list_fit <- base::list.files(
      path_fit,
      pattern = ".RDS",
      all.files = T,
      full.names = FALSE
    )
    
    save_path <- paste0(path_fit, name_fit, ".RDS")
    init_path <- str_replace(save_path, "fit\\_", "summary\\_")
    path_info <- str_replace(save_path, "/fit", "/info")
    path_prior <- str_replace(save_path, "/fit", "/prior")
    print(save_path)
    print(init_path)
    print(path_info)
    print(path_prior)
    if (str_detect(init_path, "Specialmax")) {
      init_path <- str_replace(init_path, "Specialmax")
      init_path <- str_remove(init_path, "Specialmean")
    }
    if (!file.exists(save_path)) {
      if (model$link == "powerlogit") {
        init_path <- str_replace(init_path, "powerlogit", "logit")
        print(init_path)
      }
      if (model$link == "skewnormal") {
        init_path <- str_replace(init_path, "skewnormal", "probit")
        print(init_path)
      }
      if (file.exists(init_path)) {
        init <- readRDS(init_path)
        init <- init$mode
        print(init)
      } else{
        init <- NULL
      }
      
      if (str_detect(model$link, "powerlogit")) {
        init <- c(0, 0, init)
      }
      if (str_detect(model$link, "skewnormal")) {
        init <- c(0, 0, init)
      }
      print(init)
      
      #Make sure its 0 for the special case
      print("Fitting model")
      print(name_fit)
      print(init)
      
      Whsv1 <- return_matrix(model$whsv1)
      Whsv2 <- return_matrix(model$whsv2)
      type <- model$type
      stmodel <- function_group(model$type)
      
      if (str_detect(save_path, "powerlogit") |
          str_detect(save_path, "skewnormal")) {
        cmp <- ~ 0 + baseline(
          main = levelintercept,
          model = "iidkd",
          constr = TRUE,
          order = 4,
          n = 4
        )
      } else{
        cmp <- ~ 0 + baseline(
          main = levelintercept,
          model = "iidkd",
          constr = FALSE,
          order = 4,
          n = 4
        )
      }
      
      
      cmp <- update(
        cmp,
        . ~ . +
          agecommon(
            main = age,
            group = females + 1,
            values = vector_age,
            ngroup = 2,
            model = "rw2",
            constr = T,
            scale.model = T,
            hyper = list(prec = list(
              prior = "pc.prec", param = c(1, 0.01)
            )),
            control.group = list(model = "exchangeable", hyper = list(rho = list(
              prior = "normal", param = c(0, 0.2)
            )))
          ) +
          agehsv2(
            main = age,
            weights = hsv2,
            group = females + 1,
            values = vector_age,
            ngroup = 2,
            model = "rw2",
            constr = T,
            scale.model = T,
            hyper = list(prec = list(
              prior = "pc.prec", param = c(1, 0.01)
            )),
            control.group = list(model = "exchangeable", hyper = list(rho = list(
              prior = "normal", param = c(0, 0.2)
            )))
          )
      )
      
      
      if (str_detect(model$type, "Type")) {
        cmp <- update(
          cmp,
          . ~ . +
            spacetimecommon(
              main = idspace,
              group = idtime,
              values = vector_space ,
              ngroup = length(vector_date),
              model = "bym2",
              graph = Whsv1,
              scale.model = TRUE,
              adjust.for.con.comp = TRUE,
              constr = TRUE,
              hyper = list(
                prec = list(prior = "pc.prec", param = c(u = 1, a = 0.01)),
                phi = list(prior = "pc", param =  c(u = 0.5, a = 0.5))
              ),
              control.group = stmodel
            ) +
            spacetimehsv2(
              main = idspace,
              weights = hsv2,
              group = idtime,
              values = vector_space,
              ngroup = length(vector_date),
              model = "bym2",
              graph = Whsv2,
              scale.model = TRUE,
              adjust.for.con.comp = TRUE,
              constr = TRUE,
              hyper = list(
                prec = list(prior = "pc.prec", param = c(u = 1, a = 0.01)),
                phi = list(prior = "pc", param =  c(u = 0.5, a = 0.5))
              ),
              control.group = stmodel
            )
        )
      } else{
        cmp <- update(
          cmp,
          . ~ . +
            spacetimecommon(
              main = list(space = geometry, time = idtime),
              model = stmodel,
              mapper = bru_mapper_multi(
                list(
                  space = bru_mapper_fmesher(mesh_space),
                  time = bru_mapper_fmesher(tmesh2)
                )
              )
            ) +
            spacetimehsv2(
              main = list(space = geometry, time = idtime),
              weights = hsv2,
              model = stmodel,
              mapper = bru_mapper_multi(
                list(
                  space = bru_mapper_fmesher(mesh_space),
                  time = bru_mapper_fmesher(tmesh2)
                )
              )
            )
        )
      }
      
      if (str_detect(model$type, "Special")) {
        cmp <- update(
          cmp,
          . ~ . +
            baselinepoisson(
              main = females + 1,
              model = "iidkd",
              constr = FALSE,
              order = 2,
              n = 2
            ) +
            agepoisson(
              main = age,
              group = females + 1,
              values = vector_age,
              ngroup = 2,
              model = "rw2",
              constr = T,
              scale.model = T,
              hyper = list(prec = list(
                prior = "pc.prec", param = c(1, 0.01)
              )),
              control.group = list(model = "exchangeable", hyper = list(rho = list(
                prior = "normal", param = c(0, 0.2)
              )))
            ) +
            spacetimepoisson(
              main = idspace,
              group = idtime,
              values = vector_space,
              ngroup = length(vector_date),
              model = "bym2",
              graph = Whsv1,
              scale.model = TRUE,
              adjust.for.con.comp = TRUE,
              constr = TRUE,
              hyper = list(
                prec = list(prior = "pc.prec", param = c(u = 1, a = 0.01)),
                phi = list(prior = "pc", param =  c(u = 0.5, a = 0.5))
              ),
              control.group = stmodel
            )
        )
      }
      
      link_list <- list(control.link = list(model = "logit"))
      cmp_to_use <- bru_used(
        labels = c(
          "baseline",
          "agecommon",
          "agehsv2",
          "spacetimecommon",
          "spacetimehsv2"
        )
      )
      
      formula_input <-  result ~ .
      
      if (model$link == "logit") {
        
      } else if (model$link == "powerlogit") {
        formula_input <-  result ~ .
        link_list <- list(control.link = list(
          model = "powerlogit",
          hyper = list(
            intercept = list(
              initial = 0,
              prior = "logitbeta",
              param = c(1, 1)
            ),
            power = list(
              initial = 0,
              param = c(0, 1),
              prior = "normal",
              fixed = FALSE
            )
          )
        ))
      } else if (model$link == "skewnormal") {
        formula_input <-  result ~ .
        link_list <- list(control.link = list(
          model = "sn",
          hyper = list(
            skew = list(prior = "pc.sn", param = 10),
            intercept = list(prior = "linksnintercept", param = c(0, 0))
          )
        ))
        
      } else if (model$link %in% c("probit", "cauchit", "cloglog", "loglog")) {
        link_list <- list(control.link = list(model = model$link))
      } else if (str_detect(model$link, "robit")) {
        link_list <- list(control.link = list(model = "robit", hyper = list(dof = list(
          initial = log(as.numeric(
            str_remove(model$link, "robit\\_")
          ) - 2), fixed = TRUE
        ))))
      }
      
      
      lik <- bru_obs(
        data = dfbinomial,
        Ntrials = N,
        control.family = link_list,
        family = "binomial",
        formula = formula_input,
        used = cmp_to_use
      )
      
      print(cmp_to_use)
      print(formula_input)
      print(link_list)
      print(cmp)
      print(lik)
      
      if (str_detect(model$type, "Special")) {
        df_robit_tmp <-  as.numeric(str_remove(model$link, "robit\\_"))
        if (str_detect(model$type, "betaweights")) {
          cmp <- update(
            cmp,
            . ~ .   +
              scale(
                1,
                model = "linear",
                prec.linear = 1,
                mean.linear = 0,
                marginal = bru_mapper_marginal(
                  qfun = qexp,
                  pfun = pexp,
                  dfun = dexp,
                  rate = -log(0.01) / 1
                )
              )
          )
          
          
          
          fml_poisson <- N ~ baselinepoisson + agepoisson + spacetimepoisson + scale * smoothedpmax(
            INLA::inla.link.invrobit(
              baseline_eval(main = 1 +  .data.[["females"]]) +
                agecommon_eval(main = .data.[["age"]], group =  .data.[["females"]] + 1) +
                spacetimecommon_eval(main = list(u = .data.[["idspace"]]), group = .data.[["idtime"]]),
              df = df_robit_tmp
            ) +
              INLA::inla.link.invrobit(
                baseline_eval(main = 3 +  .data.[["females"]]) +
                  agecommon_eval(main = .data.[["age"]], group =  .data.[["females"]] + 1) +
                  spacetimecommon_eval(main = list(u = .data.[["idspace"]]), group = .data.[["idtime"]]) +
                  agehsv2_eval(main = .data.[["age"]], group =  .data.[["females"]] + 1) +
                  spacetimehsv2_eval(main = list(u =  .data.[["idspace"]]), group = .data.[["idtime"]]),
                df = df_robit_tmp
              ),
            type = 2
          )
          
          if (str_detect(model$type, "mean")) {
            cmp <- update(
              cmp,
              . ~ .  +
                weightssum(
                  1,
                  model = "linear",
                  prec.linear = 1,
                  mean.linear = 0,
                  marginal = bru_mapper_marginal(
                    qfun = qbeta,
                    pfun = pbeta,
                    dfun = dbeta,
                    shape1 = 1,
                    shape2 = 1
                  )
                )
            )
            
            
            fml_poisson <- N ~ baselinepoisson + agepoisson + spacetimepoisson + scale * (
              weightssum * INLA::inla.link.invrobit(
                baseline_eval(main = 1 +  .data.[["females"]]) +
                  agecommon_eval(main = .data.[["age"]], group =  .data.[["females"]] + 1) +
                  spacetimecommon_eval(main = list(u = .data.[["idspace"]]), group = .data.[["idtime"]]),
                df = df_robit_tmp
              ) + (1 - weightssum) *
                INLA::inla.link.invrobit(
                  baseline_eval(main = 3 +  .data.[["females"]]) +
                    agecommon_eval(main = .data.[["age"]], group =  .data.[["females"]] + 1) +
                    spacetimecommon_eval(main = list(u = .data.[["idspace"]]), group = .data.[["idtime"]]) +
                    agehsv2_eval(main = .data.[["age"]], group =  .data.[["females"]] + 1) +
                    spacetimehsv2_eval(main = list(u =  .data.[["idspace"]]), group = .data.[["idtime"]]),
                  df = df_robit_tmp
                )
            )
          }
        }
        else if (str_detect(model$type, "correlatedscaling")) {
          cmp <- update(
            cmp,
            . ~ .   +
              scaling(
                main = females + 1,
                model = "iidkd",
                constr = FALSE,
                order = 2,
                n = 2
              )
          )
          
          fml_poisson <- N ~ baselinepoisson + agepoisson + spacetimepoisson + exp(scaling) * smoothedpmax(
            x = INLA::inla.link.invrobit(
              baseline_eval(main = 1 +  .data.[["females"]]) +
                agecommon_eval(main = .data.[["age"]], group =  .data.[["females"]] + 1) +
                spacetimecommon_eval(main = list(u = .data.[["idspace"]]), group = .data.[["idtime"]]),
              df = df_robit_tmp
            ),
            y = INLA::inla.link.invrobit(
              baseline_eval(main = 3 +  .data.[["females"]]) +
                agecommon_eval(main = .data.[["age"]], group =  .data.[["females"]] + 1) +
                spacetimecommon_eval(main = list(u = .data.[["idspace"]]), group = .data.[["idtime"]]) +
                agehsv2_eval(main = .data.[["age"]], group =  .data.[["females"]] + 1) +
                spacetimehsv2_eval(main = list(u =  .data.[["idspace"]]), group = .data.[["idtime"]]),
              df = df_robit_tmp
            ),
            type = 2
          )
          if (str_detect(model$type, "mean")) {
            cmp <- update(
              cmp,
              . ~ .   +
                weightssum(
                  1,
                  model = "linear",
                  prec.linear = 1,
                  mean.linear = 0,
                  marginal = bru_mapper_marginal(
                    qfun = qbeta,
                    pfun = pbeta,
                    dfun = dbeta,
                    shape1 = 1,
                    shape2 = 1
                  )
                )
            )
            
            fml_poisson <- N ~ baselinepoisson + agepoisson + spacetimepoisson + exp(scaling) * (
              weightssum * INLA::inla.link.invrobit(
                baseline_eval(main = 1 +  .data.[["females"]]) +
                  agecommon_eval(main = .data.[["age"]], group =  .data.[["females"]] + 1) +
                  spacetimecommon_eval(main = list(u = .data.[["idspace"]]), group = .data.[["idtime"]]),
                df = df_robit_tmp
              ) + (1 - weightssum) *
                INLA::inla.link.invrobit(
                  baseline_eval(main = 3 +  .data.[["females"]]) +
                    agecommon_eval(main = .data.[["age"]], group =  .data.[["females"]] + 1) +
                    spacetimecommon_eval(main = list(u = .data.[["idspace"]]), group = .data.[["idtime"]]) +
                    agehsv2_eval(main = .data.[["age"]], group =  .data.[["females"]] + 1) +
                    spacetimehsv2_eval(main = list(u =  .data.[["idspace"]]), group = .data.[["idtime"]]),
                  df = df_robit_tmp
                )
            )
          }
        }
      }
      
      
      
      opt <- bru_options(
        bru_verbose = 4,
        verbose = T,
        bru_max_iter = 100,
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
          restart = 0,
          control.vb = list(strategy = "variance")
        ),
        control.mode = list(theta = init, restart = TRUE),
        num.threads = CPUX,
        safe = TRUE
      )
      
      print(opt$control.inla$strategy)
      print(opt$control.inla$int.strategy)
      
      if (str_detect(model$type, "Special")) {
        liknzpoisson <- bru_obs(data = dfpoisson,
                                family = "nzPoisson",
                                formula = fml_poisson)
        print(fml_poisson)
        print(liknzpoisson)
        
        fit <- bru(components = cmp,
                   lik,
                   liknzpoisson,
                   options = opt)
        gc()
      }
      else{
        print("Fit regular model")
        fit <- bru(components = cmp, lik, options = opt)
      }
      
      saveRDS(extracting_model_part(fit), init_path)
      saveRDS(fit$bru_info$model, path_info)
      saveRDS(
        list(
          "type" = model$type,
          "Whsv1" = Whsv1,
          "Whsv2" = Whsv2,
          "stmodel" = stmodel
        ),
        path_prior
      )
      saveRDS(fit, save_path)
      
      rm("fit")
      rm("init")
      gc()
      if (hpc == T) {
        setwd(path1hpc)
      }
    }
    gc()
  }
}

print("Computing GCPO")
if (hpc == T) {
  pset <- c(path1hpc, path2hpc)
} else{
  pset <- getwd()
}

for (l in pset) {
  print(l)
  setwd(l)
  list_fit <- list.files(path_fit, full.names = T, pattern = "fit_")
  list_fit <- list_fit[which(str_detect(list_fit, "220"))]
  print(list_fit)
  for (p in list_fit) {
    print(p)
    gc()
    if (hpc == T) {
      setwd(path1hpc)
    }
    if (file.exists(str_replace(p, "fit_", "gcpo_"))) {
      print("Import GCPO list")
      gcpo_list <-  readRDS(str_replace(p, "fit_", "gcpo_"))
    } else{
      print("Creating GCPO list")
      gcpo_list <- list()
    }
    
    rm("fit")
    rm("Whsv1")
    rm("Whsv2")
    rm("stmodel")
    rm("prior")
    for (x in c(-1, 5, 10, 15, 20, 25, 30)) {
      print(x)
      if (x > 0) {
        namex <- paste0("lgocv.m", x)
      } else{
        namex <- "loo"
      }
      if (!rlang::has_name(gcpo_list, namex)) {
        if (!exists("fit")) {
          setwd(l)
          print(p)
          fit <- readRDS(p)
          prior <- readRDS(str_replace(p, "fit_", "prior_"))
          list2env(prior, .GlobalEnv)
        }
        
        print("Performing LGOCV")
        
        gcpo_list[[namex]] <- inla.group.cv(
          fit,
          num.level.sets = x,
          strategy = "posterior",
          size.max = 32
        )
        print("Save")
        if (hpc == T) {
          setwd(path1hpc)
        }
        saveRDS(gcpo_list, str_replace(p, "fit_", "gcpo_"))
      }
      else{
        print("Already available")
      }
    }
  }
  gc()
}

stop()

# Exporting information about version -------------------------------------
for (l in pset) {
  print(l)
  setwd(l)
  list_fit <- list.files(path_fit, full.names = T, pattern = "fit_")
  print(list_fit)
  for (p in list_fit) {
    print(p)
    gc()
    path_version <- str_replace(p, "/fit_", "/version_")
    print(path_version)
    if (!file.exists(path_version)) {
      fit <- readRDS(p)
      version <- tibble(
        "bru_version" = fit$bru_info$inlabru_version,
        "inla_version" = fit$bru_info$INLA_version
      )
      print(version)
      saveRDS(version, path_version)
    }
  }
  gc()
}

if (hpc == T) {
  setwd(path1hpc)
}


# Assessment metrics ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if (!file.exists("pcr_hsv/clean_data/assessment_metric.RDS")) {
  list_fit <- base::list.files(
    path_fit,
    pattern = "gcpo",
    all.files = T,
    full.names = T
  ) %>%
    tibble(x = .)
  
  print(list_fit)
  
  summary <- list_fit %>%
    split(~ x) %>%
    imap(function(.x, .y) {
      print(.x$x)
      model <- readRDS(.x$x)
      
      metric <- tibble(
        fitname = str_replace(.x$x, "gcpo", "fit"),
        loocv = -mean(log(model$loo$cv[seq(1, nrow(dfbinomial))])),
        lgocv.m5 = -mean(log(model$lgocv.m5$cv[seq(1, nrow(dfbinomial))])),
        lgocv.m10 = -mean(log(model$lgocv.m10$cv[seq(1, nrow(dfbinomial))])),
        lgocv.m15 = -mean(log(model$lgocv.m15$cv[seq(1, nrow(dfbinomial))])),
        lgocv.m20 = -mean(log(model$lgocv.m20$cv[seq(1, nrow(dfbinomial))])),
        lgocv.m25 = -mean(log(model$lgocv.m25$cv[seq(1, nrow(dfbinomial))])),
        lgocv.m30 = -mean(log(model$lgocv.m30$cv[seq(1, nrow(dfbinomial))]))
      )
      
      rm("model")
      rm("fit")
      gc()
      
      return(metric)
    })
  
  gc(full = T)
  assessment_metric <- summary %>%
    imap(function(.x, .y) {
      .x
    }) %>%
    do.call(rbind, .)
  
  #https:/github.com/inlabru-org/inlabru/discussions/144
  
  
  print(assessment_metric)
  assessment_metric %>% filter(lgocv.m30 == min(lgocv.m30)) %>% .$fitname
  assessment_metric %>% arrange(lgocv.m30) %>% .$fitname
  #Saving assessment
  saveRDS(assessment_metric,
          "pcr_hsv/clean_data/assessment_metric.RDS")
} else{
  assessment_metric <- readRDS("pcr_hsv/clean_data/assessment_metric.RDS")
}



if (!file.exists("pcr_hsv/clean_data/list_gcpo.RDS")) {
  #Importing df of GCPO
  list_gcpo <- base::list.files(
    path_fit,
    pattern = "gcpo",
    all.files = T,
    full.names = T
  ) %>%
    tibble(x = .) %>%
    split( ~ x) %>%
    imap(function(.x, .y) {
      readRDS(.x$x)
    })
  
  list_gcpo <- list_gcpo %>%
    imap(function(.x, .y) {
      m <- .y
      .x %>%
        imap(function(.x, .y) {
          .x$cv[seq(1, nrow(dfbinomial))] %>%
            tibble('cv' = ., 'name' = .y) %>%
            mutate(id = row_number())
        }) %>%
        do.call(rbind, .) %>%
        mutate(fitname = m)
    }) %>%
    do.call(rbind, .) %>%
    mutate(
      fitname_init = fitname,
      fitname = str_remove(fitname, "pcr_hsv/clean_data/inla_fit/gcpo_"),
      fitname = str_remove(fitname, "_queen_queen_spatiotemporal_"),
      fitname = str_remove(fitname, ".RDS"),
      spacetime = case_when(
        str_detect(fitname, "Type1") ~ "Spatiotemporal (RW(1))",
        str_detect(fitname, "Type2") ~ "Spatiotemporal (RW(2))",
        str_detect(fitname, "Type3") ~ "Spatiotemporal (AR(1))",
        str_detect(fitname, "Type4") ~ "Spatial-Temporal (RW(1))",
        str_detect(fitname, "Type5") ~ "Spatial-Temporal (RW(2))",
        str_detect(fitname, "Type6") ~ "Spatial-Temporal (AR(1))",
        TRUE ~ paste0("Spatiotemporal (DEMF(", fitname, "))")
      ),
      special = case_when(
        str_detect(fitname, "Specialcorrelatedweights") ~ "Weighted average (plogis)",
        str_detect(fitname, "Specialbetaweights") ~ "Weighted average (Beta)",
        TRUE ~ NA
      ),
      link = case_when(
        str_detect(fitname, "logit") &
          !str_detect(fitname, "power") ~ "Logit",
        str_detect(fitname, "probit") ~ "Probit",
        str_detect(fitname, "cauchit") ~ "Cauchit",
        str_detect(fitname, "powerlogit") ~ "Powerlogit",
        str_detect(fitname, "skewnormal") ~ "Skew-Gaussian",
        str_detect(fitname, "cloglog") ~ "Cloglog",
        TRUE ~ fitname
      ),
      name = ifelse(name == "loo", 1, as.numeric(str_remove(name, "lgocv.m"))),
      name = paste0("GCPO(", name, ")"),
      name = factor(name, levels = paste0("GCPO(", c(1, seq(
        5, 30, 5
      )), ')'))
    )
  
  
  saveRDS(list_gcpo, "pcr_hsv/clean_data/list_gcpo.RDS")
} else{
  list_gcpo <- readRDS("pcr_hsv/clean_data/list_gcpo.RDS")
}
if (!file.exists("pcr_hsv/clean_data/list_gcpo_with_id.RDS")) {
  list_gcpo <- list_gcpo %>%
    left_join(
      dfbinomial %>% st_drop_geometry() %>% select(-hsv1, -hsv2, -N) %>% mutate(id = row_number()),
      by = "id"
    )
  saveRDS(list_gcpo, "pcr_hsv/clean_data/list_gcpo_with_id.RDS")
}
# Stacking weights --------------------------------------------------------
function_weights <- function(gcpo_choice = "GCPO(30)") {
  #Matrix of cpo
  gcpo_matrix <- list_gcpo %>%
    filter(as.character(name) == gcpo_choice) %>%
    select(id, cv, fitname_init) %>%
    pivot_wider(names_from = "fitname_init", values_from = "cv") %>%
    select(-id) %>%
    as.matrix
  #Stacking weights
  stacking_weights <- loo::stacking_weights(gcpo_matrix)
  #Formating
  stacking_weights <- stacking_weights %>%
    as.data.frame %>%
    tibble::rownames_to_column(var = "model") %>%
    rename(weights = x) %>%
    mutate(index  = as.numeric(str_remove(model, "model"))) %>%
    mutate(weights = round(weights, 2))
  
  stacking_weights$weights <- stacking_weights$weights / sum(stacking_weights$weights)
  
  stacking_weights <- stacking_weights  %>%
    add_column(tibble("fit" = str_replace(colnames(gcpo_matrix), "/gcpo", "/fit")))
  
  #Output
  out <- list("gcpo_matrix" = gcpo_matrix, "stacking_weights" = stacking_weights)
  out
}
if (!file.exists("pcr_hsv/clean_data/weights.RDS")) {
  weights <- list()
  for (x in c(1, seq(5, 30, 5))) {
    x <- paste0("GCPO(", x, ")")
    print(x)
    weights[[x]] <- function_weights(gcpo_choice = x)
  }
  
  saveRDS(weights, "pcr_hsv/clean_data/weights.RDS")
} else{
  weights <- readRDS("pcr_hsv/clean_data/weights.RDS")
}
NdrawsTot <- 5000
NdrawsAME <- 3000
weights <- weights %>%
  imap( ~ .x %>% imap(~ if (.y == "stacking_weights") {
    .x %>%
      mutate(
        draws = weights * NdrawsTot,
        draws = round(draws, 0),
        draws = as.integer(draws),
        draws_ame = weights * NdrawsAME,
        draws_ame = round(draws_ame, 0),
        draws_ame = as.integer(draws_ame)
      )
  } else{
    .x
  }))


saveRDS(weights, "pcr_hsv/clean_data/weights.RDS")

max_draw <- weights %>%
  imap(~ .x$stacking_weights) %>% do.call(rbind, .)

passpp <- T
passsummary <- T
passsex <- T
passtime <- T
passage <- F
passdistrict <- T
passdistricttime <- T
set_gcp_analysis <- c(30)

print(set_gcp_analysis)
print("Start postfit analysis")
# ANALYSES ----------------------------------------------------------------
for (gcpo_values in set_gcp_analysis) {
  if (gcpo_values != set_gcp_analysis[1]) {
    passpp <- T
    passsummary <- T
  }
  print(gcpo_values)
  if (gcpo_values %in% as.numeric(str_remove(str_remove(names(weights), 'GCPO\\('), '\\)'))) {
    path_weights <- paste0("pcr_hsv/clean_data/weights_used_gcpo_",
                           gcpo_values,
                           ".RDS")
    print(path_weights)
    path_pp_draws <- paste0("pcr_hsv/clean_data/draws_pp_gcpo_", gcpo_values, ".RDS")
    print(path_pp_draws)
    path_pp <- paste0("pcr_hsv/clean_data/results_pp_gcpo_",
                      gcpo_values,
                      ".RDS")
    print(path_pp)
    path_draws_age <- paste0("pcr_hsv/clean_data/draws_age_gcpo_",
                             gcpo_values,
                             ".RDS")
    print(path_draws_age)
    path_age_results <- paste0("pcr_hsv/clean_data/age_result_gcpo_",
                               gcpo_values,
                               ".RDS")
    print(path_age_results)
    path_ame_sex <- paste0("pcr_hsv/clean_data/ame_result_sex_refmales_gcpo_",
                           gcpo_values,
                           ".RDS")
    print(path_ame_sex)
    path_ame_time <- paste0("pcr_hsv/clean_data/ame_result_time_ref1stweek_gcpo_",
                            gcpo_values,
                            ".RDS")
    print(path_ame_time)
    path_ame_age <- paste0("pcr_hsv/clean_data/ame_result_age_ref0year_gcpo_",
                           gcpo_values,
                           ".RDS")
    print(path_ame_age)
    path_ame_space <- paste0("pcr_hsv/clean_data/ame_result_space_refparis_gcpo_",
                             gcpo_values,
                             ".RDS")
    print(path_ame_space)
    path_ame_spacetime <- paste0(
      "pcr_hsv/clean_data/ame_result_spacetime_refdistricttime1_gcpo_",
      gcpo_values,
      ".RDS"
    )
    print(path_ame_spacetime)
    print("NDRAWS")
    # N Draws -----------------------------------------------------------------
    if (!file.exists(path_weights)) {
      set_model_to_stack <- weights[[paste0('GCPO(', gcpo_values, ')')]][["stacking_weights"]]
      set_model_to_stack <- set_model_to_stack  %>%
        select(fit, model, weights, index, draws, draws_ame)
      saveRDS(set_model_to_stack, path_weights)
    } else{
      print("Already availabe")
      set_model_to_stack <- readRDS(path_weights)
    }
    
    print(set_model_to_stack)
    formula_function <- function(path_input) {
      if (str_detect(path_input, "logit") &
          !str_detect(path_input, "powerlogit")) {
        fml <- ~ tibble(p = plogis(
          baseline + agecommon  + agehsv2 + spacetimecommon + spacetimehsv2
        ))
      } else if (str_detect(path_input, "powerlogit")) {
        #not used, not working
        #see https://inla.r-inla-download.org/r-inla.org/doc/link/powerlogit.pdf for required function
        fml <- ~ muimap(
          power.link(Link_power_logit_power)$cdf(
            Link_power_logit_intercept + baseline + agecommon + agehsv2 + spacetimecommon + spacetimehsv2
          )
        )
      } else if (str_detect(path_input, "probit")) {
        fml <- ~ tibble(p = pnorm(
          baseline + agecommon  + agehsv2 + spacetimecommon + spacetimehsv2
        ))
      } else if (str_detect(path_input, "skewnormal")) {
        fml <- ~ tibble(
          p = inla.link.invsn(
            x = baseline + agecommon  + agehsv2 + spacetimecommon + spacetimehsv2,
            skew = Link_sn_skew,
            intercept = Link_sn_intercept
          )
        )
      } else if (str_detect(path_input, "robit")) {
        dof <- case_when(
          str_detect(path_input, "robit_3") ~ 3,
          str_detect(path_input, "robit_4") ~ 4,
          str_detect(path_input, "robit_5") ~ 5,
          str_detect(path_input, "robit_6") ~ 6
        )
        assign("dof", dof, envir = .GlobalEnv)
        print(dof)
        fml <- ~ tibble(p = inla.link.invrobit(
          baseline + agecommon + agehsv2 + spacetimecommon + spacetimehsv2,
          df = dof
        ))
      } else if (str_detect(path_input, "loglog") &
                 !str_detect(path_input, "cloglog")) {
        fml <- ~ tibble(p = inla.link.invloglog(
          baseline + agecommon + agehsv2 + spacetimecommon + spacetimehsv2
        ))
      } else if (str_detect(path_input, "cloglog")) {
        fml <- ~ tibble(p = inla.link.invcloglog(
          baseline + agecommon + agehsv2 + spacetimecommon + spacetimehsv2
        ))
      } else if (str_detect(path_input, "cauchit")) {
        fml <- ~ tibble(p = inla.link.invcauchit(
          baseline + agecommon + agehsv2 + spacetimecommon + spacetimehsv2
        ))
      }
      assign("fml", fml, envir = .GlobalEnv)
      print(fml)
    }
    setting_paths <- function(i) {
      print(i)
      #Removing object
      ifrm <- function(obj, env = globalenv()) {
        obj <- deparse(substitute(obj))
        if (exists(obj, envir = env)) {
          rm(list = obj, envir = env)
        }
      }
      ifrm("fit")
      ifrm("path_fit")
      ifrm("path_state")
      ifrm("path_info")
      ifrm("path_prior")
      ifrm("ndraws")
      ifrm("ndrawsame")
      ifrm("fml")
      ifrm("stmodel")
      ifrm("Whsv1")
      ifrm("Whsv2")
      
      path_fit <- set_model_to_stack$fit[i]
      print(path_fit)
      path_state <- set_model_to_stack$state[i]
      print(path_state)
      path_info <- set_model_to_stack$info[i]
      print(path_info)
      path_prior <- set_model_to_stack$prior[i]
      print(path_prior)
      ndraws <- set_model_to_stack$draws[i]
      print(ndraws)
      ndrawsame <- set_model_to_stack$draws_ame[i]
      print(ndrawsame)
      
      wd <- set_model_to_stack$pathwd[i]
      print(wd)
      
      assign("path_fit", path_fit, envir = .GlobalEnv)
      assign("path_state", path_state, envir = .GlobalEnv)
      assign("path_info", path_info, envir = .GlobalEnv)
      assign("path_prior", path_prior, envir = .GlobalEnv)
      assign("ndraws", ndraws, envir = .GlobalEnv)
      assign("ndrawsame", ndrawsame, envir = .GlobalEnv)
      assign("wd", wd, envir = .GlobalEnv)
      
      formula_function(path_fit)
      
      if (hpc == T & wd == "pathwd2") {
        setwd(path2hpc)
        print(getwd())
      } else if (hpc == T & wd == "pathwd1") {
        setwd(path1hpc)
        print(getwd())
      } else{
        print(getwd())
      }
      
      
      if (file.exists(path_info)) {
        ifrm("info")
        info <- readRDS(path_info)
        assign("info", info, envir = .GlobalEnv)
      }
      
      if (file.exists(path_state)) {
        ifrm("state")
        state <- readRDS(path_state)
        assign("state", state, envir = .GlobalEnv)
      }
      
      if (file.exists(path_prior)) {
        ifrm("prior")
        prior <- readRDS(path_prior)
        list2env(prior, .GlobalEnv)
      }
      
    }
    
    # Saving info to avoid charging again and again each model ----------------
    set_model_to_stack <- set_model_to_stack %>%
      mutate(
        state = NA,
        info = NA,
        prior = NA,
        pathwd = NA
      )
    
    print(sum(set_model_to_stack$draws))
    print(sum(set_model_to_stack$draws_ame))
    print(set_model_to_stack)
    order_model <- seq(1, nrow(set_model_to_stack))
    for (i in order_model) {
      print(i)
      set_model_to_stack$state[i] <- str_replace(set_model_to_stack$fit[i], "/fit", paste0("/state"))
      set_model_to_stack$info[i] <- str_replace(set_model_to_stack$fit[i], "/fit", paste0("/info"))
      set_model_to_stack$prior[i] <- str_replace(set_model_to_stack$fit[i], "/fit", "/prior")
      set_model_to_stack$pathwd[i] <- w_set %>% mutate(
        name_fit = paste0("fit_", paste(
          ifelse(str_detect(type, "Type"), type, ""), link, whsv1, whsv2, sep = "_"
        )) %>%
          str_replace(., "__", "_"),
        name_fit = ifelse(
          type %in% c("202", "102", "220", "121"),
          paste0(name_fit, "_spatiotemporal_", type),
          name_fit
        )
      ) %>%
        filter(str_remove(name_fit, "fit\\_") == str_remove(gsub(
          ".*fit_", "", set_model_to_stack$fit[i]
        ), ".RDS")) %>%
        .$pathwd
      
      
      print(set_model_to_stack[i, ])
      
      # Generate draws ----------------------------------------------------------
      setting_paths(i)
      if (!file.exists(path_state)) {
        data_grid <- expand.grid(
          age = vector_age,
          females = c(0, 1),
          hsv2 = c(0, 1),
          idtime = vector_date,
          idspace = vector_space
        ) %>%
          mutate(
            levelintercept = case_when(
              hsv2 == 0 & females == 0 ~ 1,
              hsv2 == 0 & females == 1 ~ 2,
              hsv2 == 1 & females == 0 ~ 3,
              hsv2 == 1 & females == 1 ~ 4
            ),
            intercept = 1
          ) %>%
          left_join(
            dfbinomial %>% select(idspace, geometry) %>% distinct(idspace, .keep_all = T) %>% arrange(idspace),
            by = 'idspace'
          ) %>%
          st_as_sf
        
        max_weight_to_draw <- max_draw %>% filter(fit == path_fit) %>% .$draws %>% max
        print("To draw:")
        print(max_weight_to_draw)
        if (max_weight_to_draw > 0) {
          fit <- readRDS(path_fit)
          
          state <- evaluate_state(
            model = info,
            result = fit,
            data = data_grid,
            property = "sample",
            n = 10000,
            seed = 123
          )
          
          saveRDS(state, path_state)
          rm("fit")
          rm("info")
          gc()
        }
      }
    }
    
    ##Keeping only models with ndraws > 0
    set_model_to_stack <- set_model_to_stack %>%
      filter(draws > 0)
    print(set_model_to_stack)
    test_hpc(hpc)
    
    # PP ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    print("PP")
    if (passpp == F) {
      print("PP draws")
      if (!file.exists(path_pp_draws)) {
        draws_pp <- list()
      }
      else{
        print("Draws PP already available")
        draws_pp <- readRDS(path_pp_draws)
      }
      print("Using draws for PP")
      for (i in 1:nrow(set_model_to_stack)) {
        if (i > length(draws_pp)) {
          setting_paths(i)
          print("Evaluate model")
          state_to_use <- state[seq(1, ndraws, by = 1L)]
          draw_tmp <- inlabru::evaluate_model(
            model = info,
            state = state_to_use,
            data = dfbinomial,
            predictor = fml
          )
          rm("state_to_use")
          rm("state")
          print(draw_tmp %>% imap(~ .x %>% filter(p > 1)))
          print(draw_tmp %>% imap(~ .x %>% filter(p < 0)))
          draws_pp[[i]] <- draw_tmp
          
          
          test_hpc(hpc)
          saveRDS(draws_pp, path_pp_draws)
          gc()
        }
      }
      if (!file.exists(path_pp)) {
        draws_pp <- readRDS(path_pp_draws)
        pp_check <- pp_check_function(draws = draws_pp)
        print("Save")
        saveRDS(pp_check, path_pp)
      }
      else{
        print("PP Already available")
      }
    }
    
    # Age ------------------------------------------------------------------
    tible_pred <- expand.grid(
      age = vector_age,
      females = c(0, 1),
      hsv2 = c(0, 1),
      idtime = c(1, max(dfbinomial$idtime)),
      idspace = vector_space
    ) %>%
      mutate(
        levelintercept = case_when(
          hsv2 == 0 & females == 0 ~ 1,
          hsv2 == 0 & females == 1 ~ 2,
          hsv2 == 1 & females == 0 ~ 3,
          hsv2 == 1 & females == 1 ~ 4
        ),
        intercept = 1
      ) %>%
      left_join(
        dfbinomial %>% select(idspace, geometry) %>% distinct(idspace, .keep_all = T) %>% arrange(idspace),
        by = 'idspace'
      ) %>%
      st_as_sf
    
    
    if (passsummary == F) {
      print("Draws age")
      # Generate ----------------------------------------------------------------
      if (!file.exists(path_draws_age)) {
        draws_age <- list()
      } else{
        print("Draws age already available")
        draws_age <- readRDS(path_draws_age)
      }
      
      print("Age analysis")
      for (i in 1:nrow(set_model_to_stack)) {
        if (i > length(draws_age)) {
          setting_paths(i)
          
          print("Using state")
          state_to_use <- state[seq(1, ndraws, by = 1L)]
          draws_age[[i]] <- inlabru::evaluate_model(
            model = info,
            state = state_to_use,
            data = tible_pred,
            predictor = fml
          )
          
          rm("state_to_use")
          rm("state")
          test_hpc(hpc)
          saveRDS(draws_age, path_draws_age)
        }
      }
      gc()
      
      if (!file.exists(path_age_results)) {
        age_result <- list()
      } else{
        age_result <- readRDS(path_age_results)
      }
      print(length(age_result))
      gc()
      library(tidytable)
      setDTthreads(threads = cpudt)
      if (length(age_result) < 17) {
        print("Modifying draws")
        gen_age_sex <- draws_age %>%
          imap(function(.x, .y) {
            m <- .y
            .x %>%
              imap(
                ~ .x %>%
                  mutate(sample = paste(m, .y, sep = "_")) %>%
                  cbind(., tible_pred %>% st_drop_geometry())
              ) %>%
              do.call(rbind, .)
          }) %>%
          do.call(rbind, .) %>%
          numeric_to_text()
        gc()
      }
      
      if (length(age_result) < 1) {
        print("Prevalence")
        # Prevalence --------------------------------------------------------------
        age_result[["prev_age_sex"]] <- gen_age_sex %>%
          group_by(year_week,
                   idtime,
                   age,
                   groupsex,
                   virus,
                   insee_dep,
                   nom,
                   insee_reg,
                   nom_region) %>%
          rename(value = p) %>%
          summary_draw()
        gc()
        saveRDS(age_result, path_age_results)
      }
      print("Quantile prevalence")
      # Prevalence --------------------------------------------------------------
      function_for_quantile <- function(quantile) {
        df1 <- gen_age_sex %>%
          group_by(year_week, idtime, age, groupsex, virus, sample) %>%
          summarise(value = quantile(p, probs = quantile)) %>%
          group_by(year_week, idtime, age, groupsex, virus) %>%
          summary_draw() %>%
          mutate(metric = "prevalence")
        
        df2 <- gen_age_sex %>%
          group_by(year_week, idtime, age, groupsex, sample) %>%
          arrange(virus) %>%
          mutate(p = p / (p[1] + p)) %>%
          filter(virus == "HSV-2") %>%
          summarise(value = quantile(p, probs = quantile)) %>%
          group_by(year_week, idtime, age, groupsex) %>%
          summary_draw() %>%
          mutate(metric = "proportion") %>%
          mutate(virus = "HSV-2")
        gc()
        rbind(df1, df2)
      }
      if (length(age_result) < 2) {
        print("Quantile 1")
        age_result[["quantile_inf_inf"]] <- function_for_quantile(0.025)
        gc()
        saveRDS(age_result, path_age_results)
        
      }
      if (length(age_result) < 3) {
        print("Quantile 2")
        age_result[["quantile_inf"]] <- function_for_quantile(0.05)
        gc()
        saveRDS(age_result, path_age_results)
        
      }
      if (length(age_result) < 4) {
        print("Quantile 3")
        age_result[["quantile_median"]] <- function_for_quantile(0.5)
        gc()
        saveRDS(age_result, path_age_results)
        
      }
      if (length(age_result) < 5) {
        print("Quantile 4")
        age_result[["quantile_sup"]] <- function_for_quantile(0.95)
        gc()
        saveRDS(age_result, path_age_results)
      }
      if (length(age_result) < 6) {
        print("Quantile 5")
        age_result[["quantile_sup_sup"]] <- function_for_quantile(0.95)
        gc()
        saveRDS(age_result, path_age_results)
      }
      
      if (length(age_result) < 7) {
        fun_age_ref <- function(age_ref) {
          gen_age_sex %>%
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
            group_by(
              year_week,
              idtime,
              age,
              groupsex,
              virus,
              insee_dep,
              nom,
              insee_reg,
              nom_region
            )
        }
        print("Change prevalence wrt 20")
        diff_age_sex <- fun_age_ref(age_ref = 20)
        age_result[["rd_age"]] <- diff_age_sex %>%
          mutate(value = rd) %>%
          summary_draw()
        gc()
        saveRDS(age_result, path_age_results)
      }
      
      if (length(age_result) < 8) {
        print("Change prevalence wrt sex")
        diff_age_sex <- gen_age_sex %>%
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
        
        age_result[["rd_sex"]] <- diff_age_sex %>%
          mutate(value = rd) %>%
          summary_draw()
        gc()
        saveRDS(age_result, path_age_results)
      }
      
      if (length(age_result) < 9) {
        print("Change prevalence wrt virus")
        diff_virus <- gen_age_sex %>%
          group_by(
            year_week,
            idtime,
            age,
            groupsex,
            insee_dep,
            nom,
            insee_reg,
            nom_region,
            sample
          ) %>%
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
        
        rd_virus <- diff_virus %>%
          mutate(value = rd) %>%
          summary_draw()
        
        age_result[["prop_virus"]] <- diff_virus %>%
          mutate(value = prop) %>%
          summary_draw
        gc()
        saveRDS(age_result, path_age_results)
      }
      
      if (length(age_result) < 10) {
        print("Change prevalence wrt districts")
        diff_district <- gen_age_sex %>%
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
        
        age_result[["rd_district"]] <- diff_district %>%
          mutate(value = rd) %>%
          summary_draw()
        gc()
        saveRDS(age_result, path_age_results)
        
      }
      
      if (length(age_result) < 13) {
        print("Time")
        diff_time <- gen_age_sex %>%
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
          mutate(
            or = (p / (1 - p)) / (p[1] / (1 - p[1])),
            rr = p / p[1],
            rd = p - p[1],
            dummy = ifelse(p - p[1] > 0 , 1, 0)
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
        
        if (length(age_result) < 11) {
          age_result[["rd_time"]] <- diff_time %>%
            mutate(value = rd) %>%
            summary_draw()
          gc()
          saveRDS(age_result, path_age_results)
        }
        if (length(age_result) < 12) {
          age_result[["count_district_positive_trend"]] <- diff_time %>%
            ungroup() %>%
            group_by(idtime, groupsex, age, virus, sample) %>%
            summarise(value = sum(dummy)) %>%
            group_by(idtime, groupsex, age, virus) %>%
            summary_draw()
          gc()
          saveRDS(age_result, path_age_results)
        }
        if (length(age_result) < 13) {
          age_result[["count_district_positive_trend_diff_virus"]] <- diff_time %>%
            ungroup() %>%
            group_by(idtime, groupsex, age, virus, sample) %>%
            summarise(value = sum(dummy)) %>%
            group_by(age, groupsex, sample) %>%
            arrange(virus) %>%
            mutate(value = value - value[1]) %>%
            filter(virus == "HSV-2") %>%
            group_by(idtime, virus, age, groupsex) %>%
            summary_draw()
          gc()
          saveRDS(age_result, path_age_results)
        }
      }
      
      if (length(age_result) < 15) {
        print("Sex")
        diff_age_sex <- gen_age_sex %>%
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
        
        if (length(age_result) < 14) {
          age_result[["count_age_sex"]] <- diff_age_sex %>%
            ungroup() %>%
            group_by(year_week, idtime, age, virus) %>%
            summary_draw
          gc()
          saveRDS(age_result, path_age_results)
        }
        if (length(age_result) < 15) {
          age_result[["count_age_sex_diff_between_virus"]] <- diff_age_sex %>%
            ungroup() %>%
            group_by(year_week, idtime, age, sample) %>%
            arrange(virus) %>%
            mutate(value = value - value[1]) %>%
            filter(virus == "HSV-2") %>%
            group_by(year_week, idtime, age, virus) %>%
            summary_draw %>%
            mutate(var = "Difference b/w HSV-2 and 1 in the number of district with share prevalence greater among females than males")
          gc()
          saveRDS(age_result, path_age_results)
        }
      }
      
      if (length(age_result) < 17) {
        print("Virus")
        diff_age_sex <- gen_age_sex %>%
          group_by(
            year_week,
            idtime,
            groupsex,
            age,
            insee_dep,
            nom,
            insee_reg,
            nom_region,
            sample
          ) %>%
          mutate(x = ifelse(virus == "HSV-2", 1, 0)) %>%
          arrange(x) %>%
          mutate(dummy = ifelse(p / (p + p[1]) > 0.5, 1, 0)) %>%
          filter(virus == "HSV-2") %>%
          ungroup() %>%
          group_by(year_week, idtime, age, groupsex, sample) %>%
          summarise(value = sum(dummy))
        
        if (length(age_result) < 16) {
          age_result[["count_age_virus"]] <- diff_age_sex %>%
            group_by(year_week, idtime, age, groupsex) %>%
            summary_draw
          gc()
          saveRDS(age_result, path_age_results)
        }
        if (length(age_result) < 17) {
          age_result[["count_age_virus_diff_between_sex"]] <- diff_age_sex %>%
            ungroup() %>%
            group_by(year_week, idtime, age, sample) %>%
            arrange(groupsex) %>%
            mutate(value = value - value[1]) %>%
            filter(groupsex == "Males") %>%
            group_by(year_week, idtime, age, groupsex) %>%
            summary_draw %>%
            mutate(var = "Difference b/w males and females in the number of district with share HSV-2>50%")
          gc()
          saveRDS(age_result, path_age_results)
        }
      }
      detach("package:tidytable", unload = TRUE)
      rm("tible_pred")
      gc()
    }
    print("AME")
    ##Keeping only models with ndraws AME> 0
    set_model_to_stack_old <- set_model_to_stack
    set_model_to_stack <- set_model_to_stack %>%
      filter(draws_ame > 0)
    print(set_model_to_stack)
    
    print(sum(set_model_to_stack$draws_ame))
    # AME sex------------------------------------------------------------------
    if (passsex == F) {
      if (!file.exists(path_ame_sex)) {
        ame_result <- list()
        ame_result <- ame_function(
          model_set = set_model_to_stack,
          type_ame = "sexe",
          ref_ame = "Males",
          ref = 0,
          new = 1
        )
        gc()
        saveRDS(ame_result, path_ame_sex)
      }
    }
    # AME time------------------------------------------------------------------
    if (passtime == F) {
      if (!file.exists(path_ame_time)) {
        ame_result <- list()
        ame_result <- ame_function(
          model_set = set_model_to_stack,
          type_ame = "time",
          ref_ame = "1st week of 2022",
          ref = 1,
          new = max(dfdate$idtime)
        )
        gc()
        saveRDS(ame_result, path_ame_time)
      }
    }
    
    # AME district------------------------------------------------------------------
    if (passdistrict == F) {
      ###May require several run if lack of RAM
      if (!file.exists(path_ame_space)) {
        ame_result <- list()
      } else{
        ame_result <- readRDS(path_ame_space)
      }
      
      for (i in spatial_objects$geometry_dep$idspace) {
        d <- spatial_objects$geometry_dep %>%
          filter(idspace == i) %>%
          .$nom
        if (i > length(ame_result)) {
          ame_result[[i]] <- ame_function(
            model_set = set_model_to_stack,
            type_ame = "district",
            ref_ame = "Paris",
            ref = spatial_objects$geometry_dep %>% filter(nom == "Paris") %>% .$idspace,
            new = i
          )
          gc()
          saveRDS(ame_result, path_ame_space)
        }
      }
    }
    
    # AME age------------------------------------------------------------------
    if (passage == F) {
      ###May require several run if lack of RAM
      if (!file.exists(path_ame_age)) {
        ame_result <- list()
      } else{
        ame_result <- readRDS(path_ame_age)
      }
      
      for (i in seq(1, 101)) {
        if (i > length(ame_result)) {
          if (i == 101) {
            i <- 0
          }
          print(i)
          
          index <- ifelse(i == 0, 101, i)
          ame_result[[index]] <- ame_function(
            model_set = set_model_to_stack,
            type_ame = "age",
            ref_ame = "age 0",
            ref = 0,
            new = i
          )
          gc()
          saveRDS(ame_result, path_ame_age)
        }
      }
    }
    
    
    # AME district/time------------------------------------------------------------------
    if (passdistricttime == F) {
      ###May require several run if lack of RAM
      if (!file.exists(path_ame_spacetime)) {
        ame_result <- list()
      } else{
        ame_result <- readRDS(path_ame_spacetime)
      }
      
      for (i in spatial_objects$geometry_dep$idspace) {
        d <- spatial_objects$geometry_dep %>%
          filter(idspace == i) %>%
          .$nom
        if (i > length(ame_result)) {
          ame_result[[i]] <- ame_function(
            model_set = set_model_to_stack,
            type_ame = "districttime",
            ref_ame = paste0('Period 1 district ', i, ' (', d, ')'),
            ref = i,
            new = i
          )
          gc()
          saveRDS(ame_result, path_ame_spacetime)
        }
      }
    }
  }
  else{
    print('GCPO not available, pass')
  }
}
