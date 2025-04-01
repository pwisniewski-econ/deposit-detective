# Title: Double ML / PLM / IRM
# Description: Run (and tune) a partially linear model along an IRM

# Library imports ----
library(parallel)
library(tidyverse)
library(arrow)
library(future)
library(DoubleML)
library(mlr3)
library(mlr3learners)
library(mlr3tuning)
library(mlr3pipelines)

set.seed(5132483)

# Data imports ----
MICE_ADD_FULL <- read_rds("results_building/bank-additional-full-mice6.rds")
MICE_FULL <- read_rds("results_building/bank-full-mice6.rds")

# Shared ressources ----
rubin_pooling <- function(estimates, ses, name) {
  #function used to pool our aggregates
  m <- length(estimates)
  
  # Pooled estimate and variance components
  pooled_estimate <- mean(estimates)
  within_var <- mean(ses^2)
  between_var <- var(estimates)
  
  # Total variance using Rubin's formula
  total_var <- within_var + (1 + 1/m) * between_var
  pooled_se <- sqrt(total_var)
  
  # Results
  z_stat <- pooled_estimate / pooled_se
  p_value <- 2 * (1 - pnorm(abs(z_stat)))
  ci_lower <- pooled_estimate - 1.96 * pooled_se
  ci_upper <- pooled_estimate + 1.96 * pooled_se
  
  data.frame(
    name = name,
    estimate = pooled_estimate, 
    std_error = pooled_se,
    p_value = p_value,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )
}

parallel_irm <- function(DATA, xcols, dcol, ycol, lrn_m, lrn_g){
  #helper function to run the IRM in parallel
  library(DoubleML)
  library(mlr3)
  library(mlr3learners)
  library(mlr3tuning)
  library(mlr3pipelines)
  
  dml_data <- DoubleMLData$new(
    as.data.table(DATA), 
    y_col = ycol, 
    d_cols = dcol, 
    x_cols = xcols
  )
  
  dml_irm <- DoubleMLIRM$new(
    dml_data, 
    ml_m = lrn_m, 
    ml_g = lrn_g, 
    score = "ATE",
    n_folds = 8
  )
  
  dml_irm$fit()
  
  return(dml_irm)
  
}

parallel_plm <- function(DATA, xcols, dcol, ycol, lrn_m, lrn_l){
  #helper function to run the PLM in parallel
  
  library(DoubleML)
  library(mlr3)
  library(mlr3learners)
  library(mlr3tuning)
  library(mlr3pipelines)
  
  dml_data <- DoubleMLData$new(
    as.data.table(DATA), 
    y_col = ycol, 
    d_cols = dcol, 
    x_cols = xcols
  )
  
  dml_irm <- DoubleMLPLR$new(
    dml_data, 
    ml_m = lrn_m, 
    ml_l = lrn_l, 
    score = "partialling out",
    n_folds = 8
  )
  
  dml_irm$fit()
  
  return(dml_irm)
  
}

x_names <- c(
  "age", "job", "marital", "education", "default", "housing", "loan", 
  "previous", "poutcome", "contact", "pdays", "euribor_12mo", "hicp", 
  "cons_confidence", "unemployment", "stoxx_return"
)

x_names_full <- x_names[x_names!="pdays"]

# Main logic ----

## Stage 1: Tuning ----

tune_dml <- function(DATA, xcols, treatment, method){
  #function to tune both IRM and PLM
  plan(multisession, workers = parallel::detectCores(logical = FALSE))
  
  dml_data <- DoubleMLData$new(
    as.data.table(DATA), 
    y_col = "y", 
    d_cols = treatment, 
    x_cols = xcols
  )
  
  if(method=="PLM"){
    dml_plr <- DoubleMLPLR$new(
      dml_data, 
      ml_m = lrn("classif.ranger"), 
      ml_l = lrn("regr.ranger"), 
      score = "partialling out",
      n_folds = 8
    )
    param_grid = list(
      "ml_l" = paradox::ps(
        num.trees = paradox::p_int(lower = 50, upper = 500),
        mtry = paradox::p_int(lower = 1, upper = 5),
        min.node.size = paradox::p_int(lower = 5, upper = 10))
    )
  }else if(method=="IRM"){
    dml_plr <- DoubleMLIRM$new(
      dml_data, 
      ml_m = lrn("classif.ranger"), 
      ml_g = lrn("regr.ranger"), 
      score = "ATE",
      n_folds = 8
    )
    param_grid = list(
      "ml_g" = paradox::ps(
        num.trees = paradox::p_int(lower = 50, upper = 500),
        mtry = paradox::p_int(lower = 1, upper = 5),
        min.node.size = paradox::p_int(lower = 5, upper = 10))
    )
  }
  
  param_grid = append(param_grid, list(
    "ml_m" = paradox::ps(
      num.trees = paradox::p_int(lower = 50, upper = 500),
      mtry = paradox::p_int(lower = 1, upper = 5),
      min.node.size = paradox::p_int(lower = 5, upper = 10))
  ))
  
  tune_settings <- list(
    rsmp_tune = rsmp("cv", folds = 5),
    terminator = trm("evals", n_evals = 25),
    algorithm = tnr("random_search")
  )
  
  dml_plr$tune(param_set = param_grid, tune_settings = tune_settings)
  
  return(dml_plr)
}

initial_tuning <- function(DATA, xcols, dcol, type){
  #helper function to start the previous function
  dml <- tune_dml(DATA, xcols, dcol, type)
  if(type=="PLM"){
    dml_2 <- dml$tuning_res[[1]]$ml_l$params[[1]]
  }else if(type=="IRM"){
    dml_2 <- dml$tuning_res[[1]]$ml_g1$params[[1]]
  }
  dml_m <- dml$tuning_res[[1]]$ml_m$params[[1]]
  tuned_dml_2 <- lrn("regr.ranger", num.trees = dml_2$num.trees, mtry = dml_2$mtry, min.node.size = dml_2$min.node.size)
  tuned_dml_m <- lrn("classif.ranger", num.trees = dml_m$num.trees, mtry = dml_m$mtry, min.node.size = dml_m$min.node.size)
  return(list(tuned_dml_m, tuned_dml_2))
}

### Tune First Dataset and re-use settings ----
tuned_week_plm <- initial_tuning(MICE_ADD_FULL[[1]], x_names, "mid_week", "PLM")
tuned_month_plm <- initial_tuning(MICE_FULL[[1]], x_names_full, "end_month", "PLM")

tuned_week_irm <- initial_tuning(MICE_ADD_FULL[[1]], x_names, "mid_week", "IRM")
tuned_month_irm <- initial_tuning(MICE_FULL[[1]], x_names_full, "end_month", "IRM")


## Stage 2: All imputations ----
extract_pooled <- function(results, name){
  #helper function to extract results from our estimates
  estimates  <- sapply(results, "[[", "coef")
  ses <- sapply(results, "[[", "se")
  rubin_pooling(estimates, ses, name)
}

## Parallel Setup ----
cl <- makeCluster(min(6, detectCores(logical = FALSE)))
clusterSetRNGStream(cl, iseed = 21548561)

## PLM for middle of the week ----
results_week <- parLapply(cl, MICE_ADD_FULL, parallel_plm, xcols = x_names, dcol = "mid_week", ycol = "y", tuned_week_plm[[1]], tuned_week_plm[[2]])
RESULTS_WEEK_DML <- extract_pooled(results_week, "MLPLR")
write_feather(RESULTS_WEEK_DML, "results_analysis/weekdays_dml.feather", compression = "zstd")

## IRM for middle of the week ----
results_week_irm <- parLapply(cl, MICE_ADD_FULL, parallel_irm, xcols = x_names, dcol = "mid_week", ycol = "y", tuned_week_irm[[1]], tuned_week_irm[[2]])
RESULTS_WEEK_IRM <- extract_pooled(results_week_irm, "MLIRM")
write_feather(RESULTS_WEEK_IRM, "results_analysis/weekdays_irm.feather", compression = "zstd")

## PLM for end of month ----
results_eom <- parLapply(cl, MICE_FULL, parallel_plm, xcols = x_names_full, dcol = "end_month", ycol = "y", tuned_month_plm[[1]], tuned_month_plm[[2]])
RESULTS_EOM_DML <- extract_pooled(results_eom, "MLPLR")
write_feather(RESULTS_EOM_DML, "results_analysis/eom_dml.feather", compression = "zstd")

## IRM for end of month ----
results_eom_irm <- parLapply(cl, MICE_FULL, parallel_plm, xcols = x_names_full, dcol = "end_month", ycol = "y", tuned_month_irm[[1]], tuned_month_irm[[2]])
RESULTS_EOM_IRM <- extract_pooled(results_eom_irm, "MLIRM")
write_feather(RESULTS_EOM_IRM, "results_analysis/eom_irm.feather", compression = "zstd")

stopCluster(cl)
