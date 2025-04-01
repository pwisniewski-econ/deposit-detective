# Title: Double ML / PLM / IRM
# Description: Run (and tune) an Ensemble Learner for a PLM

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

# Data Imports ----
MICE_ADD_FULL <- read_rds("results_building/bank-additional-full-mice6.rds")
MICE_FULL <- read_rds("results_building/bank-full-mice6.rds")

# Shared ressources ----
rubin_pooling <- function(estimates, ses, name) {
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


parallel_plm <- function(DATA, xcols, dcol, ycol, lrn_m, lrn_l){
  
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

extract_pooled <- function(results, name){
  estimates  <- sapply(results, "[[", "coef")
  ses <- sapply(results, "[[", "se")
  rubin_pooling(estimates, ses, name)
}

x_names <- c(
  "age", "job", "marital", "education", "default", "housing", "loan", 
  "previous", "poutcome", "contact", "pdays", "euribor_12mo", "hicp", 
  "cons_confidence", "unemployment", "stoxx_return"
)

x_names_full <- x_names[x_names!="pdays"]

# Main logic ----

## Stage 1: Tuning ----

### Define tuning settings for hyper parameter search
tune_settings <- list(
  
  terminator = trm("evals", n_evals = 10),
  algorithm = tnr("random_search")
)

tune_settings = list(
  rsmp_tune = rsmp("cv", folds = 3),
  terminator = trm("evals", n_evals = 10),
  tuner = tnr("random_search")
)

# Helper function to wrap a learner in an Auto-Tuner based on task type
auto_tuner <- function(learner, param_set, id) {
  # Set performance measure based on task type (classification or regression)
  if (grepl("^classif", learner$id)) {
    measure <- msr("classif.ce")
  } else {
    measure <- msr("regr.mse")
  }
  
  AutoTuner$new(
    learner = learner,
    resampling = rsmp("cv", folds = 5),
    measure = measure,
    search_space = param_set,
    terminator = tune_settings$terminator,
    tuner = tune_settings$tuner,
    store_models = TRUE,
    id = id
  )
}

### Ensemble for Classification ----
# Define individual learners
at_nnet <- auto_tuner(
  lrn("classif.nnet"),
  ps(
    size = p_int(1, 10),
    decay = p_dbl(0.0, 0.1)
  ),
  "nnet"
)

at_logreg <- lrn("classif.log_reg")  # Logistic regression (no tuning)

at_rf <- auto_tuner(
  lrn("classif.ranger"),
  ps(
    num.trees = p_int(50, 500),
    mtry = p_int(1, 5),
    min.node.size = p_int(5, 10)
  ),
  "ranger"
)

# Combine learners into an ensemble using a union and averaging operator
graph_ensemble_classif <- gunion(list(
  po("learner", at_nnet),
  po("learner", at_logreg),
  po("learner", at_rf)
)) %>>% po("classifavg")

# Create the ensemble learner for classification
ensemble_learner <- GraphLearner$new(graph_ensemble_classif)
ensemble_pipe_classif <- as_learner(ensemble_learner)

### Ensemble for Regression ---- 

# Define individual learners
at_lm <- lrn("regr.lm")  # Linear regression (no tuning)

at_nnet_reg <- auto_tuner(
  lrn("regr.nnet"),
  ps(
    size = p_int(1, 3),
    decay = p_dbl(0.0, 0.1)
  ),
  "nnet_reg"
)

at_rf_reg <- auto_tuner(
  lrn("regr.ranger"),
  ps(
    num.trees = p_int(75, 250),
    mtry = p_int(1, 5),
    min.node.size = p_int(4, 10)
  ),
  "rf_reg"
)

at_xgb <- auto_tuner(
  lrn("regr.xgboost"),
  ps(
    nrounds = p_int(50, 150),
    eta = p_dbl(0.01, 0.3),
    max_depth = p_int(2, 6)
  ),
  "xgb"
)

# Combine learners into an ensemble using a union and averaging operator
graph_ensemble_reg <- gunion(list(
  po("learner", at_lm),
  po("learner", at_nnet_reg),
  po("learner", at_rf_reg),
  po("learner", at_xgb)
)) %>>% po("regravg")

# Create the ensemble learner for regression
ensemble_learner_reg <- GraphLearner$new(graph_ensemble_reg)
ensemble_pipe_reg <- as_learner(ensemble_learner_reg)

## Stage 2: Double Machine Learning Setup and Execution ----
ens_wrapper <- function(DATA, xcols, treatment){
  #helper function to run the ensemble learner
  library(future)
  plan(multisession, workers = 3)
  DATA = as.data.frame(model.matrix(~ . - 1, data = DATA))
  parallel_plm(
    DATA, names(DATA)[stringr::str_detect(names(DATA), paste(xcols, collapse = "|"))],
    treatment, "y", ensemble_pipe_classif, ensemble_pipe_reg
  )
}

### Day of the week ----

#### Setup Parallel ----
cl <- makeCluster(min(5, detectCores(logical = FALSE)))
clusterSetRNGStream(cl, iseed = 81342556)
clusterExport(cl, varlist = list("ensemble_pipe_classif", "ensemble_pipe_reg", "parallel_plm"))

results_week_ens_ls <- parLapply(cl, MICE_ADD_FULL[1:5], ens_wrapper, x_names, "mid_week")
RESULTS_WEEK_ENS <- extract_pooled(results_week_ens_ls, "MLPLR-ENS") 
RESULTS_WEEK_ENS
write_feather(RESULTS_WEEK_ENS, "results_analysis/weekdays_dml_ens.feather", compression = "zstd")

stopCluster(cl)

### End of the month ----

#### Setup Parallel ----
cl <- makeCluster(min(5, detectCores(logical = FALSE)))
clusterSetRNGStream(cl, iseed = 5421257)
clusterExport(cl, varlist = list("ensemble_pipe_classif", "ensemble_pipe_reg", "parallel_plm"))

results_eom_ens_ls <- parLapply(cl, MICE_FULL[1:5], ens_wrapper, x_names_full, "end_month")
RESULTS_EOM_ENS <- extract_pooled(results_eom_ens_ls, "MLPLR-ENS") 
RESULTS_EOM_ENS
write_feather(RESULTS_EOM_ENS, "results_analysis/eom_dml_ens.feather", compression = "zstd")

stopCluster(cl)
