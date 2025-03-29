library(arrow)
library(tidyverse)
library(data.table)
library(lfe)
library(future)
library(DoubleML)
library(mlr3)
library(mlr3learners)
library(paradox)
library(mlr3tuning)
library(margins)
library(future)
library(gmm)
set.seed(21457912)

plan(multisession, workers = 8)


## Data prep ----
DATA_ADD_FULL <- read_feather("results_building/bank-additional-full.feather")

LAG_EURIBOR <- DATA_ADD_FULL |> 
  select(date_month, euribor_12mo) |> 
  unique() |> 
  mutate(euribor_12mo_l1 = lag(euribor_12mo))

REG_DATA <- DATA_ADD_FULL |> 
  mutate(y=if_else(y=="no", 0, 1)) |>
  left_join(LAG_EURIBOR, by= c("date_month", "euribor_12mo")) |>
  mutate(across(c(poutcome, marital, job, housing, default, education), ~as.factor(.x))) |> 
  filter(!is.na(euribor_12mo_l1)&campaign==1)

## 2SLS ----

start_2SLS <- Sys.time()
IV_2SLS <- felm(y ~ unemployment + housing + default +  poutcome + job + marital + education | 0 | (euribor_12mo ~ euribor_12mo_l1), data=REG_DATA)
compute_2SLS <- Sys.time() - start_2SLS

summary(IV_2SLS)

## Double ML

dml_data <- DoubleMLData$new(
  setDT(reg_data),
  y_col = "y",
  d_cols = "euribor_12mo",
  z_cols = "euribor_12mo_l1",
  x_cols = c("unemployment", "poutcome", "job", "marital", "housing", "default", "education")
)

dml_pliv <- DoubleMLPLIV$new(
  data = dml_data,
  ml_g = lrn("regr.ranger"),
  ml_m = lrn("regr.ranger"),
  ml_l = lrn("regr.ranger"),
  ml_r = lrn("regr.ranger"),
  n_folds = 8,
  n_rep = 4,
  score = "IV-type"                 
)

param_set <- list(
  ml_l = ps(num.trees = p_int(lower = 50, upper = 200)),
  ml_m = ps(num.trees = p_int(lower = 50, upper = 200)),
  ml_r = ps(num.trees = p_int(lower = 50, upper = 200))
)

tune_settings <- list(
  rsmp_tune = rsmp("cv", folds = 5),
  terminator = trm("evals", n_evals = 25),
  algorithm = tnr("random_search"),
  measure = list(ml_l = msr("regr.mse"), ml_m = msr("regr.mse"), ml_r = msr("regr.mse"))
)

dml_pliv$tune(param_set = param_set, tune_settings = tune_settings)

dml_pliv$fit()
dml_pliv$pval

dml_pliv$fit(store_predictions = TRUE)
dml_pliv$learner$ml_l

dml_pliv$summary()

print(dml_pliv)

dml_pliv$summary()
plot(dml_pliv, type = "overlap")
## Probit-IV

reg_data2 <- reg_data
cov(dml_pliv$psi_a, dml_pliv$psi_b)

# 2nd stage: include residuals in a standard logit/probit
start_PIV <- Sys.time()
first <- lm(euribor_12mo ~ euribor_12mo_l1 + unemployment + housing + default + poutcome + job + marital+education, data = reg_data2)
reg_data2$resid1 <- resid(first)
probit_iv <- glm(y ~ euribor_12mo + resid1 + unemployment  + housing + default + poutcome + job + marital+education, 
                 family = binomial(link = "probit"), data = reg_data2)
compute_PIV <- Sys.time() - start_PIV
summary(probit_iv)


mfx <- margins::margins(probit_iv, variables = "euribor_12mo")
summary(mfx)
