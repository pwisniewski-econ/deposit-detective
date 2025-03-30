library(parallel)
library(tidyverse)
library(arrow)

parallel_irm <- function(DATA, xcols, dcol, ycol){
  
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
    ml_m = lrn("classif.ranger"), 
    ml_g = lrn("regr.ranger"), 
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



MICE_ADD_FULL <- read_rds("results_building/bank-additional-full-mice6.rds")


cl <- makeCluster(min(6, detectCores(logical = FALSE)))
clusterSetRNGStream(cl, iseed = 21548561)

results <- parLapply(cl, MICE_ADD_FULL, parallel_irm, xcols = x_names, dcol = "mid_week", ycol = "y")



dml_irm <- DoubleMLIRM$new(
  dml_data, 
  ml_m = lrn("classif.ranger"), 
  ml_g = lrn("regr.ranger"), 
  score = "ATE",  
  n_folds = 8
)


dml_irm$fit()
print(dml_irm$summary())



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






estimates  <- sapply(results, "[[", "coef")
ses <- sapply(results, "[[", "se")

rubin_pooling(estimates, ses, "MLRIM")



m <- length(estimates)

beta <- mean(estimates)
U_bar <- mean(ses^2)
B <- var(estimates)
T <- U_bar + (1 + 1/m) * B
se_pooled <- sqrt(T)
ci_lower <- beta_bar - 1.96 * SE_pooled
ci_upper <- beta_bar + 1.96 * SE_pooled

data.frame(name, beta, se_pooled, ci_lower, ci_upper)

name = "MLIRM"

z <- beta_bar / SE_pooled
pval <- 2 * (1 - pnorm(abs(z)))







unlist(results[[1]]$se)

results[[1]]$coef


estimates <- c(...)  # e.g., 10 coefficient estimates
ses <- c(...)        # corresponding standard errors

# Number of imputations
m <- length(estimates)

# Pooled estimate
beta_bar <- mean(estimates)

# Within-imputation variance
U_bar <- mean(ses^2)

# Between-imputation variance
B <- var(estimates)

# Total variance
T <- U_bar + (1 + 1/m) * B

# Pooled standard error
SE_pooled <- sqrt(T)

# 95% CI (normal approx)
ci_lower <- beta_bar - 1.96 * SE_pooled
ci_upper <- beta_bar + 1.96 * SE_pooled

# p-value
z <- beta_bar / SE_pooled
pval <- 2 * (1 - pnorm(abs(z)))

# Output
cat("Pooled estimate:", beta_bar, "\n")
cat("95% CI:", ci_lower, "to", ci_upper, "\n")
cat("p-value:", pval, "\n")















for (j in i){
  imputed_data_additional_temp <- as.data.table(complete(imputed_data_additional, action = i)) 
  dml_data <- DoubleMLData$new(imputed_data_additional_temp, y_col = y, d_cols = d, x_cols = x)
  dml_irm <- DoubleMLIRM$new(dml_data,
                             ml_m = learner_classif_m,
                             ml_g = learner_g,
                             score = "ATE", #or "ATTE"
                             n_folds = 10, n_rep = 1)
  dml_irm$fit()
  print(dml_irm$summary())
} 