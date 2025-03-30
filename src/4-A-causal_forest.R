library(grf)
library(data.table)

x_names <- c(
  "age", "job", "marital", "education", "default", "housing", "loan", 
  "previous", "poutcome", "contact", "pdays", "euribor_12mo", "hicp", 
  "cons_confidence", "unemployment", "stoxx_return"
)

MICE_ADD_FULL <- read_rds("results_building/bank-additional-full-mice6.rds")
MICE_FULL <- read_rds("results_building/bank-full-mice6.rds")

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

job_ate <- function(job_type, data, target, cf){
  if(job_type=="all"){
    job_ate <- average_treatment_effect(cf, target.sample = target)
  } else{
    job_ate <- average_treatment_effect(cf, subset = (data$job == job_type), target.sample = target)
  }
  data.frame(job_type, beta = job_ate[1], se = job_ate[2], row.names = NULL)
}

compute_cf <- function(data, x_names, y_name, w_name, job_list, target){
  mm <- model.matrix(~ ., data = data |> select(all_of(x_names)))
  y <- data |> select(all_of(y_name)) |> unlist()
  w <- data |> select(all_of(w_name)) |> unlist()
  cf <- causal_forest(
    mm, y, w,  
    num.threads = detectCores(logical=F), sample.fraction = 0.296923, 
    mtry = 12, min.node.size = 2, honesty.fraction = 0.799754, 
    honesty.prune.leaves = 1, imbalance.penalty = 1.015844,
    alpha = 0.01950416, num.trees = 1000
  ) 
  lapply(job_list, job_ate, data, target, cf) |> rbindlist()
}

job_rubin <- function(job, data){
  data <- data |> filter(job_type == job)
  rubin_pooling(data$beta, data$se, paste0("causal-forest-", job))
}


job_list <- c(levels(MICE_ADD_FULL[[1]]$job), "all")

WEEKDAYS <- lapply(MICE_ADD_FULL, compute_cf, x_names, "y", "mid_week", job_list, "overlap") |>
  rbindlist()

WEEKDAYS_results <- lapply(job_list, job_rubin, WEEKDAYS) |> rbindlist()

EOM <- lapply(MICE_FULL, compute_cf, c(x_names, "balance"), "y", "end_month", job_list, "overlap") |>
  rbindlist()

EOM_results <- lapply(job_list, job_rubin, EOM) |> rbindlist()

write_feather(WEEKDAYS_results, "results_analysis/weekdays_causal_forest.feather")
write_feather(EOM_results, "results_analysis/eom_causal_forest.feather")

