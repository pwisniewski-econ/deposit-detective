# Title: Logit AME (Average Marginal Effects)
# Description: Run a basic logit on imputed data

# Library imports ----
library(tidyverse)
library(data.table)
library(arrow)
library(parallel)
requireNamespace("margins")

# importing data ----
MICE_ADD_FULL <- read_rds("results_building/bank-additional-full-mice6.rds")
MICE_FULL <- read_rds("results_building/bank-full-mice6.rds")

# Shared Resources ----
job_rubin <- function(job, data, model_name){
  #function to extract heterogeneity results
  data <- data |> filter(job_type == job)
  rubin_pooling(data$beta, data$se, paste0(model_name, job))
}

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

compute_logit <- function(DATA, xcols, dcol, job_list){
  library(tidyverse)
  library(data.table)
  
  job_ame <- function(job_type, variable, DATA, model){
    if(job_type=="all"){
      ame_out <- margins::margins(model, variables = variable) |>
        summary()
    } else{
      ame_out <- margins::margins(model, variables = variable,
                                  data = DATA |> filter(job==job_type)) |>
        summary()
    }
    data.frame(job_type, beta = ame_out$AME, se = ame_out$SE, row.names = NULL)
  }
  
  DATA <- DATA |> select(y, all_of(c(dcol,xcols)))
  logit_model <- glm(y ~ ., data = DATA, family = binomial)
  lapply(job_list, job_ame, dcol, DATA, logit_model) |> 
    rbindlist() #returns results from logit estimation (marginal effects); 
}

x_names <- c(
  "age", "job", "marital", "education", "default", "housing", "loan", 
  "previous", "poutcome", "contact", "pdays", "euribor_12mo", "hicp", 
  "cons_confidence", "unemployment", "stoxx_return"
)

job_list <- c(levels(MICE_ADD_FULL[[1]]$job), "all") 


# Main computations -----

## Parallel Setup -----
cl <- makeCluster(min(6,detectCores(logical = FALSE)))

## Compute Results for Weekdays ----
WEEKDAY_LOGIT <- parLapply(cl, MICE_ADD_FULL, compute_logit, x_names, "mid_week", job_list) |> 
  rbindlist()

## Compute Results for End of Month ----
EOM_LOGIT <- parLapply(cl, MICE_FULL, compute_logit, c(x_names, "balance"), "end_month", job_list) |> 
  rbindlist()

stopCluster(cl)

## Agregate with Rubin ----
RESULTS_WEEKDAY_LOGIT <- lapply(job_list, job_rubin, WEEKDAY_LOGIT, "logit-") |> rbindlist()
RESULTS_EOM_LOGIT <- lapply(job_list, job_rubin, EOM_LOGIT, "logit-") |> rbindlist()

## Export ----
write_feather(RESULTS_WEEKDAY_LOGIT, "results_analysis/weekdays_logit.feather", compression = "zstd")
write_feather(RESULTS_EOM_LOGIT, "results_analysis/eom_logit.feather", compression = "zstd")
