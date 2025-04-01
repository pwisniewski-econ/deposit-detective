# Title: Data Imputations
# Description: Check for "unknowns" and imputate values with MICE

# Library imports ----
library(parallel)
library(tidyverse)
library(arrow)
requireNamespace("mice") 
#requires package but avoids namespace collisions by attaching it

# Shared functions ----
par_mice <- function(DATA){
  DATA <- mice::mice(DATA, m=1) |>
    mice::complete(action = 1)
}

unknown_imputations <- function(DATA, cl){
  char_u <- DATA |> 
    summarise(across(where(is.character), ~ any(grepl("unknown", . , ignore.case = TRUE)))) 
  
  char_u <- names(char_u)[char_u==T]
  # For DATA_ADD_FULL unknowns are:
  # [1] "job" "marital" "education" "default" "housing" "loan"  
  
  DATA <- DATA |>
    mutate(across(all_of(char_u), ~if_else(.x == "unknown", NA_character_, .x))) |>
    mutate(across(where(is.character), ~as.factor(.x)))
  
  data_ls <- parLapply(cl, rep(list(DATA), 6), par_mice)
  
  return(data_ls) #this returns data after imputation
}

# Import Data -----
DATA_FULL <- read_feather("results_building/bank-full.feather") |> filter(campaign==1)
DATA_ADD_FULL <- read_feather("results_building/bank-additional-full.feather") |> filter(campaign==1)

# Main imputations ------

## Parallel Setup (with RNG management) ----- 
cl <- makeCluster(min(6, detectCores(logical = FALSE)))
clusterSetRNGStream(cl, iseed = 584461256)

## Run Imputations and Export ----
MICE_FULL <- unknown_imputations(DATA_FULL, cl)
write_rds(MICE_FULL, "results_building/bank-full-mice6.rds", compress = "xz")

MICE_ADD_FULL <- unknown_imputations(DATA_ADD_FULL, cl)
write_rds(MICE_ADD_FULL, "results_building/bank-additional-full-mice6.rds", compress = "xz")

stopCluster(cl)
