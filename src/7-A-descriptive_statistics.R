# Title: Macro Econ
# Description: Aim to estimate effects of the interest rate on take-up

# Library imports ----
library(arrow)
library(tidyverse)
library(data.table)

# Data imports ----
DATA_ADD_FULL <- read_feather("results_building/bank-additional-full.feather")
DATA_FULL <- read_feather("results_building/bank-full.feather")

