# Title: Data Import
# Description: This script imports, clean the raw datasets and appends external data.

#Library imports ----
library(tidyverse)
library(arrow)
library(data.table)

#Shared resources ----
feather_export <- function(df, filename){
  write_feather(as.data.frame(df), paste0("results_building/", filename), compression = "zstd")
}

## Conversion Table ----
MONTH2NUM <- data.frame(
  month = tolower(month.abb),  
  month_num = 1:12
)

# Clean and prepare external data ------

## ECB Inflation and EURIBOR ----
import_ecb <- function(filename, colname){
  ECB_DATA <- fread(
    file = paste0("data/external/", filename), 
    select = c(1,3),
    col.names = c("date_month", colname)
  ) 
  ECB_DATA <- ECB_DATA |> mutate(date_month = as.Date(date_month + 1))
  return(ECB_DATA)
}

ECB_EURIBOR <- import_ecb("ecb-12mo_euribor.csv", "euribor_12mo")
ECB_HICP <- import_ecb("ecb-hicp.csv", "hicp")

## EUROSTAT Confidence ----
EUROSTAT_CC <- fread(
  file = "data/external/eurostat-consumer_confidence.csv", 
  select = c(8,9),
  col.names = c("date_month", "cons_confidence")
) 

EUROSTAT_CC <- EUROSTAT_CC |> 
  mutate(
    date_month = as.Date(paste0(date_month, "-01"))
  )

## ILO Unemployment ----
ILO_UNEMP <- fread(
  file = "data/external/ilo-unemployment.csv", 
  select = c(4:7),
  col.names = c("sex", "age", "date_month", "unemployment")
) 

ILO_UNEMP <- ILO_UNEMP |>
  filter(sex=="Sex: Total", age=="Age (Youth, adults): 25+") |>
  mutate(date_month = as.Date(paste0(
   substr(date_month, 1, 4), "-", 
   substr(date_month, 6, 7), "-01"
  ))) |>
  select(date_month, unemployment)

## EUROSTOXX 600 ----
EUROSTOXX <- fread(
  file = "data/external/investing-eurostoxx600.csv", 
  select = c(1, 7),
  col.names = c("date_month", "stoxx_return")
) 

EUROSTOXX <- EUROSTOXX |>
  mutate(
    month = tolower(substr(date_month, 1, 3)), 
    stoxx_return = str_remove(stoxx_return, "%") |> as.numeric()
  ) |>
  left_join(MONTH2NUM, by = "month") |>
  mutate(date_month = as.Date(paste0(
    substr(date_month, 8, 11), "-", 
    str_pad(month_num, 2, pad = "0"), "-01"
  ))) |>
  select(date_month, stoxx_return)

supplementary_ls <- list(
  ECB_EURIBOR, ECB_HICP, EUROSTAT_CC, ILO_UNEMP, EUROSTOXX
)

# Main Data Import ----

main_import <- function(filename, supplementary_ls, add=F){
  
  cols <- c("y", "age", "date_month", "day")
  
  if(add==T){ cols <- cols[-4] }
  
  DATA <- fread(paste0("data/marketing/", filename))
  
  DATA <- DATA |>
  left_join(MONTH2NUM, by = "month") |>
    mutate(
      y = if_else(y=="yes", 1, 0),
      new_year = if_else(month_num < lag(month_num), 1, 0, missing = 0),
      year = 2008 + cumsum(new_year), 
      date_month = as.Date(paste0(
        year, "-",
        str_pad(month_num, width = 2, pad="0"), "-01"))
    ) |>
    select(all_of(cols), everything(), -c(month, year, new_year, month_num))
  
  for(DF in supplementary_ls){
    DATA <- DATA |>
      left_join(DF, by = "date_month")
  }
  
  return(DATA)

}

DATA_FULL <- main_import("bank-full.csv", supplementary_ls) |>
  mutate(end_month = if_else(day>16, 1, 0)) #Treatment variable

DATA_ADD_FULL <- main_import("bank-additional-full.csv", supplementary_ls, add=T) |>
  # We simplify pdays to have more meaningful groups 
  mutate(
    mid_week = if_else(day_of_week%in%c("mon", "fri"), 0, 1), #Treatment variable
    pdays = case_when( 
      pdays == 999 ~ "nc",
      pdays < 7 ~ "<7",
      pdays > 6 ~ "7+",
      TRUE ~ NA_character_
    )
  ) 

# Export ----
feather_export(DATA_FULL, "bank-full.feather")
feather_export(DATA_ADD_FULL, "bank-additional-full.feather")
