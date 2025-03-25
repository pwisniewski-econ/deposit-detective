library(dplyr)
library(Synth)
library(SCtools)
library(tidyverse)
library(skimr)
library(arrow)
library(ggplot2)

# import
datafull = read_feather("results_building/bank-additional-full.feather")

# Count the number of rows for each date
date_counts <- datafull %>%
  group_by(date_month) %>%
  summarise(Count = n())


# Plot the histogram using ggplot2
ggplot(date_counts, aes(x = date_month, y = Count)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(title = "Number of Rows for Each Date",
       x = "Date",
       y = "Count") +
  theme_minimal()

# count date if take up
yes_counts <- datafull %>%
  filter(y == "yes") %>% 
  group_by(date_month) %>%
  summarise(Count = n())

# keep only first calls
firstcall_yes <- datafull %>% 
  filter(campaign == 1) %>% 
  group_by(date_month) %>%
  summarise(Count = n())

# keep only second calls
secondcall_yes <- datafull %>%
  filter(campaign == 2) %>% 
  group_by(date_month) %>%
  summarise(Count = n())

# dataframe for first and second calls
date_counts_complete <- firstcall_yes %>%
  left_join(secondcall_yes, by = "date_month") %>%
  mutate(Count.y = ifelse(is.na(Count.y), 0, Count.y))  # Replace NAs with 0

# append to date_counts
date_counts$yes = yes_counts$Count # nb of yes
date_counts$success = date_counts$yes/date_counts$Count # compute share of yes
date_counts$first = firstcall_yes$Count/date_counts$Count # compute share of first calls
date_counts$second = date_counts_complete$Count.y/date_counts$Count # share of second calls


# Plot the histogram using ggplot2
ggplot() +
  geom_bar(data = date_counts, aes(x = date_month, y = Count), 
           stat = "identity", fill = "blue", alpha = 0.6) +
  geom_point(data = date_counts, aes(x = date_month, y = success*10000), 
             color = "red", size = 3) +
  geom_point(data = date_counts, aes(x = date_month, y = first*10000), 
             color = "green", size = 3) +
  scale_y_continuous(
    name = "Numbre of calls",           # Primary y-axis label
    sec.axis = sec_axis(~ ./10000, name = "Take-up rate and share of 1st calls %")  # Secondary y-axis label
  ) +
  # Show all months on x-axis
  scale_x_date(
    name = "Date",
    date_breaks = "2 month",   # Show all months
    date_labels = "%b %y"         # Display month names (Jan, Feb, Mar, etc.)
  ) +
  labs(title = "Left axis: Calls (bars) Right axis: Take-up rate (red), share of 1st calls (green) and 2nd Calls (orange)",
       x = "Date",
       y = "Count") +
  theme_minimal()


# check high effectiveness
middle = datafull %>% 
  filter(day_of_week %in% c("tue", "wed", "thu")) %>% 
  group_by(date_month) %>% 
  summarise(Count = n())

date_counts$middle = middle$Count/date_counts$Count # share of calls in middle of the week

# interest rate
euribor <- datafull %>%
  distinct(date_month, euribor_12mo)



# Plot the histogram using ggplot2
ggplot() +
  geom_bar(data = date_counts, aes(x = date_month, y = middle), 
           stat = "identity", fill = "blue", alpha = 0.6) +
  geom_point(data = date_counts, aes(x = date_month, y = success), 
             color = "red", size = 3) +
  geom_point(data = euribor, aes(x = date_month, y = euribor_12mo/10), 
             color = "orange", size = 3) +
  scale_y_continuous(
    name = "Share of calls and take-up rate %",           # Primary y-axis label
    sec.axis = sec_axis(~ .*10, name = "Interest rate %")  # Secondary y-axis label
  ) +
  # Show all months on x-axis
  scale_x_date(
    name = "Date",
    date_breaks = "2 month",   # Show all months
    date_labels = "%b %y"         # Display month names (Jan, Feb, Mar, etc.)
  ) +
  labs(title = "Share of calls in the middle of the week (bars) Take-up (red) Euribor 12m (orange)",
       x = "Date",
       y = "Count") +
  theme_minimal()


# checking percentage of yes
nb_yes = datafull %>% 
  filter(y == "yes") %>% 
  summarise(Count = n())
nb_call = datafull %>% 
  filter(y %in% c("yes", "no")) %>% 
  summarise(Count = n())

nb_yes/nb_call

counts <- datafull %>%
  group_by(date_month) %>%
  summarise(Count = n())

sum(counts$Count) == nb_call

month_count = data.frame(counts$Count)
success = data.frame(date_counts$success)
check = success * (month_count / as.integer(nb_call))
check = sum(check)
check == nb_yes/nb_call
# check done

