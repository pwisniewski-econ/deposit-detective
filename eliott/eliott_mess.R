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
  # geom_point(data = date_counts, aes(x = date_month, y = first*10000), 
             # color = "green", size = 3) +
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

consconf <- datafull %>%
  distinct(date_month, cons.conf.idx)

inf <- datafull %>%
  distinct(date_month, cons.price.idx)

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


# checking macro trend to see change in socio-economic characteristics


consconf <- datafull %>%
  distinct(date_month, cons.conf.idx)

inf <- datafull %>%
  distinct(date_month, cons.price.idx)

unemp <- datafull %>%
  distinct(date_month, unemployment)

emp <- datafull %>%
  distinct(date_month, nr.employed)


# Plot the histogram using ggplot2
ggplot() +
  geom_line(data = unemp, aes(x = date_month, y = unemployment), 
            stat = "identity", fill = "red", alpha = 0.6) +
  scale_y_continuous(
    name = "Unemployment rate %",           # Primary y-axis label
  ) +
  # Show all months on x-axis
  scale_x_date(
    name = "Date",
    date_breaks = "2 month",   # Show all months
    date_labels = "%b %y"         # Display month names (Jan, Feb, Mar, etc.)
  ) +
  labs(title = "Unemployment rate",
       x = "Date",
       y = "Count") +
  theme_minimal()



# Plot the histogram using ggplot2
ggplot() +
  geom_line(data = consconf, aes(x = date_month, y = cons.conf.idx), 
            stat = "identity", fill = "blue", alpha = 0.6) +
  scale_y_continuous(
    name = "Consumer confidence %",           # Primary y-axis label
  ) +
  # Show all months on x-axis
  scale_x_date(
    name = "Date",
    date_breaks = "2 month",   # Show all months
    date_labels = "%b %y"         # Display month names (Jan, Feb, Mar, etc.)
  ) +
  labs(title = "Consumer confidence",
       x = "Date",
       y = "Count") +
  theme_minimal()



# Plot the histogram using ggplot2
ggplot() +
  geom_line(data = emp, aes(x = date_month, y = nr.employed), 
            stat = "identity", fill = "blue", alpha = 0.6) +
  scale_y_continuous(
    name = "Consumer confidence %",           # Primary y-axis label
  ) +
  # Show all months on x-axis
  scale_x_date(
    name = "Date",
    date_breaks = "2 month",   # Show all months
    date_labels = "%b %y"         # Display month names (Jan, Feb, Mar, etc.)
  ) +
  labs(title = "Consumer confidence",
       x = "Date",
       y = "Count") +
  theme_minimal()



# understanding successful previous campaign
psuccess <- datafull %>% 
  filter(poutcome == "success") %>% 
  group_by(date_month)

# checking percentage of yes
nb_yes = psuccess %>% 
  filter(y == "yes") %>% 
  summarise(Count = n())
nb_call = psuccess %>% 
  filter(y %in% c("yes", "no")) %>% 
  summarise(Count = n())

result_inner <- nb_yes %>%
  inner_join(nb_call, by = "date_month")

result_inner$yes_prop = result_inner$Count.x/result_inner$Count.y

# failure

punsuccess <- datafull %>% 
  filter(poutcome == "failure") %>% 
  group_by(date_month)

nb_yes = punsuccess %>% 
  filter(y == "yes") %>% 
  summarise(Count = n())
nb_call = punsuccess %>% 
  filter(y %in% c("yes", "no")) %>% 
  summarise(Count = n())

result_inner_fail <- nb_yes %>%
  inner_join(nb_call, by = "date_month")

result_inner_fail$yes_prop = result_inner_fail$Count.x/result_inner_fail$Count.y

# never called

pnever <- datafull %>% 
  filter(poutcome == "nonexistent") %>% 
  group_by(date_month)

nb_yes = pnever %>% 
  filter(y == "yes") %>% 
  summarise(Count = n())
nb_call = pnever %>% 
  filter(y %in% c("yes", "no")) %>% 
  summarise(Count = n())

result_inner_never <- nb_yes %>%
  inner_join(nb_call, by = "date_month")

result_inner_never$yes_prop = result_inner_never$Count.x/result_inner_never$Count.y

# Plot the histogram using ggplot2
ggplot() +
  geom_line(data = result_inner, aes(x = date_month, y = yes_prop), 
            stat = "identity", color = "green", alpha = 0.6) +
  geom_line(data = result_inner_fail, aes(x = date_month, y = yes_prop), 
            stat = "identity", color = "red", alpha = 0.6) +
  geom_line(data = result_inner_never, aes(x = date_month, y = yes_prop), 
            stat = "identity", color = "black", alpha = 0.6) +
  scale_y_continuous(
    name = "Take-up rate %",           # Primary y-axis label
  ) +
  # Show all months on x-axis
  scale_x_date(
    name = "Date",
    date_breaks = "2 month",   # Show all months
    date_labels = "%b %y"         # Display month names (Jan, Feb, Mar, etc.)
  ) +
  labs(title = "Take-up for people who accepted (green), declined (red) or were never proposed (black) the product of a previous campaign",
       x = "Date",
       y = "Count") +
  theme_minimal()


# duration bias
yes_only = datafull %>% 
  filter(y == "yes") 

no_only = datafull %>% 
  filter(y == "no") 

# plot 
sample_no = no_only %>% slice_sample(n=5000)

ggplot() +
  # geom_histogram(data= sample_no, aes(x = duration), binwidth = 100, alpha = 0.6, fill = "red", color = "black") +
  geom_histogram(data= yes_only, aes(x = duration), binwidth = 100, alpha = 0.6, fill = "green", color = "black") +
  geom_vline(aes(xintercept = median(yes_only$duration)), 
             color = "green", 
             linetype = "dashed", 
             size = 1) +
  geom_vline(aes(xintercept = median(sample_no$duration)), 
             color = "red", 
             linetype = "dashed", 
             size = 1) +
  labs(title = "Distribution of call duration",
       x = "duration",
       y = "Frequency") +
  theme_minimal()

# reschedule
asked_reschedule = secondcall_yes$Count/(firstcall_yes$Count+secondcall_yes$Count)
firstcall_yes$reschedule = data.frame(asked_reschedule)
# Plot the histogram using ggplot2
ggplot() +
  geom_bar(data = date_counts, aes(x = date_month, y = Count/10000), 
           stat = "identity", fill = "blue", alpha = 0.6) +
  geom_point(data = firstcall_yes, aes(x = date_month, y = reschedule$asked_reschedule), 
            stat = "identity", fill = "red", alpha = 0.6) +
  scale_y_continuous(
    name = "Proportion who rescheduled a second call %",           # Primary y-axis label
  ) +
  # Show all months on x-axis
  scale_x_date(
    name = "Date",
    date_breaks = "2 month",   # Show all months
    date_labels = "%b %y"         # Display month names (Jan, Feb, Mar, etc.)
  ) +
  labs(title = "People rescheduling a second call",
       x = "Date",
       y = "Count") +
  theme_minimal()


