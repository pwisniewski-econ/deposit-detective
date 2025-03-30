DATA_ADD_FULL <- read_feather("results_building/bank-additional-full.feather")
DATA_FULL <- read_feather("results_building/bank-full.feather")




DATA_FULL |>
  mutate(
    Share = if_else(y==1, "Sucess", "Failure"),
    `Last observed call` = if_else(campaign>8, "9+", as.character(campaign)), 
  ) |>
  group_by( `Last observed call`, Share)|>
  summarise(count = n()) |>
  ggplot(aes(x=`Last observed call`, y=count, fill = `Share`))+
  geom_bar(stat="identity") +
  theme_minimal()+
  theme(text = element_text(family = "Times"), 
        legend.position = "bottom", 
        legend.margin = margin(-10,0,0,0),
        legend.title = element_text(size = 9), legend.text = element_text(size = 8))+
  labs(x= "", fill="", y= "Number of Calls", caption = "Bank Marketing, PT, 2008-2010")+
  guides(fill = guide_legend(override.aes = list(size = 0)))+
  scale_fill_brewer(palette = "Blues")


DATA_ADD_FULL |>
  group_by(Outcome = outcome) |>
  summarise(`Count` = n()) |>
  mutate(`Share (%)` = `Count`/sum(`Count`) * 100) |>
  knitr::kable(digits = 2, align = "c", caption = "Observed Outcomes")












DATA_FULL |>
  mutate(
    `First Call` = if_else(campaign>1, T, F),
    date = as.Date(paste(year(date_month), month(date_month), day, sep="-"))
  ) |>
  group_by(date_month, `First Call`)|>
  summarise(count = n()) |>
  ggplot(aes(x=date_month, y=count, fill = `First Call`))+
  geom_bar(stat="identity") +
  theme_minimal()+
  theme(text = element_text(family = "Times"), 
        legend.position = "bottom", 
        legend.margin = margin(-10,0,0,0),
        legend.title = element_text(size = 9), legend.text = element_text(size = 8))+
  labs(x= "", y= "Number of Calls", caption = "Bank Marketing, PT, 2008-2010")+
  guides(fill = guide_legend(override.aes = list(size = 0)))+
  scale_fill_brewer(palette = "Blues")

DATA_FULL |>
  mutate(
    `First Call` = if_else(campaign>1, T, F),
    date = as.Date(paste(year(date_month), month(date_month), day, sep="-"))
  ) |>
  group_by(date, `First Call`)|>
  summarise(count = n()) |>
  ggplot(aes(x=date, y=count, fill = `First Call`))+
  geom_bar(stat="identity") +
  theme_minimal()+
  theme(text = element_text(family = "Times", size=12), legend.position = "bottom")+
  labs(x= "", y= "Num", caption = "Bank Marketing, PT, 2008-2010")+
  scale_fill_brewer(palette = "Set1")

DATA_C <- DATA_ADD |>
  left_join(MONTH2NUM, by = "month") |>
  mutate(
    new_year = if_else(month_num < lag(month_num), 1, 0, missing = 0),
    year = 2008 + cumsum(new_year)
  )

LAST_CONTACTS |>
  ggplot(aes(x=date, y=count))+
  geom_line(linewidth = 1) +
  theme_minimal()+
  theme(text = element_text(family = "mono", size=18))+
  labs(x= "Date", y= "Num", 
       title= "Number of last contacts per date", 
       subtitle = "Full dataset")





DATA_ADD_FULL |>
  mutate(Share = if_else(y==1, "Sucess", "Failure"), day_of_week = factor(
    day_of_week, 
    levels = c("mon", "tue", "wed", "thu", "fri"), 
    labels = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday")
  )) |>
  group_by(day_of_week, Share) |>
  summarise(cnt = n()) |>
  mutate(prop = cnt/sum(cnt)) |>
  ungroup() |>
  select(-cnt) |>
  pivot_wider(names_from = "day_of_week", values_from = "prop") |>
  knitr::kable(digits = 3, align = "c", caption = "Share of sucessful calls")
