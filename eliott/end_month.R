
### R Code for Exploring Selection into Treatment (Motivating Matching)

# Load necessary packages
library(glmnet)
library(ggplot2)
library(dplyr)
library(MASS)
library(arrow)
library(hdm)
library(Matching)
library(pROC)
library(caret)  # For confusionMatrix function
library(tableone)
library(AIPW)
library(SuperLearner)
library(margins)
library(emmeans)

# importing data
dataset = read_feather("results_building/bank-full.feather")

################################################ Preparing the data

# Dropping rows with missing data
dataset <- dataset %>%
  filter(loan != "unknown" &
           housing != "unknown" &
           default != "unknown" &
           job != "unknown" &
           marital != "unknown" &
           education != "unknown" &
           poutcome != "unknown" &
           contact != "unknown"
  )

# removing some macro variables because we added by their equivalent
dataset <- dataset[, !names(dataset) %in% c("euribor3m", "cons.price.idx", "cons.conf.idx", "emp.var.rate", "nr.employed")]

# removing post-treatment variables
dataset <- subset(dataset, select = -duration)

# Convert the false numeric-like variables to categorical
str(dataset)
# pdays (nb days since called from a previous campaign)
pdays = dataset$pdays
pdays = pdays[pdays != 999]

# Plot histogram of pdays
hist(pdays, breaks = seq(min(pdays), max(pdays)),
     main = "Histogram of pdays (Bin size = 5)",
     xlab = "Values", 
     col = "skyblue", 
     border = "black")

summary(pdays) 

#  we transform pdays into a categorical variable, with 3 categories:
# never called (previously 999)
# called less than or equal to 6 days ago (the median)
# called more than 6 days ago
dataset$pdays <- as.factor(ifelse(dataset$pdays == 999, 
                                  "Never Contacted",
                                  ifelse(dataset$pdays <= 6, 
                                         "<=6 days ago", 
                                         ">6 days ago")))
str(dataset$pdays)

# keep only observation called once to ensure treatment is random given covariates
dataset = dataset %>% filter(campaign == 1)
dataset <- subset(dataset, select = -campaign)

# convert date_month to categorical
dataset$date_month = as.factor(dataset$date_month) 

# generate treatment

dataset <- dataset %>%
  mutate(
    # Convert date_month to a proper Date object if not already
    date_month = as.Date(date_month),
    
    # Create treatment variable: 1 if it's the last day of the month, else 0
    treatment = if_else(day > 16, 1, 0)
  )

# binary outcome
dataset$outcome <- ifelse(dataset$y == "yes", 1, 0)
dataset <- subset(dataset, select = -y)

# Check the structure of your data
str(dataset)

# imbalanced outcome
prop.table(table(dataset$outcome))

# Covariate Balance Table: checking
covariates = colnames(dataset)[colnames(dataset) != "treatment"]
covariates = covariates[covariates != "outcome"]

table1 <- CreateTableOne(vars = covariates, strata = "treatment", data = dataset, test = FALSE)
print(table1, smd = TRUE)

# distribution of treatment is quite balanced across control and treatment groups
# we can proceed safely with logit





###################################################### LOGIT

# strategy: incorporate month fixed effects (as dummies) and remove all macro variables
# to prevent instabability of coefficients due to correlation
data_logit <- subset(dataset, select = - c(euribor_12mo, hicp, cons_confidence, unemployment, stoxx_return))
str(data_logit)

# Run logistic regression
# we put date_month as factor to act as month fixed effects, thereby capturing
# economic recovery, and any effect related to the current state of the world in 
# Portugal at the time of the call
logit_model <- glm(outcome ~ treatment + .,
                   data = data_logit, family = binomial)

# Print summary of the model
summary(logit_model)
exp(coef(logit_model))

# Analysis of results
pred_logit <- predict(logit_model, data = dataset, type = "response")

# Generate ROC Curve and calculate AUC
roc_obj <- roc(dataset$outcome, pred_logit)
auc_value <- auc(roc_obj)

# Generate Confidence Interval
ci_obj <- ci.auc(roc_obj)
print(paste("AUC Confidence Interval:", round(ci_obj[1], 3), "-", round(ci_obj[3], 3)))

# Plotting the ROC Curve with ggplot
ggplot() +
  geom_line(aes(x = 1-roc_obj$specificities, y = roc_obj$sensitivities), color = "blue", size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(title = paste("ROC Curve (AUC =", round(auc_value, 3), 
                     "with CI:", round(ci_obj[1], 3), "-", round(ci_obj[3], 3), ")"),
       x = "1 - Specificity",
       y = "Sensitivity") +
  theme_minimal()

# Comment on the results
# The coefficient on treatment is significant at the 0.001% confidence level
# 
# The area under the curve is around 0.8 which reveals that the model is good and 
# is able to discriminate between individuals subscribing upon treatment and those who do not



######################### ATE ########################################


# Average Marginal Effect
ame <- margins(logit_model, variables = "treatment")
summary(ame)

# Heterogeneous
marital_ame <- margins(logit_model, variables = "treatment", at = list(marital = c("married", "single")))
summary(marital_ame)

me_by_job <- margins(logit_model, variables = "treatment", at = list(job = unique(data_logit$job)))
me_by_job = summary(me_by_job)

# plot 

me_df <- as.data.frame(me_by_job)
me_df$lower <- me_df$AME - 1.96 * me_df$SE
me_df$upper <- me_df$AME + 1.96 * me_df$SE

ggplot(me_df, aes(x = reorder(job, AME), y = AME)) +
  geom_point(size = 2) +  # Point estimate
  geom_linerange(aes(ymin = lower, ymax = upper), size = 0.6) +  # CI bounds
  coord_flip() +  # Horizontal layout
  labs(
    x = "Job Category",
    y = "Marginal Effect of Treatment (AME)",
    title = "Heterogeneous Treatment Effect by Job"
  ) +
  theme_minimal(base_size = 13)



################################################### AIPW
str(data_logit)

data_ipw = data_logit
W = subset(data_logit, select = -c(outcome, treatment))
X_numeric <- model.matrix( ~ . - 1, data = W)  # "-1" removes the intercept column
str(X_numeric)


#create an object
aipw_sl <- AIPW$new(Y=data_logit$outcome, 
                    A=data_logit$treatment,
                    W=X_numeric,
                    Q.SL.library = c("SL.glm", "SL.mean"),
                    g.SL.library = c("SL.glm", "SL.mean"),
                    k_split=10,verbose=TRUE)
#fit the object
aipw_sl$stratified_fit()
#calculate the results
aipw_sl$summary(g.bound = 0.025)
#check the propensity scores by exposure status after truncation
aipw_sl$plot.p_score()
print(aipw_sl$result, digits = 5)

# recompute LOGIT coefficient
# Values from your output
p0 <- aipw_sl$result["Risk of control", "Estimate"]  # Risk under control
att_rd <- aipw_sl$result["ATT Risk Difference", "Estimate"]
p1 <- p0 + att_rd

# Compute odds
odds0 <- p0 / (1 - p0)
odds1 <- p1 / (1 - p1)

# Compute odds ratio and logit coefficient
or <- odds1 / odds0
logit_coef <- log(or)

# Output
cat("Implied OR:", round(or, 4), "\n")
cat("Implied Logit Coefficient (log OR):", round(logit_coef, 4), "\n")


## plot 
# From logistic regression
logit_est <- coef(logit_model)["treatment"]
logit_se <- summary(logit_model)$coefficients["treatment", "Std. Error"]
logit_lower <- logit_est - 1.96 * logit_se
logit_upper <- logit_est + 1.96 * logit_se

# From AIPW - let's say you’re comparing log OR
# Use your actual values here
aipw_log_or <- log(aipw_sl$result["Odds Ratio", "Estimate"])  # OR from AIPW output
aipw_se <- aipw_sl$result["Odds Ratio", "SE"] / aipw_sl$result["Odds Ratio", "Estimate"]  # Delta method: SE(log OR) ≈ SE(OR) / OR
aipw_lower <- aipw_log_or - 1.96 * aipw_se
aipw_upper <- aipw_log_or + 1.96 * aipw_se

# Create comparison data frame
df <- tibble(
  method = c("Logit", "AIPW"),
  estimate = c(logit_est, aipw_log_or),
  lower = c(logit_lower, aipw_lower),
  upper = c(logit_upper, aipw_upper)
)

# Plot
ggplot(df, aes(x = method, y = estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  labs(title = "Treatment Effect (log OR): Logit vs AIPW",
       y = "Log Odds Ratio", x = "") +
  theme_minimal()

## Marginal effects
# 1. Extract marginal effect for "treatment" from logit model
ame_summary = summary(ame)
logit_marginal <- ame_summary[ame_summary$factor == "treatment", ]
logit_est <- logit_marginal$AME
logit_se <- logit_marginal$SE
logit_lower <- logit_est - 1.96 * logit_se
logit_upper <- logit_est + 1.96 * logit_se

# 2. Extract AIPW Risk Difference (ATE) from your object
aipw_result <- aipw_sl$result
aipw_rd_row <- aipw_result["Risk Difference", ]
aipw_est <- aipw_rd_row["Estimate"]
aipw_se <- aipw_rd_row["SE"]
aipw_lower <- aipw_est - 1.96 * aipw_se
aipw_upper <- aipw_est + 1.96 * aipw_se

# 3. Combine into a tidy data frame
df_plot <- tibble(
  Method = c("Logit AME", "AIPW Risk Difference"),
  Estimate = c(logit_est, aipw_est),
  Lower = c(logit_lower, aipw_lower),
  Upper = c(logit_upper, aipw_upper)
)

# 4. Plot
ggplot(df_plot, aes(x = Method, y = Estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.15) +
  labs(title = "Average Treatment Effect",
       y = "ATE (Probability Scale)", x = "") +
  theme_minimal()






### ##################################################
#ROBUSTNESS
##################################################


###################################################### Logit without month dummies but with macro 

data_logit <- subset(dataset, select = - c(date_month))
str(data_logit)

# Run logistic regression
# we put date_month as factor to act as month fixed effects, thereby capturing
# economic recovery, and any effect related to the current state of the world in 
# Portugal at the time of the call
logit_model <- glm(outcome ~ treatment + .,
                   data = data_logit, family = binomial)

# Print summary of the model
summary(logit_model)
exp(coef(logit_model))

# Analysis of results
pred_logit <- predict(logit_model, data = dataset, type = "response")

# Generate ROC Curve and calculate AUC
roc_obj <- roc(dataset$outcome, pred_logit)
auc_value <- auc(roc_obj)

# Generate Confidence Interval
ci_obj <- ci.auc(roc_obj)
print(paste("AUC Confidence Interval:", round(ci_obj[1], 3), "-", round(ci_obj[3], 3)))

# Plotting the ROC Curve with ggplot
ggplot() +
  geom_line(aes(x = 1-roc_obj$specificities, y = roc_obj$sensitivities), color = "blue", size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(title = paste("ROC Curve (AUC =", round(auc_value, 3), 
                     "with CI:", round(ci_obj[1], 3), "-", round(ci_obj[3], 3), ")"),
       x = "1 - Specificity",
       y = "Sensitivity") +
  theme_minimal()

############################################# logit with all covariates
data_logit <- dataset
logit_model <- glm(outcome ~ treatment + .,
                   data = data_logit, family = binomial)

# Print summary of the model
summary(logit_model)
exp(coef(logit_model))


####################################################  AIPW without month dummies with macro
data_logit <- subset(dataset, select = - c(date_month))
str(data_logit)

data_ipw = data_logit
W = subset(data_logit, select = -c(outcome, treatment))
X_numeric <- model.matrix( ~ . - 1, data = W)  # "-1" removes the intercept column
str(X_numeric)


#create an object
aipw_sl <- AIPW$new(Y=data_logit$outcome, 
                    A=data_logit$treatment,
                    W=X_numeric,
                    Q.SL.library = c("SL.glm", "SL.mean"),
                    g.SL.library = c("SL.glm", "SL.mean"),
                    k_split=10,verbose=TRUE)
#fit the object
aipw_sl$stratified_fit()
#calculate the results
aipw_sl$summary(g.bound = 0.025)
#check the propensity scores by exposure status after truncation
aipw_sl$plot.p_score()
print(aipw_sl$result, digits = 5)

# recompute LOGIT coefficient
# Values from your output
p0 <- aipw_sl$result["Risk of control", "Estimate"]  # Risk under control
att_rd <- aipw_sl$result["ATT Risk Difference", "Estimate"]
p1 <- p0 + att_rd

# Compute odds
odds0 <- p0 / (1 - p0)
odds1 <- p1 / (1 - p1)

# Compute odds ratio and logit coefficient
or <- odds1 / odds0
logit_coef <- log(or)

# Output
cat("Implied OR:", round(or, 4), "\n")
cat("Implied Logit Coefficient (log OR):", round(logit_coef, 4), "\n")


## plot 
# From logistic regression
logit_est <- coef(logit_model)["treatment"]
logit_se <- summary(logit_model)$coefficients["treatment", "Std. Error"]
logit_lower <- logit_est - 1.96 * logit_se
logit_upper <- logit_est + 1.96 * logit_se

# From AIPW - let's say you’re comparing log OR
# Use your actual values here
aipw_log_or <- log(aipw_sl$result["Odds Ratio", "Estimate"])  # OR from AIPW output
aipw_se <- aipw_sl$result["Odds Ratio", "SE"] / aipw_sl$result["Odds Ratio", "Estimate"]  # Delta method: SE(log OR) ≈ SE(OR) / OR
aipw_lower <- aipw_log_or - 1.96 * aipw_se
aipw_upper <- aipw_log_or + 1.96 * aipw_se

# Create comparison data frame
df <- tibble(
  method = c("Logit", "AIPW"),
  estimate = c(logit_est, aipw_log_or),
  lower = c(logit_lower, aipw_lower),
  upper = c(logit_upper, aipw_upper)
)

# Plot
ggplot(df, aes(x = method, y = estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  labs(title = "Treatment Effect (log OR): Logit vs AIPW",
       y = "Log Odds Ratio", x = "") +
  theme_minimal()

## Marginal effects
# 1. Extract marginal effect for "treatment" from logit model
ame_summary = summary(ame)
logit_marginal <- ame_summary[ame_summary$factor == "treatment", ]
logit_est <- logit_marginal$AME
logit_se <- logit_marginal$SE
logit_lower <- logit_est - 1.96 * logit_se
logit_upper <- logit_est + 1.96 * logit_se

# 2. Extract AIPW Risk Difference (ATE) from your object
aipw_result <- aipw_sl$result
aipw_rd_row <- aipw_result["Risk Difference", ]
aipw_est <- aipw_rd_row["Estimate"]
aipw_se <- aipw_rd_row["SE"]
aipw_lower <- aipw_est - 1.96 * aipw_se
aipw_upper <- aipw_est + 1.96 * aipw_se

# 3. Combine into a tidy data frame
df_plot <- tibble(
  Method = c("Logit AME", "AIPW Risk Difference"),
  Estimate = c(logit_est, aipw_est),
  Lower = c(logit_lower, aipw_lower),
  Upper = c(logit_upper, aipw_upper)
)

# 4. Plot
ggplot(df_plot, aes(x = Method, y = Estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.15) +
  labs(title = "Average Treatment Effect",
       y = "ATE (Probability Scale)", x = "") +
  theme_minimal()


