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


# importing data
dataset = read_feather("results_building/bank-additional-full.feather")

### Preparing the data

# Dropping rows with missing data
dataset <- dataset %>%
  filter(loan != "unknown" &
           housing != "unknown" &
           default != "unknown" &
           job != "unknown" &
           marital != "unknown" &
           education != "unknown"
         )

# removing some variables added by us
dataset <- dataset[, !names(dataset) %in% c("euribor_12mo", "hicp", "cons_confidence", "unemployment", "stoxx_return")]

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

# generate treatment
dataset$treatment <- ifelse(dataset$day_of_week  %in% c("tue", "wed", "thu"), 1, 0)
dataset <- subset(dataset, select = -day_of_week)

# binary outcome
dataset$outcome <- ifelse(dataset$y == "yes", 1, 0)
dataset <- subset(dataset, select = -y)

# Check the structure of your data
str(dataset)

# Fit a logistic regression to estimate propensity scores
propensity_model <- glm(treatment ~ .,
                        data = dataset, family = binomial)

# Summary of the model
summary(propensity_model)

# Predicted propensity scores
propensity_scores <- predict(propensity_model, type = 'response')
dataset$propensity_score <- propensity_scores

# Plot distribution of propensity scores by treatment status
p1 <- ggplot(dataset, aes(x = propensity_score, fill = factor(treatment))) +
  geom_density(alpha = 0.5) +
  labs(title = 'Distribution of Propensity Scores by Treatment Status', x = 'Propensity Score', y = 'Density')

print(p1)

# The distributions show a good overlap.
# The treated group's propensity scores are comparable to that of the control group.
# This indicates that the two groups are comparable, and that the treatment
# assignment can be considered random once controlling for covariates. 
# When estimating a PSM estimator, treated individuals will have comparable controls,
# this is necessary condition to estimate unbiased treatment effects. 
# For this reason, we conclude that we can safely use PSM.

# remove propensity score
dataset <- subset(dataset, select = -propensity_score)

# imbalanced outcome
prop.table(table(dataset$outcome))

# Weights to mitigate imbalance in outcome 
weights <- ifelse(dataset$outcome == "Yes", 3.72, 1)

### LOGIT

# ALL NUMERIC

# turning all data to numeric
# Create the numeric matrix (model matrix) for predictors
X_numeric <- model.matrix(outcome ~ . - 1, data = data_smote)  # "-1" removes the intercept column
X_numeric = data.frame(X_numeric)

data_numeric = X_numeric
data_numeric$outcome = dataset$outcome

logit_model <- glm(outcome ~ treatment + .
                   - euribor3m
                   - emp.var.rate - cons.price.idx - cons.conf.idx - nr.employed
                   ,
                   data = data_numeric, family = binomial, weights = weights)

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



#################### NESTED CROSS VALIDATION #######################################


y <- data_numeric$outcome
x_matrix = X_numeric

# Define cross-validation folds
outer_folds <- createFolds(y, k = 5, list = TRUE, returnTrain = TRUE)


# Storing results
results <- data.frame(Fold = integer(), AUC = double(), Best_Threshold = double(), Accuracy = double(), F1_Score = double())
best_models <- list()
inner_results <- data.frame(Fold = integer(), Weight = double(), Threshold = double(), Mean_AUC = double())

for (i in 1:length(outer_folds)) {
  
  # Outer training and testing data
  train_index <- outer_folds[[i]]
  x_train <- x_matrix[train_index, ]
  y_train <- y[train_index]
  x_test <- x_matrix[-train_index, ]
  y_test <- y[-train_index]
  
  # Inner Cross-Validation Setup
  inner_folds <- createFolds(y_train, k = 5, list = TRUE)  # Inner 5-fold CV
  
  weight_grid <- seq(1, 10, by = 0.5)  # Weights to test
  threshold_grid <- seq(0.1, 0.9, by = 0.05)  # Thresholds to test
  best_auc <- 0
  best_model <- NULL
  best_weight <- 1
  best_threshold <- 0.5
  
  # Randomized Search CV
  search_iterations <- 30  # Number of random combinations to test
  
  for (k in 1:search_iterations) {
    
    # Randomly sample a weight and threshold combination
    weight <- sample(weight_grid, 1)
    threshold <- sample(threshold_grid, 1)
    auc_values <- c()  
    
    # Inner Cross-Validation Loop
    for (j in 1:length(inner_folds)) {
      inner_train_index <- inner_folds[[j]]
      
      # Shuffle the inner training data before each fold
      shuffle_index <- sample(seq_len(length(inner_train_index)))
      inner_train_index <- inner_train_index[shuffle_index]
      
      x_inner_train <- x_train[inner_train_index, ]
      y_inner_train <- y_train[inner_train_index]
      x_inner_test <- x_train[-inner_train_index, ]
      y_inner_test <- y_train[-inner_train_index]
      
      # Define weights (higher for positive class to handle imbalance)
      weights <- ifelse(y_inner_train == 1, weight, 1)
      
      # Fit logistic regression model
      cv_model <- glm(y_inner_train ~ ., family = binomial, weights = weights, data = as.data.frame(x_inner_train))
      
      # Make predictions
      predicted_probs_inner <- predict(cv_model, newdata = as.data.frame(x_inner_test), type = "response")
      auc_value_inner <- auc(y_inner_test, predicted_probs_inner)
      auc_values <- c(auc_values, auc_value_inner)
    }
    
    # Average AUC over inner folds
    mean_auc <- mean(auc_values)
    
    # Save inner results for plotting
    inner_results <- rbind(inner_results, data.frame(
      Fold = i,
      Weight = weight,
      Threshold = threshold,
      Mean_AUC = mean_auc
    ))
    
    # Update best model if average AUC improves
    if (mean_auc > best_auc) {
      best_auc <- mean_auc
      best_model <- cv_model
      best_weight <- weight
      best_threshold <- threshold
    }
  }
  
  # Test the best model on outer test set
  predicted_probs <- predict(best_model, newdata = as.data.frame(x_test), type = "response")
  predicted_class <- ifelse(predicted_probs > best_threshold, 1, 0)
  
  # Confusion Matrix
  confusion <- table(Predicted = predicted_class, Actual = y_test)
  precision <- confusion[2, 2] / sum(confusion[2, ])
  recall <- confusion[2, 2] / sum(confusion[, 2])
  f1_score <- 2 * (precision * recall) / (precision + recall)
  
  if (is.na(f1_score)) f1_score <- 0
  
  # Store results
  results <- rbind(results, data.frame(
    Fold = i,
    AUC = best_auc,
    Best_Threshold = best_threshold,
    Accuracy = sum(predicted_class == y_test) / length(y_test),
    F1_Score = f1_score
  ))
  
  best_models[[i]] <- list(model = best_model, weight = best_weight, threshold = best_threshold)
}

# Display results
print(results)

# Display average performance across folds
average_results <- colMeans(results[, -1])
print(average_results)

# Plot the 5 best combinations from inner loop
inner_results_top10 <- inner_results[order(-inner_results$Mean_AUC), ][1:100, ]

inner_plot <- ggplot(inner_results_top10, aes(x = Mean_AUC, y = paste0("(", Weight, ", ", Threshold, ")"))) +
  geom_point(color = 'red', size = 3) +
  labs(title = "Top 5 Best AUC Combinations from Inner Loop",
       x = "Mean AUC",
       y = "(Weight, Threshold)") +
  theme_minimal()

print(inner_plot)



####################################################

# Create a sequence of thresholds from 0.1 to 0.9
thresholds <- seq(0.05, 0.9, by = 0.05)

# Initialize a data frame to store results
results <- data.frame(Threshold = numeric(), 
                      Accuracy = numeric(), 
                      Sensitivity = numeric(), 
                      Specificity = numeric(), 
                      Precision = numeric(), 
                      F1_Score = numeric(), 
                      AUC = numeric())

# Loop over thresholds and calculate metrics
for (threshold in thresholds) {
  # Convert probabilities to predicted classes based on threshold
  predicted_classes <- ifelse(pred_logit > threshold, "1", "0")
  
  # Create a confusion matrix
  cm <- confusionMatrix(as.factor(predicted_classes), as.factor(dataset$outcome), positive = "1")
  
  # Extract metrics
  accuracy <- cm$overall["Accuracy"]
  sensitivity <- cm$byClass["Sensitivity"]
  specificity <- cm$byClass["Specificity"]
  precision <- cm$byClass["Precision"]
  f1_score <- cm$byClass["F1"]
  auc_value <- auc(roc(dataset$outcome, pred_logit))
  
  # Store metrics in the results data frame
  results <- rbind(results, data.frame(Threshold = threshold,
                                       Accuracy = accuracy,
                                       Sensitivity = sensitivity,
                                       Specificity = specificity,
                                       Precision = precision,
                                       F1_Score = f1_score,
                                       AUC = auc_value))
}

print(results)


# Find the threshold with the highest F1-Score
best_threshold <- results$Threshold[which.max(results$F1_Score)]
print(paste("Best Threshold (F1-Score):", best_threshold))

# Plotting Sensitivity, Specificity, and F1-Score against Thresholds
ggplot(results, aes(x = Threshold)) + 
  geom_line(aes(y = Sensitivity, color = "Sensitivity")) +
  geom_line(aes(y = Specificity, color = "Specificity")) +
  geom_line(aes(y = F1_Score, color = "F1-Score")) +
  labs(title = "Performance Metrics vs. Threshold",
       x = "Threshold",
       y = "Metric Value") +
  scale_color_manual(name = "Metrics", values = c("blue", "red", "green")) +
  theme_minimal()


# No weights 

# Fit a logistic regression model WITHOUT weights
logit_model_unweighted <- glm(outcome ~ ., data = data_numeric, family = binomial)

# Generate predicted probabilities without weights
predicted_probs_unweighted <- predict(logit_model_unweighted, newdata = data_numeric, type = "response")

# Compare ROC Curves (Weighted vs. Unweighted)
roc_weighted <- roc(dataset$outcome, pred_logit)
roc_unweighted <- roc(dataset$outcome, predicted_probs_unweighted)

# Plot the ROC Curves
plot(roc_weighted, col = "blue", main = "ROC Curves: Weighted vs. Unweighted Models", legacy.axes = TRUE)
plot(roc_unweighted, col = "red", add = TRUE)
legend("bottomright", legend = c("Weighted Model", "Unweighted Model"), col = c("blue", "red"), lty = 1)






#####################################################


### Double Post-LASSO Estimation (Without Matching)

# Define the outcome, treatment, and covariates
Y <- dataset$outcome
T <- dataset$treatment

# Convert 'month' to a factor variable (if not already)
# dataset$date_month <- as.factor(dataset$date_month)

# Create a matrix of covariates including 'date_month'
X <- model.matrix(~ . + factor(date_month) - date_month , data = dataset)[, -1]  # Removes intercept, one-hot-encode date_month

# Cross-validated LASSO for Propensity Score Estimation
cv_lasso_treatment <- cv.glmnet(X, dataset$treatment, family = "binomial", alpha = 1)

# Extract the best model based on cross-validation
lasso_treatment_model <- glmnet(X, dataset$treatment, family = "binomial", alpha = 1, lambda = cv_lasso_treatment$lambda.min)

# Cross-validated LASSO for Outcome Model Estimation
cv_lasso_outcome <- cv.glmnet(X, Y, family = "binomial", alpha = 1)

# Extract the best model based on cross-validation
lasso_outcome_model <- glmnet(X, Y, family = "binomial", alpha = 1, lambda = cv_lasso_outcome$lambda.min)

# Selected Covariates (Non-zero coefficients)
selected_covariates <- which(coef(lasso_treatment_model) != 0 | coef(lasso_outcome_model) != 0)
X_selected <- X[, selected_covariates, drop = FALSE]

# Post-LASSO Estimation
post_lasso_model <- glm(Y ~ T + X_selected, family = binomial)
summary(post_lasso_model)


### Double LASSO with Propensity Score Matching

# Propensity Score Matching
matching_result <- Match(Y = Y, Tr = T, X = propensity_scores, M = 1)

# Extract matched data
matched_data <- dataset[unique(c(matching_result$index.treated, matching_result$index.control)), ]

# Extract matched X matrix using only selected covariates
X_matched <- model.matrix(~ . + factor(date_month), data = matched_data[, !(names(matched_data) %in% c('Outcome', 'Treatment'))])
X_matched_selected <- X_matched[, selected_covariates, drop = FALSE]

# Fit a logistic regression with matched data
post_lasso_matched_model <- glm(matched_data$outcome ~ matched_data$treatment + X_matched_selected, family = binomial)
summary(post_lasso_matched_model)


# Double LASSO without Matching
lasso_outcome <- rlasso(Y ~ X)
lasso_treatment <- rlasso(T ~ X)

# Identify selected covariates
selected_covariates <- which(coef(lasso_outcome) != 0 | coef(lasso_treatment) != 0)
X_selected <- X[, selected_covariates, drop = FALSE]

# Post-LASSO estimation
post_lasso_model <- glm(Y ~ T + X_selected, family = binomial)
  summary(post_lasso_model)

# Double LASSO with Propensity Score Matching
matching_result <- Match(Y = Y, Tr = T, X = propensity_scores, M = 1)

# Outcome regression on matched sample
matched_data <- dataset[unique(c(matching_result$index.treated, matching_result$index.control)), ]
X_matched <- model.matrix(~ . + factor(date_month), data = matched_data)[, -1]
post_lasso_matched_model <- glm(matched_data$outcome ~ matched_data$treatment + X_matched, family = binomial)
summary(post_lasso_matched_model)



### ROBUSTNESS
library(smotefamily)
library(ROSE)

# Run logistic regression
# we put date_month as factor to act as month fixed effects, thereby capturing
# economic recovery, and any effect related to the current state of the world in 
# Portugal at the time of the call
logit_model <- glm(outcome ~ treatment + . + factor(date_month) - date_month
                   - euribor3m
                   - emp.var.rate - cons.price.idx - cons.conf.idx - nr.employed
                   ,
                   data = dataset, family = binomial, weights = weights)

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


# SMOTE

# create new dataset for SMOTE
data_smote = dataset

# Convert the outcome variable to a factor if it isn't already
data_smote$outcome <- as.factor(data_smote$outcome)
data_smote$date_month <- as.factor(data_smote$date_month)

# Create the numeric matrix (model matrix) for predictors
X_numeric <- model.matrix(outcome ~ . - 1, data = data_smote)  # "-1" removes the intercept column
X_numeric = data.frame(X_numeric)

# Applying SMOTE
smote_data <- SMOTE(X = X_numeric, 
                    target = data_smote$outcome,
                    K = 5,          # Number of nearest neighbors to consider
                    dup_size = 4)   # Amount of oversampling (2x = 200%)

# Extracting the new dataset
smote_data <- smote_data$data
smote_data$class <- as.numeric(smote_data$class)
smote_data$outcome = smote_data$class
smote_data <- subset(smote_data, select = -class)

# Checking the distribution of the outcome variable
prop.table(table(dataset$outcome)) # original dataset 
table(dataset$outcome) # original dataset 
table(smote_data$outcome) # SMOTE, oversamples minority class

# LOGIT SMOTE
logit_model <- glm(outcome ~ treatment + .,
                   data = smote_data, family = binomial)

# Print summary of the model
summary(logit_model)
exp(coef(logit_model))

# Analysis of results
pred_logit <- predict(logit_model, data = smote_data, type = "response")

# Generate ROC Curve and calculate AUC
roc_obj <- roc(smote_data$outcome, pred_logit)
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

# ROSE

# create new dataset for ROSE (take formatted for smote)
data_rose = X_numeric
data_rose$outcome = dataset$outcome

# perform ROSE
rose_data <- ROSE(outcome ~ ., data = data_rose, seed = 1)$data

# Checking the distribution of the outcome variable
table(dataset$outcome) # original dataset 
table(rose_data$outcome) # ROSE, oversamples minority and undersamples majority

logit_model <- glm(outcome ~ treatment + .,
                   data = rose_data, family = binomial)

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







#  CHECK 1: increasing dup_size in SMOTE, that is oversampling even more minority class (yes)
# create new dataset for SMOTE
data_smote = dataset

# Convert the outcome variable to a factor if it isn't already
data_smote$outcome <- as.factor(data_smote$outcome)
data_smote$date_month <- as.factor(data_smote$date_month)

# Create the numeric matrix (model matrix) for predictors
X_numeric <- model.matrix(outcome ~ . - 1, data = data_smote)  # "-1" removes the intercept column
X_numeric = data.frame(X_numeric)

# Applying SMOTE
smote_data <- SMOTE(X = X_numeric, 
                    target = data_smote$outcome,
                    K = 5,          # Number of nearest neighbors to consider
                    dup_size = 5)   # Amount of oversampling (2x = 200%)

# Extracting the new dataset
smote_data <- smote_data$data
smote_data$class <- as.numeric(smote_data$class)
smote_data$outcome = smote_data$class
smote_data <- subset(smote_data, select = -class)

# Checking the distribution of the outcome variable
prop.table(table(dataset$outcome)) # original dataset 
table(dataset$outcome) # original dataset 
table(smote_data$outcome) # SMOTE, oversamples minority class

# LOGIT SMOTE
logit_model <- glm(outcome ~ treatment + .,
                   data = smote_data, family = binomial)

# Print summary of the model
summary(logit_model)
exp(coef(logit_model))

# Analysis of results
pred_logit <- predict(logit_model, data = smote_data, type = "response")

# Generate ROC Curve and calculate AUC
roc_obj <- roc(smote_data$outcome, pred_logit)
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

# Comments
# AUC remains the same
# but intercept of odd-ratio now is more in favour of yes

#  CHECK 2
#  decreasing dup_size in SMOTE, that is oversampling less minority class (yes)
# create new dataset for SMOTE
data_smote = dataset

# Convert the outcome variable to a factor if it isn't already
data_smote$outcome <- as.factor(data_smote$outcome)
data_smote$date_month <- as.factor(data_smote$date_month)

# Create the numeric matrix (model matrix) for predictors
X_numeric <- model.matrix(outcome ~ . - 1, data = data_smote)  # "-1" removes the intercept column
X_numeric = data.frame(X_numeric)

# Applying SMOTE
smote_data <- SMOTE(X = X_numeric, 
                    target = data_smote$outcome,
                    K = 5,          # Number of nearest neighbors to consider
                    dup_size = 3)   # Amount of oversampling (2x = 200%)

# Extracting the new dataset
smote_data <- smote_data$data
smote_data$class <- as.numeric(smote_data$class)
smote_data$outcome = smote_data$class
smote_data <- subset(smote_data, select = -class)

# Checking the distribution of the outcome variable
prop.table(table(dataset$outcome)) # original dataset 
table(dataset$outcome) # original dataset 
table(smote_data$outcome) # SMOTE, oversamples minority class

# LOGIT SMOTE
logit_model <- glm(outcome ~ treatment + .,
                   data = smote_data, family = binomial)

# Print summary of the model
summary(logit_model)
exp(coef(logit_model))

# Analysis of results
pred_logit <- predict(logit_model, data = smote_data, type = "response")

# Generate ROC Curve and calculate AUC
roc_obj <- roc(smote_data$outcome, pred_logit)
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

# Comments
# some coefficients are turning into NA, they cannot be estimated

