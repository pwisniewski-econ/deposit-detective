


### ##################################################
#ROBUSTNESS
##################################################

library(smotefamily)
library(ROSE)

###################################################### Logit without month fixed effects but with macro 

data_logit <- subset(dataset, select = - c(date_month))
str(data_logit)

# Run logistic regression
# we put date_month as factor to act as month fixed effects, thereby capturing
# economic recovery, and any effect related to the current state of the world in 
# Portugal at the time of the call
logit_model <- glm(outcome ~ treatment + .,
                   data = data_logit, family = binomial, weights = weights)

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
                   data = data_logit, family = binomial, weights = weights)

# Print summary of the model
summary(logit_model)
exp(coef(logit_model))




##################################################### Boostrap Logit


set.seed(123)

# Number of bootstrap samples
B <- 10
boot_coefs <- numeric(B)

n <- nrow(data_logit)  # your dataframe
boot_ames <- numeric(0)


for (i in 1:B) {
  index <- sample(1:n, replace = TRUE)
  boot_data <- data_logit[index, ]
  
  # Check that treatment has at least two levels
  if (length(unique(boot_data$treatment)) < 2) {
    next  # skip this iteration
  }
  
  # Try to fit the model and catch errors (e.g., separation)
  boot_model <- tryCatch({
    glm(outcome ~ treatment + ., data = boot_data, family = binomial())
  }, error = function(e) NULL)
  
  # If model failed, skip
  if (!is.null(boot_model)) {
    coef_treat <- tryCatch({
      coef(boot_model)["treatment"]
    }, error = function(e) NA)
    
    if (!is.na(coef_treat)) {
      boot_coefs <- c(boot_coefs, coef_treat)
    }
  }
  if (!is.null(boot_model)) {
    ame <- tryCatch({
      summary(margins(boot_model))[which(summary(margins(boot_model))$factor == "treatment"), "AME"]
    }, error = function(e) NA)
    
    if (!is.na(ame)) {
      boot_ames <- c(boot_ames, ame)
    }
  }
}

# Compute bootstrap mean, SE, and 95% CI
boot_mean <- mean(boot_coefs)
boot_se <- sd(boot_coefs)
boot_ci <- quantile(boot_coefs, probs = c(0.025, 0.975))

# Summary of bootstrap marginal effects
mean_ame <- mean(boot_ames)
se_ame <- sd(boot_ames)
ci_ame <- quantile(boot_ames, c(0.025, 0.975))

###################################################### IPW
prop.table(table(dataset$outcome))

# 
data_ipw = data_logit
data_ipw$date_month <- as.factor(data_ipw$date_month)
str(data_ipw)

# Fit a logistic regression to estimate propensity scores
propensity_model <- glm(treatment ~ .,
                        data = data_ipw, family = binomial)

# Summary of the model
summary(propensity_model)

# Predicted propensity scores
propensity_scores <- predict(propensity_model, type = 'response')
data_ipw$propensity_score <- propensity_scores

# Plot distribution of propensity scores by treatment status
p1 <- ggplot(data_ipw, aes(x = propensity_score, fill = factor(treatment))) +
  geom_density(alpha = 0.5) +
  labs(title = 'Distribution of Propensity Scores by Treatment Status', x = 'Propensity Score', y = 'Density')

print(p1)


# The distributions show a good overlap.
# The treated group's propensity scores are comparable to that of the control group.
# This indicates that the two groups are comparable, and that the treatment
# assignment can be considered random once controlling for covariates. 
# When estimating a PSM estimator, treated individuals will have comparable controls,
# this is necessary condition to estimate unbiased treatment effects. 
# For this reason, we conclude that we can safely use PSM or logit

# Compute IPW Weights
data_ipw$ipw <- ifelse(data_ipw$treatment == 1, 
                       1 / data_ipw$propensity_score, 
                       1 / (1 - data_ipw$propensity_score))


# remove propensity score
data_ipw <- subset(data_ipw, select = -propensity_score)

# Basic IPW logit
model_ipw <- glm(outcome ~ treatment, data = data_ipw, 
                 family = binomial(), weights = ipw)
summary(model_ipw)


library(survey)

# Create survey design using IPW
design <- svydesign(ids = ~1, weights = ~ipw, data = data_ipw)

# Estimate treatment effect
model_svy <- svyglm(outcome ~ treatment, design = design, family = quasibinomial())
summary(model_svy)



######################################## Logit LASSO

# Design matrix for covariates
X <- model.matrix(~ . -date_month, data = dataset)[, -1]  # drop intercept

# Outcome and treatment
Y <- dataset$outcome
D <- dataset$treatment

cvfit <- cv.glmnet(X, Y, family = "binomial", alpha = 1)

# Plot CV curve
plot(cvfit)

# Get optimal lambda
lambda_min <- cvfit$lambda.min
lambda_1se <- cvfit$lambda.1se
coef(cvfit, s = lambda_min)  # model coefficients at lambda.min


########################################### Post double Logit LASSO

# Design matrix for covariates
X <- model.matrix(~ ., data = dataset)[, -1]  # drop intercept

# Outcome and treatment
Y <- dataset$outcome
D <- dataset$treatment

lasso_treat <- cv.glmnet(X, D, family = "binomial", alpha = 1)
coef_treat <- coef(lasso_treat, s = "lambda.min")
vars_treat <- rownames(coef_treat)[which(coef_treat != 0)]

lasso_outcome <- cv.glmnet(X, Y, family = "binomial", alpha = 1)
coef_outcome <- coef(lasso_outcome, s = "lambda.min")
vars_outcome <- rownames(coef_outcome)[which(coef_outcome != 0)]

selected_vars <- union(vars_treat, vars_outcome)
selected_vars <- setdiff(selected_vars, "(Intercept)")  # remove intercept if present

X_selected <- as.data.frame(X)[, selected_vars, drop = FALSE]

final_data <- cbind(treatment = D, X_selected, outcome = Y)

# Fit logistic regression
model_post_lasso <- glm(outcome ~ treatment + ., data = final_data, family = binomial())
summary(model_post_lasso)


################################################################# SMOTE

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

################################################################ ROSE

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


###################################################### SMOTE/ROSE robustness

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

