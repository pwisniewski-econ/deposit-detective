### 1- Trimming

# Trimming dataset to keep observations within common support region
trim_threshold <- 0.6  # Based on visual inspection, you can adjust this if needed
trimmed_data <- dataset[dataset$propensity_score <= trim_threshold, ]

# Check overlap after trimming
p_trimmed <- ggplot(trimmed_data, aes(x = propensity_score, fill = factor(treatment))) +
  geom_density(alpha = 0.5) +
  labs(title = 'Trimmed Propensity Score Distribution by Treatment Status',
       x = 'Propensity Score',
       y = 'Density')

print(p_trimmed)


### 2- Caliper Matching 
library(Matching)

# Define propensity scores and treatment indicator
propensity_scores <- dataset$propensity_score
T <- dataset$treatment
Y <- dataset$outcome

# Caliper width (0.1 * SD of logit of propensity score)
logit_ps <- log(propensity_scores / (1 - propensity_scores))
caliper_width <- 0.1 * sd(logit_ps)

# Perform matching with caliper
matching_result <- Match(Y = Y, Tr = T, X = propensity_scores, M = 1, caliper = caliper_width)

# Extract matched data
matched_data <- dataset[unique(c(matching_result$index.treated, matching_result$index.control)), ]

# Check overlap after caliper matching
p_matched <- ggplot(matched_data, aes(x = propensity_score, fill = factor(treatment))) +
  geom_density(alpha = 0.5) +
  labs(title = 'Matched Propensity Score Distribution by Treatment Status',
       x = 'Propensity Score',
       y = 'Density')

print(p_matched)

### 3- IPW with Truncation

# One-hot-encoding variables
covariates <- fastDummies::dummy_cols(covariates, remove_first_dummy = TRUE, remove_selected_columns = TRUE)


# Set truncation threshold (e.g., discarding extreme weights)
truncate_threshold <- 0.01  # You can adjust this to be more or less restrictive

# Calculate weights
weights <- ifelse(T == 1, 1 / propensity_scores, 1 / (1 - propensity_scores))

# Truncate extreme weights
weights[weights > quantile(weights, 1 - truncate_threshold)] <- quantile(weights, 1 - truncate_threshold)
weights[weights < quantile(weights, truncate_threshold)] <- quantile(weights, truncate_threshold)

# Apply weighted regression
ipw_model <- glm(Y ~ T, weights = weights, family = binomial)
summary(ipw_model)

# Required libraries
library(tableone)  # For SMD calculation
library(fastDummies) # for one-hot-encoding

# 1. Checking Effective Sample Size (ESS)
ess <- sum(weights)^2 / sum(weights^2)
cat("Effective Sample Size (ESS):", ess, "\n")

# 2. Assessing Covariate Balance Before and After Weighting

calculate_weighted_mean <- function(var, weights) {
  weighted.mean(var, w = weights)
}

calculate_weighted_smd <- function(var, treat, weights) {
  treated_mean <- calculate_weighted_mean(var[treat == 1], weights[treat == 1])
  control_mean <- calculate_weighted_mean(var[treat == 0], weights[treat == 0])
  pooled_sd <- sqrt((var(treat == 1) + var(treat == 0)) / 2)
  smd <- abs(treated_mean - control_mean) / pooled_sd
  return(smd)
}

# Get covariates
covariates <- dataset[, !(colnames(dataset) %in% c('Outcome', 'Treatment', 'propensity_score', 'weights'))]

# Calculate SMDs before and after weighting
smds_before <- sapply(covariates, function(x) calculate_weighted_smd(x, dataset$treatment, rep(1, length(x))))
smds_after <- sapply(covariates, function(x) calculate_weighted_smd(x, dataset$treatment, dataset$weights))

# Display the SMDs
smd_table <- data.frame(Variable = names(smds_before),
                        SMD_Before = smds_before,
                        SMD_After = smds_after)

print(smd_table)

# 3. Plotting Weighted Propensity Score Distribution
dataset$weights <- weights
p_weighted <- ggplot(dataset, aes(x = propensity_score, weight = weights, fill = factor(treatment))) +
  geom_density(alpha = 0.5) +
  labs(title = 'Weighted Propensity Score Distribution by Treatment Status',
       x = 'Propensity Score',
       y = 'Weighted Density')
print(p_weighted)

# 4. Visualize Weight Distribution
p_weights <- ggplot(dataset, aes(x = weights)) +
  geom_histogram(bins = 50, fill = 'blue', color = 'black') +
  labs(title = 'Distribution of Weights After Truncation',
       x = 'Weights',
       y = 'Frequency')
print(p_weights)

# checking different truncation thresholds
check_truncation <- function(threshold) {
  # Truncate weights
  weights_truncated <- weights
  weights_truncated[weights_truncated > quantile(weights_truncated, 1 - threshold)] <- quantile(weights_truncated, 1 - threshold)
  weights_truncated[weights_truncated < quantile(weights_truncated, threshold)] <- quantile(weights_truncated, threshold)
  
  # Calculate Effective Sample Size (ESS)
  ess <- sum(weights_truncated)^2 / sum(weights_truncated^2)
  
  # Calculate SMDs after truncation
  smds_truncated <- sapply(covariates, function(x) calculate_weighted_smd(x, dataset$Treatment, weights_truncated))
  
  # Store results
  list(
    ESS = ess,
    SMDs = smds_truncated
  )
}

# Test different thresholds
thresholds <- c(0.001, 0.005, 0.01, 0.02, 0.05)
results <- lapply(thresholds, check_truncation)

# Print results
for (i in 1:length(thresholds)) {
  cat("\nTruncation Threshold:", thresholds[i], "\n")
  cat("Effective Sample Size (ESS):", results[[i]]$ESS, "\n")
  print(data.frame(Variable = names(results[[i]]$SMDs), SMD_After_Truncation = results[[i]]$SMDs))
}
