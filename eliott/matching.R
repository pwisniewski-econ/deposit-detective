# Load libraries
library(MatchIt)
library(cobalt)
library(dplyr)
library(survey)


# Assume 'data' exists with Treatment (binary), Y (binary), and covariates
# Ensure Treatment is binary 0/1 and Y is factor with valid levels
data = data_numeric

# -----------------------------
# 1. PROPENSITY SCORE MATCHING
# -----------------------------

# Estimate propensity scores and perform matching
psm_model <- matchit(treatment ~ . -outcome, data = data, method = "nearest", distance = "logit")

# Check covariate balance
summary(psm_model)
# Optional visual
love.plot(psm_model, binary = "std")

# Extract matched data
matched_data <- match.data(psm_model)

# Estimate causal effect (e.g., difference in outcome)
psm_effect <- glm(outcome ~ treatment, data = matched_data, family = binomial())
summary(psm_effect)


# -----------------------------
# 2. INVERSE PROBABILITY WEIGHTING (IPW)
# -----------------------------

# Estimate propensity scores
ps_model <- glm(treatment ~ . -outcome, data = data, family = binomial())
data$pscore <- predict(ps_model, type = "response")

# Compute IPW weights
# Stabilized weights (optional, more stable)
data$weight_ipw <- with(data, ifelse(treatment == 1, 
                                     mean(treatment) / pscore, 
                                     (1 - mean(treatment)) / (1 - pscore)))

# Set up survey design for weighted estimation
ipw_design <- svydesign(ids = ~1, weights = ~weight_ipw, data = data)

# Estimate causal effect using weighted logistic regression
ipw_effect <- svyglm(outcome ~ treatment, design = ipw_design, family = quasibinomial())
summary(ipw_effect)
