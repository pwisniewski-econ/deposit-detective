# Logit, propensity score matching/IPW, Lasso
# imbalanced data ->  sensitivity analysis
# month fixed effects
# Collider post-treatment variables: duration
# Backdoor/confounders:
# variables on causal path:
# one hot encode variables for LASSO
# double post lasso
# proposensity matching - double post lasso
# ipw - double post lasso

# cv, train test, average treatment effect
# lasso otherwise overfit, lanbda use cv
# neyman ortho

# Install if not already installed
library(ROSE)

# keep only called twice
data_run = datafull %>% 
  filter(campaign %in% c(1, 2))

# generate treatment
data_run$treatment <- ifelse(data_run$day_of_week  %in% c("tue", "wed", "thu"), 1, 0)

# binary outcome
data_run$y <- ifelse(data_run$y == "yes", 1, 0)



# Applying ROSE to generate a balanced dataset
df_balanced <- ROSE(y ~ treatment +
                      campaign + age + job + day_of_week +
                      + duration + contact + loan + housing +
                      default + education +
                      previous + poutcome + 
                      emp.var.rate + cons.price.idx + cons.conf.idx + euribor3m +
                      date_month,
                      data = data_run, seed = 123)$data

# Checking new class distribution
table(df_balanced$Pass)

# PICO formulation:

* Population : Who are we interested in? - Customers from a Portuguese bank, who are called for the first time, as part of a telemarketing campaign to subscribe to a term deposit product
* Intervention : What treatment/intervention do we study? - Receiving a call during the middle of the week (Tue, Wed, Thur) campaign
* Comparison : What are we comparing it to? - Customers called on Monday or Friday 
* Outcome : What are we interested in? - Subscribing or not to the term deposit product
* Time : Outcome is determined after the call ends (directly after treatment)



# Run logistic regression
logit_model <- glm(y ~ treatment +
                     campaign + age + job + day_of_week +
                     + duration + contact + loan + housing +
                     default + education +
                     previous + poutcome + 
                     emp.var.rate + cons.price.idx + cons.conf.idx + euribor3m +
                     factor(date_month),
                   data = data_run, family = binomial)

# Print summary of the model
summary(logit_model)
exp(coef(logit_model))

pred_run <- predict(logit_model, data=data_run, type="response")
pred_run = data.frame(pred_run)

# Generate ROC Curve and calculate AUC
roc_obj <- roc(dataset$outcome, pred)
auc_value <- auc(roc_obj)

# Plotting the ROC Curve
ggroc(roc_obj) + 
  ggtitle(paste("ROC Curve (AUC =", round(auc_value, 3), ")")) +
  theme_minimal()
roc.curve(data_run$y, pred_run, 
          main="ROC curve \n (Half circle depleted data balanced by ROSE)")

