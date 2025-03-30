library(caret)
library(pROC)
library(dplyr)
library(ggplot2)

data = data_numeric

# Parameters for tuning
num_trials <- 50  # Reduced to 20 trials for testing
inner_folds <- 5  # Inner cross-validation folds
outer_folds <- 5  # Outer cross-validation folds

# Store results
results <- data.frame(Weight = numeric(), Threshold = numeric(),
                      AUC = numeric(), Accuracy = numeric(),
                      F1 = numeric(), Precision = numeric())

# Randomly sample weights and thresholds
set.seed(42)
sampled_params <- data.frame(
  Weight = runif(num_trials, 2, 8),  # Changed weight range to [2, 8]
  Threshold = runif(num_trials, 0.2, 0.8)
)

# Inner Cross-Validation Tuning
for (i in 1:num_trials) {
  weight <- sampled_params$Weight[i]
  threshold <- sampled_params$Threshold[i]
  
  # Train model using glm
  model <- glm(outcome ~ treatment + ., data = data, family = binomial(),
               weights = ifelse(outcome == "Yes", weight, 1))
  
  # Get predictions
  pred_probs <- predict(model, data, type = "response")
  predictions <- ifelse(pred_probs >= threshold, "1", "0")
  
  # Calculate metrics
  auc <- roc(data$outcome, pred_probs)$auc
  cm <- confusionMatrix(as.factor(predictions), as.factor(data$outcome), positive = "1")
  accuracy <- cm$overall['Accuracy']
  precision <- cm$byClass['Pos Pred Value']
  f1 <- cm$byClass['F1']
  
  # Store results
  results <- rbind(results, data.frame(Weight = weight, Threshold = threshold, AUC = auc,
                                       Accuracy = accuracy, F1 = f1, Precision = precision))
}

# Display results table
print(results)


# Select top 5 models based on AUC
best_models <- results %>% arrange(desc(AUC)) %>% head(5)

# Outer Cross-Validation
outer_cv_results <- data.frame()
folds <- createFolds(data$outcome, k = outer_folds, list = TRUE)

for (j in 1:outer_folds) {
  training_data <- data[-folds[[j]], ]
  testing_data <- data[folds[[j]], ]
  
  for (i in 1:nrow(best_models)) {
    weight <- best_models$Weight[i]
    threshold <- best_models$Threshold[i]
    
    # Train model using glm
    model <- glm(outcome ~ treatment + ., data = training_data, family = binomial(),
                 weights = ifelse(training_data$outcome == "Yes", weight, 1))
    
    # Get predictions
    pred_probs <- predict(model, testing_data, type = "response")
    predictions <- ifelse(pred_probs >= threshold, "1", "0")
    
    # Calculate metrics
    auc <- roc(as.numeric(testing_data$outcome) - 1, pred_probs)$auc
    cm <- confusionMatrix(as.factor(predictions), as.factor(testing_data$outcome), positive = "1")
    accuracy <- cm$overall['Accuracy']
    
    # Store results
    outer_cv_results <- rbind(outer_cv_results, data.frame(Model = paste0("(Weight=", round(weight, 2), ", Threshold=", round(threshold, 2), ")"),
                                                           Fold = j, AUC = auc, Accuracy = accuracy))
  }
}

# Plot AUC Boxplot
outer_cv_results$Model <- factor(outer_cv_results$Model)  # Convert Model to factor

ggplot(outer_cv_results, aes(x = Model, y = AUC)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "AUC Boxplot for Top 5 Models", x = "Model (Weight, Threshold)", y = "AUC") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot Accuracy Boxplot
ggplot(outer_cv_results, aes(x = Model, y = Accuracy)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Accuracy Boxplot for Top 5 Models", x = "Model (Weight, Threshold)", y = "Accuracy") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
