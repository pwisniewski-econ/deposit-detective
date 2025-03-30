library(grf)

Y <- dataset$outcome
D <- dataset$treatment
X <- dataset 

str(X)
str(D)
str(Y)

# Create the numeric matrix (model matrix) for predictors
str(X)

X_numeric <- model.matrix(~ . - 1, data = X)  # "-1" removes the intercept column
X_numeric = data.frame(X_numeric)
X = X_numeric
str(X)

cf <- causal_forest(X, Y, D)
ate_overlap <- average_treatment_effect(cf, target.sample = "overlap")
print(ate_overlap)
