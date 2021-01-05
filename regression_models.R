

library(tidyverse)
library(xgboost)
library(e1071)
library(glmnet)

to_use <- NULL # set to NULL to use all features
clono_level <- F
unique_sequences <- T # nt level
encoding <- "onehot"

# Reading data ------------------------------------------------------------

if (clono_level) level <- "^clono" else level <- "^cell"

temp_files <- list()
for (file in list.files("features", level, full.names = T)) {
  temp_files[[sub(".*/", "", file)]] <- read.csv(file)
}

feature_matrix <- bind_rows(temp_files) %>% filter(!is.na(octet)) %>% select(-ELISA_bind)

if (!is.null(to_use)) {feature_matrix <- feature_matrix %>% select(all_of(append(to_use, "octet")))}
if (unique_sequences == T) {feature_matrix <- feature_matrix %>% distinct(sequence_HC, sequence_LC, .keep_all=T)}

# Label encoding ----------------------------------------------------------

if (encoding == "onehot") {
  f_temp <- feature_matrix %>% select_if(is.character) %>% select(-matches("clon|barcode|octet|chain_"))
  f_temp <- lapply(f_temp, function(x)
    data.frame(str_split_fixed(x, "", max(sapply(x, str_length))), stringsAsFactors = T))

  f_encoded <- lapply(f_temp, function(x) data.frame(sapply(x, as.numeric)))
  f_encoded <- cbind.data.frame(f_encoded)
  f_encoded <- f_encoded[vapply(f_encoded, function(x) length(unique(x)) > 1, logical(1L))]

  target <- feature_matrix %>% pull(octet)
}

if (encoding == "kmer") {

  f_temp <- feature_matrix %>% select_if(is.character) %>% select(-matches("clon|barcode|octet|chain_"))
  getKmers <- function(sequence, size=5) {for (i in (str_split(sequence, ""))) {i:i+size}}

  f_temp <- lapply(f_temp, getKmers)

  f_encoded <- lapply(f_temp, function(x) data.frame(sapply(x, as.numeric)))
  f_encoded <- cbind.data.frame(f_encoded)
  f_encoded <- f_encoded[vapply(f_encoded, function(x) length(unique(x)) > 1, logical(1L))]

  target <- feature_matrix %>% pull(octet)
}

# train test split --------------------------------------------------------

sample <- sample.int(n = nrow(f_encoded), size = floor(.66*nrow(f_encoded)), replace = F, useHash = F)

x_train <- f_encoded[sample, ]
x_test  <- f_encoded[-sample, ]
y_train <- target[sample]
y_test <- target[-sample]

# Model generation and training --------------------------------------------------------

# Linear regression
# lm <- lm(y_train ~., x_train)

# Ridge regression
lambdas <- 10^seq(4, -3, by = -.1)
cv_ridge = cv.glmnet(as.matrix(x_train), y_train, nlambda=25, alpha=0, family='gaussian', lambda=lambdas)
best_lambda <- cv_ridge$lambda.min
ridge <- glmnet(as.matrix(x_train), y_train, nlambda=25, alpha=0, family='gaussian', lambda=best_lambda)

# Setting alpha = 1 implements lasso regression
cv_lasso <- cv.glmnet(as.matrix(x_train), y_train, alpha = 1, lambda = lambdas, standardize = T, nfolds = 5)
best_lambda <- cv_lasso$lambda.min
lasso <- glmnet(as.matrix(x_train), y_train, alpha = 1, lambda = best_lambda, standardize = T)

# XGBRegressor
xgb <- xgboost(as.matrix(x_train), y_train, nrounds = 10, max.depth=3, eta=1)

importance_matrix <- xgb.importance(model = xgb)
print(importance_matrix)
xgb.plot.importance(importance_matrix = importance_matrix)

# Model comparison --------------------------------------------------------

# Compute R^2 from true and predicted values
eval_results <- function(true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))

  data.frame(RMSE = RMSE, Rsquare = R_square)
}

pred_ridge <- predict(ridge, as.matrix(x_test))[,1]
pred_lasso <- predict(lasso, as.matrix(x_test))[,1]
pred_xgb <- predict(xgb, as.matrix(x_test))

rmse_ridge <- eval_results(y_test, pred_ridge, x_test)
rmse_lasso <- eval_results(y_test, pred_lasso, x_test)
rmse_xgb <- eval_results(y_test, pred_xgb, x_test)

# Plots -------------------------------------------------------------------

plot_data <- gather(data.frame(ridge = pred_ridge, lasso = pred_lasso, xgb = pred_xgb), key = model, value = pred) %>%
  mutate(true = rep(y_test, 3))

ggplot(plot_data, aes(x=true, y=pred, color=model)) + geom_point() + geom_smooth(method="lm") +
  scale_y_log10(limits = c(1, ceiling(max(y_test)))) + scale_x_log10(limits = c(1, ceiling(max(y_test))))


