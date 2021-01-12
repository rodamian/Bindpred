#' Assings to the repertoire sequencing data a scalar value representing the predicted affinity for each cell/clonotype (the lower the value the higher the affinity). This can be done both at the clonotype level or at the single cell level.
#' @param features List of dataframes containing the extracted features. This is the output of load_data function.
#' @param to.use Character vector indicating which features to use. If not supplied all the features will be used
#' @param clono.level Logical indicating if the classified sequences are at the single cell level or at the clonotype level. Default is set to FALSE
#' @param unique.sequences Logical indicating if only unique aa sequences are to be considered. If all the sequences in one clonotype are identical then all of them will be predicted to have the same affinity.
#' @param encoding Character indicating which encoding strategy to use. Options are "onehot", "kmer". Default is set to "onehot".
#' @return This function plots the predicted affinity scores against the actual affinities for each model and plots feature importance for XGBoost.
#' @export
#' @examples
#' \dontrun{
#' check_predict_affinity <- predict_affinity(features = output.load_data, to.use = FALSE, clono.level = FALSE, unique.sequences = FALSE, encoding = "onehot")
#' }
#'

predict_affinity <- function(features, to.use, clono.level, unique.sequences, encoding) {

  require(tidyverse)
  require(xgboost)
  require(e1071)
  require(glmnet)

  if (missing(to.use)) to.use <- NULL # set to NULL to use all features
  if (missing(clono.level)) clono.level <- F
  if (missing(unique.sequences)) unique.sequences <- T # aa level
  if (missing(encoding)) encoding <- "onehot"

  # Reading data ------------------------------------------------------------

  features <- features %>% bind_rows %>% select(-ELISA_bind) %>% filter(!is.na(octet))

  if (!is.null(to.use)) {features <- features %>% select(all_of(append(to_use, "octet")))}
  if (unique.sequences == T) {features <- features %>% distinct(cdr3s_aa, .keep_all = T)}

  # Label encoding ----------------------------------------------------------

  if (encoding == "onehot") {
    f_temp <- features %>% select_if(is.character) %>% select(-matches("clon|barcode|ELISA_bind|chain_"))
    f_temp <- lapply(f_temp, function(x)
      data.frame(str_split_fixed(x, "", max(sapply(x, str_length))), stringsAsFactors = T))

    f_encoded <- lapply(f_temp, function(x) data.frame(sapply(x, as.numeric)))
    f_encoded <- cbind.data.frame(f_encoded)
    f_encoded <- f_encoded[vapply(f_encoded, function(x) length(unique(x)) > 1, logical(1L))]

    target <- features %>% pull(octet)
  }

  if (encoding == "kmer") {

    f_temp <- features %>% select_if(is.character) %>% select(-matches("clon|barcode|ELISA_bind|chain_"))
    getKmers <- function(sequence, size=5) {for (i in (str_split(sequence, ""))) {i:i+size}}

    f_temp <- lapply(f_temp, getKmers)

    f_encoded <- lapply(f_temp, function(x) data.frame(sapply(x, as.numeric)))
    f_encoded <- cbind.data.frame(f_encoded)
    f_encoded <- f_encoded[vapply(f_encoded, function(x) length(unique(x)) > 1, logical(1L))]

    target <- features %>% pull(octet)
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

  plot_data <- gather(data.frame(ridge = pred_ridge, lasso = pred_lasso, xgb = pred_xgb), key = model, value = predicted) %>%
    mutate(true = rep(y_test, 3))

  ggplot(plot_data, aes(x=true, y=predicted, color=model)) + geom_point() + geom_smooth(method="lm") +
    scale_y_log10(limits = c(1, ceiling(max(y_test)))) + scale_x_log10(limits = c(1, ceiling(max(y_test)))) +
      geom_abline(slope=1, linetype = "dashed", color = "grey")

}

predict_affinity(features, unique.sequences = T)
