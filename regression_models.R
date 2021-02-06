#' Assings to the repertoire sequencing data a scalar value representing the predicted affinity for each cell/clonotype (the lower the value the higher the affinity). This can be done both at the clonotype level or at the single cell level.
#' @param features List of dataframes containing the extracted features. This is the output of load_data function.
#' @param unique.sequences Names of the sequences to be kept. Default is c("aa_sequence_HC", "aa_sequence_LC") which keeps every cell with a unique combination of heavy and light chain sequences. Other options include "clonotype_id", "cdr3s_aa", "aa_sequence_HC".
#' @param encoding Character indicating which encoding strategy to use. Options are "onehot", "kmer". Default is set to "onehot".
#' @param to.use Character vector indicating which features to use. If not supplied all the features will be used
#' @return This function plots the predicted affinity scores against the actual affinities for each model and plots feature importance for XGBoost.
#' @export
#' @examples
#' \dontrun{
#' check_predict_affinity <- predict_affinity(features = output.load_data, unique.sequences = c("aa_sequence_HC", "aa_sequence_LC"), encoding = "onehot", to.use = NULL)
#' }
#'

predict_affinity <- function(features, unique.sequences, encoding, to.use) {

  require(tidyverse)
  require(xgboost)
  require(e1071)
  require(protr)
  require(glmnet)

  if (missing(unique.sequences)) unique.sequences <- c("aa_sequence_HC", "aa_sequence_LC")
  if (missing(encoding)) encoding <- "onehot"

  f_encoded <- Bindpred::encode_features(features, encoding, unique.sequences, to.use) %>% drop_na(octet)
  target <- f_encoded %>% pull(octet)
  f_encoded <- f_encoded %>% select(-c(ELISA_bind, octet))

  # train test split --------------------------------------------------------
  sample <- sample.int(n = nrow(f_encoded), size = floor(0.7 * nrow(f_encoded)), replace = F, useHash = F)

  x_train <- f_encoded[sample, ]
  x_test  <- f_encoded[-sample, ]
  y_train <- target[sample]
  y_test <- target[-sample]

  # Model generation and training --------------------------------------------------------

  # Ridge regression
  lambdas = 10^seq(4, -3, by = -.1)
  cv_ridge = cv.glmnet(as.matrix(x_train), y_train, nlambda=25, alpha=0, family="gaussian", lambda=lambdas)
  best_lambda <- cv_ridge$lambda.min
  ridge <- glmnet(as.matrix(x_train), y_train, nlambda=25, alpha=0, family="gaussian", lambda=best_lambda)

  # Setting alpha = 1 implements lasso regression
  cv_lasso <- cv.glmnet(as.matrix(x_train), y_train, alpha = 1, lambda = lambdas, standardize = T, nfolds = 5)
  best_lambda <- cv_lasso$lambda.min
  lasso <- glmnet(as.matrix(x_train), y_train, alpha = 1, lambda = best_lambda, standardize = T)

  # XGBRegressor
  xgb <- xgboost(as.matrix(x_train), as.numeric(as.factor(y_train))-1, nrounds = 10, max.depth = 8, eta=1)

  importance_matrix <- xgb.importance(model = xgb)
  importance_plot <- xgb.ggplot.importance(importance_matrix = importance_matrix, top_n = 25)

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
    mutate(true = rep(y_test, 3)) %>%
      mutate(RMSE = c(rep(rmse_ridge$RMSE, length(y_test)),
                      rep(rmse_lasso$RMSE, length(y_test)),
                      rep(rmse_xgb$RMSE, length(y_test))))

  # Log scale
  output.plot <- list()

  output.plot[["reg"]] <- ggplot(plot_data, aes(x=true, y=predicted, color=model)) +
      geom_point() + geom_smooth(method="lm", aes(fill=model), alpha=0.2) +
         scale_y_log10(limits = c(1, ceiling(max(y_test)))) + scale_x_log10(limits = c(1, ceiling(max(y_test)))) +
            geom_abline(slope=1, linetype = "dashed", color = "grey") +
                ggtitle(paste0("encoding: ", encoding)) + labs(linetype = "model")

  output.plot[["importance"]] <- importance_plot

  return(output.plot)

}

