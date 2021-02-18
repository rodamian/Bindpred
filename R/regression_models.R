#' Assings to the repertoire sequencing data a scalar value representing the predicted affinity for each cell/clonotype (the lower the value the higher the affinity). This can be done both at the clonotype level or at the single cell level.
#' @param features List of dataframes containing the extracted features. This is the output of load_data function.
#' @param unique.sequences Names of the sequences to be kept. Default is c("aa_sequence_HC", "aa_sequence_LC") which keeps every cell with a unique combination of heavy and light chain sequences. Other options include "clonotype_id", "cdr3s_aa", "aa_sequence_HC".
#' @param encoding Character indicating which encoding strategy to use. Options are "onehot", "kmer". Default is set to "onehot".
#' @param to.use Character vector indicating which features to use. If not supplied all the features will be used
#' @param cv Numeric indicating the number of folds used in cross validation. Default is 5.
#' @return This function plots the predicted affinity scores against the actual affinities for each model and plots feature importance for XGBoost.
#' @export
#' @examples
#' \dontrun{
#' check_predict_affinity <- predict_affinity(features = output.load_data, unique.sequences = c("aa_sequence_HC", "aa_sequence_LC"), encoding = "onehot", to.use = NULL)
#' }
#'

predict_affinity <- function(features,
                             unique.sequences = c("aa_sequence_HC", "aa_sequence_LC"),
                             encoding = "onehot",
                             to.use = c("cdr3s_aa", "cdr3s_nt", "aa_sequence_HC", "aa_sequence_LC"),
                             cv = 5) {

  require(tidyverse)
  require(xgboost)
  require(protr)
  require(glmnet)
  theme_set(theme_bw())
  update_geom_defaults("point", list(size = 0.7))
  update_geom_defaults("smooth", list(alpha = 0.25, size = 0.7))

  f_encoded <- encode_features(features, encoding, unique.sequences, to.use)
  f_encoded <- f_encoded %>% drop_na(octet) %>% dplyr::select(-ELISA_bind)

  # Cross validation
  folds <- createFolds(f_encoded$octet, k = cv)

  crossv <- lapply(folds, function(n) {

    # train test split --------------------------------------------------------
    x_train <- as.matrix(f_encoded[-n,] %>% dplyr::select(-octet))
    x_test  <- as.matrix(f_encoded[n,] %>% dplyr::select(-octet))
    y_train <- f_encoded[-n,]$octet
    y_test <- f_encoded[n,]$octet

    # Model generation and training -------------------------------------------
    models <- list()

    # Ridge regression
    lambdas = 10^seq(4, -3, by = -.1)
    cv_ridge = cv.glmnet(as.matrix(x_train), y_train, nlambda=25, alpha=0, family="gaussian", lambda=lambdas)
    best_lambda <- cv_ridge$lambda.min
    models[["ridge"]] <- glmnet(as.matrix(x_train), y_train, nlambda=25, alpha=0, family="gaussian", lambda=best_lambda)

    # Setting alpha = 1 implements lasso regression
    cv_lasso <- cv.glmnet(as.matrix(x_train), y_train, alpha = 1, lambda = lambdas, standardize = T, nfolds = 5)
    best_lambda <- cv_lasso$lambda.min
    models[["lasso"]] <- glmnet(as.matrix(x_train), y_train, alpha = 1, lambda = best_lambda, standardize = T)

    # XGBRegressor
    models[["xgb"]] <- xgboost(as.matrix(x_train), y_train,
                               booster = "gbtree",
                               lambda = 0.5,
                               objective = "reg:tweedie",
                               eval_metric = "rmse",
                               nrounds = 14,
                               eta=1)

    # Prediction
    preds <- data.frame(sapply(models, predict, x_test))
    colnames(preds) <- names(models)
    rownames(preds) <- n

    return(preds)
  })

  # Model evaluation --------------------------------------------------------
  preds <- crossv %>% purrr::reduce(rbind)
  preds <- preds[order(as.numeric(rownames(preds))),]

  # Model comparison --------------------------------------------------------

  # Compute R squared from true and predicted values
  RMSE <- function (pred, obs) round(sqrt(mean((pred - obs)^2)), 3)

  rmse <- lapply(preds, RMSE, f_encoded$octet)

  # Plots -------------------------------------------------------------------

  plot_data <- as.data.frame(preds) %>% gather(key = model, value = predicted) %>%
    mutate(true = rep(f_encoded$octet, ncol(preds))) %>% mutate(model = paste(model, " RMSE: ", as.character(rmse[model])))

  # Log scale
  output.plot <- list()

  output.plot[["reg"]] <- ggplot(plot_data, aes(x=true, y=predicted, color=model)) +
      geom_point() + geom_smooth(method="lm", aes(fill=model), alpha=0.2) +
         scale_y_log10(limits = c(1, ceiling(max(f_encoded$octet)))) + scale_x_log10(limits = c(1, ceiling(max(f_encoded$octet)))) +
            geom_abline(slope=1, linetype = "dashed", color = "grey") +
                ggtitle(paste0("encoding: ", encoding)) + theme(legend.position="bottom")

  # importance_matrix <- xgb.importance(model = models[["xgb"]])
  # output.plot[["importance"]] <- xgb.ggplot.importance(importance_matrix = importance_matrix, top_n = 25)

  return(output.plot)
}

