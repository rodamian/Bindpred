#' Classifies the repertoire sequencing data binarily (binding / not binding). This can be done both at the clonotype level or at the single cell level.
#' @param features List of dataframes containing the extracted features. This is the output of load_data function.
#' @param unique.sequences Names of the sequences to be kept. Default is c("aa_sequence_HC", "aa_sequence_LC") which keeps every cell with a unique combination of heavy and light chain sequences. Other options include "clonotype_id", "cdr3s_aa", "aa_sequence_HC".
#' @param encoding Character indicating which encoding strategy to use. Options are "onehot", "kmer", "protr". To set the kmer size set encoding = "5mer" for size 5. If only "kmer" is given the default size is 3. The default overall is set to "onehot".
#' @param to.use Character vector indicating which features to use. If not supplied all the features will be used
#' @param cv Numeric indicating the number of folds used in cross validation. Default is 5.
#' @return This function plots AUC scores for each model, feature importance for XGBoost and return the predicted labels (Binding / not Binding)
#' @export
#' @examples
#' \dontrun{
#' check_classify_data <- classify_data(features = output.load_data, to.use = NULL, unique.sequences = c("aa_sequence_HC", "aa_sequence_LC"), encoding = "onehot")
#' }
#'

classify_data <- function(features,
                          unique.sequences = "cdr3s_aa",
                          encoding = "onehot",
                          to.use = c("cdr3s_aa", "cdr3s_nt", "aa_sequence_HC", "aa_sequence_LC"),
                          cv = 5) {

  require(tidyverse)
  require(fastNaiveBayes)
  require(xgboost)
  require(e1071)
  require(pROC)
  require(caret)
  require(mboost)
  require(Ckmeans.1d.dp)
  theme_set(theme_bw())

  f_encoded <- encode_features(features, encoding, unique.sequences, to.use)
  f_encoded <- f_encoded %>% filter(ELISA_bind %in% c("yes", "no")) %>% dplyr::select(-octet)

  # Setting target variable to 0,1
  f_encoded$ELISA_bind <- as.numeric(as.factor(f_encoded$ELISA_bind)) - 1

  folds <- createFolds(f_encoded$ELISA_bind, k = cv)

  # Cross validation
  crossv <- lapply(folds, function(n) {

    # train test split --------------------------------------------------------
    x_train <- as.matrix(f_encoded[-n,] %>% dplyr::select(-ELISA_bind))
    x_test  <- as.matrix(f_encoded[n,] %>% dplyr::select(-ELISA_bind))
    y_train <- f_encoded[-n,]$ELISA_bind
    y_test <- f_encoded[n,]$ELISA_bind

    # Model generation and training -------------------------------------------
    models <- list()

    # models[["gnb"]] <- fastNaiveBayes(
    #                         x_train,
    #                         y_train)

    models[["svm"]] <- e1071::svm(formula = y_train ~ .,
                            data = x_train,
                            kernel = 'radial',
                            probability = TRUE)

    models[["xgb"]] <- xgb.train(params = list(
                            booster = "gbtree",
                            objective = "binary:logistic",
                            eval_metric = "auc",
                            max_depth = 8),
                            data = xgb.DMatrix(
                                label = y_train,
                                data = x_train),
                            nrounds = 25)

    models[["glm"]] <- glmboost(x_train,
                            y_train,
                            family = Gaussian(),
                            control = boost_control())


    # Calculating probabilities
    probs <- data.frame(sapply(models, predict, newdata = x_test))
    colnames(probs) <- names(models)
    rownames(probs) <- n

    return(probs)
  })

  # Model evaluation --------------------------------------------------------
  probs <- crossv %>% purrr::reduce(rbind)
  probs <- probs[order(as.numeric(rownames(probs))),]
  rocs <- lapply(probs, function(x) roc(f_encoded$ELISA_bind, x, ci = T))

  # Confidence intervals
  # ci.sp <- lapply(rocs, ci.sp, sensitivities=seq(0, 1, .01), boot.n=100)
  conf <- lapply(rocs, ci)

  # Model comparison --------------------------------------------------------
  output.plot <- list()

  output.plot[["roc"]] <- ggroc(rocs) +
    geom_abline(slope=1, linetype = "dashed", color = "grey", intercept = 1) +
      ggtitle(paste0("encoding: ", encoding)) + labs(color = "model") +
        annotate(geom = "text", x = 0.25, y = seq(0.1, 0.4, length.out = length(rocs)),
          label = paste("auc", names(rocs), ":", sapply(rocs, function(x) as.character(round(x$auc, 3)))))

  # importance_matrix <- xgb.importance(model = models[["xgb"]])
  # output.plot[["importance"]] <- xgb.ggplot.importance(importance_matrix = importance_matrix, top_n = 25, measure = "Gain")

  return(output.plot)
}
