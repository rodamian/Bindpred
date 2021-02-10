#' Classifies the repertoire sequencing data binarily (binding / not binding). This can be done both at the clonotype level or at the single cell level.
#' @param features List of dataframes containing the extracted features. This is the output of load_data function.
#' @param to.use Character vector indicating which features to use. If not supplied all the features will be used
#' @param unique.sequences Names of the sequences to be kept. Default is c("aa_sequence_HC", "aa_sequence_LC") which keeps every cell with a unique combination of heavy and light chain sequences. Other options include "clonotype_id", "cdr3s_aa", "aa_sequence_HC".
#' @param encoding Character indicating which encoding strategy to use. Options are "onehot", "kmer", "protr". To set the kmer size set encoding = "5mer" for size 5. If only "kmer" is given the default size is 3. The default overall is set to "onehot".
#' @return This function plots AUC scores for each model, feature importance for XGBoost and return the predicted labels (Binding / not Binding)
#' @export
#' @examples
#' \dontrun{
#' check_classify_data <- classify_data(features = output.load_data, to.use = NULL, unique.sequences = c("aa_sequence_HC", "aa_sequence_LC"), encoding = "onehot")
#' }
#'

classify_data <- function(features, to.use, unique.sequences, encoding) {

  require(tidyverse)
  require(naivebayes)
  require(xgboost)
  require(e1071)
  require(pROC)
  library(kernlab)
  library(speedglm)
  library(caret)
  require(Ckmeans.1d.dp)

  theme_set(theme_bw())
  if (missing(unique.sequences)) unique.sequences <- c("aa_sequence_HC", "aa_sequence_LC")
  if (missing(encoding)) encoding <- "onehot"

  f_encoded <- encode_features(features, encoding, unique.sequences, to.use)
  target <- f_encoded %>% pull(ELISA_bind) %in% c("yes", "no") %>% as.factor
  f_encoded <- f_encoded %>% dplyr::select(-c(ELISA_bind, octet))

  # train test split --------------------------------------------------------
  train_id <- sample(1:nrow(f_encoded), size = floor(0.75 * nrow(f_encoded)), replace=F)

  x_train <- f_encoded[train_id, ]
  x_test  <- f_encoded[-train_id, ]
  y_train <- target[train_id]
  y_test <- target[-train_id]

  # Model generation and training -------------------------------------------
  models <- list()

  models[["gnb"]] <- gaussian_naive_bayes(as.matrix(x_train), y_train)

  # models[["svm"]] <- ksvm(y_train, data=x_train,
  #             kernel = "rbfdot",
  #             kpar = list(sigma=0.015),
  #             C = 70,
  #             cross = 5,
  #             prob.model = TRUE)
  #
  # models[["svm"]] <- svm_model <- e1071::svm(
  #             x_train,
  #             as.numeric(y_train),
  #             kernel = "linear",
  #             cost = 10,
  #             probability = TRUE)

  # models[["svm"]] <- svm(x_train, y_train, kernel = "radial", probability = TRUE)

  models[["xgb"]] <- xgboost(as.matrix(x_train), label = as.numeric(y_train)-1, nround = 10,
              params = list(
              booster = "gbtree",
              max_depth = 8,
              lambda = 0.5,
              objective = "binary:logistic",
              eval_metric = "auc"))

  # models[["glm"]] <- speedglm(y_train ~., family = binomial(link='logit'), data = x_train, control = list(maxit = 5))

  # Calculating probabilities
  probs <- lapply(models, predict, as.matrix(x_test), type = "prob")

  # Model evaluation --------------------------------------------------------
  ROC <- lapply(probs, function(x) if (NCOL(x) == 2) roc(y_test, x[,2]) else roc(y_test, x))

  # Model comparison --------------------------------------------------------
  output.plot <- list()

  output.plot[["roc"]] <- ggroc(ROC) +
    geom_abline(slope=1, linetype = "dashed", color = "grey", intercept = 1) +
      ggtitle(paste0("encoding: ", encoding)) + labs(color = "model") +
        annotate(geom = "text", x = 0.25, y = seq(0.1, 0.4, length.out = length(ROC)),
          label = paste(names(ROC), "AUC: ", sapply(ROC, function(x) as.character(round(x$auc, 3)))))

  importance_matrix <- xgb.importance(model = models[["xgb"]])
  output.plot[["importance"]] <- xgb.ggplot.importance(importance_matrix = importance_matrix, top_n = 25, measure = "Gain")


  return(output.plot)
}
