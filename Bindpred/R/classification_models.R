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
  require(Ckmeans.1d.dp)
  require(pROC)

  if (missing(unique.sequences)) unique.sequences <- c("aa_sequence_HC", "aa_sequence_LC")
  if (missing(encoding)) encoding <- "protr.cdr3"

  f_encoded <- encode_features(features, encoding, unique.sequences)
  target <- f_encoded %>% pull(ELISA_bind) %in% c("yes", "no") %>% as.factor
  f_encoded <- f_encoded %>% select(-c(ELISA_bind, octet))

  # train test split --------------------------------------------------------
  train_id <- sample(1:nrow(f_encoded), size = floor(0.75 * nrow(f_encoded)), replace=F)

  x_train <- f_encoded[train_id, ]
  x_test  <- f_encoded[-train_id, ]
  y_train <- target[train_id]
  y_test <- target[-train_id]

  # Model generation and training -------------------------------------------
  gnb <- gaussian_naive_bayes(as.matrix(x_train), y_train)

  svm <- svm(x_train, y_train, kernel = "radial", probability = T,
             scale = vapply(f_encoded, function(x) length(unique(x)) > 5, logical(1L)))

  xgb <- xgboost(as.matrix(x_train), label = y_train, nround = 10, params = list(
    booster = "gbtree",
    max_depth = 5,
    lambda = 0.5,
    objective = "binary:logistic",
    eval_metric = "auc"))

  # glm <- glm(y_train ~., family = binomial(link='logit'), data = as.data.frame(x_train), control = list(maxit = 50))

  # Calculating probabilities
  probs_gnb <- gnb %prob% as.matrix(x_test)
  probs_xgb <- predict(xgb, as.matrix(x_test))
  probs_svm <- predict(svm, as.matrix(x_test), probability = T)
  # probs_glm <- predict(glm, as.data.frame(x_test), probability = T)

  # Model evaluation --------------------------------------------------------
  roc_xgb <- roc(y_test, probs_xgb)
  roc_gnb <- roc(y_test, probs_gnb[,2])
  roc_svm <- roc(y_test, attr(probs_svm, "probabilities")[,1], levels=c("no", "yes"))
  # roc_glm <- roc(y_test, probs_glm, levels=c("no", "yes"))

  # Model comparison --------------------------------------------------------
  roc_plot <- plot(roc_xgb, col = "blue", print.auc=T, print.auc.y=0.9, print.auc.x=0.1)
  roc_plot <- plot(roc_gnb, col = "green", add=T, print.auc=T, print.auc.y=.8, print.auc.x=0.1)
  roc_plot <- plot(roc_svm, col = "red", add=T, print.auc=T, print.auc.y=.7, print.auc.x=0.1)
  # roc_plot <- plot(roc_glm, col = "grey", add=T, print.auc=T, print.auc.y=.6, print.auc.x=0.1)

  roc_plot <- legend("bottomright", legend=c("XGB", "NB", "SVM", "LogReg"), col=c("blue", "green", "red", "grey"), lwd=2)

  importance_matrix <- xgb.importance(model = xgb)
  importance_plot <- xgb.ggplot.importance(importance_matrix = importance_matrix, top_n = 25)

  return(list(roc_plot, importance_plot))
}
