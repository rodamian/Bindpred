#' Classifies the repertoire sequencing data binarily (binding / not binding). This can be done both at the clonotype level or at the single cell level.
#' @param features List of dataframes containing the extracted features. This is the output of load_data function.
#' @param to.use Character vector indicating which features to use. If not supplied all the features will be used
#' @param clono.level Logical indicating if the classified sequences are at the single cell level or at the clonotype level. Default is set to FALSE
#' @param unique.sequences Logical indicating if only unique aa sequences are to be considered. If all the sequences in one clonotype are identical then all of them will be predicted with the same label.
#' @param encoding Character indicating which encoding strategy to use. Options are "onehot", "kmer". Default is set to "onehot".
#' @return This function plots AUC scores for each model, feature importance for XGBoost and return the predicted labels (Binding / not Binding)
#' @export
#' @examples
#' \dontrun{
#' check_classify_data <- classify_data(features = output.load_data, to.use = FALSE, clono.level = FALSE, unique.sequences = FALSE, encoding = "onehot")
#' }
#'

classify_data <- function(features, to.use, clono.level, unique.sequences, encoding) {

  require(tidyverse)
  require(naivebayes)
  require(xgboost)
  require(e1071)
  require(pROC)

  if (missing(to.use)) to.use <- NULL # set to NULL to use all features
  if (missing(clono.level)) clono.level <- F
  if (missing(unique.sequences)) unique.sequences <- T # aa level
  if (missing(encoding)) encoding <- "onehot"

  # Reading data ------------------------------------------------------------

  features <- features %>% bind_rows %>% filter(ELISA_bind %in% c("yes", "no")) %>% select(-octet)

  if (!is.null(to.use)) {features <- features %>% select(all_of(append(to_use, "ELISA_bind")))}
  if (unique.sequences == T) {features <- features %>% distinct(cdr3s_aa, .keep_all = T)}

  # Label encoding ----------------------------------------------------------

  if (encoding == "onehot") {
    f_temp <- features %>% select_if(is.character) %>% select(-matches("clon|barcode|ELISA|chain_"))
    f_temp <- lapply(f_temp, function(x)
      data.frame(str_split_fixed(x, "", max(sapply(x, str_length))), stringsAsFactors = T))

    f_encoded <- lapply(f_temp, function(x) data.frame(sapply(x, as.numeric)))
    f_encoded <- cbind.data.frame(f_encoded)
    f_encoded <- f_encoded[vapply(f_encoded, function(x) length(unique(x)) > 1, logical(1L))]

    target <- as.factor(features %>% pull(ELISA_bind))
  }

  if (encoding == "kmer") {
    f_temp <- features %>% select_if(is.character) %>% select(-matches("clon|barcode|ELISA|chain_"))
    f_temp <- lapply(f_temp, function(x)
      data.frame(str_split_fixed(x, "", max(sapply(x, str_length))), stringsAsFactors = T))

    f_encoded <- lapply(f_temp, function(x) data.frame(sapply(x, as.numeric)))
    f_encoded <- cbind.data.frame(f_encoded)
    f_encoded <- f_encoded[vapply(f_encoded, function(x) length(unique(x)) > 1, logical(1L))]

    target <- as.factor(features %>% pull(ELISA_bind))
  }

  # train test split --------------------------------------------------------

  sample <- sample.int(n = nrow(f_encoded), size = floor(.75*nrow(f_encoded)), replace = F, useHash = F)

  x_train <- f_encoded[sample, ]
  x_test  <- f_encoded[-sample, ]
  y_train <- target[sample]
  y_test <- target[-sample]

  # Model generation and training -------------------------------------------

  gnb <- gaussian_naive_bayes(as.matrix(x_train), y_train)

  svm <- svm(x_train, y_train, kernel = "radial", probability = T,
             scale = vapply(f_encoded, function(x) length(unique(x)) > 5, logical(1L)))

  xgb <- xgboost(as.matrix(x_train), label = as.numeric(factor(y_train))-1, nround = 10, params = list(
             booster = "gbtree", max_depth = 3, lambda = 0.5, objective = "binary:logistic", eval_metric = "auc"))

  glm <- glm(y_train ~., family = binomial(link='logit'), data = x_train)

  # probs_knn <- knn(x_train, x_test, y_train, prob=T)
  probs_gnb <- gnb %prob% as.matrix(x_test)
  probs_xgb <- predict(xgb, as.matrix(x_test))
  probs_svm <- predict(svm, as.matrix(x_test), probability = T)
  probs_glm <- predict(glm, x_test, probability = T)

  # Model evaluation --------------------------------------------------------

  roc_xgb <- roc(y_test, probs_xgb, levels=c("no", "yes"))
  roc_gnb <- roc(y_test, probs_gnb[,2], levels=c("no", "yes"))
  roc_svm <- roc(y_test, attr(probs_svm, "probabilities")[,1], levels=c("no", "yes"))
  roc_glm <- roc(y_test, probs_glm, levels=c("no", "yes"))
  # roc_knn <- roc(y_test, attr(probs_knn, "probabilities")[,1])

  # Model comparison --------------------------------------------------------

  roc_plot <- plot(roc_xgb, col = "blue", print.auc=T, print.auc.y=0.9, print.auc.x=0.1)
  roc_plot <- plot(roc_gnb, col = "green", add=T, print.auc=T, print.auc.y=.8, print.auc.x=0.1)
  roc_plot <- plot(roc_svm, col = "red", add=T, print.auc=T, print.auc.y=.7, print.auc.x=0.1)
  roc_plot <- plot(roc_glm, col = "grey", add=T, print.auc=T, print.auc.y=.6, print.auc.x=0.1)

  legend("bottomright", legend=c("XGB", "NB", "SVM", "GLM"), col=c("blue", "green", "red", "grey"), lwd=2)

  importance_matrix <- xgb.importance(model = xgb)
  print(head(importance_matrix))
  xgb.plot.importance(importance_matrix = importance_matrix)

}


classify_data(features)
