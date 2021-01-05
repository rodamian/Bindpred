

library(tidyverse)
library(naivebayes)
library(xgboost)
library(e1071)
library(pROC)
library(parallel)
if(.Platform$OS.type == "unix")
  {options(mc.cores=detectCores()-1, print("Using parallel package"))} else {options(mc.cores=1)}

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

feature_matrix <- bind_rows(temp_files) %>% filter(!is.na(ELISA_bind)) %>%
  filter(ELISA_bind %in% c("yes", "no")) %>% select(-octet)

if (!is.null(to_use)) {feature_matrix <- feature_matrix %>% select(all_of(append(to_use, "ELISA_bind")))}
if (unique_sequences == T) {feature_matrix <- feature_matrix %>% distinct(clonotype_id, .keep_all = T)}

# Label encoding ----------------------------------------------------------

if (encoding == "onehot") {
  f_temp <- feature_matrix %>% select_if(is.character) %>% select(-matches("clon|barcode|ELISA|chain_"))
  f_temp <- lapply(f_temp, function(x)
    data.frame(str_split_fixed(x, "", max(sapply(x, str_length))), stringsAsFactors = T))

  f_encoded <- lapply(f_temp, function(x) data.frame(sapply(x, as.numeric)))
  f_encoded <- cbind.data.frame(f_encoded)
  f_encoded <- f_encoded[vapply(f_encoded, function(x) length(unique(x)) > 1, logical(1L))]

  target <- as.factor(feature_matrix %>% pull(ELISA_bind))
}

if (encoding == "kmer") {
  f_temp <- feature_matrix %>% select_if(is.character) %>% select(-matches("clon|barcode|ELISA|chain_"))
  f_temp <- lapply(f_temp, function(x)
    data.frame(str_split_fixed(x, "", max(sapply(x, str_length))), stringsAsFactors = T))

  f_encoded <- lapply(f_temp, function(x) data.frame(sapply(x, as.numeric)))
  f_encoded <- cbind.data.frame(f_encoded)
  f_encoded <- f_encoded[vapply(f_encoded, function(x) length(unique(x)) > 1, logical(1L))]

  target <- as.factor(feature_matrix %>% pull(ELISA_bind))
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


