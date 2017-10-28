options(stringsAsFactors = F)


#' Caret ensemble list prediction
#'
#' This function takes the model list and runs the prediction for one model
#'
#' @param i model index
#' @param model_list caret list containing list of models
#' @param data data on which prediction should be performed
#'
#' @return Vector of predictions for p(y=1)
#'
#' @examples
#' pred_new(i = 1, model_list = caretlist, data = train_data)
#'
#' @export

pred_new <- function(i, model_list, data) {
  model <- model_list[[i]]
  train_pred <- predict(model, data, type = "prob")
  if(names(model_list)[i]=="bartMachine")
    train_pred <- 1 - train_pred
  return(train_pred[,2])
}

#' Binary log loss
#'
#' This function calculates binary logloss given the actual classes and
#'  predictions
#'
#' @param data Vector of observed classes
#' @param lev (Optional) Levels of the outcome variable
#' @param model Model predictions for p(y=1)
#'
#' @return List of accuracy, kappa and binary logloss
#'
#' @examples
#' LogLossSummary(data = c(0,1), model = c(0.1, 0.9))
#'
#' @export

LogLosSummary <- function (data, lev = NULL, model = NULL) {
  LogLoss <- function(actual, pred, eps = 1e-15) {
    stopifnot(all(dim(actual) == dim(pred)))
    pred[pred < eps] <- eps
    pred[pred > 1 - eps] <- 1 - eps
    return(-sum(actual * log(pred))/nrow(pred))
  }
  if (is.character(data$obs)) data$obs <- factor(data$obs, levels = lev)
  pred <- data[, "pred"]
  obs <- data[, "obs"]
  isNA <- is.na(pred)
  pred <- pred[!isNA]
  obs <- obs[!isNA]
  data <- data[!isNA, ]
  cls <- levels(obs)

  if (length(obs) + length(pred) == 0) {
    out <- rep(NA, 2)
  } else {
    pred <- factor(pred, levels = levels(obs))
    require("e1071")
    out <- unlist(e1071::classAgreement(table(obs, pred)))[c("diag", "kappa")]

    probs <- data[, cls]
    actual <- model.matrix(~ obs - 1)
    out2 <- LogLoss(actual = actual, pred = probs)
  }
  out <- c(out, out2)
  names(out) <- c("Accuracy", "Kappa", "LogLoss")

  if (any(is.nan(out))) out[is.nan(out)] <- NA

  return(out)
}

#' Threshold calculation
#'
#' This function calculates threshold based on given condition
#'
#' @param predict Vector of predicted p(y=1)
#' @param response Vector of actual class labels
#' @param type (Optional) Threshold type. Currently "ppv" and "sens+spec"
#'  (default) are supported
#'
#' @return Threshold value
#'
#' @examples
#' threshold2(response = c(0,1), predict = c(0.1, 0.9))
#'
#' @export

threshold2 <- function(predict, response, type = "sens+spec") {
  r <- pROC::roc(response, predict)
  threshold <- c()
  if(type == "ppv") {
    for(i in 1:length(r$thresholds)) {
      tryCatch({
        cm <- confusionMatrix(data = ifelse(predict>r$thresholds[i], "Yes", "No"),
                              response)
        ppv <- cm$table[2,2]/sum(cm$table[2,])
        if(ppv>=0.7)
          threshold <- c(threshold, r$thresholds[i])
      }, error = function(e) NULL)
    }
    if(length(threshold)>0) {
      return(min(threshold))
    } else {
      return(1)
    }
  } else {
    return(min(r$thresholds[which.max(r$sensitivities + r$specificities)]))
  }
}

#' Accuracy curve
#'
#' This function plots the accuracy of model at different thresholds
#'
#' @param pred Vector of predicted p(y=1)
#' @param actual Vector of actual class labels
#'
#' @return None
#'
#' @examples
#' accuracy_curve(actual = c(0,1), pred = c(0.1, 0.9))
#'
#' @export

accuracy_curve <- function(pred, actual) {
  pred <- ROCR::prediction(pred, actual)
  perf <- ROCR::performance(pred, "acc")
  plot(perf@x.values[[1]], perf@y.values[[1]], type='l')
}

#' TPR - PPV curve
#'
#' This function plots the true positive rate and positive predictive value
#'  of model at different thresholds in the same line chart
#'
#' @param pred Vector of predicted p(y=1)
#' @param actual Vector of actual class labels
#'
#' @return None
#'
#' @examples
#' tpr_ppv_curve(actual = c(0,1), pred = c(0.1, 0.9))
#'
#' @export

tpr_ppv_curve <- function(pred, actual) {
  pred <- ROCR::prediction(pred, actual)
  perf <- ROCR::performance(pred, "ppv")
  plot(perf@x.values[[1]], perf@y.values[[1]], type='l')
  perf <- ROCR::performance(pred, "tpr")
  lines(perf@x.values[[1]], perf@y.values[[1]], type='l')
}

#' Model Summaries
#'
#' This function builds the list of models mentioned and saves the model summaries
#'  in files for train, test, blind (optional) and validation set.
#'
#' @param form1 Formula for outcome as function of independent variables
#' @param models Vector of model names
#' @param i Model index
#' @param train1 Training set
#' @param test Testing set
#' @param val Number of variables
#' @param cl CaretList containing previously built models
#' @param vars Vector of variables used for modeling
#' @param wd Project working directory
#' @param resamples Bootstrap resamples generated from analytical dataset
#' @param type Type of prediction. Currently supports only "prob"
#' @param blind_data Blind data set
#' @param outcome Outcome variable name
#'
#' @return None
#'
#' @examples
#' get_model_summaries(form1 = Species~., models = c("glm", "rf"), i = 1,
#'  train1 = train_set, test = test_set, val = validation_set, cl = caretList(),
#'  vars = num_vars, wd = "C:\\", resamples = createResample(), type = "prob",
#'  blind_data = c(), outcome = "Species")
#'
#' @export

get_model_summaries <- function(form1, models, i, train1, test, val, cl, vars, wd,
                                resamples, type, blind_data = c(),
                                outcome = "DSAT_Flag") {
  if(models[i]=="glmStepAIC") {
    model <- train(form = form1, data = train1, method = "glmStepAIC",
                   trControl = trainControl(
                     savePredictions = T, classProbs = T,
                     index = resamples), tuneLength = 5, trace = F)
  } else if(models[i]=="glmnet") {
    model <- train(form = form1, data = train1, method = "glmnet",
                   trControl = trainControl(
                     savePredictions = T, classProbs = T,
                     index = resamples, search = "random"), tuneLength = 5,
                   tuneGrid = expand.grid(.alpha = seq(.05, 1, length = 15),
                                          .lambda = c((1:5)/10)))
  } else if(models[i]=="rf") {
    model <- train(form = form1, data = train1, method = "rf",
                   tuneGrid = expand.grid(mtry = seq(2, min(12, vars), 2)),
                   trControl = trainControl(
                     savePredictions = T, classProbs = T,
                     index = resamples), tuneLength = 5,
                   family = "binomial", metric = "ROC")
    saveRDS(model, "rf_model.Rds")
  } else if(models[i]=="bartMachine") {
    model <- train(form = form1, data = train1, method = "bartMachine",
                   trControl = trainControl(
                     savePredictions = T, classProbs = T,
                     index = resamples), tuneLength = 5, verbose = FALSE,
                   seed = 1)
  } else if(models[i]=="C5.0") {
    c50_grid <- expand.grid(.winnow = c(TRUE,FALSE),
                            .trials=c(1,5,10,15,20),
                            .model="tree")
    model <- train(form = form1, data = train1, method = "C5.0",
                   trControl = trainControl(
                     savePredictions = T, classProbs = T, index = resamples),
                   verbose = F, tuneLength = 5)
  } else if(models[i]=="earth") {
    egrid <- data.frame(.degree = 1, .nprune = (2:4)*2)
    model <- train(form = form1, data = train1, method = "earth",
                   trControl = trainControl(
                     savePredictions = T, classProbs = T, index = resamples),
                   tuneGrid = egrid)
  } else if(models[i]=="xgbTree") {
    xgb_grid <- expand.grid(
      eta = c(0, 0.1, 0.3),
      max_depth = c(2, 4, 6, 8, 10),
      nrounds = 400,
      gamma = 0,               #default=0
      colsample_bytree = 1,    #default=1
      min_child_weight = 1,     #default=1
      subsample = 0.5
    )
    xgb_trcontrol <- trainControl(index = resamples, classProbs = TRUE,
                                  allowParallel = TRUE, savePredictions = "all")
    model <- train(form = form1, data = train1, method = "xgbTree",
                   metric = "logloss", tuneGrid = xgb_grid, verbose = F,
                   trControl = xgb_trcontrol, tuneLength = 5)
  } else if(models[i]=="xgbLinear") {
    xgb_grid <- expand.grid(
      eta = c(0, 0.1, 0.3),
      lambda = c(0, 0.1, 0.3),
      nrounds = 400,
      alpha = c(0, 0.1, 0.3)
    )
    xgb_trcontrol <- trainControl(index = resamples, classProbs = TRUE,
                                  allowParallel = TRUE, savePredictions = "all")
    model <- train(form = form1, data = train1, method = "xgbLinear",
                   metric = "logloss", tuneGrid = xgb_grid, verbose = F,
                   trControl = xgb_trcontrol, tuneLength = 5)
  } else {
    model <- train(form = form1, data = train1, method = models[i],
                   trControl = trainControl(
                     savePredictions = T, classProbs = T,
                     index = resamples), tuneLength = 5)
  }
  saveRDS(model, paste0(models[i], ".Rds"))
  cl[[i]] <- model
  train_pred <- predict(model, train1, type = "prob")[,2]
  test_pred <- predict(model, test, type = "prob")[,2]
  val_pred <- predict(model, val, type = "prob")[,2]
  if(!is.null(blind_data))
    blind_pred <- predict(model, blind_data, type = "prob")[,2]
  if(models[i]=="bartMachine") {
    train_pred <- 1-train_pred
    test_pred <- 1-test_pred
    val_pred <- 1-val_pred
    blind_pred <- 1-blind_pred
  }
  dir.create("Train", showWarnings = F)
  setwd("Train")
  png(paste0(type, "_", models[i], "_train_accuracy_curve.png"),
      width = 600, height = 400)
  tryCatch({
    accuracy_curve(pred = train_pred, actual = train1[[outcome]])
  }, error = function(e) NULL)
  dev.off()
  png(paste0(type, "_", models[i], "_train_tpr_fnr.png"),
      width = 600, height = 400)
  tryCatch({
    tpr_ppv_curve(pred = train_pred, actual = train1[[outcome]])
  }, error = function(e) NULL)
  dev.off()
  setwd(wd)
  dir.create("Test", showWarnings = F)
  setwd("Test")
  png(paste0(type, "_", models[i], "_test_accuracy_curve.png"),
      width = 600, height = 400)
  tryCatch({
    accuracy_curve(pred = test_pred, actual = test[[outcome]])
  }, error = function(e) NULL)
  dev.off()
  png(paste0(type, "_", models[i], "_test_tpr_fnr.png"),
      width = 600, height = 400)
  tryCatch({
    tpr_ppv_curve(pred = test_pred, actual = test[[outcome]])
  }, error = function(e) NULL)
  dev.off()
  setwd(wd)
  dir.create("Val", showWarnings = F)
  setwd("Val")
  png(paste0(type, "_", models[i], "_val_accuracy_curve.png"),
      width = 600, height = 400)
  tryCatch({
    accuracy_curve(pred = val_pred, actual = val[[outcome]])
  }, error = function(e) NULL)
  dev.off()
  png(paste0(type, "_", models[i], "_val_tpr_fnr.png"),
      width = 600, height = 400)
  tryCatch({
    tpr_ppv_curve(pred = val_pred, actual = val[[outcome]])
  }, error = function(e) NULL)
  dev.off()

  setwd(wd)
  dir.create("Blind", showWarnings = F)
  setwd("Blind")
  if(!is.null(blind_data)) {
    blind_data[[outcome]] <- "No"
    blind_data[blind_data$SCORE<2, outcome] <- "Yes"
    png(paste0(type, "_", models[i], "_blind_accuracy_curve.png"),
        width = 600, height = 400)
    tryCatch({
      accuracy_curve(pred = blind_pred, actual = blind_data[[outcome]])
    }, error = function(e) NULL)
    dev.off()
    png(paste0(type, "_", models[i], "_test_tpr_fnr.png"),
        width = 600, height = 400)
    tryCatch({
      tpr_ppv_curve(pred = blind_pred, actual = blind_data[[outcome]])
    }, error = function(e) NULL)
    dev.off()
  }
  thr <- threshold2(predict = test_pred, response = test[[outcome]])
  setwd("../Train")
  write.csv(table(train1[[outcome]], train_pred>thr, dnn = c("Actual", "Pred")),
            paste0(type, "_", models[i], "_train_Conf.csv"))
  setwd("../Test")
  write.csv(table(test[[outcome]], test_pred>thr, dnn = c("Actual", "Pred")),
            paste0(type, "_", models[i], "_test_Conf.csv"))
  setwd("../Val")
  write.csv(table(val[[outcome]], val_pred>thr, dnn = c("Actual", "Pred")),
            paste0(type, "_", models[i], "_val_Conf.csv"))
  if(!is.null(blind_data)) {
    setwd("../Blind")
    write.csv(
      table(blind_data[[outcome]], blind_pred>thr, dnn = c("Actual", "Pred")),
              paste0(type, "_", models[i], "_blind_Conf.csv"))
  }
  setwd(wd)
  return(cl)
}

#' Flooring Capping
#'
#' This function performs outlier treatment of a single column by floor/cap for
#'  given threshold.
#'
#' @param col_vec Vector of continuous variables
#' @param floor_thr Threshold for performing floor/cap
#'
#' @return Vector of values after floor/cap
#'
#' @examples
#' floor_cap_column(col_vec = c(1,2,3), floor_thr = 0.05)
#'
#' @export

floor_cap_column <- function(col_vec, floor_thr) {
  lower <- quantile(x = col_vec, floor_thr/2)
  upper <- quantile(x = col_vec, 1-(floor_thr/2))
  col_vec[col_vec<lower] <- lower
  col_vec[col_vec>upper] <- upper
  return(col_vec)
}

do_floor_cap <- function(full_data, cont_vars, floor_thr) {
  for(col in cont_vars) {
    full_data[[col]] <- floor_cap_column(col_vec = full_data[[col]],
                                         floor_thr = floor_thr)
  }
  return(full_data)
}

#' PCA
#'
#' This function performs PCA on continuous independent variables and returns
#'  principal component projections for different thresholds of variance explained
#'
#' @param full_data Data set
#' @param cont_vars Set of continuous variables
#' @param threshold Threshold for converting probability to class
#' @param variables Complete list of all variables
#'
#' @return Data set with principal components as columns
#'
#' @examples
#' do_pca(full_data = data.frame(a=c(1,2), b=c(3,4)), threshold = 0.95,
#'  cont_vars = c("a", "b"), variables = c("a", "b"))
#'
#' @export

do_pca <- function(full_data, cont_vars, threshold, variables) {
  pca_results <- princomp(full_data[, cont_vars], scores = T)
  comps <- predict(pca_results, full_data[, cont_vars])
  # pca_results$scores <- data.frame(sapply(pca_results$scores, as.numeric))
  vars <- apply(pca_results$scores, 2, var)
  props <- vars / sum(vars)
  props <- cumsum(props)
  index <- which(props>threshold)
  index <- index[-c(1)]
  exclude <- names(index)
  for(k in 1:length(cont_vars)) {
    full_data[[cont_vars[k]]] <- NULL
  }
  full_data <- data.frame(cbind(full_data, comps))
  for(k in 1:length(exclude)) {
    full_data[[exclude[k]]] <- NULL
  }
  variables <- setdiff(c(variables, colnames(pca_results$scores)),
                       c(cont_vars, exclude))
  cont_vars <- setdiff(colnames(pca_results$scores), exclude)
  saveRDS(pca_results, "pca_results.Rds")
  return(list(variables = variables, cont_vars = cont_vars,
              full_data = full_data, exclude = exclude))
}

#' Scaling
#'
#' This function performs scaling on continuous independent variables and returns
#'  principal component projections for different thresholds of variance explained
#'
#' @param full_data Data set
#' @param cont_vars Set of continuous variables
#' @param variables Complete list of all variables
#'
#' @return Data set with columns standardized (mean = 0, sd = 1)
#'
#' @examples
#' do_scale(full_data = data.frame(a=c(1,2), b=c(3,4)), variables = c("a", "b"),
#'  cont_vars = c("a", "b"))
#'
#' @export

do_scale <- function(full_data, cont_vars, variables) {
  full_data[, cont_vars] <- scale(full_data[, cont_vars])
  nas <- sapply(full_data, function(col) sum(is.na(col)))
  nas <- names(nas)[nas>0]
  if(length(nas)>0) {
    for(col in nas) {
      full_data[[col]] <- NA
    }
  }
  cont_vars <- setdiff(cont_vars, nas)
  variables <- setdiff(variables, nas)
  lis <- list(full_data = full_data, cont_vars = cont_vars, variables = variables)
  return(lis)
}

#=================================================================================
# What it does: Performs multivariate missing value imputation
#               using chained equations
# Input: Input full data, continous variables and variable list
# Output: Transformed dataset after imputing missing values

impute_nas <- function(full_data, variables, cont_vars) {
  nas <- sapply(full_data[, variables], function(col) sum(is.na(col)))
  nas <- names(nas)[nas>0]
  imputation <- mice(data = full_data[, variables], m = 50, method = "pmm",
                     diagnostics = F, seed = 1, printFlag = F)
  discrete <- setdiff(nas, cont_vars)
  continuous <- intersect(nas, cont_vars)
  removes <- c()
  if(length(continuous)>0) {
    if(length(continuous)>1) {
      for(i in 1:length(continuous)) {
        if(!is.null(imputation$imp[[continuous[i]]])) {
          full_data[is.na(full_data[[continuous[i]]]), continuous[i]] <-
            rowMeans(imputation$imp[[continuous[i]]])
        } else {
          removes <- c(removes, continuous[[i]])
        }
      }
    } else {
      full_data[is.na(full_data[,continuous]), continuous] <-
        mean(imputation$imp[[continuous]])
    }
  }
  if(length(discrete)>0) {

  }
  variables <- setdiff(variables, removes)
  cont_vars <- setdiff(cont_vars, removes)
  lis <- list(full_data = full_data, variables = variables, cont_vars = cont_vars)
  return(lis)
}

#=================================================================================
# What it does: Builds multiple models and gives extensive summary of models
# Input: Full data, blind data, variable list, continuous variable list,
#        scaling flag, PCA flag, PCA variance explained threshold,
#        model list, flooring capping flag, flooring capping threshold,
#        stratified sampling flag, number of clusters you want over the optimum
#        cluster. Whether you want the cluster number in the predictor variable
#        or not, impute missing value flag
# Output: If single model is given with PCA = True, then model object is saved.
#         PCA model object based on variance threshold, ROC curve, accuracy curve
#         and confusion matrix are saved on train, test and validation dataset
#         If single model is given witout  PCA, then model object is saved.
#         ROC curve, accuracy curve and confusion matrix are saved on train, test
#         and validation dataset

model_and_summarize <- function(full_data, variables, cont_vars, scale = F, type,
                                pca = F, threshold = 0.95, models = "rf", wd,
                                floor_cap = F, floor_thr = 0.05, cluster = F,
                                over_k = 10, include_cluster = F, impute_na = F,
                                blind_data = NA, outcome = "DSAT_Flag") {
  setwd(wd)
  for(i in 1:length(cont_vars)) {
    full_data[[cont_vars[i]]] <- as.numeric(full_data[[cont_vars[i]]])
    if(!is.null(blind_data))
      blind_data[[cont_vars[i]]] <- as.numeric(blind_data[[cont_vars[i]]])
  }
  if(impute_na) {
    tryCatch({
      full_data <- impute_nas(full_data, variables, cont_vars)
      variables <- full_data[[2]]
      cont_vars <- full_data[[3]]
      full_data <- full_data[[1]]
    }, error = function(e) NULL)
    if(!is.null(blind_data))
      tryCatch({
        blind_data <- impute_nas(blind_data, variables, cont_vars)
        variables <- blind_data[[2]]
        cont_vars <- blind_data[[3]]
        blind_data <- blind_data[[1]]
      }, error = function(e) NULL)
  } else {
    na_rows <- apply(full_data, 1, function(row) sum(is.na(row)))
    full_data <- full_data[na_rows==0, ]
    na_rows <- apply(blind_data, 1, function(row) sum(is.na(row)))
    if(!is.null(blind_data))
      blind_data <- blind_data[na_rows==0, ]
  }
  if(scale) {
    full_data <- do_scale(full_data = full_data, cont_vars = cont_vars,
                          variables = variables)
    if(!is.null(blind_data))
      blind_data <- do_scale(full_data = blind_data, cont_vars = cont_vars,
                             variables = variables)
    cont_vars <- full_data[[2]]
    variables <- full_data[[3]]
    full_data <- full_data[[1]]
    if(!is.null(blind_data)) {
      cont_vars <- intersect(cont_vars, blind_data[[2]])
      variables <- intersect(variables, blind_data[[3]])
      blind_data <- blind_data[[1]]
    }
  }
  if(floor_cap) {
    full_data <- do_floor_cap(full_data = full_data, cont_vars = cont_vars,
                              floor_thr = floor_thr)
    if(!is.null(blind_data))
      blind_data <- do_floor_cap(full_data = blind_data, cont_vars = cont_vars,
                                 floor_thr = floor_thr)
  }
  if(pca) {
    full_data <- do_pca(full_data = full_data, cont_vars = cont_vars,
                        threshold = threshold, variables = variables)
    pca_model <- readRDS("pca_results.Rds")
    variables <- full_data[[1]]
    cont_vars <- full_data[[2]]
    exclude <- full_data[[4]]
    full_data <- full_data[[3]]
    for(i in 1:length(cont_vars))
      full_data[, cont_vars[i]] <- as.numeric(full_data[, cont_vars[i]])
    print(cont_vars)
    print(sapply(full_data[,cont_vars], class))
    if(!is.null(blind_data)) {
      blind_data_pc <- predict(pca_model, blind_data[, cont_vars])
      blind_data[, cont_vars] <- blind_data_pc[, cont_vars]
    }
  }
  set.seed(1)
  if(cluster) {
    full_data <- do_scale(full_data = full_data, cont_vars = cont_vars,
                          variables = variables)
    cont_vars <- full_data[[2]]
    variables <- full_data[[3]]
    full_data <- full_data[[1]]
    pamk_best <- pamk(full_data[, cont_vars])
    k_best <- pamk_best$nc
    clustering <- kcca(full_data[, cont_vars], k = k_best + over_k,
                       kccaFamily("kmeans"))
    cluster <- predict(clustering)
    if(include_cluster) {
      full_data$cluster <- as.factor(cluster)
      variables <- c(variables, "cluster")
      if(!is.null(blind_data))
        blind_data$cluster <- as.factor(
          predict(clustering, blind_data[, cont_vars]))
    }
    full_data <- split(full_data, cluster)
    val_test_train <- lapply(full_data, function(df) {
      in_train <- createFolds(y = df[[outcome]], k = 3, list = F, returnTrain = F)
      train1 <- df[in_train==1, ]
      test <- df[in_train==2, ]
      val <- df[in_train==3, ]
      lis <- list(train1 = train1, test = test, val = val)
      return(lis)
    })
    train1 <- lapply(val_test_train, function(df_list) {
      return(df_list$train1)
    })
    test <- lapply(val_test_train, function(df_list) {
      return(df_list$test)
    })
    val <- lapply(val_test_train, function(df_list) {
      return(df_list$val)
    })
    train1 <- do.call("rbind.fill", train1)
    test <- do.call("rbind.fill", test)
    val <- do.call("rbind.fill", val)
    rownames(train1) <- rownames(test) <- rownames(val) <- NULL
  } else {
    in_train <- createFolds(y = full_data[[outcome]], k = 3, list = F,
                            returnTrain = F)
    train1 <- df[in_train==1, ]
    test <- df[in_train==2, ]
    val <- df[in_train==3, ]
    rownames(train1) <- rownames(test) <- rownames(val) <- NULL
  }
  if("CATEGORY_LEVEL_1" %in% variables) {
    common_cats <- intersect(train1$CATEGORY_LEVEL_1, test$CATEGORY_LEVEL_1)
    common_cats <- intersect(common_cats, val$CATEGORY_LEVEL_1)
    test <- test[test$CATEGORY_LEVEL_1 %in% common_cats, ]
    val <- val[val$CATEGORY_LEVEL_1 %in% common_cats, ]
    train1 <- train1[train1$CATEGORY_LEVEL_1 %in% common_cats, ]
    if(!is.null(blind_data))
      blind_data <- blind_data[blind_data$CATEGORY_LEVEL_1 %in% common_cats, ]
  }
  if("DISPOSITION_LEVEL_2" %in% variables) {
    common_disps <- intersect(train1$DISPOSITION_LEVEL_2, test$DISPOSITION_LEVEL_2)
    common_disps <- intersect(common_disps, val$DISPOSITION_LEVEL_2)
    test <- test[test$DISPOSITION_LEVEL_2 %in% common_disps, ]
    val <- val[val$DISPOSITION_LEVEL_2 %in% common_disps, ]
    train1 <- train1[train1$DISPOSITION_LEVEL_2 %in% common_disps, ]
    if(!is.null(blind_data))
      blind_data <- blind_data[blind_data$DISPOSITION_LEVEL_2 %in% common_disps, ]
  }
  if("DISPOSITION_LEVEL_1" %in% variables) {
    common_disps <- intersect(train1$DISPOSITION_LEVEL_1, test$DISPOSITION_LEVEL_1)
    common_disps <- intersect(common_disps, val$DISPOSITION_LEVEL_1)
    test <- test[test$DISPOSITION_LEVEL_1 %in% common_disps, ]
    val <- val[val$DISPOSITION_LEVEL_1 %in% common_disps, ]
    train1 <- train1[train1$DISPOSITION_LEVEL_1 %in% common_disps, ]
    if(!is.null(blind_data))
      blind_data <- blind_data[blind_data$DISPOSITION_LEVEL_1 %in% common_disps, ]
  }
  if("Problem" %in% variables) {
    common_probs <- intersect(train1$Problem, test$Problem)
    common_probs <- intersect(common_probs, val$Problem)
    test <- test[test$Problem %in% common_probs, ]
    val <- val[val$Problem %in% common_probs, ]
    train1 <- train1[train1$Problem %in% common_probs, ]
    if(!is.null(blind_data))
      blind_data <- blind_data[blind_data$Problem %in% common_probs, ]
  }
  form1 <- as.formula(paste0(outcome, " ~ ", paste0(variables, collapse = " + ")))
  resamples <- createResample(train1[[outcome]], 5)
  print(form1)
  if(length(models)==1) {
    cl <- list()
    cl <- get_model_summaries(form1 = form1, models = models, i = 1,
                              train1 = train1, test = test, val = val, cl = cl,
                              vars = length(variables), wd = wd,
                              resamples = resamples, type = type,
                              blind_data = blind_data, outcome = outcome)
    train_pred <- predict(cl[[1]], train1, type = "prob")
    test_pred <- predict(cl[[1]], test, type = "prob")
    val_pred <- predict(cl[[1]], val, type = "prob")
    if(!is.null(train_pred)) {
      train_pred <- train_pred[,2]
      test_pred <- test_pred[,2]
      val_pred <- val_pred[,2]
    }
    if(models[i]=="bartMachine") {
      train_pred <- 1-train_pred
      test_pred <- 1-test_pred
      val_pred <- 1-val_pred
    }
    setwd(wd)
  } else {
    cl <- list()
    for(k in 1:length(models)) {
      cl <- get_model_summaries(form1 = form1, models = models, train1 = train1,
                                test = test, val = val, cl = cl, i = k, wd = wd,
                                vars = length(variables), resamples = resamples,
                                type = type, blind_data = blind_data,
                                outcome = outcome)
      setwd(wd)
    }
    names(cl) <- models
    class(cl) <- "caretList"
    # stack <- caretStack(cl, method = "glm", metric = "LogLoss",
    #                     trControl = trainControl(
    #                       method = "repeatedcv", number = 20,
    #                       savePredictions = "final", classProbs = TRUE,
    #                       summaryFunction = LogLosSummary))
    # print(stack$error)
    train_pred <- data.frame(sapply(1:length(cl), function(i)
      pred_new(i, cl, train1)))
    test_pred <- data.frame(sapply(1:length(cl), function(i)
      pred_new(i, cl, test)))
    val_pred <- data.frame(sapply(1:length(cl), function(i)
      pred_new(i, cl, val)))
    colnames(test_pred) <- colnames(train_pred) <-
      colnames(val_pred) <- names(cl)
    # train_pred$ensemble <- predict(stack, train1, type = "prob")
    # test_pred$ensemble <- predict(stack, test, type = "prob")
    # val_pred$ensemble <- predict(stack, val, type = "prob")
    val_pred[[outcome]] <- val[[outcome]]
    test_pred[[outcome]] <- test[[outcome]]
    train_pred[[outcome]] <- train1[[outcome]]
    form <- as.formula(paste0(outcome, "~."))
    stack <- glm(form, data = test_pred, family = binomial)
    val_pred$ensemble <- predict(stack, val_pred, type = "response")
    train_pred$ensemble <- predict(stack, train_pred, type = "response")
    test_pred$ensemble <- predict(stack, test_pred, type = "response")
    setwd(wd)
    setwd("Train")
    thr <- threshold2(predict = test_pred$ensemble, response = test[[outcome]])
    png(paste0(type, "_ensemble_train_accuracy_curve.png"),
        width = 600, height = 400)
    tryCatch({
      accuracy_curve(pred = train_pred$ensemble, actual = train1[[outcome]])
    }, error = function(e) NULL)
    dev.off()
    png(paste0(type, "_ensemble_train_tpr_fnr.png"),
        width = 600, height = 400)
    tryCatch({
      tpr_ppv_curve(pred = train_pred$ensemble, actual = train1[[outcome]])
    }, error = function(e) NULL)
    dev.off()
    write.csv(table(train1[[outcome]], train_pred$ensemble>thr,
                    dnn = c("Actual", "Pred")),
              paste0(type, "_ensemble_train_Conf.csv"))
    setwd("../Test")
    png(paste0(type, "_ensemble_test_accuracy_curve.png"),
        width = 600, height = 400)
    tryCatch({
      accuracy_curve(pred = test_pred$ensemble, actual = test[[outcome]])
    }, error = function(e) NULL)
    dev.off()
    png(paste0(type, "_ensemble_test_tpr_fnr.png"),
        width = 600, height = 400)
    tryCatch({
      tpr_ppv_curve(pred = test_pred$ensemble, actual = test[[outcome]])
    }, error = function(e) NULL)
    dev.off()
    write.csv(table(test[[outcome]], test_pred$ensemble>thr,
                    dnn = c("Actual", "Pred")),
              paste0(type, "_ensemble_test_Conf.csv"))
    setwd("../Val")
    png(paste0(type, "_ensemble_val_accuracy_curve.png"),
        width = 600, height = 400)
    tryCatch({
      accuracy_curve(pred = val_pred, actual = val[[outcome]])
    }, error = function(e) NULL)
    dev.off()
    png(paste0(type, "_ensemble_val_tpr_fnr.png"),
        width = 600, height = 400)
    tryCatch({
      tpr_ppv_curve(pred = val_pred, actual = val[[outcome]])
    }, error = function(e) NULL)
    dev.off()
    write.csv(table(val[[outcome]], val_pred$ensemble>thr,
                    dnn = c("Actual", "Pred")),
              paste0(type, "_ensemble_val_Conf.csv"))
    setwd(wd)
  }
  if(!is.null(blind_data)) {
    blind_pred <- data.frame(sapply(1:length(cl), function(i)
      pred_new(i, cl, blind_data)))
    colnames(blind_pred) <- names(cl)
    blind_pred$ensemble <- predict(stack, blind_pred, type = "response")
    blind_data[[outcome]] <- "No"
    blind_data[blind_data$SCORE<=2, outcome] <- "Yes"
    write.csv(table(blind_data[[outcome]], blind_pred$ensemble>thr,
                    dnn = c("Actual", "Pred")),
              paste0(type, "_ensemble_blind_Conf.csv"))
  }
  # write.csv(x = caTools::colAUC(X = val_pred,
  #                               y = val[[outcome]]),
  #           file = paste0(type, "_val_AUC.csv"))
  # write.csv(x = caTools::colAUC(X = test_pred,
  #                               y = test[[outcome]]),
  #           file = paste0(type, "_test_AUC.csv"))
  # write.csv(x = caTools::colAUC(X = train_pred,
  #                               y = train1[[outcome]]),
  #           file = paste0(type, "_train_AUC.csv"))
  setwd(wd)
}

