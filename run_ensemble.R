library(caret)
library(caretEnsemble)
library(caTools)
library(fpc)
library(plyr)
library(mice)
library(flexclust)
library(readr)
library(EnsembleModel)
options(stringsAsFactors = F)

#=================================================================================
# Reading arguments from command line

args <- commandArgs(trailingOnly = TRUE)
csat <- 0
compare_full_data_with_blind_data <- F
if(length(args)==2) {
  if(args[1]=="csat_cleaning_modeling.R")
    csat <- 1
  source(args[1])
} else if(length(args)==0) {
  source("csat_cleaning_modeling.R")
  csat <- 1
  compare_full_data_with_blind_data <- T
}

#=================================================================================
# Build models and summarize for different thresholds of variance explained in PCA

for(threshold in seq(from = 0.95, to = 0.95, by = 0.05)) {
  if(csat==0) {
    source(args[2])
  } else {
    build_csat_model(main_wd, threshold, voice_only_data, voice_variables,
                     voice_cont_vars, voice_only_blind_data, outcome,
                     email_only_data, chat_only_data, email_chat_variables,
                     email_chat_cont_vars, email_only_blind_data,
                     chat_only_blind_data)
  }
}

#=================================================================================
# Comparison full data joint distribution with blind data joint distribution

# if(!is.null(blind_data) & compare_full_data_with_blind_data) {
#   vars <- union(email_chat_variables, voice_variables)
#   cont_var <- union(email_chat_cont_vars, voice_cont_vars)
#   for(var in cont_var) {
#     blind_data[[var]] <- as.integer(blind_data[[var]])
#     blind_data[is.na(blind_data[[var]]), var] <- 0
#   }
#
#   #===============================================================================
#   # Check for significant difference in univariate distributions of predictors
#   # between training data set and blind data set
#
#   diff_vars <- c()
#
#   for(var in intersect(email_chat_variables, voice_variables)) {
#     if(var %in% email_chat_cont_vars) {
#       diff <- ks.test(x = full_data[[var]], y = blind_data[[var]])
#       if(diff$p.value < 0.05)
#         diff_vars <- c(diff_vars, var)
#     } else {
#       t1 <- table(full_data[[var]])
#       t2 <- table(blind_data[[var]])
#       t1 <- t1[intersect(names(t2), names(t1))]
#       t2 <- t2[intersect(names(t2), names(t1))]
#       ex <- t1/sum(t1)
#       ex <- sum(t2)*ex
#       X2 <- sum((t2-ex)^2/ex)
#       p_val <- 1-pchisq(q = X2, df = length(t1)-1)
#       if(p_val < 0.05)
#         diff_vars <- c(diff_vars, var)
#     }
#   }
#
#   #===============================================================================
#   # Test for joint distribution using generative model (naive Bayes)
#   form1 <- as.formula(paste0(outcome, " ~ ", paste0(vars, collapse = " + ")))
#   model <- naiveBayes(formula = form1, data = full_data)
#   full_pred <- predict(model, full_data[, vars], type = "raw")[,2]
#   blind_pred <- predict(model, blind_data, type = "raw")[,2]
#   df <- data.frame(blind_pred = blind_pred, full_pred = full_pred)
#   df1 <- data.frame(IRT = full_data$issue_resolve_time, type = "full")
#   df2 <- data.frame(IRT = blind_data$issue_resolve_time, type = "blind")
#   df <- rbind(df1, df2)
#   ggplot(df, aes(IRT))+
#     geom_histogram(data = subset(df, type == 'blind'),
#                    fill = "red", alpha = 0.4, bins = 60,
#                    aes(y = ..count../sum(..count..))) +
#     geom_histogram(data = subset(df, type == 'full'),
#                    fill = "blue", alpha = 0.4, bins = 60,
#                    aes(y = ..count../sum(..count..))) +
#     xlab("Issue Resolve Time") + ylab("Relative Frequency")
# }
