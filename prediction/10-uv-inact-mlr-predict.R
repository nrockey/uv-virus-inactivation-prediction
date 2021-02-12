# 10 - Multiple linear regression model
# Predicts UV virus inactivation rates for prediction set using models developed 
# from training set. 

# 2020.10.17
# Nicole Rockey

# LIBRARIES---------------------------------------------------------------------
library(tidyverse); library(glmnet)

# SET WORKING DIRECTORY---------------------------------------------------------
setwd("~/Documents/GitHub/uv-inactivation-model")

# DATA INPUT--------------------------------------------------------------------

# Loads in 'id_vars, ind_vars, dep_vars, class_vars, pred_data.'
# load("data/uv-inact-model-predict-inputs.RData")
bar = load("data/uv-inact-model-predict-inputs-test.RData")

# ACTION ITEM: Update which of these lines are and are not commented to predict
# inactivation using different developed models.
# Loads in mlr_opt models depending on the virus subset used for training.
# load("development/development-results/mlr/mlr-opt-rep-TRUE-host-TRUE-all.RData")
# load("development/development-results/mlr/mlr-opt-rep-FALSE-host-FALSE-all.RData")
foo = load("development/development-results/mlr/mlr-opt-rep-TRUE-host-TRUE-dsdna.RData")
# load("development/development-results/mlr/mlr-opt-rep-FALSE-host-FALSE-dsdna.RData")
# load("development/development-results/mlr/mlr-opt-rep-FALSE-host-FALSE-plusssrna.RData")

# VARIABLE SET-UP---------------------------------------------------------------

# Creates data set for virus id.
pred_id = pred_data[, id_vars_pred]

# Removes all class_vars from prediction set in prediction set without class_vars_pred.
pred_X_no_cat = model.matrix( ~ ., pred_data[, setdiff(ind_vars_pred, cat_vars_pred)] 
)[, -1]

# N: Work here!!
# Creates the prediction matrix including only class_vars_pred.
#! Set the levels for the prediction variables:
pred_data$repair = factor(pred_data$repair, levels = c('0', '1'))
pred_data$host = factor(pred_data$host, levels = c('0', '2'))

# N: Work here!!
#! Main error is that "type" is included in cat_vars_pred but shouldn't
#  be b/c it is not part of the mlr_opt model. 
cat_vars_pred = setdiff(cat_vars_pred, 'type')

if ( !is.null(cat_vars_pred) ) {
  pred_X_cat = 
    model.matrix( ~ ., pred_data[, cat_vars_pred, drop = FALSE] )[, -1] 
} else {
  pred_X_cat = NULL
}

# PREDICTION--------------------------------------------------------------------

# Generates prediction matrix for each virus.
pred_mat = vector(length = length(uniq_class), mode = 'list')
names(pred_mat) = names(mlr_opt)

# Fills out prediction matrix with predicted, experimental, and error for each
# virus in subset.
for (class in names(pred_mat)) {
  # Subset prediction data to only have the same categories as what is in the
  # optimal model.
  pred_X_subset = pred_X_no_cat[ , mlr_opt[[class]]$vars]
  
  # Standardize data here (or can this be done in the function?)
  pred_X_std = t( {t(pred_X_subset) - mlr_opt[[class]]$mean} / mlr_opt[[class]]$sd)
  
  # Calculate pca values for independent vars.
  pred_X_pca = pred_X_std %*% mlr_opt[[class]]$loadX
  
  # N: 10/13/20 - How to make this so that if there are no variables in here that are different,
  # it won't set off an error...
  # Make the class prediction set the same as what is in the model development set.
  pred_X_cat_subset = pred_X_cat[ , mlr_opt[[class]]$cat_vars]
  
  # Calculates predicted inactivation.
  kpred = cbind(1, pred_X_pca, pred_X_cat_subset) %*% mlr_opt[[class]]$beta

  # Prediction matrix to summarize info about mlr predictions.
  pred_mat[[class]] = kpred
  row.names(pred_mat[[class]]) = pred_id
  colnames(pred_mat[[class]]) = 'kpred'
}

# OUTPUT------------------------------------------------------------------------

# Creates an R Data file name for prediction data.
file_pred = sprintf('prediction/prediction-results/mlr/mlr-predict-rep-%s-host-%s-%s.RData',
                    include_repair, include_host, data_subset)

# Saves an R Data file with the predicted data from the mlr model.
save(pred_mat, data_subset, file = file_pred)

for (class in names(pred_mat)) {
  file = sprintf('prediction/csv-outputs/mlr-predict-summ-rep-%s-host-%s-%s-%s.csv',
                 include_repair, include_host, data_subset, class)
  write.csv(pred_mat[[class]], file = file)
}

