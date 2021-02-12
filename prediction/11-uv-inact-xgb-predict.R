# 11 - Boosted trees model prediction
# Uses xgboost to make a model for predicting UV virus inactivation rates.
# Prediction models were assessed using all data in the training set with the
# optimized hyperparameters.

# 2020.10.17
# Nicole Rockey

# LIBRARIES---------------------------------------------------------------------
library(tidyverse); library(xgboost)

# SET WORKING DIRECTORY---------------------------------------------------------
setwd('~/uv-virus-inactivation-prediction/')

#DATA INPUT---------------------------------------------------------------------

# Loads in 'id_vars, ind_vars, dep_vars, class_vars, pred_data.'
load("data/uv-inact-model-predict-inputs.RData")

# ACTION ITEM: Update which of these lines are and are not commented to predict
# inactivation using different developed models.
# Loads in xgb_opt models depending on the virus subset used for training.
load("development/development-results/xgb/xgb-opt-dsdna.RData")

# VARIABLE SET-UP---------------------------------------------------------------

# Creates data set for virus id.
pred_id = pred_data[ , id_vars_pred]

# Creates independent data set for prediction data.
pred_X = model.matrix( ~ ., pred_data[, ind_vars_pred])[ ,-1]

# PREDICTION--------------------------------------------------------------------

# Generates prediction matrix for each virus.
pred_mat = vector(length = length(uniq_class), mode = 'list')
names(pred_mat) = names(xgb_opt)

# Fills out prediction matrix with predicted, experimental, and error for each
# virus in subset.
for (class in names(pred_mat)) {
  pred_X_subset = pred_X[ , xgb_opt[[class]][['feature_names']]]
  # Calculates predicted inactivation.
  kpred = predict(xgb_opt[[class]], pred_X_subset)
  # Prediction matrix to summarize info about xgb predictions.
  pred_mat[[class]] = as.matrix(kpred)
  row.names(pred_mat[[class]]) = pred_id
  colnames(pred_mat[[class]]) = 'kpred'
}

# OUTPUT------------------------------------------------------------------------

# Creates an R Data file name for prediction data.
file_pred = sprintf('prediction/prediction-results/xgb-predict-%s.RData', data_subset)

# Saves an R Data file with the predicted data from the elnt model.
save(pred_mat, data_subset, file = file_pred)

for (class in names(pred_mat)) {
  file = sprintf('prediction/csv-outputs/xgb-predict-summ-%s%s.csv',
                 data_subset,
                 ifelse(data_subset == '-all', class, '')
                 )
  write.csv(pred_mat[[class]], file = file)
}

