# 07 - Random forest model
# Objectives of this rf file:
#      1) estimate loocv mse's for all rf models.
#      2) summarize mse by virus type.
#      3) generate optimal models for each virus type.
#      4) save all mse summaries and optimal models in R Data files for
#         downstream use.

# 2020.10.19
# Nicole Rockey

# LIBRARIES---------------------------------------------------------------------

library(tidyverse); library(randomForest); library(xgboost)

# SET WORKING DIRECTORY---------------------------------------------------------

setwd("~/uv-virus-inactivation-prediction/")

#DATA INPUT---------------------------------------------------------------------

# ACTION ITEM: Update which of these lines are and are not commented to develop
# model with different data.
# Loads in 'id_vars, dep_vars, ind_vars, train_data, n_virus, n_exp.'
# load("data/uv-inact-model-data-inputs-all.RData")
# load("data/uv-inact-model-data-inputs-dsdna.RData")
load("data/uv-inact-model-data-inputs-plusssrna.RData")

# FUNCTIONS---------------------------------------------------------------------

# Function to define rf.
rf = function(virus_id, X, y, wt) {
  
  # Defines # of viruses in data set.
  n_virus = length(virus_id)
    
  # Generates k matrix for each virus.
  k_mat = matrix(NA, nrow = nrow(X), ncol = 1)
    
  # Loop through all virus folds for LOOCV.
  for (fold in 1:n_virus) {
    # Fold virus id.
    id_fold = virus_id[fold]
      
    # Virus id for all experiments except the fold virus.
    id_set = virus_id[-which(virus_id == id_fold)]
      
    # Independent variables for the fold virus.
    xfold = X[which(virus_id == id_fold), , drop = FALSE]
      
    # Inactivation rate for the fold virus.
    yfold = y[which(virus_id == id_fold)]
      
    # Independent variables for all viruses except the fold virus.
    xset = X[-which(virus_id == id_fold), ]
      
    # Inactivation rate constant for all viruses except the fold virus.
    yset = y[-which(virus_id == id_fold)]
    
    # Weights for all viruses except the fold virus.
    wset = wt[-which(virus_id == id_fold)]
    # Re-normalize weights to sum to n_virus in training set for each fold.
    wset = length(wset) * wset / sum(wset)
    
    # Hyperparameters to include in random forests model. Attempting to match
    # as closely to randomForest defaults as possible, with exception that
    # min_child_weight is set to 1.
    params = list(colsample_bytree = 1 / 3, # will randomly select 1/3 of variables for each tree.
                  learning_rate = 1, # in random forests models, this is only 1.
                  min_child_weight = 1,
                  num_parallel_tree = 500, # sets number of trees
                  subsample = 0.632 # fraction of data set to sample for each tree.
                  )
    
    dtrain = xgb.DMatrix(data = xset,
                         label = yset,
                         weight = wset)
    
    fit = xgb.train(params = params,
                    data = dtrain,
                    nrounds = 1
                    )

    # Predicts inactivation rate.
    yval = predict(fit, xfold)
      
    # Calculate MSE for this fold.
    k_mat[fold, ] = yval
  }
  
  # Create list to output round and k data.
  rf_info = vector(mode = 'list', length = 1)
  
  # Define information from rf.
  rf_info[[1]] = k_mat

  rf_info
}

# MODEL OPTIMIZATION------------------------------------------------------------

# Defines vector with virus names for all experiments.
train_id_data = train_data[, id_vars]

# Creates independent data set for all training data.
train_ind_data = model.matrix( ~ ., train_data[ , ind_vars])[ ,-1]

# Creates dependent data set for all training data.
train_dep_data = train_data[ , dep_vars]

# Creates weight vector for all training data.
w = train_data[ , weight]
names(w) = rownames(train_data)

# Normalize weights so w sums to number of viruses
w_norm = n_virus * w / sum(w)

# Run random forest function.
rf_info = rf(X = train_ind_data,
             y = train_dep_data,
             virus_id = train_id_data,
             wt = w_norm
            )

k_rf_all = rf_info[[1]]

# Defines the matrix rows with virus names.
row.names(k_rf_all) = train_id_data

class_dict = c('0' = 'dsdna',
              '1' = 'plusssrna',
              '2' = 'ssdna',
              '3' = 'negssrna',
              '4' = 'dsrna' 
              )

class_dict_rev = names(class_dict)
names(class_dict_rev) = class_dict

if ( 'class' %in% names(train_data) ) {
  # Generates list of unique virus types.
  uniq_class = class_dict[unique(train_data[, 'class'])]
  
  # Creates a list to contain kpred, kest, and wt data for each virus.
  k_rf_list = w_list = k_est_list = vector( length = length(uniq_class), mode = 'list')
  names(k_rf_list) = names(w_list) = names(k_est_list) = uniq_class
  
  for ( class in names(k_rf_list) ) {
# N: Run this code to make sure it works!!!
    idx = which(train_data[ , 'class'] == class_dict_rev[class])
    
    # Hold mse results in a single list
    k_rf_list[[class]] = k_rf_all[idx, , drop = FALSE] 
    
    # Hold weight info in a single list.
    w_list[[class]] = w[idx] / sum(w[idx])
    
    # Hold estimated k info in a single list.
    k_est_list[[class]] = train_dep_data[idx]
  }
} else {
  uniq_class = data_subset 
  k_rf_list = list(k_rf_all) # length one
  names(k_rf_list) = uniq_class # a single class
  
  w_list = list(w / sum(w)) # length one
  names(w_list) = uniq_class # a single class
  
  k_est_list = list(train_dep_data) # length one
  names(k_est_list) = uniq_class # a single class
}

# Creates lists to place mse values and summary stats, and names those
# lists.
mse_rf_list = mse_wt_rf_list = mse_sum_list = mse_var_list = mse_se_list = 
  rmse_sum_list = rmse_se_list = mse_rf_summ_list = 
  vector(length = length(uniq_class), mode = 'list')
names(mse_rf_list) = names(mse_wt_rf_list) = names(mse_sum_list) =
  names(mse_var_list) = names(mse_se_list) = names(rmse_sum_list) =
  names(rmse_se_list) = names(mse_rf_summ_list) = uniq_class

for ( class in names(mse_rf_list) ) {
  # Calculates mse for each virus.
  mse_rf_list[[class]] = (k_rf_list[[class]] - k_est_list[[class]])^2
  # Calculate weighted mse values for each virus.
  mse_wt_rf_list[[class]] = mse_rf_list[[class]] * w_list[[class]]
  # Hold averaged mse results in a single list.
  mse_sum_list[[class]] = sum(mse_wt_rf_list[[class]] / k_est_list[[class]]^2)
  # Hold var of mse results in a single list. This is the bias corrected
  # weighted sample variance.
  mse_var_list[[class]] = 
    sum(w_list[[class]] * (mse_rf_list[[class]]) / k_est_list[[class]]^2) / (1 - sum(w_list[[class]]^2))
  # Hold se of mse results in a single list.
  mse_se_list[[class]] = sqrt(mse_var_list[[class]] / nrow(mse_rf_list[[class]]))
  # Hold averaged rmse results in a single list.
  rmse_sum_list[[class]] = sqrt(mse_sum_list[[class]])
  # Hold se of rmse results in a single list.
  rmse_se_list[[class]] = sqrt(mse_se_list[[class]])
  # Summarize the mean and sd in the same list.
  mse_rf_summ_list[[class]] = cbind(mse_sum_list[[class]],
                                    mse_se_list[[class]],
                                    rmse_sum_list[[class]],
                                    rmse_se_list[[class]])
  # Name columns of each list with mean and sd.
  colnames(mse_rf_summ_list[[class]]) = c('mse_sum',
                                          'mse_se',
                                          'rmse_sum',
                                          'rmse_se'
                                          )
}

# PREDICTION MODELS-------------------------------------------------------------

# Creates the model to be used in prediction for all viruses. Use entire
# training set as data for model.
# Makes a list to contain the predictive model for each virus subset.
rf_opt = vector( length = length(uniq_class), mode = 'list')
names(rf_opt) = uniq_class

# parameters for random forests model.
params_opt = list(colsample_bytree = 1 / 3, # will randomly select 1/3 of variables for each tree.
                  learning_rate = 1, # in random forests models, this is only 1.
                  min_child_weight = 1,
                  num_parallel_tree = 500, # sets number of trees
                  subsample = 0.632 # fraction of data set to sample for each tree.
                  )
  
d_opt = xgb.DMatrix(data = train_ind_data,
                    label = train_dep_data,
                    weight = w_norm
                    )

# Creates the predictive model for each virus subset.
for (class in names(rf_opt)) {
  rf_opt[[class]] = xgb.train(params = params_opt,
                              data = d_opt,
                              nrounds = 1,
                              weight = w_norm
                              )
} 

# OUTPUT------------------------------------------------------------------------

# Generates an R Data file name to include the predictive model for the data
# subset.
file_model = sprintf('development/development-results/rf/rf-opt-%s.RData', data_subset)

# Saves the predictive model list to an RData file.
save(rf_opt, data_subset, uniq_class, cat_vars, file = file_model)

# Generate an R Data file name for mse training info.
file_mse = sprintf('development/development-results/rf/rf-summ-loocv-%s.RData', data_subset)

# Saves the loocv mse information for each virus subset and the
# summary information for each virus subset to an RData file.
save(k_rf_all, k_rf_list, mse_rf_list, mse_rf_summ_list, uniq_class,
     file = file_mse)

# VIRUS-SPECIFIC LOOCV SUMMARY--------------------------------------------------

# Creates a list to place the summary info for each virus class for
# the top model - but for all virus experiments.
summ_rf_all = vector(length = length(uniq_class), mode = 'list')
names(summ_rf_all) = uniq_class

# Provides summary information for each class of viruses from
# the top model.
for ( class in names(summ_rf_all) ) {
  # Hold all the predicted and mean virus rate constants from the top model in
  # one list.
  summ_rf_all[[class]] =
    tibble(
      virus = rownames(k_rf_list[[class]]), 
      kpred = k_rf_list[[class]], 
      kest = k_est_list[[class]],
      mspe = mse_rf_list[[class]]
    ) %>%
    mutate( rmspe = sqrt(mspe) )
  
  # Generate an R Data file name for mse training info.
  file_summ = sprintf('development/csv-outputs/rf-summ-kval-%s%s.csv',
                      class,
                      ifelse(data_subset == 'all', '-all', '')
                      )
  
  # Outputs rmspe summary info from top models to a csv - for si table s4,
  # fig 1, si fig s2.
  write.csv(summ_rf_all[[class]], file_summ)
}


