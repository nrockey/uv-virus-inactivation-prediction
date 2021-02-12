# 06 - Boosted trees model
# Objectives of this boosted trees file:
#      1) estimate loocv mse's for all boosted trees models.
#      2) summarize mse by virus type.
#      3) generate optimal models for each virus type.
#      4) save all mse summaries and optimal models in R Data files for
#         downstream use.

# 2020.10.19
# Nicole Rockey

# LIBRARIES---------------------------------------------------------------------
library(tidyverse); library(xgboost); library(tibble)

# SET WORKING DIRECTORY---------------------------------------------------------

setwd("~/uv-virus-inactivation-prediction/")

#DATA INPUT---------------------------------------------------------------------
# ACTION ITEM: Update which of these lines are and are not commented to develop
# model with different data.
# Loads in 'id_vars, dep_vars, ind_vars, train_data, n_virus, n_exp.'
# load("data/uv-inact-model-data-inputs-all.RData")
# foo = load("data/uv-inact-model-data-inputs-dsdna.RData")
load("data/uv-inact-model-data-inputs-plusssrna.RData")

# FUNCTIONS---------------------------------------------------------------------

# Function to define xgb.
xgb = function(virus_id, X, y, wt) {
  
  # Defines # of viruses in data set.
  n_virus = length(virus_id)

  # Rounds of boosting iterations to undertake
  n_rounds = 100
  
  # Generates k matrix for each virus.
  k_mat = matrix(NA, nrow = nrow(X), ncol = 1)
  
  # Generates matrix to contain training error from each boosting iteration.
  train_error = matrix(NA, nrow = nrow(X), ncol = n_rounds)
  
  # Generates matrix to contain validation error from each boosting iteration.
  val_error = matrix(NA, nrow = nrow(X), ncol = n_rounds)
  
  # Loop through all virus folds for LOOCV.
  for (fold in 1:n_virus) {
    
    # Fold virus id.
    id_fold = virus_id[fold]
      
    # Virus id for all experiments except the fold virus.
    id_set = virus_id[-which(virus_id %in% id_fold)]
      
    # Independent variables for the fold virus.
    xfold = X[which(virus_id %in% id_fold), , drop = FALSE]
      
    # Inactivation rate for the fold virus.
    yfold = y[which(virus_id %in% id_fold)]
    
    # Weights for the fold virus.
    wfold = wt[which(virus_id %in% id_fold)]
      
    # Independent variables for all experiments except the fold virus.
    xset = X[-which(virus_id %in% id_fold), ]
      
    # Inactivation rate for all experiments except the fold virus.
    yset = y[-which(virus_id %in% id_fold)]
    
    # Weights for all viruses except the fold virus.
    wset = wt[-which(virus_id %in% id_fold)]
    # Re-normalize weights to sum to n_virus in training set for each fold.
    wset = length(wset) * wset / sum(wset)
    
    # Define variables to use in training and validation.
    vars = colnames(xset)
    
    # Create training and validating sets for xgb.train function.
    dtrain = xgb.DMatrix(xset[ , vars, drop = FALSE],
                         label = yset,
                         weight = wset
                         )
    dvalid = xgb.DMatrix(xfold[ , vars, drop = FALSE],
                         label = yfold,
                         weight = wfold)

    watchlist = list(train = dtrain, validation = dvalid)
    
    # Run xgb.train to evaluate training and validation error for each boosting
    # iteration.
    fit_eval = xgb.train(data = dtrain,
                         nround = n_rounds,
                         eta = 0.5,
                         watchlist = watchlist
                        )
    
    # Predicts inactivation rate.
    yval = predict(fit_eval, xfold)
    
    # Training error for each boosting round.
    train_error[fold, ] =
      rep(fit_eval$evaluation_log$train_rmse)
    
    # Validation error for each boosting round.
    val_error[fold, ] =
      rep(fit_eval$evaluation_log$validation_rmse)
    
  }
  
  # Take an average of val_error and train_error by column to determine the 
  # average RMSE for all folds and combine in one matrix.
  train_error_ave = apply(train_error^2, 2, mean)
  train_error_ave = sqrt(train_error_ave)
  val_error_ave = apply(val_error^2, 2, mean)
  val_error_ave = sqrt(val_error_ave)
  
  # Makes a matrix with the training and validation error for each round of
  # boosting.
  rounds_error = cbind(train_error_ave, val_error_ave)
  
  # Run LOOCV for the optimal number of rounds of boosting.
  # Identifies optimal number of rounds of boosting.
  opt_rounds = which(rounds_error[ , 2] == min(rounds_error[ , 2]))
  opt_rounds = opt_rounds[1]
  
  # Loop through all virus folds for LOOCV.
  for (fold in 1:n_virus) {
    cat(fold, '\n')
    # Fold virus id.
    id_fold = virus_id[fold]
    
    # Virus id for all experiments except the fold virus.
    id_set = virus_id[-which(virus_id == id_fold)]
    
    # Number of experiments for the fold virus.
    exp_per_fold = length(which(virus_id == id_fold))
    
    # Independent variables for the fold virus.
    xfold = X[which(virus_id == id_fold), , drop = FALSE]
    
    # Inactivation rate for the fold virus.
    yfold = y[which(virus_id == id_fold)]
    
    # Calculates the mean inactivation rate constant for the fold virus.
    ybar = mean(yfold)
    
    # Independent variables for all experiments except the fold virus.
    xset = X[-which(virus_id == id_fold), ]
    
    # Inactivation rate for all experiments except the fold virus.
    yset = y[-which(virus_id == id_fold)]
    
    # Weights for all viruses except the fold virus.
    wset = wt[-which(virus_id %in% id_fold)]
    # Re-normalize weights to sum to n_virus in training set for each fold.
    wset = length(wset) * wset / sum(wset)
    
    # Run xgb to evaluate training and validation error for each boosting
    # iteration for the optimal number of rounds.
    fit_eval = xgboost(data = xset,
                       label = yset,
                       nrounds = opt_rounds,
                       eta = 0.5,
                       weight = wset
                       )
  
    # Predicts inactivation rate constant.
    yval = predict(fit_eval, xfold)
    
    # Includes yval for this fold.
    k_mat[fold, ] = yval
  }
  
  # Create list to output round and mse data.
  xgb_info = vector(mode = 'list', length = 2)

  # Define information from xgboosting rounds.
  xgb_info[[1]] = k_mat
  xgb_info[[2]] = rounds_error
  
  xgb_info

}


# MODEL OPTIMIZATION------------------------------------------------------------

# Defines vector with virus names for all experiments.
train_id_data = train_data[, id_vars]

# the -1 index removes the intercept, which glmnet should add by default.
train_ind_data = model.matrix( ~ ., train_data[ , ind_vars])[ ,-1]

# Creates dependent data set for all training data.
train_dep_data = train_data[, dep_vars]

# Creates weighted data set for all training data.
w = train_data[ , weight]
# Normalize weights so w sums to number of viruses
w_norm = n_virus * w / sum(w)

# Run xgboost function.
xgb_info = xgb(X = train_ind_data,
               y = train_dep_data,
               virus_id = train_id_data,
               wt = w_norm
               )

# ORGANIZE MODEL OUTPUTS--------------------------------------------------------

# Defines mse for each loocv round using xgb.
k_xgb_all = xgb_info[[1]]
# Defines error from 100 rounds of xgb.
rounds_error = xgb_info[[2]]

# Defines the matrix rows with virus names.
row.names(k_xgb_all) = train_id_data

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
  k_xgb_list = w_list = k_est_list = vector( length = length(uniq_class), mode = 'list')
  names(k_xgb_list) = names(w_list) = names(k_est_list) = uniq_class
  
  for ( class in names(k_xgb_list) ) {
    # N: Check this to make sure it works correctly!!
    idx = which(train_data[ , 'class'] == class_dict_rev[class])
    
    # Hold kpred results in a single list.
    k_xgb_list[[class]] = k_xgb_all[idx, , drop = FALSE] 
    
    # Hold weight info in a single list.
    w_list[[class]] = w[idx] / sum(w[idx])
    
    # Hold estimated k info in a single list.
    k_est_list[[class]] = train_dep_data[idx]
  }
} else {
  uniq_class = data_subset 
  k_xgb_list = list(k_xgb_all) # length one
  names(k_xgb_list) = uniq_class # a single class
  
  w_list = list(w / sum(w)) # length one
  names(w_list) = uniq_class # a single class
  
  k_est_list = list(train_dep_data) # length one
  names(k_est_list) = uniq_class # a single class
}

# Creates a list where I will place the sum, variance, etc. of mse values.
mse_xgb_list = mse_wt_xgb_list = mse_sum_list = mse_var_list = mse_se_list = 
  rmse_sum_list = rmse_se_list = mse_xgb_summ_list = 
  vector(length = length(uniq_class), mode = 'list')
names(mse_xgb_list) = names(mse_wt_xgb_list) = names(mse_sum_list) =
  names(mse_var_list) = names(mse_se_list) = names(rmse_sum_list) = 
  names(rmse_se_list) = names(mse_xgb_summ_list) = uniq_class

for ( class in names(k_xgb_list) ) {
  # Calculates mse for each virus.
  mse_xgb_list[[class]] = (k_xgb_list[[class]] - k_est_list[[class]])^2
  # Calculate weighted mse values for each virus.
  mse_wt_xgb_list[[class]] = mse_xgb_list[[class]] * w_list[[class]]
  # Hold averaged mse results in a single list.
  mse_sum_list[[class]] = sum(mse_wt_xgb_list[[class]] / k_est_list[[class]]^2)
  # Hold var of mse results in a single list. This is called the bias corrected 
  # weighted sample variance.
  mse_var_list[[class]] = 
    sum(w_list[[class]] * (mse_xgb_list[[class]]) / k_est_list[[class]]^2) / (1 - sum(w_list[[class]]^2))
  # Hold se of mse results in a single list.
  mse_se_list[[class]] = sqrt(mse_var_list[[class]] / nrow(mse_xgb_list[[class]]))
  # Hold averaged rmse results in a single list.
  rmse_sum_list[[class]] = sqrt(mse_sum_list[[class]])
  # Hold se of rmse results in a single list.
  rmse_se_list[[class]] = sqrt(mse_se_list[[class]])
  # Summarize the squared and rooted mean, se in the same list.
  mse_xgb_summ_list[[class]] = cbind(mse_sum_list[[class]], mse_se_list[[class]],
                                 rmse_sum_list[[class]], rmse_se_list[[class]])
  # Name columns of each list with mse and rmse mean and se.
  colnames(mse_xgb_summ_list[[class]]) = c('mse_sum',
                                           'mse_se',
                                           'rmse_sum',
                                           'rmse_se'
                                           )
}

# Redefine opt_rounds found from rounds_error
opt_rounds = which(rounds_error[ , 2] == min(rounds_error[ , 2]))
opt_rounds = opt_rounds[1]

# PREDICTION MODELS-------------------------------------------------------------

# Makes a list to contain the predictive model for each virus subset.
xgb_opt = vector( length = length(uniq_class), mode = 'list')
names(xgb_opt) = uniq_class

# Creates the predictive model for each virus subset.
for (class in names(xgb_opt)) {
  xgb_opt[[class]] = xgboost(data = train_ind_data,
                            # label = log(train_dep_data),
                            label = train_dep_data,
                            nrounds = opt_rounds,
                            weight = w_norm
  )
}

# OUTPUT------------------------------------------------------------------------

# Generates an R Data file name to include the predictive model for the data
# subset.
file_model = sprintf('development/development-results/xgb/xgb-opt-%s.RData', data_subset)

# Saves the predictive model list to an RData file.
save(xgb_opt, data_subset, uniq_class, cat_vars, file = file_model)

# Generate an R Data file name for mse training info.
file_mse = sprintf('development/development-results/xgb/xgb-summ-loocv-%s.RData', data_subset)

# Saves the loocv mse information for each virus subset and the
# summary information for each virus subset to an RData file.
save(k_xgb_all, k_xgb_list, mse_xgb_list, mse_xgb_summ_list, uniq_class,
     file = file_mse
)
# Export training/validation error in each boosting round to a csv file.
# file_round_error = sprintf('development/csv-outputs/xgb-rounds-error-%s.csv', data_subset)
# write.csv(rounds_error, file = file_round_error)

# VIRUS-SPECIFIC LOOCV SUMMARY--------------------------------------------------

# Creates a list to place the summary info for each virus class for
# the top model - but for all virus experiments.
summ_xgb_all = vector(length = length(uniq_class), mode = 'list')
names(summ_xgb_all) = uniq_class

# Provides summary information for each class of viruses from
# the top model.
for ( class in names(summ_xgb_all) ) {
  # Hold all the predicted and mean virus rate constants from the top model in
  # one list.
  summ_xgb_all[[class]] =
    tibble(
      virus = rownames(k_xgb_list[[class]]), 
      kpred = k_xgb_list[[class]], 
      kest = k_est_list[[class]],
      mspe = mse_xgb_list[[class]]
    ) %>%
    mutate( rmspe = sqrt(mspe) )
  
  # Generate an R Data file name for mse training info.
  file_summ = sprintf('development/csv-outputs/xgb-summ-kval-%s%s.csv',
                      class,
                      ifelse(data_subset == 'all', '-all', '')
                      )
  
  # Outputs rmspe summary info from top models to a csv - for si table s4,
  # fig 1,si fig s2.
  write.csv(summ_xgb_all[[class]], file_summ)
}


