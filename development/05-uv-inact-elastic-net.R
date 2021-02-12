# 05 - Elastic net model
# Objectives of this elnt file:
#      1) estimate loocv mse's for all elastic net models.
#      2) summarize mse by virus type.
#      3) generate optimal models for each virus type.
#      4) save all mse summaries and optimal models in R Data files for
#         downstream use.

# 2020.10.19
# Nicole Rockey

# LIBRARIES---------------------------------------------------------------------
library(tidyverse); library(glmnet)

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

# Function to order a matrix based on column.
ord_mat = function(x, col) {
  # inputs: x - a matrix
  #         col - the column to order the matrix by
  # Orders rows in matrix x by column col.
  x[order(x[ , col]), ]
}

# Function to define elnt.
elnt = function(alphas, virus_id, X, y, wt) {
  # Inputs:
  # alphas - hyperparamater to use in model developmenbt.
  # virus_id - virus id for determining folds
  # X - independent training data for loocv.
  # y - target inactivation rates for training

  # Creates a vector with only unique viruses listed.
  uniq_id = unique(virus_id)
  
  # Defines # of viruses in data set.
  n_virus = length(uniq_id)
  
  # Defines # of experiments in data set.
  n_exp = length(virus_id)
  
  # Creates a list to store MSE values for all alpha/lambdas.
  k_elnt_list = vector(mode = 'list',
                       length = length(alphas)
                       )
  
  lambda_list = vector(length = length(alphas), mode = 'list')
  n_lambda = vector(length = length(alphas), mode = 'list')
  
  # Goes through each alpha assigned.
  # Need this to test out different alphas to see what the optimal alpha value
  # to use will be. 
  for (i in 1:length(alphas)) {

    # Sets alpha for this loop at the index of i.
    a = alphas[i]
    
    # Currently including weights even in the alpha selection.
    model_prep = glmnet(x = X,
                        y = y,
                        alpha = a,
                        weights = wt
                        )
    
    # Define the set of lambda values to be used in hyperparameter
    # optimization.
    lambda = model_prep$lambda
      
    # Generates k matrix for all lambdas for this alpha for each virus.
    k_mat = matrix(NA, nrow = nrow(X), ncol = length(lambda))
    # colnames(mse_mat) = modelnames
    
    # Loop through all virus folds for LOOCV.
    for (fold in 1:n_virus) {
      # Fold virus id.
      id_fold = uniq_id[fold]
      
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
      
      # Weights for all viruses excpet the fold virus.
      wset = wt[-which(virus_id == id_fold)]
      # Re-normalize weights to sum to n_virus in training set for each fold.
      wset = length(wset) * wset / sum(wset)
      
      # Create an elastic net model with one alpha value and a vector of 100
      # different lambda values.
      fit = glmnet(x = xset,
                   y = yset,
                   alpha = a,
                   lambda = lambda,
                   weights = wset
                  )
      
      # Predicts inactivation rate.
      yval = predict(fit, xfold)
      
      # Output predicted rate constant for this virus fold and alpha/lambda.
      k_mat[fold, ] = yval
      
    } # ends folds loop
    
    # Save the lambda information in a list.
    lambda_list[[i]] = lambda
    n_lambda[[i]] = length(lambda)
    # Save the set of k_mat in the first level of the list for alpha.
    k_elnt_list[[i]] = k_mat
    
  } # ends alpha loop
  
  elnt_info = vector(length = 3, mode = 'list')
  
  #! J - I'd suggest using a named list, e.g. below, to make this easier
  #      to remember
  elnt_info[['lambda_list']] = lambda_list
  elnt_info[['n_lambda']] = n_lambda
  elnt_info[['k_elnt_list']] = k_elnt_list
  
  # returns
  elnt_info
}

# VARIABLE SET-UP & MODEL OPTIMIZATION------------------------------------------

# Creates a numeric that increments from 0 to 1 with a 0.1 step.
alphas = seq(0, 1, 0.1)

# Defines the number of lambda values desired.
n_lambda = 100

# Defines vector with virus names for all experiments.
train_id_data = train_data[, id_vars]

# Creates independent data set for all training data.
train_ind_data = model.matrix( ~ ., train_data[ , ind_vars])[ ,-1]

# Creates dependent data set for all training data.
train_dep_data = train_data[ , dep_vars]

# Creates weights for all training data.
w = train_data[ , weight]
names(w) = rownames(train_data)

# Normalize the weights so w sums to number of viruses.
w_norm = n_virus * w / sum(w)

# Run elnt function.
elnt_info = elnt(alphas = alphas,
                 X = train_ind_data,
                 y = train_dep_data,
                 virus_id = train_id_data,
                 wt = w_norm
                 )

# ORGANIZE MODEL OUTPUTS--------------------------------------------------------

# Combines all models made for each alpha/lambda combination into a matrix.
k_elnt_all = do.call("cbind", elnt_info[['k_elnt_list']])
row.names(k_elnt_all) = train_id_data
# Gets all lambdas from all models developed.
lambda_list = elnt_info[['lambda_list']]
# Combines all lambda values for model into a vector.
lambda_all = do.call("c", lambda_list)

# Make a matrix with all alpha values.
alphas_each = matrix(ncol = length(alphas), nrow = n_lambda)

for (i in 1:length(alphas)) {
  alphas_each[1:(length(lambda_list[[i]])), i] = alphas[i]
}

# Make a vector with all alpha values.
alpha_all = vector(length = ncol(k_elnt_all), mode = 'numeric')
alpha_all = c(alphas_each)
alpha_all = na.omit(alpha_all)

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
  k_elnt_list = w_list = k_est_list = vector( length = length(uniq_class), mode = 'list')
  names(k_elnt_list) = names(w_list) = names(k_est_list) = uniq_class
  
  for ( class in names(k_elnt_list) ) {
    idx = which(train_data[, 'class'] == class_dict_rev[class])
    
    # Hold mse results in a single list
    k_elnt_list[[class]] = k_elnt_all[idx, ] 
    
    # Hold weight info in a single list.
    w_list[[class]] = w[idx] / sum(w[idx])
    
    # Hold estimated k info in a single list.
    k_est_list[[class]] = train_dep_data[idx]
  }
} else {
  uniq_class = data_subset 
  k_elnt_list = list(k_elnt_all) # length one
  names(k_elnt_list) = uniq_class # a single class
  
  w_list = list(w / sum(w)) # length one
  names(w_list) = uniq_class # a single class
  
  k_est_list = list(train_dep_data) # length one
  names(k_est_list) = uniq_class
}

# Creates lists where I will place mse values and summary stats, and names those
# lists.
mse_elnt_list = mse_wt_elnt_list = mse_sum_list = mse_var_list = mse_se_list = 
  rmse_sum_list = rmse_se_list = mse_elnt_summ_list = mse_elnt_order_list = 
  vector(length = length(uniq_class), mode = 'list')
names(mse_elnt_list) = names(mse_wt_elnt_list) = names(mse_sum_list) =
  names(mse_var_list) = names(mse_se_list) = names(rmse_sum_list) =
  names(rmse_se_list) = names(mse_elnt_summ_list) =
  names(mse_elnt_order_list) = uniq_class

for ( class in names(mse_elnt_list) ) {
  # Calculates mse for each virus.
  mse_elnt_list[[class]] = (k_elnt_list[[class]] - k_est_list[[class]])^2
  # Calculate weighted mse values for each virus.
  mse_wt_elnt_list[[class]] = (mse_elnt_list[[class]] * w_list[[class]])
  # Hold summed mse results in a single list.
  mse_sum_list[[class]] = apply(mse_wt_elnt_list[[class]] / k_est_list[[class]]^2, 2, sum)
  # Hold var of mse results in a single list. This is the biased corrected
  # weighted variance.
  mse_var_list[[class]] =
    apply((w_list[[class]] * mse_elnt_list[[class]] / k_est_list[[class]]^2), 2, sum) / (1 - sum(w_list[[class]]^2))
  # Hold se of mse results in a single list.
  mse_se_list[[class]] = sqrt(mse_var_list[[class]] / nrow(mse_elnt_list[[class]]))
  # Hold rmse of sum mse results in a single list.
  rmse_sum_list[[class]] = sqrt(mse_sum_list[[class]])
  # Hold sd of rmse results in a single list.
  rmse_se_list[[class]] = sqrt(mse_se_list[[class]])
  # Summarize the mse, rmse, error, alpha, and lambda in the same list.
  mse_elnt_summ_list[[class]] = cbind(mse_sum_list[[class]], mse_se_list[[class]],
                                 rmse_sum_list[[class]], rmse_se_list[[class]],
                                 alpha_all, lambda_all)
  # Name columns of each list with mse and rmse mean and sd.
  colnames(mse_elnt_summ_list[[class]]) = c('mse_mean',
                                            'mse_se',
                                            'rmse_mean',
                                            'rmse_se',
                                            'alpha',
                                            'lambda'
                                       )
  # Orders values in the matrix so model summary information is shown from
  # best-performing (i.e., lowest average mse across loocv rounds) to worst-
  # performing (ie.., highest average mse across loocv rounds) for each virus
  # subset.
  mse_elnt_order_list[[class]] = ord_mat(x = mse_elnt_summ_list[[class]], col = 1)
}

# PREDICTION MODELS-------------------------------------------------------------

# Makes a list to contain the predictive model for each virus subset.
elnt_opt = vector( length = length(uniq_class), mode = 'list')
names(elnt_opt) = uniq_class

# Creates the predictive model for each virus subset.
for (class in names(elnt_opt)) {
  elnt_opt[[class]] = glmnet(x = train_ind_data,
                            y = train_dep_data,
                            alpha = mse_elnt_order_list[[class]][1, 'alpha'],
                            lambda = mse_elnt_order_list[[class]][1, 'lambda'],
                            weights = w_norm
                            )
}

# DATA OUTPUT-------------------------------------------------------------------

# Generates an R Data file name to include the predictive model for the data
# subset.
file_model = sprintf('development/development-results/elnt/elnt-opt-%s.RData', data_subset)

# Saves the predictive model list to an RData file.
save(elnt_opt, data_subset, uniq_class, class_vars, file = file_model)

# Generate an R Data file name for mse training info.
file_mse = sprintf('development/development-results/elnt/elnt-summ-loocv-%s.RData', data_subset)

# Saves the loocv mse information for each virus subset and the
# summary information for each virus subset to an RData file.
save(k_elnt_all, k_elnt_list, mse_elnt_list, mse_wt_elnt_list, mse_elnt_summ_list,
     mse_elnt_order_list, uniq_class, file = file_mse
     )

# VIRUS-SPECIFIC LOOCV SUMMARY--------------------------------------------------

# Creates a list to place the summary info for each virus class for
# the top model - but for all virus experiments.
summ_elnt_all = vector(length = length(uniq_class), mode = 'list')
names(summ_elnt_all) = uniq_class

# Provides summary information for each class of viruses from
# the top model.
for ( class in names(summ_elnt_all) ) {
  # Find the top model in the elnt data set for this class.
  top_model = which(mse_elnt_summ_list[[class]][ , 'alpha'] ==
                      mse_elnt_order_list[[class]][1, 'alpha'] &
                      mse_elnt_summ_list[[class]][ , 'lambda'] ==
                      mse_elnt_order_list[[class]][1, 'lambda'])
  
  # Hold all the predicted and mean virus rate constants from the top model in
  # one list.
  summ_elnt_all[[class]] =
    tibble(
      virus = rownames(k_elnt_list[[class]]), 
      kpred = k_elnt_list[[class]][ , top_model], 
      kest = k_est_list[[class]],
      mspe = mse_elnt_list[[class]][ , top_model]
    ) %>%
    mutate( rmspe = sqrt(mspe) )
  
  # Generate an R Data file name for mse training info.
  file_summ = sprintf('development/csv-outputs/elnt-summ-kval-%s%s.csv',
                      class,
                      ifelse(data_subset == 'all', '-all', '')
                      )
  
  # Outputs rmspe summary info from top models to a csv - for si table s4,
  # fig 1,si fig s2.
  write.csv(summ_elnt_all[[class]], file_summ)
}

