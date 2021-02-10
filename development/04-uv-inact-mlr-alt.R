# 04 - MLR Model-alt
# Objectives of this mlr file:
#      1) estimate loocv mse's for all mlr models.
#      2) summarize mse by virus type.
#      3) generate optimal models for each virus type.
#      4) save all mse summaries and optimal models in R Data files for
#         downstream use.

# 2020.10.25
# Nicole Rockey

# LIBRARIES---------------------------------------------------------------------

library(tidyverse); library(glmnet)

# SET WORKING DIRECTORY---------------------------------------------------------

setwd("~/uv-virus-inactivation-prediction/")

# FUNCTIONS---------------------------------------------------------------------

# Sources the function file.
source("alternate/mlr-funcs.R")

ord_mat = function(x, col) {
  # inputs: x - a matrix
  #         col - the column to order the matrix by
  # Orders rows in matrix x by column col.
  x[order(x[, col]), ]
}

# DATA INPUTS-------------------------------------------------------------------

# N: Currently manually toggling through raw data to load. Could make this 
# automated once everything is good to go.

# ACTION ITEM: Update which of these lines are and are not commented to develop
# model with different data.
# Loads in 'id_vars, dep_vars, ind_vars, train_data, n_virus, n_exp.'
# foo_bar = load("data/uv-inact-model-data-inputs-all.RData")
# load("data/uv-inact-model-data-inputs-dsdna.RData")
foo = load("data/uv-inact-model-data-inputs-plusssrna.RData")

#J: If you want to loop over all of these ...
#J: Let me know if you want to run these in parallel instead.
#all_data_subsets = c('all', 'dsna', 'ssdna', ...)
#for ( data_subset in all_data_subsets ) {
#  load( sprintf("data/uv-inact-model-data-inputs-%s.RData", data_subset) )
#} # move this to the very end of the script and then indent everything.

# DATA MANIPULATION-------------------------------------------------------------

# Manually toggle whether these categorical variables should be included in
# the model or removed.
include_repair = FALSE
include_host = FALSE

if (include_repair == FALSE) {
  # Drops repair from ind_vars.
  ind_vars = setdiff(ind_vars, 'repair')
  cat_vars = setdiff(cat_vars, 'repair')
  if (is_empty(cat_vars) == TRUE) {
    cat_vars = NULL
  }
}

if (include_host == FALSE) {
  # Drops host cell from ind_vars.
  ind_vars = setdiff(ind_vars, 'host')
  cat_vars = setdiff(cat_vars, 'host')
  if (is_empty(cat_vars) == TRUE) {
    cat_vars = NULL
  }
}

# Number of independent variables to be used in modeling.
tot_vars = length(ind_vars)

# Counts number of class variables.
n_cat_vars = length(cat_vars)

# Number of models per n variables.
models_pernvar = npc_func(1:(tot_vars - n_cat_vars), n_virus) *
  choose((tot_vars - n_cat_vars), 1:(tot_vars - n_cat_vars))

# Number of models.
n_models = sum(npc_func(1:(tot_vars - n_cat_vars), n_virus) *
                 choose((tot_vars - n_cat_vars),
                        1:(tot_vars - n_cat_vars)))

# Creates an empty vector with mse_allmodels 
k_models_list = vector(mode = 'list',
                       length = (tot_vars - n_cat_vars))

# Creates an empty vector with allmodel_names.
allmodel_names = c()

# Defines vector with virus names for all training_data.
train_id_data = train_data[ , id_vars]

# Creates dependent data set for all training data.
train_dep_data = train_data[, dep_vars]
names(train_dep_data) = rownames(train_data)

# Creates a weights data set for all training data.
w = train_data[ , weight]
names(w) = rownames(train_data)
# Normalize weights so w sums to number of viruses
w_norm = n_virus * w / sum(w)

# Removes all cat_vars from training set in training set without class_vars.
train_ind_no_cat = 
  model.matrix( ~ ., train_data[ , setdiff(ind_vars, cat_vars)] 
  )[ , -1]

# Creates the training matrix including only class_vars
if ( !is.null(cat_vars) ) {
  train_ind_cat = 
    model.matrix( ~., train_data[ , cat_vars, drop = FALSE] )[
      , -1, drop = FALSE] 
} else {
  train_ind_cat = NULL
}

# Loop over # of variables in model
for ( n_vars in 1:(tot_vars - n_cat_vars) ) {
  # Defines total number of possible variable combinations.
  sets = combn(1:(tot_vars - n_cat_vars), n_vars)
  
  # Makes list for each variable combination, defined by "sets."
  k_sets = vector(mode = 'list', length = ncol(sets))
  
  # Loop over combinations of the "sets" variables.
  for ( set in 1:ncol(sets) )  {
    X_ind = train_ind_no_cat[ , sets[ , set], drop = FALSE]
    
    # Run lmpc function.
    k_sets[[set]] = lmpc(y = train_dep_data, 
                         virus_id = train_id_data,
                         X = X_ind, 
                         X_cat = train_ind_cat, 
                         wt = w_norm,
                         maxpcs = 3)
  } # ends sets loop
  
  # Defines set_names
  set_names = apply(sets, 2, paste, collapse = "-")
  
  # Defines the list of model names.
  model_names = sprintf("n%i_v%s_p%%s", n_vars, set_names)
  
  # Updates the list of model names to include a copy of the model name for
  # the total number of PCs for that model.
  model_names = rep(model_names, each = npc_func(n_vars, n_virus - 1) )
  
  # Updates the list of model names to include the # of PCs in the particular
  # model.
  model_names = sprintf(model_names, 1:npc_func(n_vars, n_virus - 1) )
  
  # Adds  model_names to allmodel_names for each n_vars.
  allmodel_names = c(allmodel_names, model_names)
  
  # Makes the k_sets list into a matrix.
  k_models_list[[n_vars]] = do.call("cbind", k_sets)
  colnames(k_models_list[[n_vars]]) = model_names
  
} # ends n_vars loop

# The list k_models_list has length 16 (# of genetic vars) and each
# component is a choose(16 [genetic vars], n_vars) by 3 (pcs) matrix.

# Binds the list of k values into one large matrix
# This matrix has dimensions n_exp x choose(16, 1) + 2 * choose(16, 2) +
#   3 * sum_{k=3}^16 (choose(16, k)

k_mlr_all = do.call("cbind", k_models_list)
colnames(k_mlr_all) = allmodel_names
row.names(k_mlr_all) = train_id_data

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
  k_mlr_list = w_list = k_est_list = vector( length = length(uniq_class), mode = 'list')
  names(k_mlr_list) = names(w_list) = names(k_est_list) = uniq_class
  
  for ( class in names(k_mlr_list) ) {
    idx = which(train_data[, 'class'] == class_dict_rev[class])
    
    # Hold k results in a single list
    k_mlr_list[[class]] = k_mlr_all[idx, ]
    
    # Hold weight info in a single list.
    w_list[[class]] = w[idx] / sum(w[idx])
    
    # Hold estimated k info in a single list.
    k_est_list[[class]] = train_dep_data[idx]
  }
} else {
  uniq_class = data_subset 
  k_mlr_list = list(k_mlr_all) # length one
  names(k_mlr_list) = uniq_class # a single class
  
  w_list = list(w / sum(w)) # length one
  names(w_list) = uniq_class # a single class
  
  k_est_list = list(train_dep_data) # length one
  names(k_est_list) = uniq_class
}

# Creates lists to place mse values and summary stats, and names those
# lists.
mse_mlr_list = mse_wt_mlr_list = mse_sum_list = mse_var_list =
  mse_se_list = rmse_sum_list = rmse_se_list = mse_mlr_summ_list =
  mse_mlr_order_list = vector(length = length(uniq_class), mode = 'list')
names(mse_mlr_list) = names(mse_sum_list) = names(mse_var_list) =
  names(mse_se_list) = names(rmse_sum_list) = names(rmse_se_list) =
  names(mse_mlr_summ_list) = names(mse_mlr_order_list) = uniq_class

for ( class in names(mse_mlr_list) ) {
  # Calculates mse for each virus.
  mse_mlr_list[[class]] = (k_mlr_list[[class]] - k_est_list[[class]])^2
  # Calculate weighted mse values for each virus.
  mse_wt_mlr_list[[class]] = (mse_mlr_list[[class]] * w_list[[class]])
  # Hold averaged mse results in a single list.
  mse_sum_list[[class]] = apply(mse_wt_mlr_list[[class]] / k_est_list[[class]]^2, 2, sum)
  # Hold variance of mse results in a single list. This is the bias
  # corrected weighted sample variance.
  mse_var_list[[class]] =
    apply((w_list[[class]] * mse_mlr_list[[class]] / k_est_list[[class]]^2), 2, sum) / (1 - sum(w_list[[class]]^2))
  # Hold standard error of mse results in a single list.
  mse_se_list[[class]] = sqrt(mse_var_list[[class]] / nrow(mse_mlr_list[[class]]))
  # Hold averaged rmse results in a single list.
  rmse_sum_list[[class]] = sqrt(mse_sum_list[[class]])
  # Hold sd of rmse results in a single list.
  rmse_se_list[[class]] = sqrt(mse_se_list[[class]])
  # Summarize the mean and sd in the same list.
  mse_mlr_summ_list[[class]] = cbind(mse_sum_list[[class]], mse_se_list[[class]],
                                     rmse_sum_list[[class]], rmse_se_list[[class]])
  # Name columns of each list with relative mean and se.
  colnames(mse_mlr_summ_list[[class]]) = c('mse_sum',
                                           'mse_se',
                                           'rmse_sum',
                                           'rmse_se'
  )
  # Orders values in the matrix so model summary information is shown from
  # best-performing (i.e., lowest average mse across loocv rounds) to
  # worst-performing (ie.., highest average mse across loocv rounds) for
  # each virus subset.
  mse_mlr_order_list[[class]] = ord_mat(x = mse_mlr_summ_list[[class]], col = 1)
}

# PREDICTION MODELS-------------------------------------------------------------

# Makes a list to contain top models for each virus type.
mlr_opt_params = vector( length = length(uniq_class), mode = 'list')
names(mlr_opt_params) = names(mse_mlr_list)

# Pull out the top model for each subset of viruses.
for ( class in names(mse_mlr_order_list) ) {
  mlr_opt_params[[class]] = decode_name(rownames(mse_mlr_order_list[[class]])[1])
} 

# Now make all the optimal models for each virus subset.
# Makes a list to contain the predictive model for each virus subset.
mlr_opt = vector( length = length(uniq_class), mode = 'list')
names(mlr_opt) = uniq_class

# Creates the predictive model for each virus subset.
for ( class in uniq_class ) {
  
  mlr_opt[[class]] = lmpc_fit(y = train_dep_data,
                              X = train_ind_no_cat,
                              X_cat = train_ind_cat,
                              wt = w_norm,
                              npc = mlr_opt_params[[class]]$npcs,
                              vars = mlr_opt_params[[class]]$vars
  )
}

# DATA OUTPUT-------------------------------------------------------------
# Generates an R Data file name to include the predictive model for the
# data subset.
file_model = sprintf('results/mlr/mlr-opt-rep-%s-host-%s-%s.RData',
                     include_repair,
                     include_host,
                     data_subset)

# Saves the predictive model list to an RData file.
save(mlr_opt, data_subset, uniq_class, cat_vars, include_repair,
     include_host,
     file = file_model
)

# Generate an R Data file name for mse training info.
file_mse = sprintf('results/mlr/mlr-summ-loocv-rep-%s-host-%s-%s.RData',
                   include_repair,
                   include_host,
                   data_subset
)

# Saves the loocv mse information for each virus subset and the
# summary information for each virus subset to an RData file.
save(k_mlr_all, k_mlr_list, mse_mlr_list, mse_mlr_summ_list, mse_mlr_order_list,
     uniq_class, file = file_mse)

# VIRUS-SPECIFIC LOOCV SUMMARY--------------------------------------------

# Creates a list to place the summary info for each virus class for
# the top model - but for all virus experiments.
summ_mlr_all = vector(length = length(uniq_class), mode = 'list')
names(summ_mlr_all) = uniq_class

# Provides summary information for each class of viruses from
# the top model.
for ( class in names(summ_mlr_all) ) {
  # Find the top model in the elnt data set for this class.
  top_model = which(row.names(mse_mlr_summ_list[[class]]) ==
                      row.names(mse_mlr_order_list[[class]])[1])
  
  # Hold all the predicted and mean virus rate constants from the top
  # model in one list.
  summ_mlr_all[[class]] =
    tibble(
      virus = rownames(k_mlr_list[[class]]), 
      kpred = k_mlr_list[[class]][ , top_model], 
      kest = k_est_list[[class]],
      spe = mse_mlr_list[[class]][ , top_model]
    ) %>%
    mutate( rspe = sqrt(spe) ) %>%
    mutate( perc_error = rspe/kest)
  
  # Generate an R Data file name for mse training info.
  file_summ = sprintf('csv-outputs/mlr-summ-kval-%s%s.csv',
                      class,
                      ifelse(data_subset == 'all', '-all', '')
  )
  
  # Outputs rmspe summary info from top models to a csv - for si table s4,
  # fig 1,si fig s2.
  write.csv(summ_mlr_all[[class]], file_summ)
} # end virus summary

