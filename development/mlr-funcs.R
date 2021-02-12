# MLR - functions
# Generates multiple linear regression models using features of virus genome
# and biological function characteristics as independent variables. 
# The dimension of genomic features is reduced using PCA.
# Model accuracy is assessed using LOOCV at *virus* level.
#
# The 'lmpc' acronym is for 'linear model w/ principal components'
# 
# 2020.10.19
# Nicole Rockey

# FUNCTIONS---------------------------------------------------------

#Function to define the number of PCs.
npc_func = function(p, n, maxpcs = 3) {
  # p the number of (genomic) columns in the training data subject to PCA
  # n the number of viruses (unique rows) available for training
  pmax(1, pmin(n, pmin(maxpcs, p - 1 ) ) )
}

# Function to define lmpc.
lmpc = function(y, virus_id, X, X_cat, wt, maxpcs = 3) {
  # Do LOOCV for response y using up to maxpcs of the PCs of X.
  # Inputs:
  #  y  - target inactivation rates for training
  #  id - virus id for determining folds
  #  X  - training model matrix with nrow(X) == length(y),
  #               columns in this matrix are subject to PCA
  #  X_cat - training model matrix for 'categorical' variables not subject to PCA
  #  w - weights for use in modeling.
  #  maxpcs - the maximum number of principal components to use
  #  
  # Returns: 
  #  a matrix of "leave one virus out" k predictions,
  #    with dimensions nrow(X) rows by npc columns where
  #    npc = max(1, min(maxpcs, ncol(X) - length(cat_vars))
  
  # Defines # of viruses in data.
  n_virus = length(virus_id)
  
  # Defines maximum # of PCs to go through in models.
  npc = npc_func(ncol(X), n_virus - 1, maxpcs)
  
  # Predicted k storage
  k_mat = matrix(NA, nrow = nrow(X), ncol = npc)
  
  # Loop over validation-sets for LOOCV (so for n viruses).
  for ( fold in 1:n_virus )  {
    # Fold virus id.
    id_fold = virus_id[fold]
    
    # Virus id for all viruses except the fold virus.
    id_set = virus_id[-which(virus_id == id_fold)]
    
    # Inactivation rate constant for the fold virus.
    yfold = y[which(virus_id == id_fold)]
    
    # Inactivation rates for all viruses except the fold virus.
    yset = y[-which(virus_id == id_fold)]
    
    # Weights for all viruses except the fold virus.
    wset = wt[-which(virus_id == id_fold)]
    # Re-normalize weights to sum to n_virus in training set for each fold.
    wset = length(wset) * wset / sum(wset)
    
    # Categorical information for the fold virus.
    cat_fold = X_cat[which(virus_id == id_fold), , drop = FALSE]
    
    # Categorical information for all experiments except the fold virus.
    cat_set = X_cat[-which(virus_id == id_fold), , drop = FALSE]
    
    # Creates training set with virus data except the fold virus from the X
    # data set for PCA, also removes the categorical variable.
    trainpcaX = X[-which(virus_id == id_fold), , drop = FALSE]
    
    # Determine the mean and standard deviation of each variable from the
    # train set before we do PCA analysis.
    mn = apply(trainpcaX, 2, mean)
    sd = apply(trainpcaX, 2, sd)
    
    # Validation set for the fold virus, removes the categorical variables.
    Xval = X[which(virus_id == id_fold), , drop = FALSE]
    
    # Standardize validation set using mean and standard deviation.
    Xval = t( {t(Xval) - mn} / sd)
    
    # Creates rotated data matrix. scale = TRUE standardizes data for
    # PCA analysis.
    pca = prcomp(trainpcaX, scale = TRUE)
    pcaX = pca$x
    
    # Creates matrix with PCA loadings for all variables.
    loadX = pca$rotation
    
    #Linear model fitting either 1, 2, or 3 PCs.
    for ( pcmax in 1:npc ) {
      #Fits a linear model using PCA data.
      pca_fit = lm.wfit(cbind(1, pcaX[ , 1:pcmax, drop = FALSE], cat_set),
                    yset, w = wset)
      
      # Multiply test variables by loadings.
      pcavalX = Xval %*% loadX[, 1:pcmax, drop = FALSE]
      
      #Predict k for the test set based on the model developed.
      yval = cbind(1, pcavalX, cat_fold) %*% pca_fit$coefficients
      
      # Calculating MSE for fold "fold" using y.
      k_mat[fold, pcmax] = yval
    }
  }
    
  # return
  k_mat
  
}

lmpc_fit = function(y, X, X_cat, npc, wt, vars){
  # Fits a least-squares model using PCs of genomic variables.
  # Inputs:
  #  y - target inactivation rates for training data
  #  X  - training model matrix with nrow(X) == length(y),
  #               columns in this matrix are subject to PCA
  #  X_cat - training model matrix for 'categorical' variables not subject to
  #               PCA
  #  npc - number of pcs to use in fit
  #  vars - the *indices* of the genetic variables from X[,-cat_vars]
  #         (would be better to pass by name in the final version)
  #
  # Returns:
  #  mean, sd - means and sds from training data for standardization
  #  loadX - the loading matrix to compute rotated PCs
  #  beta - lm coefficients
  #  vars - the *names* of the vars to use in the PCA
  #  cat_vars - same as input
  
  # Defines training set as only the variables defined in vars.
  trainpcaX = X[ , vars, drop = FALSE]
  
  # Determine the mean and standard deviation of each variable from the
  # train set before we do PCA analysis.
  mn = apply(trainpcaX, 2, mean)
  sd = apply(trainpcaX, 2, sd)
  
  #Creates rotated data matrix.
  pca = prcomp(trainpcaX, scale = TRUE)
  pcaX = pca$x
  
  # Creates matrix with PCA loadings for all variables.
  loadX = pca$rotation

  # lm coefficients for these values with y.
  pca_fit = lm.wfit(cbind(1, pcaX[ , 1:npc, drop = FALSE], X_cat),
                   y, w = wt)
  
  # return
  list(
    mean = mn,
    sd = sd,
    loadX = loadX[ , 1:npc], # subset here as we only need npc components,
    beta = coef(pca_fit),
    vars = colnames(trainpcaX),
    cat_vars = colnames(X_cat)
  )
}

# Here is a function to pull out the relevant information to determine the
# best model. 
decode_name = function(x){
  # a column name identifying a model from mse_* above
  params = stringr::str_split(x, '_')[[1]]
  params = stringr::str_sub(params, 2, -1)
  list( n_vars = as.numeric(params[1]),
        vars = as.numeric(stringr::str_split(params[2], '-')[[1]]),
        npcs = as.numeric(params[3])
  )
}
