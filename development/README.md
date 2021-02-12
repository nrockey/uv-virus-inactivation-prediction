# uv-virus-inactivation-prediction

# 'development' folder

# Files in this folder were used to set up model inputs, run model training/validation for four different model classes, and determine model performance of each model.

# This code can be used to reproduce the model training/validation findings we share in the publication: Rockey, Nicole C., Henderson, James B., Chin, Kaitlyn, Raskin, Lutgarde, Wigginton, Krista R., "Predictive Modeling of Virus Inactivation by UV". Environmental Science & Technology, 2021.

# Files

# 01-raw-data-analysis.R
Reads in rate constants, standard errors (if available), viruses, and study ids from the systematic review. Removes outliers from the data set, determines k-bar and weights to use in model training/validation. Defines a final model data set that removes viruses with no standard error values. Finally, removes viruses with no available full-length genome information, and removes viruses that have unique categorical variable values that would not work in loocv.

# 02-uv-inact-input-vars.R
Reads in training/validation virus set and defines independent variables using information from 00 (in 'sequence-set-up' folder) and 01.

# 03-uv-inact-model-data-inputs.R
Creates training data sets to include in model training/validation with virus
information from 02.

# 04-uv-inact-mlr.R
Generates multiple linear regression models using virus genome and biological
function characteristics as predictors. Model accuracy is assessed using LOOCV
on the training/validation set.

# 05-uv-inact-elastic-net.R
Uses glmnet to make a model for predicting UV virus inactivation rate constants. Model accuracy is assessed using LOOCV on the training/validation set.

# 06-uv-inact-boosted-trees.R
Uses xgboost to make a model for predicting UV virus inactivation rate constants. Models are assessed using LOOCV on the training/validation set.

# 07-uv-inact-rand-forest.R
Uses xgboost (with settings similar to those in randomForest) to make a model
for predicting UV virus inactivation rate constants. Models are assessed using LOOCV on the training/validation set.

# mlr-funcs.R
Numerous functions used for the multiple linear regression model training/validation (file 04) - these run loocv with principal components from the predictors used in modeling. Functions include 'npc_func', 'lmpc', and 'lmpc_fit'.