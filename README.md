# uv-virus-inactivation-prediction

#The objective of this work is to identify a model that most accurately predicts the UV inactivation of different viruses.

# Files

# 00-uv-inact-sequence-attributes.R
Defines virus genome sequence information.

# 01-raw-data-analysis.R
Reads in rate constants, standard errors (if available), viruses, and study ids
from the systematic review. Removes outliers from the data set, determines k-bar
and weights to use in model development. Defines a final model data set that
removes viruses with no standard error values. Finally, removes viruses with no
available full-length genome information, and removes viruses that have unique
categorical variable values that would not work in loocv.

# 02-uv-inact-input-vars.R
Reads in training/validation virus set and defines independent variables using
information from 00 and 01.

# 03-uv-inact-model-data-inputs.R
Creates training data sets to include in model training/validation with virus
information from 02.

# 04-uv-inact-mlr.R
Generates multiple linear regression models using virus genome and biological
function characteristics as predictors. Model accuracy is assessed using LOOCV
on the training/validation set.

# 05-uv-inact-elastic-net.R
Uses glmnet to make a model for predicting UV virus inactivation rates. Model
accuracy is assessed using LOOCV on the training/validation set.

# 06-uv-inact-boosted-trees.R
Uses xgboost to make a model for predicting UV virus inactivation rates.
Models are assessed using LOOCV on the training/validation set.

# 07-uv-inact-rand-forest.R
Uses xgboost (with settings similar to those in randomForest) to make a model
for predicting UV virus inactivation rates. Models are assessed using LOOCV on the training/validation set.

# 08-uv-inact-predict-vars.R
Reads in prediction virus set and defines independent variables.

# 09-uv-inact-model-predict-inputs.R
Creates a prediction data set with virus information from 08.

# 10-uv-inact-mlr-predict.R
Uses the optimized mlr model from 04 to predict virus inactivation of the
viruses read in from 08.

# 12-uv-inact-elnt-predict.R
Uses the optimized elastic net model from 05 to predict virus inactivation of the viruses read in from 08.

# 12-uv-inact-xgb-predict.R
Uses the optimized boosted trees model from 06 to predict virus inactivation of the viruses read in from 08.

# 13-uv-inact-rf-predict.R
Uses the optimized random forests model from 07 to predict virus inactivation of the viruses read in from 08.