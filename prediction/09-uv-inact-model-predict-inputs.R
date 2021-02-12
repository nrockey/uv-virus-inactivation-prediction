# 09 - Prediction data inputs - UV inactivation models
# Creates optimized models from all training data and generates prediction data
# sets to include in final model.

# 2020.10.19
# Nicole Rockey

# SET WORKING DIRECTORY---------------------------------------------------------

setwd("~/uv-virus-inactivation-prediction/")

#DATA INPUT---------------------------------------------------------------------

# Loads in 'pred_data.'
load("data/virus-predict-data.RData")

# DATA MANIPULATION-------------------------------------------------------------

# Creates a vector that identifies the columns with identifiers.
id_vars_pred = 'virus'

# Defines vector with virus names for all predictions.
pred_id = pred_data[ , id_vars_pred]

# Creates a vector with the virus classification.
class_vars_pred = 'class'
names(class_vars_pred) = class_vars_pred

# Creates a vector with the virus categorical variables.
cat_vars_pred = c('type', 'host', 'repair')
names(cat_vars_pred) = cat_vars_pred

# Creates a vector with independent variables.
ind_vars_pred = setdiff(names(pred_data), c(id_vars_pred, class_vars_pred))

# Replaces the text for the virus classifcation variable with binaries 
# (dsdna = 0, plusssrna = 1, ssdna = 2, negssrna = 3, dsrna = 4)
pred_data[which(pred_data[ , class_vars_pred[1]] == "dsdna"), class_vars_pred[1]] = 0
pred_data[which(pred_data[ , class_vars_pred[1]] == "plusssrna"), class_vars_pred[1]] = 1
pred_data[ which(pred_data[ , class_vars_pred[1]] == "ssdna"), class_vars_pred[1]] = 2
pred_data[ which(pred_data[ , class_vars_pred[1]] == "negssrna"), class_vars_pred[1]] = 3
pred_data[ which(pred_data[ , class_vars_pred[1]] == "dsrna"), class_vars_pred[1]] = 4

# Replaces the text for the virus nucleic acid type variable with numerics.
# (ds = 0, ss = 1)
pred_data[ which(pred_data[ , cat_vars_pred['type']] == "ds"), cat_vars_pred['type']] = 0
pred_data[ which(pred_data[ , cat_vars_pred['type']] == "ss"), cat_vars_pred['type']] = 1

# OUTPUTS-----------------------------------------------------------------------

# Outputs prediction data as an R Data file.
save(id_vars_pred,
     ind_vars_pred,
     cat_vars_pred,
     class_vars_pred,
     pred_data,
     file = "data/uv-inact-model-predict-inputs.RData"
    )
