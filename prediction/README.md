# uv-virus-inactivation-prediction

# 'prediction' folder

# Files in this folder were used to set up prediction model inputs and run predictions of virus inactivation rate constants.

# This code can be used to reproduce the predictions we share in the publication: Rockey, Nicole C., Henderson, James B., Chin, Kaitlyn, Raskin, Lutgarde, Wigginton, Krista R., "Predictive Modeling of Virus Inactivation by UV". Environmental Science & Technology, 2021.

# Files

# 08-uv-inact-predict-vars.R
Reads in prediction virus set and defines predictors.

# 09-uv-inact-model-predict-inputs.R
Creates a prediction data set with virus information from 08.

# 10-uv-inact-mlr-predict.R
Uses the optimized mlr model from 04 to predict virus inactivation of the
viruses read in from 08.

# 11-uv-inact-xgb-predict.R
Uses the optimized boosted trees model from 06 to predict virus inactivation of the viruses read in from 08.