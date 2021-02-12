# uv-virus-inactivation-prediction

# The objective of this work is to identify a model that most accurately predicts the UV inactivation of different viruses.

# The folders included in this repository contain raw data, R files, and results from model development and prediction that are described in our study: Rockey, Nicole C., Henderson, James B., Chin, Kaitlyn, Raskin, Lutgarde, Wigginton, Krista R., "Predictive Modeling of Virus Inactivation by UV". Environmental Science & Technology, 2021.

# Four of these folders ('data', 'sequence-set-up', 'development', 'prediction') contain sets of data, R files, and results that were used to obtain publication findings. These files have been revised from their original versions to improve readability and workflow.

# A fifth folder ('prediction-examples') contains example code and associated data files that can be used to predict the UV rate constant(s) of (a) virus(es) of interest.

# Folders

# 'data'
This folder contains raw data that was used in model development and prediction, including txt files with genome sequence information and csv files with additional virus attributes such as virus class (e.g., dsdna, ssrna, etc.).

# 'sequence-set-up'
This folder contains an R file for reading in genome sequences from specified txt files included in the 'data' folder and collecting genome sequence information to be used in modeling work.

# 'development'
This folder contains R files with code used to set up virus data sets, develop models using the data, and evaluate model performance. The folder also contains results (csv and RData files) generated from the code to easily access findings.

# 'prediction'
This folder contains R files with code used to set up virus prediction data sets and predict inactivation of those viruses. The folder also contains results (csv and RData files) generated from the code to easily access findings.

# 'prediction-examples'
This folder contains example R files and input data for predicting the UV inactivation rate constants for viruses of interest.
