# 03 - Data inputs - UV Inactivation models
# Uses data set generated in 02 to set up training set, independent
# variables, dependent variables, categorical variables, and class variables to
# be used in model training/validation.

# 2020.10.25
# Nicole Rockey

# SET WORKING DIRECTORY---------------------------------------------------------

setwd("~/uv-virus-inactivation-prediction/")

# DATA INPUT--------------------------------------------------------------------

# Loads in 'tot_data.'
load("data/virus-inact-data.RData")

# DATA MANIPULATION-------------------------------------------------------------

# Total number of viruses.
n_virus = nrow(tot_data)

# Creates a vector that identifies the columns with identifiers.
id_vars = 'virus'

# Creates a vector that identifies the columns with the dependent variable.
dep_vars = 'k_bar'

# Creates a vector that identifies the columns with the error for the
# dependent variable.
dep_error = 'k_var'

# Creates a vector with the categorical virus variables.
cat_vars = c('type', 'repair', 'host')
names(cat_vars) = cat_vars

# Creates a vector with the virus classification variable.
class_vars = c('class')
names(class_vars) = class_vars

# Creates a vector with the weight.
weight = c('w')

# Creates a vector with independent variables.
ind_vars = setdiff(names(tot_data), c(id_vars, class_vars, dep_vars, dep_error, weight))

# Replaces the text for the virus classification variable with numerics. 
# (dsdna = 0, plusssrna = 1, ssdna = 2, negssrna = 3, dsrna = 4)
tot_data[ which(tot_data[ , class_vars[1]] == "dsdna"), class_vars[1]] = 0
tot_data[ which(tot_data[ , class_vars[1]] == "plusssrna"), class_vars[1]] = 1
tot_data[ which(tot_data[ , class_vars[1]] == "ssdna"), class_vars[1]] = 2
tot_data[ which(tot_data[ , class_vars[1]] == "negssrna"), class_vars[1]] = 3
tot_data[ which(tot_data[ , class_vars[1]] == "dsrna"), class_vars[1]] = 4

# Replaces the text for the virus classification variable with numerics.
# (ds = 0, ss = 1)
tot_data[ which(tot_data[ , cat_vars['type']] == "ds"), cat_vars['type']] = 0
tot_data[ which(tot_data[ , cat_vars['type']] == "ss"), cat_vars['type']] = 1

# MAKING SUBSETS OF DATA AND OUTPUTS--------------------------------------------

# All viruses.
# Creates the training set.
train_data = tot_data

# Defines the subset of data included in modeling.
data_subset = 'all'

# Number of viruses.
n_virus = length(unique(train_data$virus))

# Generates file with all viruses for model training.
save(id_vars, dep_vars, dep_error, weight, ind_vars, class_vars, cat_vars,
     train_data, n_virus, data_subset, 
     file = "data/uv-inact-model-data-inputs-all.RData")

# dsDNA viruses.
#Creates the training set.
train_data = subset(tot_data, class == 0)
# Removes any columns that do not have multiple unique values.
n_unique = apply(train_data, 2, function(x) length( unique(x) ) ) # count unique values for each column
train_data = train_data[, n_unique > 1]
  
# Define class_vars for dsDNA subset.
class_vars = NULL

# Define cat_vars for dsDNA subset.
cat_vars = c('repair', 'host')
names(cat_vars) = cat_vars

# Creates a vector with independent variables.
ind_vars = setdiff(names(train_data), c(id_vars, class_vars, dep_vars, dep_error, weight))
  
# Number of viruses.
n_virus = length(unique(train_data$virus))

# Defines the subset of data included in modeling.
data_subset = 'dsdna'

# Generates file with dsdna viruses for model training.
save(id_vars, dep_vars, dep_error, weight, ind_vars, class_vars, cat_vars,
     train_data, n_virus, data_subset,
     file = "data/uv-inact-model-data-inputs-dsdna.RData")

# ssDNA viruses.
#Creates the training set.
train_data = subset(tot_data, class == 2)
# Removes any columns that do not have multiple unique values.
n_unique = apply(train_data, 2, function(x) length( unique(x) ) ) # count unique values for each column
train_data = train_data[ , n_unique > 1]

# Define class_vars for ssDNA subset.
class_vars = NULL

# Define cat_vars for ssDNA subset.
cat_vars = NULL

# Creates a vector with independent variables.
ind_vars = setdiff(names(train_data), c(id_vars, dep_vars, dep_error, weight))

# Number of viruses.
n_virus = length(unique(train_data[ , id_vars]))

# Defines the subset of data included in modeling.
data_subset = 'ssdna'

# Generates file with ssdna viruses for model training.
save(id_vars, dep_vars, dep_error, weight, ind_vars, class_vars, cat_vars,
     train_data, n_virus, data_subset,
     file = "data/uv-inact-model-data-inputs-ssdna.RData")

# dsRNA viruses.
#Creates the training set.
train_data = subset(tot_data, class == 4)
# Removes any columns that do not have multiple unique values.
n_unique = apply(train_data, 2, function(x) length( unique(x) ) ) # count unique values for each column
train_data = train_data[, n_unique > 1]

# Define class_vars for dsRNA subset.
class_vars = NULL

# Define cat_vars for dsRNA subset.
cat_vars = NULL

# Creates a vector with independent variables.
ind_vars = setdiff(names(train_data), c(id_vars, dep_vars, dep_error, weight))

# Number of viruses.
n_virus = length(unique(train_data[ , id_vars]))

# Defines the subset of data included in modeling.
data_subset = 'dsrna'

# Generates file with dsdna viruses for model training.
save(id_vars, dep_vars, dep_error, weight, ind_vars, class_vars, cat_vars,
     train_data, n_virus, data_subset,
     file = "data/uv-inact-model-data-inputs-dsrna.RData")

# (+) ssRNA viruses.
#Creates the training set.
train_data = subset(tot_data, class == 1)
# Removes any columns that do not have multiple unique values.
n_unique = apply(train_data, 2, function(x) length( unique(x) ) ) # count unique values for each column
train_data = train_data[, n_unique > 1]
# Shouldn't need this if we set the host to FALSE in our files.
train_data = train_data[ , -2]

# Define class_vars for (+) ssRNA subset.
class_vars = NULL

# Define cat_vars for (+) ssRNA subset.
cat_vars = NULL

# Creates a vector with independent variables.
ind_vars = setdiff(names(train_data), c(id_vars, dep_vars, dep_error, weight))

# Number of viruses.
n_virus = length(unique(train_data[ , id_vars]))

# Defines the subset of data included in modeling.
data_subset = 'plusssrna'

# Generates file with (+) ssrna viruses for model training.
save(id_vars, dep_vars, dep_error, weight, ind_vars, class_vars, cat_vars,
     train_data, n_virus, data_subset,
     file = "data/uv-inact-model-data-inputs-plusssrna.RData")

# (-) ssRNA viruses.
#Creates the training set.
train_data = subset(tot_data, class == 3)
# Removes any columns that do not have multiple unique values.
n_unique = apply(train_data, 2, function(x) length( unique(x) ) ) # count unique values for each column
train_data = train_data[, n_unique > 1]

# Define class_vars for (-) ssRNA subset.
class_vars = NULL

# Define cat_vars for (-) ssRNA subset.
cat_vars = NULL

# Creates a vector with independent variables.
ind_vars = setdiff(names(train_data), c(id_vars, dep_vars, dep_error, weight))

# Number of viruses.
n_virus = length(unique(train_data[ , id_vars]))

# Defines the subset of data included in modeling.
data_subset = 'negssrna'

# Generates file with (-) ssrna viruses for model training.
save(id_vars, dep_vars, dep_error, weight, ind_vars, class_vars, cat_vars,
     train_data, n_virus, data_subset,
     file = "data/uv-inact-model-data-inputs-negssrna.RData")

# DNA (ssDNA and dsDNA) viruses.
#Creates the training set.
train_data = subset(tot_data, class == 0 | class == 2)
# Removes any columns that do not have multiple unique values.
n_unique = apply(train_data, 2, function(x) length( unique(x) ) ) # count unique values for each column
train_data = train_data[, n_unique > 1]

# Define class_vars for DNA subset.
class_vars = 'class'
names(class_vars) = class_vars

# Define cat_vars for DNA subset.
cat_vars = c('type', 'repair', 'host')
names(cat_vars) = cat_vars

# Creates a vector with independent variables.
ind_vars = setdiff(names(train_data), c(id_vars, class_vars, dep_vars, dep_error, weight))

# Number of viruses.
n_virus = length(unique(train_data[ , id_vars]))

# Defines the subset of data included in modeling.
data_subset = 'dna'

# Generates file with dna viruses for model training.
save(id_vars, dep_vars, dep_error, weight, ind_vars, class_vars, cat_vars,
     train_data, n_virus, data_subset,
     file = "data/uv-inact-model-data-inputs-dna.RData")

# RNA (dsRNA, (+) ssRNA, and (-) ssRNA) viruses.
#Creates the training set.
train_data = subset(tot_data, class == 1 | class == 3 | class == 4)
# Removes any columns that do not have multiple unique values.
n_unique = apply(train_data, 2, function(x) length( unique(x) ) ) # count unique values for each column
train_data = train_data[ , n_unique > 1]

# Define class_vars for RNA subset.
class_vars = 'class'
names(class_vars) = class_vars

# Define cat_vars for RNA subset.
cat_vars = 'type'
names(cat_vars) = cat_vars

# Creates a vector with independent variables.
ind_vars = setdiff(names(train_data), c(id_vars, dep_vars, dep_error, weight))

# Number of viruses.
n_virus = length(unique(train_data[ , id_vars]))

# Defines the subset of data included in modeling.
data_subset = 'rna'

# Generates file with rna viruses for model training.
save(id_vars, dep_vars, dep_error, weight, ind_vars, class_vars, cat_vars,
     train_data, n_virus, data_subset, 
     file = "data/uv-inact-model-data-inputs-rna.RData")

# ssRNA ((+) ssRNA and (-) ssRNA) viruses.
#Creates the training set.
train_data = subset(tot_data, class == 1 | class == 3)
# Removes any columns that do not have multiple unique values.
n_unique = apply(train_data, 2, function(x) length( unique(x) ) ) # count unique values for each column
train_data = train_data[, n_unique > 1]

# Define class_vars for ssRNA subset.
class_vars = 'class'
names(class_vars) = class_vars

# Define cat_vars for ssRNA subset.
cat_vars = NULL

# Creates a vector with independent variables.
ind_vars = setdiff(names(train_data), c(id_vars, class_vars, dep_vars, dep_error, weight))

# Number of viruses.
n_virus = length(unique(train_data[ , id_vars]))

# Defines the subset of data included in modeling.
data_subset = 'ssrna'

# Generates file with ssrna viruses for model training.
save(id_vars, dep_vars, dep_error, weight, ind_vars, class_vars, cat_vars,
     train_data, n_virus, data_subset,
     file = "data/uv-inact-model-data-inputs-ssrna.RData")

