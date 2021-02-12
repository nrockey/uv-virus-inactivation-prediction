# 02 - UV inactivation virus data set for training/testing viruses.
# Takes raw sequence data from 00 and rate constant data from 01 to create a
# matrix containing all pertinent virus genome attributes and additional
# characteristics to be used for predictive modeling training/validation.

# 2020.10.25
# Nicole Rockey

# LIBRARIES---------------------------------------------------------------------
 
library(stringr); library(tibble); library(mefa4)

# SET WORKING DIRECTORY---------------------------------------------------------

setwd("~/uv-virus-inactivation-prediction/")

# DATA INPUT--------------------------------------------------------------------

# Loads in 'seq_extract' matrix containing all virus sequence information.
load("data/virus-seq-attributes.RData")

# Reads in csv file with viruses to use in modeling work, along with
# inactivation rates (dependent variable) and additional independent variables
# (genome repair ability, virus genome class) beyond sequence info.
adnl_virus_vars = read.csv("data/adnl-virus-vars.csv", header = TRUE, sep = ",",
                           stringsAsFactors = FALSE
                           )

# Loads in 'k_bar' and 'k_var' which contain the estimate and variance of the
# virus inactivation rate constants based on the data collected from the 
# systematic review.
foo = load('data/uv-inact-meta-k-bar.RData') 

# DATA MANIPULATION-------------------------------------------------------------

# Number of virus inputs to include in final virus data set.
num_virus = length(k_bar)

# Defines a matrix to update the set of viruses so they match with the viruses
# for which we have 'k_est' and 'k_var'.
virus_vars = matrix(data = NA, nrow = num_virus, ncol = 4)
rownames(virus_vars) = names(k_bar)
colnames(virus_vars) = c('virus', 'class', 'repair', 'host')
virus_vars[ , 'virus'] = rownames(virus_vars)

for (i in 1:num_virus) {
  virus_vars[i, 'class'] = adnl_virus_vars[which(adnl_virus_vars[ , 'virus'] ==
                                                   rownames(virus_vars)[i]),
                                           which(colnames(adnl_virus_vars) == 'class')]
  virus_vars[i, 'repair'] = adnl_virus_vars[which(adnl_virus_vars[, 'virus'] ==
                                                    rownames(virus_vars)[i]),
                                            which(colnames(adnl_virus_vars) == 'repair')]
  virus_vars[i, 'host'] = adnl_virus_vars[which(adnl_virus_vars[, 'virus'] ==
                                                    rownames(virus_vars)[i]),
                                            which(colnames(adnl_virus_vars) == 'host')]
}

# Independent variables we include besides the genomic variables.
num_adnl_vars = ncol(virus_vars)

# Defines all independent variables to include from sequence information.
tot_seq_vars = c("length","T", "TT", "TTT", "TTTT", "TTTTT", "C", "CT", "TC",
                 "U", "UU","UUU","UUUU","UUUUU", "CU" , "UC"
                )

# Define empty matrix tot_seq_data to fill in for all viruses in data set.
tot_seq_data = matrix(NA, nrow = num_virus, ncol = length(tot_seq_vars))

colnames(tot_seq_data) = tot_seq_vars
rownames(tot_seq_data) = names(virus_vars)

type = numeric(length = num_virus)

# Loops through each virus in the virus variables data set.
for (i in 1:num_virus) {
  # Name of the virus for this loop.
  vir_name = virus_vars[i, 1]
  
  # Returns the row number of the virus name so it can be used below.
  seq_extr_row = which(row.names(seq_extract) == vir_name)
  
  # Places data in correct locations in matrix.
  # For dsDNA viruses.
  if (virus_vars[i, 'class'] == "dsdna") {
    # Adds information on whether the virus has ds or ss nucleic acid.
    type[i] = 'ds'
    # Records sequence information for virus i in seq_info matrix if it is dsDNA.
    tot_seq_data[i, ] =  c (seq_extract[seq_extr_row, 1],
                            seq_extract[seq_extr_row, 2],
                            seq_extract[seq_extr_row, 3],
                            seq_extract[seq_extr_row, 4],
                            seq_extract[seq_extr_row, 5],
                            seq_extract[seq_extr_row, 6],
                            seq_extract[seq_extr_row, 7],
                            seq_extract[seq_extr_row, 8],
                            seq_extract[seq_extr_row, 9],
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0)
    
  }
  # For ssDNA viruses.
  if (virus_vars[i, 'class'] == "ssdna") {
    # Adds information on whether the virus has ds or ss nucleic acid.
    type[i] = 'ss'
    # Records sequence information for virus i in seq_info matrix if it is ssDNA.
    tot_seq_data[i, ] =  c (seq_extract[seq_extr_row, 1],
                            seq_extract[seq_extr_row, 25],
                            seq_extract[seq_extr_row, 26],
                            seq_extract[seq_extr_row, 27],
                            seq_extract[seq_extr_row, 28],
                            seq_extract[seq_extr_row, 29],
                            seq_extract[seq_extr_row, 15],
                            seq_extract[seq_extr_row, 30],
                            seq_extract[seq_extr_row, 31],
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0
                            )
  }
  # For (-) ssRNA viruses.
  if (virus_vars[i, 'class'] == "negssrna") {
    # Adds information on whether the virus has ds or ss nucleic acid.
    type[i] = 'ss'
    # Records sequence information for virus i in seq_info matrix if it is
    # negative sense ssRNA.
    tot_seq_data[i, ] =  c (seq_extract[seq_extr_row, 1],
                            0,
                            0,
                            0,
                            0,
                            0,
                            seq_extract[seq_extr_row, 37],
                            0,
                            0,
                            seq_extract[seq_extr_row, 32],
                            seq_extract[seq_extr_row, 33],
                            seq_extract[seq_extr_row, 34],
                            seq_extract[seq_extr_row, 35],
                            seq_extract[seq_extr_row, 36],
                            seq_extract[seq_extr_row, 38],
                            seq_extract[seq_extr_row, 39]
                           )
  }
  
  # For (+) ssRNA viruses.
  if (virus_vars[i, 'class'] == "plusssrna") {
    # Adds information on whether the virus has ds or ss nucleic acid.
    type[i] = 'ss'
    # Records sequence information for virus i in seq_info matrix if it is ssRNA.
    tot_seq_data[i, ] =  c (seq_extract[seq_extr_row, 1],
                            0,
                            0,
                            0,
                            0,
                            0,
                            seq_extract[seq_extr_row, 15],
                            0,
                            0,
                            seq_extract[seq_extr_row, 10],
                            seq_extract[seq_extr_row, 11],
                            seq_extract[seq_extr_row, 12],
                            seq_extract[seq_extr_row, 13],
                            seq_extract[seq_extr_row, 14],
                            seq_extract[seq_extr_row, 16],
                            seq_extract[seq_extr_row, 17]
                            )
  }
  # For dsRNA viruses.
  if (virus_vars[i, 'class'] == "dsrna") {
    # Adds information on whether the virus has ds or ss nucleic acid.
    type[i] = 'ds'
    
    # Records sequence information for virus i in seq_info matrix if it is dsRNA.
    tot_seq_data[i, ] =  c (seq_extract[seq_extr_row, 1],
                            0,
                            0,
                            0,
                            0,
                            0,
                            seq_extract[seq_extr_row, 7],
                            0,
                            0,
                            seq_extract[seq_extr_row, 18],
                            seq_extract[seq_extr_row, 19],
                            seq_extract[seq_extr_row, 20],
                            seq_extract[seq_extr_row, 21],
                            seq_extract[seq_extr_row, 22],
                            seq_extract[seq_extr_row, 23],
                            seq_extract[seq_extr_row, 24]
                            )

  }
  
}

# Calculate the modified inverse variance weights for these viruses to use during
# training and validation of models.This is the squared inverse coefficient of
# variation, e.g., inverse relative variance. 
## We will use k_bar^2/k_var as weight because we do not want viruses with
# smaller k values to have higher weights.
w = k_bar^2 / k_var

# Trim weights (Winsorize) to make values beyond the 90th percentile have a
# maximum value.
w_90 = quantile(w, 0.9)
w[w > w_90] = w_90

# Combines sequence attributes with additional virus variables.
tot_data = as.data.frame(cbind(
                k_bar,
                k_var,
                w,
                tot_seq_data
                ))

# Add in virus_vars data (categorical data).
tot_data = add_column(tot_data, virus_vars[, 4], .before = 1)
colnames(tot_data)[1] = 'host'
tot_data = add_column(tot_data, virus_vars[, 3], .before = 1)
colnames(tot_data)[1] = 'repair'
tot_data = add_column(tot_data, virus_vars[, 2], .before = 1)
colnames(tot_data)[1] = 'class'

# Add single-stranded or double-stranded info to the training data.
tot_data = add_column(tot_data, type, .after = 'class')

# Add in column for virus names.
tot_data = add_column(tot_data, virus_vars[, 1], .before = k_bar)
colnames(tot_data)[1] = 'virus'

# Make 'repair' and 'host' categorical.
tot_data$repair = as.character(tot_data$repair)
tot_data$host = as.character(tot_data$host)

# DATA OUTPUT-------------------------------------------------------------------

# Outputs data as an R file.
save(tot_data, file = "data/virus-inact-data.RData")

