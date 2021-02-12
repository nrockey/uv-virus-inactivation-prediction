# 08   - UV inactivation virus data set for predicting viruses.
# Takes raw sequence data from 00 and creates matrix containing all pertinent
# virus genome attributes and additional characteristics to be used for
# predictive modeling prediction using optimal models.

# 2020.10.19
# Nicole Rockey

# LIBRARIES---------------------------------------------------------------------

library(stringr); library(tibble)

# SET WORKING DIRECTORY---------------------------------------------------------

setwd("~/uv-virus-inactivation-prediction/")

# DATA INPUT--------------------------------------------------------------------

# Loads in 'seq_extract' matrix containing all virus sequence information.
load("data/virus-seq-attributes.RData")

# Reads in csv file with viruses to use in predictions and additional
# independent variables (genome repair ability, virus genome type) beyond
# sequence info.
adnl_predict_virus_vars = read.csv("data/adnl-predict-virus-vars.csv",
                                   header = TRUE, sep = ",",
                                   stringsAsFactors = FALSE
                                   )

# DATA MANIPULATION-------------------------------------------------------------

# Number of virus inputs to include in final virus data set.
num_virus = nrow(adnl_predict_virus_vars)

num_adnl_vars = ncol(adnl_predict_virus_vars)

class = "class"

# Defines all independent variables to include from sequence information.
predict_seq_vars = c("length","T", "TT", "TTT", "TTTT", "TTTTT", "C", "CT", "TC",
                 "U", "UU","UUU","UUUU","UUUUU", "CU" , "UC"
                )

# Define empty matrix tot_seq_data to fill in for all viruses in data set.
predict_seq_data = matrix(NA, nrow = num_virus, ncol = length(predict_seq_vars))

colnames(predict_seq_data) = predict_seq_vars

type = numeric(length = num_virus)

# Loops through each virus in the virus variables data set.
for (i in 1:num_virus) {

  # Name of the virus for this loop.
  vir_name = adnl_predict_virus_vars[i, 1]
  
  # Returns the row number of the virus name so it can be used below.
  seq_extr_row = which(row.names(seq_extract) == vir_name)
  
  # Places data in correct locations in matrix. 
  if (adnl_predict_virus_vars[i, class] == "dsdna") {
    # Adds information on whether the virus has ds or ss nucleic acid.
    type[i] = 'ds'
    # Records sequence information for virus i in seq_info matrix if it is DNA.
    predict_seq_data[i, ] =  c (seq_extract[seq_extr_row, 1],
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
  if (adnl_predict_virus_vars[i, class] == "ssdna") {
    # Adds information on whether the virus has ds or ss nucleic acid.
    type[i] = 'ss'
    # Records sequence information for virus i in seq_info matrix if it is ssDNA.
    predict_seq_data[i, ] =  c (seq_extract[seq_extr_row, 1],
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
  if (adnl_predict_virus_vars[i, class] == "negssrna") {
    # Adds information on whether the virus has ds or ss nucleic acid.
    type[i] = 'ss'
    # Records sequence information for virus i in seq_info matrix if it is
    # negative sense ssRNA.
    predict_seq_data[i, ] =  c (seq_extract[seq_extr_row, 1],
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
  if (adnl_predict_virus_vars[i, class] == "plusssrna") {
    # Adds information on whether the virus has ds or ss nucleic acid.
    type[i] = 'ss'
    # Records sequence information for virus i in seq_info matrix if it is RNA.
    predict_seq_data[i, ] =  c (seq_extract[seq_extr_row, 1],
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
  if (adnl_predict_virus_vars[i, class] == "dsrna") {
    # Adds information on whether the virus has ds or ss nucleic acid.
    type[i] = 'ds'
    # Records sequence information for virus i in seq_info matrix if it is dsRNA.
    predict_seq_data[i, ] =  c (seq_extract[seq_extr_row, 1],
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

# Combines sequence attributes with additional virus variables.
pred_data = cbind(adnl_predict_virus_vars, predict_seq_data)

pred_data = add_column(pred_data, type, .after = 'class')

# Make 'repair' and 'host' categorical.
pred_data$repair = as.character(pred_data$repair)
pred_data$host = as.character(pred_data$host)

# DATA OUTPUT-------------------------------------------------------------------

# Outputs data as an R file.
save(pred_data, file = "data/virus-predict-data.RData")

