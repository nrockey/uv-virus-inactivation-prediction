# 01 - Raw systematic review data analysis
# Takes the raw systematic review data and organizes/summarizes it so that
# the final modeling data set includes no outliers, no study rate constants
# without standard error estimates of k, no viruses for which there is no
# complete genome information, and no viruses that have a unique categorical
# variable type. Estimates of k and the weights for each virus are determined.

# 2020.10.19
# Nicole Rockey

# LIBRARIES---------------------------------------------------------------------

library(tidyverse)

# SET WORKING DIRECTORY---------------------------------------------------------

setwd('~/uv-virus-inactivation-prediction/')

# DATA INPUT--------------------------------------------------------------------

sys_review_data = read.csv('data/sys-review-data.csv', header = TRUE, sep = ",",
                    stringsAsFactors = FALSE)

# FUNCTIONS---------------------------------------------------------------------

# Function to remove outliers from the  virus set.
is_outlier = function(k) {
  # If k is length 3 or greater, returns TRUE whenever k is more than 1.5 IQRs
  #  above the upper quartile or below the lower quartile
  if ( length(k) < 3 ) {
    return( rep(FALSE, length(k)) )
  } else {
    return(
      k > quantile(k, .75) + 1.5 * IQR(k) |
        k < quantile(k, .25) - 1.5 * IQR(k)
    )
  }
}

# DATA MANIPULATION-------------------------------------------------------------

colnames(sys_review_data) = c('study', 'virus', 'k', 'se', 'class')

# Creates a vector that identifies the columns with identifiers.
id_vars = 'virus'

# Creates a vector that identifies the columns with the dependent variable.
dep_vars = 'k'

# Systematic review data without outliers.
sys_review_data_no_outliers = sys_review_data %>%
  #flag and remove outlying estimates
  group_by(virus) %>%
  mutate( outlier = is_outlier(k) ) %>%
  filter( !outlier )

# Manually remove FR from Study 24. Although not removed using the outlier
# approach above (because there are not a lot of rate constants for this virus),
# this rate constant is quite far from the other two rate constants, we know it
# is not accurate. Remove.
sys_review_data_no_outliers = sys_review_data_no_outliers[-which(sys_review_data_no_outliers$study == 24 & sys_review_data_no_outliers$virus == "FR" ), ]

# REMOVE VIRUSES WITHOUT ERROR ESTIMATES----------------------------------------

# Only systematic review data that includes error estimates.
sys_review_data_se = sys_review_data_no_outliers %>%
  # flag and remove outlying estimates
  group_by(virus) %>%
  filter( !is.na(se) )

# INVERSE VARIANCE WEIGHTING BY VIRUS-------------------------------------------

# Inverse variance weights.
model_data = sys_review_data_se %>%
  mutate( w = 1 / se^2 )

# This is the most classical approach to determining the estimated k. It assumes
# all studies are noisy estimaes of the same effect and no study heterogeneity.
model_data = model_data %>%
  group_by(virus) %>%
  mutate(
    kbar_w = sum(w * k)  / sum(w),
    kbar_w_se = sqrt(1 / sum(w)),
    kbar_w_lwr = kbar_w - qnorm(.975) * kbar_w_se,
    kbar_w_upr = kbar_w + qnorm(.975) * kbar_w_se
  ) %>%
  # order by k
  arrange(kbar_w) 

model_data$virus_o = with(model_data, factor(virus, unique(virus)))

model_data %>%
  group_by(virus) %>%
  mutate( 
    study_weight = w / sum(w)  # shade according to weight of experiment for virus
  ) %>% 
  ungroup() %>% 
  ggplot(aes( y = virus_o, x = kbar_w) ) +
  geom_point( pch = 15, color = 'darkgrey', size = 1.5 ) +
  geom_errorbarh( aes(xmin = kbar_w_lwr, xmax = kbar_w_upr), 
                  color = 'darkgrey') +
  geom_point( aes(x = k, color = virus, alpha = study_weight), pch = 16) +
  ylab('') + 
  xlab('k') +
  scale_x_log10() + 
  theme_bw() +
  guides(color = 'none', alpha = guide_legend(title = 'Relative Study Weight'))

# REMOVE VIRUSES WTIH NO COMPLETE GENOME SEQUENCES OR UNIQUE CAT VARS-----------

# For current model, we cannot include the following viruses. Will remove.
model_data_final = model_data %>%
  # flag and remove outlying estimates
  group_by(virus) %>%  
  # Manually remove T4, PMV.
  # Need to remove because these viruses have unique categorical values. LOOCV
  # will not work with these viruses in the data set.
  filter( !{ virus %in% c('T4', 'PMV_JC')}) %>%
  # Manually remove BCV, F2, PMV_Toronto, SA_A994, SP8.
  # Need to remove because they do not have genome sequences available.
  filter( !{ virus %in% c('BCV', 'F2', 'PMV_Toronto', 'SA_A994', 'SP8', 'T1UV')})

# SUMMARIZE---------------------------------------------------------------------
# Systematic review - all data.
sys_review_data_summ_virus = sys_review_data %>%
  group_by(virus) %>%
  summarize( n = n(), .groups = 'drop')
sys_review_data_summ_class = sys_review_data %>%
  group_by(class) %>%
  summarize( n = n(), .groups = 'drop')

# Systematic review without outliers (used in Figure 1).
outliers = sys_review_data %>%
  group_by(virus) %>%
  mutate( outlier = is_outlier(k) ) %>%
  filter( outlier )

sys_review_data_no_outliers_summ_virus = sys_review_data_no_outliers %>%
  group_by(virus) %>%
  summarize(k_mean = mean(k), n = n(), .groups = 'drop')
sys_review_data_no_outliers_summ_class = sys_review_data_no_outliers %>%
  group_by(class) %>%
  summarize(k_mean = mean(k), n = n(), .groups = 'drop')

# Modeling data - only rate constants with se included.
model_data_summ_virus = model_data %>%
  group_by(virus) %>%
  summarize( across(starts_with("k"), mean), n = n(), .groups = 'drop') 
model_data_summ_class = model_data %>%
  group_by(class) %>%
  summarize( n = n(), .groups = 'drop')

# Final modeling data set - remove viruses with unique cat vars or no genome
# sequence (used in model training and validation).
# summarize final model data set.
model_data_final_summ_virus = model_data_final %>%
  group_by(virus) %>%
  summarize( across(starts_with("k"), mean), n = n(), .groups = 'drop') 
model_data_final_summ_class = model_data_final %>%
  group_by(class) %>%
  summarize(  n = n(), .groups = 'drop') 

k_bar = model_data_final_summ_virus$kbar_w
k_var = model_data_final_summ_virus$kbar_w_se^2
names(k_bar) = names(k_var) = model_data_final_summ_virus$virus 

# OUTPUT------------------------------------------------------------------------
save(id_vars, dep_vars, k_var, k_bar, model_data_final, sys_review_data_se, sys_review_data_no_outliers,
     file = 'data/uv-inact-meta-k-bar.RData')

# Data to include in Figure 1 (this includes all experimental rate constants
# from systematic review except for the outliers).
write.csv(sys_review_data_no_outliers, 'development/csv-outputs/fig-1.csv')

