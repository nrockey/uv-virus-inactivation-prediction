# 00 - UV inactivation virus sequence attributes
# Creates matrix with raw sequence information of all viruses for which we have
# genomes using ncbi fasta files.

# 2020.10.19
# Nicole Rockey

# LIBRARIES---------------------------------------------------------------------

library(stringr)

# SET WORKING DIRECTORY---------------------------------------------------------

# Set the working directory.
setwd("~/uv-virus-inactivation-prediction/data")

# DATA INPUT--------------------------------------------------------------------

# Input all the file names that should be assessed for sequence information.
# Make sure to include the full file name EXACTLY the way it appears in the
# 'data' folder.
files = c("AdV1.txt", "AdV2.txt", "AdV5.txt", "AdV6_Tonsil99.txt", "AdV15.txt",
          "AdV40_Dugan.txt", "AdV41_Tak.txt", "B40-8.txt", "lambda.txt",
          "PMV_JC.txt", "PRD-1.txt", "T1.txt", "T2.txt", "T4.txt", "T5.txt",
          "T6.txt", "T7.txt", "T7M.txt", "AHNV-RNA1.txt",
          "AHNV-RNA2.txt","CCV.txt", "CVB3.txt", "CVB5_Faulkner.txt",
          "CVB6_Schmitt.txt", "EV1_Farouk.txt", "EV11.txt", "EV12.txt",
          "EMC.txt", "FCV.txt", "FR.txt", "GA.txt", "HAV.txt", "HEV.txt",
          "MS2.txt", "MNV_CW3.txt",  "PV1_Sabin.txt", "PV1_Mahoney.txt",
          "QB.txt", "VEEV.txt", "GII.4_Sydney.txt", "MERS_CoV.txt",
          "SARS_CoV_1.txt", "SARS_CoV_2.txt", "mhv-a59-1000.txt",
          "phiX174.txt", "parvoH1.txt", "fd.txt", "phi6-s.txt", "phi6-m.txt",
          "phi6-l.txt", "reo1-l1.txt", "reo1-l2.txt", "reo1-l3.txt",
          "reo1-m1.txt", "reo1-m2.txt", "reo1-m3.txt", "reo1-s1.txt",
          "reo1-s2.txt", "reo1-s3.txt", "reo1-s4.txt", "reo3-l1.txt",
          "reo3-l2.txt", "reo3-l3.txt", "reo3-m1.txt", "reo3-m2.txt",
          "reo3-m3.txt", "reo3-s1.txt", "reo3-s2.txt", "reo3-s3.txt",
          "reo3-s4.txt", "RVSA11-1.txt", "RVSA11-2.txt", "RVSA11-3.txt",
          "RVSA11-4.txt", "RVSA11-5.txt", "RVSA11-6.txt", "RVSA11-7.txt",
          "RVSA11-8.txt", "RVSA11-9.txt", "RVSA11-10.txt", "RVSA11-11.txt",
          "IPNV-A.txt", "IPNV-B.txt", "ISAV-1.txt", "ISAV-2.txt", "ISAV-3.txt",
          "ISAV-4.txt", "ISAV-5.txt", "ISAV-6.txt", "ISAV-7.txt", "ISAV-8.txt",
          "VHSV.txt", "dengue.txt", "zika.txt", "hs2.txt", "ebv.txt",
          "hcmv.txt", "hsv-1.txt", "hrv.txt", "variola.txt"
          )

# DATA MANIPULATION-------------------------------------------------------------

# Number of different sequences to assess.
num_files = length(files)

# Sequence information to extract from text files.
seq_vars = c("Length", "T_both", "TT_both", "TTT_both", "TTTT_both", "TTTTT_both",
             "C_both", "CT_both", "TC_both",
             "U", "UU", "UUU", "UUUU", "UUUUU", "C", "CU", "UC", "U_both", "UU_both",
             "UUU_both", "UUUU_both", "UUUUU_both", "CU_both", "UC_both", "T", "TT",
             "TTT", "TTTT", "TTTTT", "CT", "TC", "negsense_U", "negsense_UU",
             "negsense_UUU", "negsense_UUUU", "negsense_UUUUU", "negsense_C",
             "negsense_CU", "negsense_UC"
             )
  
# Define empty matrix seq_data to fill in with sequence information.
seq_extract = matrix(NA, nrow = num_files, ncol = length(seq_vars))

# Define row names
rownames(seq_extract) = str_remove(files,".txt")

# Define column names
colnames(seq_extract) = seq_vars

# Loop through each file and count and record sequence information in each file.
for (i in 1:num_files) {
  # Read in each file, separated by line.
  eachfile = scan(files[i], character(0), sep = "\n")
  
  # Count sequence information.
  length = as.numeric( nchar( eachfile[2] ))
  num_T_total = as.numeric ( str_count ( eachfile[2],"T" ))
  num_TT_total = as.numeric ( str_count ( eachfile[2],"TT" ))
  num_TTT_total = as.numeric ( str_count ( eachfile[2],"TTT" ))
  num_TTTT_total = as.numeric ( str_count ( eachfile[2],"TTTT" ))
  num_TTTTT_total = as.numeric ( str_count ( eachfile[2],"TTTTT" ))
  num_C_total = as.numeric ( str_count ( eachfile[2],"C" ))
  num_CT_total = as.numeric ( str_count ( eachfile[2],"CT" ))
  num_TC_total = as.numeric ( str_count ( eachfile[2],"TC" ))
  
  num_U = num_T_total
  num_UU = num_TT_total - num_TTT_total - num_TTTT_total
  num_UUU = num_TTT_total - num_TTTT_total
  num_UUUU = num_TTTT_total - num_TTTTT_total
  num_UUUUU = num_TTTTT_total
  num_C = num_C_total
  num_CU = num_CT_total
  num_UC = num_TC_total
  
  # Make sure to account for the other strand.
  num_A_total = as.numeric ( str_count ( eachfile[2],"A" ))
  num_AA_total = as.numeric ( str_count ( eachfile[2],"AA" ))
  num_AAA_total = as.numeric ( str_count ( eachfile[2],"AAA" ))
  num_AAAA_total = as.numeric ( str_count ( eachfile[2],"AAAA" ))
  num_AAAAA_total = as.numeric ( str_count ( eachfile[2],"AAAAA" ))
  num_G_total = as.numeric ( str_count ( eachfile[2],"G" ))
  num_GA_total = as.numeric ( str_count ( eachfile[2],"GA" ))
  num_AG_total = as.numeric ( str_count ( eachfile[2],"AG" ))
  
  # Sum up both strands and record sequence information for dsDNA.
  num_T_both = num_T_total + num_A_total
  num_TT_both = (num_TT_total - num_TTT_total - num_TTTT_total) + 
    (num_AA_total - num_AAA_total - num_AAAA_total)
  num_TTT_both = (num_TTT_total - num_TTTT_total) +
    (num_AAA_total - num_AAAA_total)
  num_TTTT_both = (num_TTTT_total - num_TTTTT_total) +
    (num_AAAA_total - num_AAAAA_total)
  num_TTTTT_both = (num_TTTTT_total) + (num_AAAAA_total)
  
  num_C_both = num_C_total + num_G_total
  
  num_CT_both = num_AG_total + num_CT_total
  num_TC_both = num_GA_total + num_TC_total
  
  # Record information needed for ssDNA.
  num_T = num_T_total
  num_TT = num_TT_total - num_TTT_total - num_TTTT_total
  num_TTT = num_TTT_total - num_TTTT_total
  num_TTTT = num_TTTT_total - num_TTTTT_total
  num_TTTTT = num_TTTTT_total

  num_CT = num_CT_total
  num_TC = num_TC_total
  
  # Record information needed for dsRNA.
  num_U_both = num_T_total + num_A_total
  num_UU_both = (num_TT_total - num_TTT_total - num_TTTT_total) + 
    (num_AA_total - num_AAA_total - num_AAAA_total)
  num_UUU_both = (num_TTT_total - num_TTTT_total) +
    (num_AAA_total - num_AAAA_total)
  num_UUUU_both = (num_TTTT_total - num_TTTTT_total) +
    (num_AAAA_total - num_AAAAA_total)
  num_UUUUU_both = (num_TTTTT_total) + (num_AAAAA_total)
  
  num_C_both = num_C_total + num_G_total
  
  num_CU_both = num_AG_total + num_CT_total
  num_UC_both = num_GA_total + num_TC_total
  
  # Record information needed for (-)ssRNA.
  num_sense_U = num_A_total
  num_sense_UU = num_AA_total - num_AAA_total - num_AAAA_total
  num_sense_UUU = num_AAA_total - num_AAAA_total
  num_sense_UUUU = num_AAAA_total - num_AAAAA_total
  num_sense_UUUUU = num_AAAAA_total
  num_sense_C = num_G_total
  num_sense_UC = num_AG_total
  num_sense_CU = num_GA_total
  
  
  # Records sequence information for virus i in seq_info matrix.
  seq_extract[i,] =  c(length, num_T_both, num_TT_both, num_TTT_both,
                       num_TTTT_both, num_TTTTT_both, num_C_both, num_CT_both,
                       num_TC_both, num_U, num_UU, num_UUU, num_UUUU, num_UUUUU,
                       num_C, num_CU, num_UC, num_U_both, num_UU_both,
                       num_UUU_both, num_UUUU_both, num_UUUUU_both, num_CU_both,
                       num_UC_both, num_T, num_TT, num_TTT, num_TTTT, num_TTTTT,
                       num_CT, num_TC, num_sense_U, num_sense_UU, num_sense_UUU,
                       num_sense_UUUU, num_sense_UUUUU, num_sense_C, num_sense_CU,
                       num_sense_UC
                      )
}
  
# Combine segmented viruses (AHNV, IPNV, ISAV, reo1, reo3, RVSA11, phi6).
# Need to combine segmented genomes - combine AHNV from two lines to one.
seq_extract[which(rownames(seq_extract) == "AHNV-RNA1"), ] =
  colSums(seq_extract[c(which(rownames(seq_extract) == "AHNV-RNA1"),
                        which(rownames(seq_extract) == "AHNV-RNA2")) , ])
# Delete the RNA2 row now.
rownames(seq_extract)[which(rownames(seq_extract) == "AHNV-RNA1")] = "AHNV"
seq_extract = seq_extract[-which(rownames(seq_extract) == "AHNV-RNA2"),]

# Need to combine segmented genomes - combine IPNV from two lines to one.
seq_extract[which(rownames(seq_extract) == "IPNV-A"), ] =
  colSums(seq_extract[c(which(rownames(seq_extract) == "IPNV-A"),
                        which(rownames(seq_extract) == "IPNV-B")) , ])
# Delete the B row now.
rownames(seq_extract)[which(rownames(seq_extract) == "IPNV-A")] = "IPNV"
seq_extract = seq_extract[-which(rownames(seq_extract) == "IPNV-B"),]

# Need to combine segmented genomes - combine ISAV lines to one.
seq_extract[which(rownames(seq_extract) == "ISAV-1"), ] =
  colSums(seq_extract[c(which(rownames(seq_extract) == "ISAV-1"),
                        which(rownames(seq_extract) == "ISAV-2"),
                        which(rownames(seq_extract) == "ISAV-3"),
                        which(rownames(seq_extract) == "ISAV-4"),
                        which(rownames(seq_extract) == "ISAV-5"),
                        which(rownames(seq_extract) == "ISAV-6"),
                        which(rownames(seq_extract) == "ISAV-7"),
                        which(rownames(seq_extract) == "ISAV-8")) , ])
# Delete rows 2 - 8 now.
rownames(seq_extract)[which(rownames(seq_extract) == "ISAV-1")] = "ISAV"
seq_extract = seq_extract[-c(which(rownames(seq_extract) == "ISAV-2"),
                             which(rownames(seq_extract) == "ISAV-3"),
                             which(rownames(seq_extract) == "ISAV-4"),
                             which(rownames(seq_extract) == "ISAV-5"),
                             which(rownames(seq_extract) == "ISAV-6"),
                             which(rownames(seq_extract) == "ISAV-7"),
                             which(rownames(seq_extract) == "ISAV-8")), ]

# Need to combine segmented genomes - combine phi6 lines to one.
seq_extract[which(rownames(seq_extract) == "phi6-s"), ] =
  colSums(seq_extract[c(which(rownames(seq_extract) == "phi6-s"),
                        which(rownames(seq_extract) == "phi6-m"),
                        which(rownames(seq_extract) == "phi6-l")) , ])
# Delete additional rows now.
rownames(seq_extract)[which(rownames(seq_extract) == "phi6-s")] = "phi6"
seq_extract = seq_extract[-c(which(rownames(seq_extract) == "phi6-m"),
                             which(rownames(seq_extract) == "phi6-l")), ]

# Need to combine segmented genomes - combine RVSA11 lines to one.
seq_extract[which(rownames(seq_extract) == "RVSA11-1"), ] =
  colSums(seq_extract[c(which(rownames(seq_extract) == "RVSA11-1"),
                        which(rownames(seq_extract) == "RVSA11-2"),
                        which(rownames(seq_extract) == "RVSA11-3"),
                        which(rownames(seq_extract) == "RVSA11-4"),
                        which(rownames(seq_extract) == "RVSA11-5"),
                        which(rownames(seq_extract) == "RVSA11-6"),
                        which(rownames(seq_extract) == "RVSA11-7"),
                        which(rownames(seq_extract) == "RVSA11-8"),
                        which(rownames(seq_extract) == "RVSA11-9"),
                        which(rownames(seq_extract) == "RVSA11-10"),
                        which(rownames(seq_extract) == "RVSA11-11")) , ])
# Delete rows 2 - 8 now.
rownames(seq_extract)[which(rownames(seq_extract) == "RVSA11-1")] = "RVSA11"
seq_extract = seq_extract[-c(which(rownames(seq_extract) == "RVSA11-2"),
                             which(rownames(seq_extract) == "RVSA11-3"),
                             which(rownames(seq_extract) == "RVSA11-4"),
                             which(rownames(seq_extract) == "RVSA11-5"),
                             which(rownames(seq_extract) == "RVSA11-6"),
                             which(rownames(seq_extract) == "RVSA11-7"),
                             which(rownames(seq_extract) == "RVSA11-8"),
                             which(rownames(seq_extract) == "RVSA11-9"),
                             which(rownames(seq_extract) == "RVSA11-10"),
                             which(rownames(seq_extract) == "RVSA11-11")), ]

# Need to combine segmented genomes - combine reo1 lines to one.
seq_extract[which(rownames(seq_extract) == "reo1-l1"), ] =
  colSums(seq_extract[c(which(rownames(seq_extract) == "reo1-l1"),
                        which(rownames(seq_extract) == "reo1-l2"),
                        which(rownames(seq_extract) == "reo1-l3"),
                        which(rownames(seq_extract) == "reo1-m1"),
                        which(rownames(seq_extract) == "reo1-m2"),
                        which(rownames(seq_extract) == "reo1-m3"),
                        which(rownames(seq_extract) == "reo1-s1"),
                        which(rownames(seq_extract) == "reo1-s2"),
                        which(rownames(seq_extract) == "reo1-s3"),
                        which(rownames(seq_extract) == "reo1-s4")) , ])
# Delete additional rows now.
rownames(seq_extract)[which(rownames(seq_extract) == "reo1-l1")] = "reo1"
seq_extract = seq_extract[-c(which(rownames(seq_extract) == "reo1-l2"),
                             which(rownames(seq_extract) == "reo1-l3"),
                             which(rownames(seq_extract) == "reo1-m1"),
                             which(rownames(seq_extract) == "reo1-m2"),
                             which(rownames(seq_extract) == "reo1-m3"),
                             which(rownames(seq_extract) == "reo1-s1"),
                             which(rownames(seq_extract) == "reo1-s2"),
                             which(rownames(seq_extract) == "reo1-s3"),
                             which(rownames(seq_extract) == "reo1-s4")), ]

# Need to combine segmented genomes - combine reo3 lines to one.
seq_extract[which(rownames(seq_extract) == "reo3-l1"), ] =
  colSums(seq_extract[c(which(rownames(seq_extract) == "reo3-l1"),
                        which(rownames(seq_extract) == "reo3-l2"),
                        which(rownames(seq_extract) == "reo3-l3"),
                        which(rownames(seq_extract) == "reo3-m1"),
                        which(rownames(seq_extract) == "reo3-m2"),
                        which(rownames(seq_extract) == "reo3-m3"),
                        which(rownames(seq_extract) == "reo3-s1"),
                        which(rownames(seq_extract) == "reo3-s2"),
                        which(rownames(seq_extract) == "reo3-s3"),
                        which(rownames(seq_extract) == "reo3-s4")) , ])
# Delete additional rows now.
rownames(seq_extract)[which(rownames(seq_extract) == "reo3-l1")] = "reo3"
seq_extract = seq_extract[-c(which(rownames(seq_extract) == "reo3-l2"),
                             which(rownames(seq_extract) == "reo3-l3"),
                             which(rownames(seq_extract) == "reo3-m1"),
                             which(rownames(seq_extract) == "reo3-m2"),
                             which(rownames(seq_extract) == "reo3-m3"),
                             which(rownames(seq_extract) == "reo3-s1"),
                             which(rownames(seq_extract) == "reo3-s2"),
                             which(rownames(seq_extract) == "reo3-s3"),
                             which(rownames(seq_extract) == "reo3-s4")), ]

# DATA OUTPUT-------------------------------------------------------------------

# Outputs data as an R file.
save(seq_extract, file = "virus-seq-attributes.RData")

