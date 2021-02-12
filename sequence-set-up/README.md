# uv-virus-inactivation-prediction

# 'sequence-set-up' folder

# Files in this folder are used to obtain genome sequence information (e.g., number of uracils in the genome, genome length) that was then used in modeling work. The genome sequence information is extracted directly from the .txt files containing virus genome sequences.

# This file must run prior to running code in the other folders to ensure that genome sequence attributes for all viruses are included in virus data sets.

# Files

# 00-uv-inact-sequence-attributes.R
Defines virus genome sequence information. NOTE: The exact txt file name must be included in the DATA INPUT section of this code to ensure genome sequence information extraction for the virus of interest.