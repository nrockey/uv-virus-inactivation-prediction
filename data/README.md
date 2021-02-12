# uv-virus-inactivation-prediction

# 'data' folder

# This folder contains all the raw data inputs (virus sequence info, systematic review data, and additional virus attributes for viruses used in model training/validation and viruses used in predictions) needed to reproduce the findings described in our study: Rockey, Nicole C., Henderson, James B., Chin, Kaitlyn, Raskin, Lutgarde, Wigginton, Krista R., "Predictive Modeling of Virus Inactivation by UV". Environmental Science & Technology, 2021.

# Files

# .txt files (numerous)
These files are in fasta file format and contain the genome sequence (including only letters A, G, T, and C) for the specified virus. There should be no spaces or carriage returns between sequence letters, otherwise, sequence info from this file will not be accurately enumerated. The first line of the file should be in fasta format (i.e., '>NAME_OF_VIRUS_HERE').

# adnl-virus-vars.csv
This file contains a list of viruses (one virus per row) and additional variables needed for each of those viruses. The viruses included in this data set were viruses returned in the systematic review for inactivation rate constants. The set of additional variables needed for each virus includes three variables ('class', 'repair', and 'host'). The exact designations given for these variables are categorical and are defined as follows:
    - 'class':      'dsdna' = virus with a double-stranded DNA genome.
                    'ssdna' = virus with a single-stranded DNA genome.
                    'dsrna' = virus with a double-stranded RNA genome.
                    'plusssrna' = virus with a postive sense, single-stranded RNA genome.
                    'negssrna' = virus with a negative sense, single-stranded RNA genome.
                    
    - 'repair':
                    0 = host cell mediated (this is used for most dsDNA viruses, which are known to have host cell mediated repair, or any dsDNA virus where additional details of repair are not known);
                    1 = virus-gene controlled using one repair system (this is used for viruses that are known to encode a gene for controlling genome repair)
                    2 = no repair (this is used for any viruses other than dsDNA viruses)
                    3 = virus-gene controlled using multiple repair systems (this is used for viruses that are known to encode multiple genes for controlling genome repair)
    - 'host':
                    0 = prokaryotic host (this is used for any virus that infects bacteria)
                    1 = eukaryotic host with reduced repair abilities (this is used for any mammalian dsDNA virus where the host is known to have reduced repair abilities. E.g., a host with xeroderma pigmentosum)
                    2 = eukaryotic host with wild type repair abilities (this is used for any mammalian dsDNA virus where the host's repair abilities are not known; This is also used this for any mammalian virus with ssDNA, ssRNA, or dsRNA)

# adnl-predict-virus-vars.csv
This file contains a list of viruses (one virus per row) and additional variables needed for each of those viruses. The viruses included in this data set were viruses returned in the systematic review for inactivation rate constants as well as viruses of interest for which rate constants were to be predicted, because rate constants were not available from the review. The set of additional variables needed for each virus includes three variables ('class', 'repair', and 'host'). The exact designations given for these variables are categorical and are defined as follows:
    - 'class':      'dsdna' = virus with a double-stranded DNA genome.
                    'ssdna' = virus with a single-stranded DNA genome.
                    'dsrna' = virus with a double-stranded RNA genome.
                    'plusssrna' = virus with a postive sense, single-stranded RNA genome.
                    'negssrna' = virus with a negative sense, single-stranded RNA genome.
                    
    - 'repair':
                    0 = host cell mediated (this is used for most dsDNA viruses, which are known to have host cell mediated repair, or any dsDNA virus where additional details of repair are not known);
                    1 = virus-gene controlled using one repair system (this is used for viruses that are known to encode a gene for controlling genome repair)
                    2 = no repair (this is used for any viruses other than dsDNA viruses)
                    3 = virus-gene controlled using multiple repair systems (this is used for viruses that are known to encode multiple genes for controlling genome repair)
    - 'host':
                    0 = prokaryotic host (this is used for any virus that infects bacteria)
                    1 = eukaryotic host with reduced repair abilities (this is used for any mammalian dsDNA virus where the host is known to have reduced repair abilities. E.g., a host with xeroderma pigmentosum)
                    2 = eukaryotic host with wild type repair abilities (this is used for any mammalian dsDNA virus where the host's repair abilities are not known; This is also used this for any mammalian virus with ssDNA, ssRNA, or dsRNA)

# sys-review-data.csv
This file contains the results of the systematic review conducted to collect UV inactivation rate constants from the literature. Each row contains the rate constant obtained from a particular study for a specific virus. Information in each row includes study ID number, virus identifier, rate constant, standard error (if it was available), and virus class (e.g., dsdna, plusssrna, etc.).

# RData files (numerous)
These files include R data (e.g., matrices, variables, etc.) that was produced from the R code in the 'sequence-set-up', 'development', and 'prediction' folders. These files are loaded into R code as needed.
