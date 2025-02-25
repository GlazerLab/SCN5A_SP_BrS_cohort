# SCN5A_SP_BrS_cohort
Analysis of SCN5A electrophysiology data and penetrance calculations for "Cohort-scale automated patch clamp data improves variant classification and penetrance stratification for SCN5A-Brugada Syndrome"

Brugada Syndrome (BrS) is an inherited disorder linked to sudden cardiac death, with about 20% of cases involving SCN5A gene variants. In a study analyzing 252 SCN5A variants from 3,335 BrS patients using a high-throughput automated patch clamp assay, 146 variants were found to be functionally abnormal, with severe loss-of-function in 100 cases. Functional evidence reclassified 110 of 225 variants of uncertain significance (VUS), identifying 104 as likely pathogenic and 6 as likely benign. Loss-of-function variants were mainly in the transmembrane domains, correlating with higher BrS risk and penetrance. This dataset refines the understanding of SCN5A variants, aiding in BrS diagnosis and management.

Our dataset is represented in 4 primary parts.
  1) SP_analysis_functions.R - this R script contains all functions used to analyze each dataset, and are called by each individual experiment and the combined analysis
  2) SP_analysis_combined.R - this R script combines data for each EP parameter across all relevant APC experiments, performs outlier removal, and merges each parameter
  3) SP_exp folder - these folders contain the raw data from the APC platform, an analysis script to process the data, and the outputted data for use in the combined R script
  4) Penetrance Math for_Github.R - this R script contains our analysis of BrS penetrance (method adopted from McGurk et al PMID 37652022)

Together, these files provide a framework for rapid analysis of these datasets. Please email andrew.m.glazer@vumc.org with questions. 

Additional raw SP data will be made available upon reasonable request. 
