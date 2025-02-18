# SCN5A_SP_BrS_cohort
Analysis of electrophysiology data for a high-throughput SCN5A Automated Patch Clamp assay of BrS cohort variants

[NOT UPDATED]
****MO UPDATE THIS *****

A major challenge of precision medicine is interpreting variant effect. The lethal condition Brugada Syndrome is associated with rare, loss-of-function variants in the gene SCN5A, which encodes for the major cardiac sodium channel. Previously, we have published a validated high-throughput Automated Patch Clamp assay using (PMID: 32533946). Now, we present a calibrated assay based on control B and P/LP variants, and compare these results with those independently obtained by collaborators at the Victor Chang Cardiovascular Research Institute. 

Our dataset is represented in 3 primary parts.
  1) SP_analysis_functions.R - this R script contains all functions used to analyze each dataset, and are called by each individual experiment and the combined analysis
  2) SP_analysis_combined.R - this R script combines data for each EP parameter across all relevant APC experiments, performs outlier removal, and merges each parameter
  3) SP_exp folder - these folders contain the raw data from the APC platform, an analysis script to process the data, and the outputted data for use in the combined R script

Together, these files provide a framework for rapid analysis of these datasets. Please email andrew.m.glazer@vumc.org with questions. These data form the basis for the paper "Multisite Validation of a Functional Assay to Adjudicate SCN5A Brugada Syndrome–Associated Variants" by Ma and O'Neill et al 2024, published in Circulation: Genomics in Precision Medicine (PMID 38953211). Additional raw SP data will be made available upon reasonable request. 
