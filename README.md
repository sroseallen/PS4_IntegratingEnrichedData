# PS4_IntegratingEnrichedData

This repository provides the scripts used to produce all outputs and for the manuscript "Integrating enriched case data from national laboratory testing with population-based case-control analyses: a novel statistical likelihood-ratio methodology for PS4 applied to 325,345 breast cancer cases and 671,006 controls".

The following scripts are included:
- `lrcalc_functions.updated-thresholds.py`: PS4-LR-Calc functions, calculates the Likelihood and Log Likelihood Ratio of the observed data being greater than the specified target odds of association (please see Rowlands et al., 2025 for further details on the PS4-LR-Calc, DOI: 10.1136/jmg-2024-110034)
- `OR_pvalue_checker.R`: functions to calculate the observed Odds Ratio plus upper and lower 95% confidence intervals, and the Fishers Exact p-value.
- `Integrating_Enriched_Data_Pipleine.R`: Data cleaning and preparation pipeline for the datasets used in the main manuscript: BRIDGES, CARRIERS, UK Biobank, UK Laboratory (NDRS), Ambry, and gnomAD v4.1. Raw data files are provided for each either as VCF or CSV files - data from gnomAD v4.1 and BRIDGES can be found at the respective links (see below), other data files are only available on application and approval from the project lead or (in the case of UK Biobank) through the project application portal. As such, due to data protection and legal requirement, raw data files cannot be publically provided for this analysis.

## Data sources
- gnomAD v4.1: https://gnomad.broadinstitute.org/data
- BRIDGES: https://www.ccge.medschl.cam.ac.uk/breast-cancer-association-consortium-bcac/data-data-access/summary-results/bridges-summary-results

