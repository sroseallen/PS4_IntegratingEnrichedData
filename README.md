# PS4_IntegratingEnrichedData

This repository provides the scripts used to produce all outputs and for the manuscript "Integrating enriched case data from national laboratory testing with population-based case-control analyses: a novel statistical likelihood-ratio methodology for PS4 applied to 325,345 breast cancer cases and 671,006 controls".

The following scripts are included:
- `lrcalc_functions.updated-thresholds.py`: PS4-LR-Calc functions, calculates the Likelihood and Log Likelihood Ratio of the observed data being greater than the specified target odds of association (please see Rowlands et al., 2025 for further details on the PS4-LR-Calc, DOI: 10.1136/jmg-2024-110034)
- `OR_pvalue_checker.R`: functions to calculate the observed Odds Ratio plus upper and lower 95% confidence intervals, and the Fishers Exact p-value.
- `Integrating_Enriched_Data_Cleaning.R`*: Data cleaning script used to reformat and prepare the datasets for combination and calculation of PS4-LLR (each file contains at minimum the variant hgvs nomenclature, per-variant case and control carrier counts, and total case and control counts). Is run first to produce the combined output file 'combined_ptvmissvar_rare_cleaned' of all rare PTV and missense variants across the input datasets.
- `Integrating_Enriched_Data_Pipeline.R`: Data analysis pipeline for the datasets used in the main manuscript: BRIDGES, CARRIERS, UK Biobank, UK Laboratory (NDRS), Ambry, and gnomAD v4.1. Takes the output file from `Integrating_Enriched_Data_Cleaning.R` and, using the lrcalc and OR_pvalue_checker functions scripts, calculates OR, LR, and PS4-LLR for each variant on both a per-dataset level and on a combined level incorporating data from across all datasets. Outputs one file of analysed PTV and missense variants, 'combined_ptvmissvar_calculated_filtered.csv'. 

*Please note: raw data files are provided for each either as VCF or CSV files - data from gnomAD v4.1 and BRIDGES can be found at the respective links (see below), other data files are only available on application and approval from the project lead or (in the case of UK Biobank) through the project application portal. As such, the raw data files and some variant counts cannot be publically shared. 

## Testing the calc_lr function
The `calc_lr` function provides the LR and PS4-LLR for variant, taking the following parameters:
- case_var_count = the number of variant carriers in cases
- total_case_count = the total number of individuals tested in cases
- control_var_count = the number of variant carriers in controls
- total_control_count = the total number of individuals tested in controls
- PATHOGENIC_OR_THRESHOLD = the target odds of association, default set at 5
- BENIGN_OR_THRESHOLD = the target odds of non-association, default set at 1

To test the calculator, we have provided a test file `test_calc_lr_function` which reads in the `lrcalc_functions.updated-thresholds.py` file and runs `test_calc_lr` using a pre-set number of variant carriers. The number of variant carriers can be adjusted to test different variant scenarios. 
Please see Rowlands et al., 2024 for further information on the PS4-LR-Calculator function, and the supplementary methods for information on the adjustments made for this analysis. 

## Data sources
- gnomAD v4.1: https://gnomad.broadinstitute.org/data
- BRIDGES: https://www.ccge.medschl.cam.ac.uk/breast-cancer-association-consortium-bcac/data-data-access/summary-results/bridges-summary-results

