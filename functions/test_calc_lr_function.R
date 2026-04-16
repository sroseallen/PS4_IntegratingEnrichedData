#use of python-based scripts
library(reticulate)
use_virtualenv(".../.virtualenvs/r-reticulate/") #set directory location for virtual environment. Environment must contain the SciPy package to enable LR calculations.
setwd(".../") #set working directory

# source python script which gives LR calculation (calc_lr)
source_python("functions/lrcalc_functions.updated-thresholds.py")

### save exact Enrichment Factors saved for each gene in each dataset here ([[1]] returns EF; [[2]] returns LCI, and [[3]] returns UCI, all outputs saved in Sup. Table 1)
ndrs_brca1 <- EF_calculator(1606,44917,685,92786)[[1]]
ndrs_brca2 <- EF_calculator(2165,44917,1141,92786)[[1]]
ambry_brca1 <- EF_calculator(3642,187642,685,92786)[[1]]
ambry_brca2 <- EF_calculator(5688,187642,1141,92786)[[1]]
ambry_palb2 <- EF_calculator(1646,187642,423,92786)[[1]]
ambry_atm <- EF_calculator(2186,187642,535,92786)[[1]]
ambry_chek2 <- EF_calculator(4539,187642,1243,92786)[[1]]

# test calc_lr for different scenarios
test_calc_lr <- function(case_var_count, total_case_count, control_var_count, total_control_count, PATHOGENIC_OR_THRESHOLD) {
  calc_lr(case_var_count = case_var_count, 
          total_case_count = total_case_count, 
          control_var_count = control_var_count, 
          total_control_count = total_control_count, 
          PATHOGENIC_OR_THRESHOLD = PATHOGENIC_OR_THRESHOLD)[[1]] -> LR
  calc_lr(case_var_count = case_var_count, 
          total_case_count = total_case_count, 
          control_var_count = control_var_count, 
          total_control_count = total_control_count, 
          PATHOGENIC_OR_THRESHOLD = PATHOGENIC_OR_THRESHOLD)[[2]] -> PS4_LLR
  new_lr <- data.frame(LR = LR, PS4_LLR = PS4_LLR)
  return(new_lr)
}

test_calc_lr(6,42062,2,44035,4) #variant counts from BRCA1 c.427G>T in the BRIDGES dataset, see Supplementary Table 2
test_calc_lr(156,187642,14,173824,(4 * ambry_brca2)) #variant counts from BRCA2 c.6275_6276del in the Ambry dataset, see Supplementary Table 2
