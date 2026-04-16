library(dplyr)
library(tidyr)
library(readr)
library(stringr)
#use of python-based scripts
library(reticulate)
use_virtualenv(".../.virtualenvs/r-reticulate/") #set directory location for virtaul environment. Environment must contain the SciPy package to enable LR calculations.
setwd(".../") #set working directory

# case:control max denominator imbalance lookup (see Sup Table 14)
read_csv("data/denom_limits.csv") -> denom_limits

# FUNCTIONS ----

# source R script for OR calculations
source("functions/OR_pvalue_checker.R")

# source python script which gives LR calculation
source_python("functions/lrcalc_functions.updated-thresholds.py")

# source enrichment factor and batch functions
source("functions/output_functions.R")

# filter for PTV variants (VEP consequences)
consequence_filter = c("stop_gained","splice_acceptor_variant","splice_donor_variant","frameshift_variant")

# MTAFs defined by VCEPs
BRCA1_MTAF <- 0.000868 #BRCA1 and BRCA2: ENIGMA v1.2.0 from the appendix (used specific minimal credible AF for BA1)
BRCA2_MTAF <- 0.000906
PALB2_MTAF <- 0.00118 #PALB2 v1.2 specifications, full specification, rounds to 0.001 (0.1%) in CSpec
CHEK2_MTAF <- 0.00625 #CHEK2: uses the same as the autosomal dominaint threshold condition for ATM (breast cancer)
ATM_MTAF <- 0.00625 #ATM v1.3 specifications, uses max credible AF for autosomal dominant disorders (0.0625%)

# Base/standard target odds of association, set at 4 for high-penetrance and 2 for moderate-penetrance by default
BRCA1_base = 4
BRCA2_base = 4
PALB2_base = 4
ATM_base = 2
CHEK2_base = 2

### save exact Enrichment Factors saved for each gene in each dataset here ([[1]] returns EF; [[2]] returns LCI, and [[3]] returns UCI, all outputs saved in Sup. Table 1)
ndrs_brca1 <- EF_calculator(1606,44917,685,92786)[[1]]
ndrs_brca2 <- EF_calculator(2165,44917,1141,92786)[[1]]
ambry_brca1 <- EF_calculator(3642,187642,685,92786)[[1]]
ambry_brca2 <- EF_calculator(5688,187642,1141,92786)[[1]]
ambry_palb2 <- EF_calculator(1646,187642,423,92786)[[1]]
ambry_atm <- EF_calculator(2186,187642,535,92786)[[1]]
ambry_chek2 <- EF_calculator(4539,187642,1243,92786)[[1]]


# ANALYSIS ----
## Odds ratios and LRs ----
### Set column names for calculated values

# combined dataset file
joined_names <- read_csv(joined_names, "combined_ptvmissvar_rare_cleaned.csv")

ukbiobank <- c("p-value_ukbiobank","OR_ukbiobank","Lower_CI_ukbiobank","Upper_CI_ukbiobank","LR_path_ukbiobank","ACMG_points_ukbiobank")
carriers <- c("p-value_carriers","OR_carriers","Lower_CI_carriers","Upper_CI_carriers","LR_path_carriers","ACMG_points_carriers")
bridges <- c("p-value_bridges","OR_bridges","Lower_CI_bridges","Upper_CI_bridges","LR_path_bridges","ACMG_points_bridges")
ndrs <- c("p-value_ndrs","OR_ndrs","Lower_CI_ndrs","Upper_CI_ndrs","LR_path_ndrs","ACMG_points_ndrs")
ambry <- c("p-value_ambry","OR_ambry","Lower_CI_ambry","Upper_CI_ambry","LR_path_ambry","ACMG_points_ambry")

# run LR and OR calculations (function = batch_stats)
joined_names %>%
  batch_stats("hgvs","UKB_case_carrier","UKB_case_total","UKB_control_carrier","UKB_control_total",ukbiobank,
              target_odds_BRCA1=BRCA1_base,target_odds_BRCA2=BRCA2_base,target_odds_PALB2=PALB2_base,target_odds_ATM=ATM_base,target_odds_CHEK2=CHEK2_base) %>%
  batch_stats("hgvs","BRIDGES_case_carrier","BRIDGES_case_total","BRIDGES_control_carrier","BRIDGES_control_total",bridges,
              target_odds_BRCA1=BRCA1_base,target_odds_BRCA2=BRCA2_base,target_odds_PALB2=PALB2_base,target_odds_ATM=ATM_base,target_odds_CHEK2=CHEK2_base) %>%
  batch_stats("hgvs","CARRIERS_case_carrier","CARRIERS_case_total","CARRIERS_control_carrier","CARRIERS_control_total",carriers,
              target_odds_BRCA1=BRCA1_base,target_odds_BRCA2=BRCA2_base,target_odds_PALB2=PALB2_base,target_odds_ATM=ATM_base,target_odds_CHEK2=CHEK2_base) %>%
  batch_stats("hgvs","NDRS_case_carrier","NDRS_case_total","NDRS_control_carrier","NDRS_control_total",ndrs,
              target_odds_BRCA1=BRCA1_base*ndrs_brca1,target_odds_BRCA2=BRCA2_base*ndrs_brca2) %>%
  batch_stats("hgvs","AMBRY_case_carrier","AMBRY_case_total","AMBRY_control_carrier","AMBRY_control_total",ambry,
              target_odds_BRCA1=BRCA1_base*ambry_brca1,target_odds_BRCA2=BRCA2_base*ambry_brca2,target_odds_PALB2=PALB2_base*ambry_palb2,target_odds_ATM=ATM_base*ambry_atm,target_odds_CHEK2=CHEK2_base*ambry_chek2) -> stats_target_4

## Reformatting
colnames(stats_target_4) -> colnames
colnames[29] <- "Target_Odds_ukbiobank"
colnames[36] <- "Target_Odds_bridges"
colnames[43] <- "Target_Odds_carriers"
colnames[50] <- "Target_Odds_ndrs"
colnames[57] <- "Target_Odds_ambry"
colnames(stats_target_4) <- colnames

stats_target_4 %>%
  relocate(ends_with("ukbiobank"), .after=UKB_control_total) %>%
  relocate(ends_with("bridges"), .after=BRIDGES_control_total) %>%
  relocate(ends_with("carriers"), .after=CARRIERS_control_total) %>%
  relocate(ends_with("ndrs"), .after=NDRS_control_total) %>%
  relocate(ends_with("ambry"), .after=AMBRY_control_total) %>%
  mutate(LR_path_ukbiobank=as.numeric(LR_path_ukbiobank),ACMG_points_ukbiobank=as.numeric(ACMG_points_ukbiobank),
         LR_path_carriers=as.numeric(LR_path_carriers),ACMG_points_carriers=as.numeric(ACMG_points_carriers),
         LR_path_bridges=as.numeric(LR_path_bridges),ACMG_points_bridges=as.numeric(ACMG_points_bridges),
         LR_path_ndrs=as.numeric(LR_path_ndrs),ACMG_points_ndrs=as.numeric(ACMG_points_ndrs),
         LR_path_ambry=as.numeric(LR_path_ambry),ACMG_points_ambry=as.numeric(ACMG_points_ambry)) -> stats_target_4

write_csv(stats_target_4, "combined_ptvmissvar_LR.csv")

## Parameter filtering ----
### add discrepancy thresholds (see Supplementary Methods)
stats_target_4 %>%
  mutate(across(starts_with("Target_Odds_"), ~ round(.x,2))) %>%
  left_join(denom_limits, by=c("UKB_control_carrier"="obs_control","UKB_case_carrier"="obs_case","Target_Odds_ukbiobank"="target_odds")) %>% rename(control_case_ratio_ukbiobank=control_case_ratio, case_control_ratio_ukbiobank=case_control_ratio) %>%
  left_join(denom_limits, by=c("BRIDGES_control_carrier"="obs_control","BRIDGES_case_carrier"="obs_case","Target_Odds_bridges"="target_odds")) %>% rename(control_case_ratio_bridges=control_case_ratio, case_control_ratio_bridges=case_control_ratio) %>%
  left_join(denom_limits, by=c("CARRIERS_control_carrier"="obs_control","CARRIERS_case_carrier"="obs_case","Target_Odds_carriers"="target_odds")) %>% rename(control_case_ratio_carriers=control_case_ratio, case_control_ratio_carriers=case_control_ratio) %>%
  left_join(denom_limits, by=c("NDRS_control_carrier"="obs_control","NDRS_case_carrier"="obs_case","Target_Odds_ndrs"="target_odds")) %>% rename(control_case_ratio_ndrs=control_case_ratio, case_control_ratio_ndrs=case_control_ratio) %>%
  left_join(denom_limits, by=c("AMBRY_control_carrier"="obs_control","AMBRY_case_carrier"="obs_case","Target_Odds_ambry"="target_odds")) %>% rename(control_case_ratio_ambry=control_case_ratio, case_control_ratio_ambry=case_control_ratio) -> stats_target_4b

## apply filtering for each dataset (see Supplementary Table 14)
stats_target_4b %>%
  parameter_filtering("ukbiobank", "UKB_case_carrier", "UKB_control_carrier", "UKB_case_total", "UKB_control_total", "ACMG_points_ukbiobank", "Upper_CI_ukbiobank", "Lower_CI_ukbiobank", "control_case_ratio_ukbiobank", "case_control_ratio_ukbiobank") %>%
  parameter_filtering("bridges", "BRIDGES_case_carrier", "BRIDGES_control_carrier", "BRIDGES_case_total", "BRIDGES_control_total", "ACMG_points_bridges", "Upper_CI_bridges", "Lower_CI_bridges", "control_case_ratio_bridges", "case_control_ratio_bridges") %>%
  parameter_filtering("carriers", "CARRIERS_case_carrier", "CARRIERS_control_carrier", "CARRIERS_case_total", "CARRIERS_control_total", "ACMG_points_carriers", "Upper_CI_carriers", "Lower_CI_carriers", "control_case_ratio_carriers", "case_control_ratio_carriers") %>%
  parameter_filtering("ndrs", "NDRS_case_carrier", "NDRS_control_carrier", "NDRS_case_total", "NDRS_control_total", "ACMG_points_ndrs", "Upper_CI_ndrs", "Lower_CI_ndrs", "control_case_ratio_ndrs","case_control_ratio_ndrs") %>%
  parameter_filtering("ambry", "AMBRY_case_carrier", "AMBRY_control_carrier", "AMBRY_case_total", "AMBRY_control_total", "ACMG_points_ambry", "Upper_CI_ambry", "Lower_CI_ambry", "control_case_ratio_ambry", "case_control_ratio_ambry") -> case_counts

## set 'inclusion' flag to 'N' for benign NDRS data, as we aren't including this in the combined LR.
case_counts %>% mutate(include_ndrs = if_else(ACMG_points_ndrs<0, "N", include_ndrs)) -> case_counts

## Reformatting
case_counts %>%
  relocate(ends_with("ukbiobank"), .after=UKB_control_total) %>%
  relocate(ends_with("bridges"), .after=BRIDGES_control_total) %>%
  relocate(ends_with("carriers"), .after=CARRIERS_control_total) %>%
  relocate(ends_with("ndrs"), .after=NDRS_control_total) %>%
  relocate(ends_with("ambry"), .after=AMBRY_control_total) -> case_counts

# single obs inclusion rule: include single obs if there is any data in another dataset (and that other dataset is NOT benign evidence from NDRS)
case_counts %>%
  mutate(include_ukbiobank_combo = case_when(include_ukbiobank=="Y" ~ "Y",
                                             # non-zero observation in another dataset and both the obs in this dataset and the other dataset are not discrepant
                                             single_obs_ukbiobank==1 & discrepant_evidence_ukbiobank==0 & ((zero_obs_bridges==0 & discrepant_evidence_bridges==0) | 
                                                                                                             (zero_obs_carriers==0 & discrepant_evidence_carriers==0) | 
                                                                                                             (zero_obs_ndrs==0 & discrepant_evidence_ndrs==0 & ACMG_points_ndrs>0) | 
                                                                                                             (zero_obs_ambry==0 & discrepant_evidence_ambry==0)) ~ "Y",
                                             .default="N")) %>%
  mutate(include_bridges_combo = case_when(include_bridges=="Y" ~ "Y",
                                           # non-zero observation in another dataset
                                           single_obs_bridges==1 & discrepant_evidence_bridges==0 & ((zero_obs_ukbiobank==0 & discrepant_evidence_ukbiobank==0) | 
                                                                                                       (zero_obs_carriers==0 & discrepant_evidence_carriers==0) | 
                                                                                                       (zero_obs_ndrs==0 & discrepant_evidence_ndrs==0 & ACMG_points_ndrs>0) | 
                                                                                                       (zero_obs_ambry==0 & discrepant_evidence_ambry==0)) ~ "Y",
                                           .default="N")) %>%
  mutate(include_carriers_combo = case_when(include_carriers=="Y" ~ "Y",
                                            # non-zero observation in another dataset
                                            single_obs_carriers==1 & discrepant_evidence_carriers==0 & ((zero_obs_ukbiobank==0 & discrepant_evidence_ukbiobank==0) | 
                                                                                                          (zero_obs_bridges==0 & discrepant_evidence_bridges==0) | 
                                                                                                          (zero_obs_ndrs==0 & discrepant_evidence_ndrs==0 & ACMG_points_ndrs>0) | 
                                                                                                          (zero_obs_ambry==0 & discrepant_evidence_ambry==0)) ~ "Y",
                                            .default="N")) %>%
  mutate(include_ndrs_combo = case_when(include_ndrs=="Y" ~ "Y",
                                        # non-zero observation in another dataset, and (for NDRS data only) the evidence is not benign as we are disregarding this 
                                        single_obs_ndrs==1 & discrepant_evidence_ndrs==0 & ACMG_points_ndrs>0 & ((zero_obs_ukbiobank==0 & discrepant_evidence_ukbiobank==0) | 
                                                                                                                   (zero_obs_bridges==0 & discrepant_evidence_bridges==0) | 
                                                                                                                   (zero_obs_carriers==0 & discrepant_evidence_carriers==0) | 
                                                                                                                   (zero_obs_ambry==0 & discrepant_evidence_ambry==0)) ~ "Y",
                                        .default="N")) %>%
  mutate(include_ambry_combo = case_when(include_ambry=="Y" ~ "Y",
                                         # non-zero observation in another dataset
                                         single_obs_ambry==1 & discrepant_evidence_ambry==0 & ((zero_obs_ukbiobank==0 & discrepant_evidence_ukbiobank==0) | 
                                                                                                 (zero_obs_bridges==0 & discrepant_evidence_bridges==0) | 
                                                                                                 (zero_obs_carriers==0 & discrepant_evidence_carriers==0) | 
                                                                                                 (zero_obs_ndrs==0 & discrepant_evidence_ndrs==0 & ACMG_points_ndrs>0)) ~ "Y",
                                         .default="N")) %>%
  
  # combines all LR marked as 'Y' for inclusion 
  # should include variants which are single observations in multiple datasets as long as the single observations all give evidence in the same direction. 
  # also adding a check here: do NOT include the NDRS data in the combined LR if the evidence is towards benignity due to under-reporting of benign variants from NDRS
  mutate(LR_combo = if_else((include_ukbiobank_combo == "Y") & (discrepant_evidence_ukbiobank==0), LR_path_ukbiobank, 1) * 
           if_else((include_bridges_combo == "Y") & (discrepant_evidence_bridges==0), LR_path_bridges, 1) *
           if_else((include_carriers_combo == "Y") & (discrepant_evidence_carriers==0), LR_path_carriers, 1) *
           if_else((include_ndrs_combo == "Y") & (discrepant_evidence_ndrs==0) & (ACMG_points_ndrs>0) & (gene %in% c("BRCA1","BRCA2")), LR_path_ndrs, 1) *
           if_else((include_ambry_combo == "Y") & (discrepant_evidence_ambry==0) * !is.na(LR_path_ambry), LR_path_ambry, 1)) %>%
  # where LR value is too low to be represented in R (<5e-324), set combined LR to the smallest number possible as a stand-in.
  mutate(LR_combo = if_else(LR_combo==0,5e-324,LR_combo)) %>%
  mutate(constituent_data = paste0(if_else(include_ukbiobank_combo == "Y", "A", ""),
                                   if_else(include_bridges_combo == "Y", "B", ""),
                                   if_else(include_carriers_combo == "Y", "C", ""),
                                   if_else(include_ndrs_combo == "Y" & gene %in% c("BRCA1","BRCA2"), "D", ""),
                                   if_else(include_ambry_combo == "Y", "E", ""))) %>%
  rowwise() %>%
  # calculate PS4-LLR (log LR to base 2.08)
  mutate(ACMG_combo=log(LR_combo,2.08)) %>%
  # remove variants which are not present in any dataset after filtering
  mutate(UKB_zero_from_part = if_else(zero_obs_ukbiobank==1 & zero_obs_bridges==1 & zero_obs_carriers==1 & zero_obs_ndrs==1 & zero_obs_ambry==1, 1, 0)) %>%
  filter(UKB_zero_from_part==0) %>%
  mutate(constituent_data = gsub("^$|^ $",NA,constituent_data)) -> case_counts_with_combo

# assign categorical evidence strengths to PS4-LLR
plot_prep(case_counts_with_combo, "UK Biobank", "ACMG_points_ukbiobank")-> case_counts_with_combo
plot_prep(case_counts_with_combo, "BRIDGES", "ACMG_points_bridges")-> case_counts_with_combo
plot_prep(case_counts_with_combo, "CARRIERS", "ACMG_points_carriers")-> case_counts_with_combo
plot_prep(case_counts_with_combo, "NDRS", "ACMG_points_ndrs")-> case_counts_with_combo
plot_prep(case_counts_with_combo, "AMBRY", "ACMG_points_ambry")-> case_counts_with_combo
plot_prep(case_counts_with_combo, "All Datasets", "ACMG_combo")-> case_counts_with_combo

# output final dataset
write_csv(case_counts_with_combo, "combined_ptvmissvar_calculated_filtered.csv")
