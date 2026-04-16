library(dplyr)
library(tidyr)
library(readr)
library(stringr)
#statistics
library(rstatix)
#plotting
library(ggplot2)
library(ggsankey)
library(RColorBrewer)
library(ggpubr)
#use of python-based scripts
library(reticulate)
use_virtualenv("C:/Users/sallen/Documents/.virtualenvs/r-reticulate/") #contains SciPy package needed for LR calculations

setwd("C:/Users/sallen/OneDrive - The Institute of Cancer Research/Documents/CanGene-CanVar Work Package 2/CGCV_WP2/ACTIVE_VARIANT_ANALYSES_PAPERS/PS4/PS4_Breast_SA")

# main colour palette
colour_palette = c("P_VSTR" = "#650b04",
                   "P_STR" = "#B2182B",
                   "P_MOD" = "#D6604D",
                   "P_SUP" = "#F4A582",
                   "None" = "#929292",
                   "B_SUP" = "#92C5DE",
                   "B_MOD" = "#4393C3",
                   "B_STR" = "#2166AC",
                   "B_VSTR" = "#130857")


# FUNCTIONS ----

# source R script for OR calculations
source("OR_pvalue_checker.R")

# source python script which gives LR calculation
source_python("lrcalc_functions.updated-thresholds.apr26.fixed.py")

# MTAFs defined by VCEPs
BRCA1_MTAF <- 0.000868 #BRCA1 and BRCA2: ENIGMA v1.2.0 from the appendix (used specific minimal credible AF for BA1)
BRCA2_MTAF <- 0.000906
PALB2_MTAF <- 0.00118 #PALB2 v1.2 specifications, full specification, rounds to 0.001 (0.1%) in CSpec
CHEK2_MTAF <- 0.00625 #CHEK2: uses the same as the autosomal dominaint threshold condition for ATM (breast cancer)
ATM_MTAF <- 0.00625 #ATM v1.3 specifications, uses max credible AF for autosomal dominant disorders (0.0625%)

# Base/standard target odds of association
BRCA1_base = 8.73
BRCA2_base = 5.68
PALB2_base = 4.30
ATM_base = 2.16
CHEK2_base = 2.40

# filter for PTV variants (VEP consequences)
consequence_filter = c("stop_gained","splice_acceptor_variant","splice_donor_variant","frameshift_variant")

### Enrichment Factor calculations (+ confidence intervals)
EF_calculator <- function(PTV_carriers_lab, tot_tested_lab, PTV_carriers_unselected, tot_tested_unselected) {
  EF <- (PTV_carriers_lab/tot_tested_lab) * (tot_tested_unselected/PTV_carriers_unselected)
  LCI <- exp(log(EF)-(1.96*sqrt(1/PTV_carriers_lab - 1/tot_tested_lab + 1/PTV_carriers_unselected - 1/tot_tested_unselected)))
  UCI <- exp(log(EF)+(1.96*sqrt(1/PTV_carriers_lab - 1/tot_tested_lab + 1/PTV_carriers_unselected - 1/tot_tested_unselected)))
  return(list(EF,LCI,UCI))
}

### exact Enrichment Factors saved for each gene in each dataset here ([[1]] returns EF; [[2]] returns LCI, and [[3]] returns UCI, outputs saved in Sup. Table 1)
ndrs_brca1 <- EF_calculator(1606,44917,685,92786)[[1]]
ndrs_brca2 <- EF_calculator(2165,44917,1141,92786)[[1]]
ambry_brca1 <- EF_calculator(3642,187642,685,92786)[[1]]
ambry_brca2 <- EF_calculator(5688,187642,1141,92786)[[1]]
ambry_palb2 <- EF_calculator(1646,187642,423,92786)[[1]]
ambry_atm <- EF_calculator(2186,187642,535,92786)[[1]]
ambry_chek2 <- EF_calculator(4539,187642,1243,92786)[[1]]

# To flag variants exceeding the BA1 MTAC threshold
## using modified Whiffin/Ware calculator find_af_filter function (since gnomAD v4.1 filtering uses the FAF for white ethnicity, this method is more comparable and allows use of the lower 95% CI of AF estimate, as calculated for gnomAD)
## Whiffin N, Minikel E, Walsh R, et al. Using high-resolution variant frequencies to empower clinical genome interpretation. Genet Med. 2017;19(10):1151-1158. doi:10.1038/gim.2017.26
find_max_ac = function(af,an,ci=.95) {
  quantile_limit = ci # ci for one-sided, 1-(1-ci)/2 for two-sided
  max_ac = qpois(quantile_limit,an*af)
  return (max_ac)
}

find_af_filter = function(ac, an, ci=.95, lower=(.1/(2*60706)), upper=2, tol=1e-7, precision=1e-6) { 
  if (is.na(ac) | is.na(an) | ac == 0 | an == 0 | ac == 1) {
    return (0.0)
  } else {
    quantile_limit = ci # ci for one-sided, 1-(1-ci)/2 for two-sided
    uniroot_result = uniroot(f = function(af,ac,an) { return (ac - 1 - qpois(p=quantile_limit,lambda=an*af)) },lower=lower,upper=upper,ac=ac,an=an,tol=tol)
    max_af = round(uniroot_result$root,-log10(precision)) # round to nearest millionth
    while(find_max_ac(af=max_af,an=an) < ac) {
      max_af = max_af + precision # increase by millionths until you find the upper bound - the highest AF for which 95%CI AC is still less than observed AC
    }
    max_af = max_af - precision # back off one unit from the AF that violated the 95%CI AC < obs AC condition
    return (max_af)
  }
}

# Split UKB controls between UKB and NDRS cases to minimise case:control ratio for each case data series
# set random seed to standard value to ensure results can be replicated
set.seed(123)

sample_controls <- function(carriers, prob=c(0.2914637,0.7085363)) {
  assign <- sample(c("A", "B"), size=carriers, replace=TRUE, prob=prob) #A = ukb, B = ndrs
  A <- sum(str_count(assign, "A"))
  B <- sum(str_count(assign, "B"))
  return (A)
}

# batch generation of OR and LR stats
batch_stats <- function(df, hgvs, col1, col2, col3, col4, new_colnames=c("p-value","OR","Lower_CI","Upper_CI","LR_path","ACMG_points","Target_Odds"), target_odds_BRCA1=4, target_odds_BRCA2=4, target_odds_PALB2=4, target_odds_ATM=2, target_odds_CHEK2=2, ci=0) {
  
  newdata_df <- data.frame(matrix(ncol=7,nrow=0))
  colnames(newdata_df) <- new_colnames
  
  for (i in 1:nrow(df)) {
    df[i,col1] -> a
    df[i,col2] -> b
    df[i,col3] -> c
    df[i,col4] -> d
    df[i, "gene"] -> gene_name
    
    a <- as.numeric(a)
    b <- as.numeric(b)
    c <- as.numeric(c)
    d <- as.numeric(d)
    
    # decide if haldane correction is necessary for OR calculation
    if (a != 0 & c == 0) {
      new_OR <- all_ORstat(a,(b-a),c,(d-c),haldane=TRUE,critical_value=1.96)
    } else {
      new_OR <- all_ORstat(a,(b-a),c,(d-c),haldane=FALSE,critical_value=1.96)
    }
    
    ## variants which have a count that prevents proper calculation of LR with PS4-LR-Calc + will need to be entered as "NC - not calculable" in the new_LR table, see Sup Methods
    brca1_highcount <- c("c.1487G>A")
    brca2_highcount <- c("c.7544C>T")
    palb2_highcount <- c()
    atm_highcount <- c("c.6067G>A","c.998C>T","c.5558A>T","c.5071A>C","c.1229T>C","c.3925G>A")
    chek2_highcount <- c("c.1100del","c.470T>C")
    
    # calculate LR and ACMG points
    if (gene_name == "BRCA1") {
      if (df[i,hgvs] %in% brca1_highcount & "LR_path_ambry" %in% new_colnames) {
        new_LR <- data.frame(LR_path = c("NC - not calculable"), ACMG_points = c("NC - not calculable"))
      } else if (a==0 & c==0) {
        new_LR <- data.frame(LR_path = c("NA"), ACMG_points = c("NA"))
      } else {
        new_LR <- as.data.frame(calc_lr(a,b,c,d, PATHOGENIC_OR_THRESHOLD=target_odds_BRCA1, CONFIDENCE_LEVEL=ci))
      }
      new_LR <- mutate(new_LR, Target_Odds := target_odds_BRCA1)
    } else if (gene_name == "BRCA2") {
      if (df[i,hgvs] %in% brca2_highcount & "LR_path_ambry" %in% new_colnames) {
        new_LR <- data.frame(LR_path = c("NC - not calculable"), ACMG_points = c("NC - not calculable"))
      } else if (a==0 & c==0) {
        new_LR <- data.frame(LR_path = c("NA"), ACMG_points = c("NA"))
      } else {
        new_LR <- as.data.frame(calc_lr(a,b,c,d, PATHOGENIC_OR_THRESHOLD=target_odds_BRCA2, CONFIDENCE_LEVEL=ci))
      }
      new_LR <- mutate(new_LR, Target_Odds := target_odds_BRCA2)
    } else if (gene_name == "PALB2") {
      if (df[i,hgvs] %in% palb2_highcount & "LR_path_ambry" %in% new_colnames) {
        new_LR <- data.frame(LR_path = c("NC - not calculable"), ACMG_points = c("NC - not calculable"))
      } else if (a==0 & c==0) {
        new_LR <- data.frame(LR_path = c("NA"), ACMG_points = c("NA"))
      } else {
        new_LR <- as.data.frame(calc_lr(a,b,c,d, PATHOGENIC_OR_THRESHOLD=target_odds_PALB2, CONFIDENCE_LEVEL=ci))
      }
      new_LR <- mutate(new_LR, Target_Odds := target_odds_PALB2)
    } else if (gene_name == "ATM") {
      if (df[i,hgvs] %in% atm_highcount & "LR_path_ambry" %in% new_colnames) {
        new_LR <- data.frame(LR_path = c("NC - not calculable"), ACMG_points = c("NC - not calculable"))
      } else if (a==0 & c==0) {
        new_LR <- data.frame(LR_path = c("NA"), ACMG_points = c("NA"))
      } else {
        new_LR <- as.data.frame(calc_lr(a,b,c,d, PATHOGENIC_OR_THRESHOLD=target_odds_ATM, CONFIDENCE_LEVEL=ci))
      }
      new_LR <- mutate(new_LR, Target_Odds := target_odds_ATM)
    } else if (gene_name == "CHEK2") {
      if (df[i,hgvs] %in% chek2_highcount & "LR_path_ambry" %in% new_colnames) {
        new_LR <- data.frame(LR_path = c("NC - not calculable"), ACMG_points = c("NC - not calculable"))
      } else if (a==0 & c==0) {
        new_LR <- data.frame(LR_path = c("NA"), ACMG_points = c("NA"))
      } else {
        new_LR <- as.data.frame(calc_lr(a,b,c,d, PATHOGENIC_OR_THRESHOLD=target_odds_CHEK2, CONFIDENCE_LEVEL=ci))
      }
      new_LR <- mutate(new_LR, Target_Odds := target_odds_CHEK2)
    }
    
    # combine OR and LR stats into one variant df
    new_stats <- cbind(new_OR, new_LR)
    colnames(new_stats) <- new_colnames
    
    # add all stats for latest variant to main df
    newdata_df <- rbind(newdata_df, new_stats)
  }
  cbind(df, newdata_df)
}

# filtering datasets based on parameters (single obs, AOC coverage, discrepancies between datasets, potential moderate penetrance)
parameter_filtering <- function(df, dataset, case, control, case_denom, control_denom, points, upper_or, lower_or, control_case_thresh, case_control_thresh, gene="gene", OR_UPPER_LIMIT_BRCA=4, OR_UPPER_LIMIT_PALB2=3, OR_UPPER_LIMIT_ATMCHEK2=2, OR_LOWER_LIMIT=1, EXCLUDE_ZERO=TRUE, EXCLUDE_SINGLE=TRUE, EXCLUDE_DISCREPANT=TRUE, EXCLUDE_MODERATE=TRUE) {
  df %>% mutate(zero_obs = 0, single_obs=0, discrepant_evidence=0, low_coverage=0, moderate=0) -> df
  
  ## Parameter 1.1: Zero observation across cases/controls
  if (EXCLUDE_ZERO==TRUE) {
    df %>% mutate(zero_obs := if_else((!!sym(case) == 0 & !!sym(control) == 0),1,0)) -> df
  }
  ## Parameter 1.2: Single observation across cases/controls
  if (EXCLUDE_SINGLE==TRUE) {
    df %>% mutate(single_obs := if_else((!!sym(case) == 1 & !!sym(control) == 0) | (!!sym(case) == 0 & !!sym(control) == 1),1,0)) -> df
  }

  # Parameter 2: Case:control ratio exceeds limit for appropriate use (i.e. would give false positive or false negative evidence)
  if (EXCLUDE_DISCREPANT==TRUE) {
    df %>% mutate(discrepant_evidence := case_when((!!sym(control_denom) / !!sym(case_denom)) > !!sym(control_case_thresh) ~ 1,
                                                   (!!sym(case_denom) / !!sym(control_denom)) > !!sym(case_control_thresh) ~ 1,
                                                   .default = 0)) -> df
  }
  
  ## Parameter 3: Variants of likely reduced penetrance (ie target OR is higher than actual OR, set to target OR of 4)
  if (EXCLUDE_MODERATE==TRUE) {
    df %>% mutate(moderate := case_when(!!sym(upper_or) <= OR_UPPER_LIMIT_BRCA & !!sym(lower_or) >= OR_LOWER_LIMIT & !!sym(gene) == "BRCA1" ~ 1,
                                        !!sym(upper_or) <= OR_UPPER_LIMIT_BRCA & !!sym(lower_or) >= OR_LOWER_LIMIT & !!sym(gene) == "BRCA2" ~ 1,
                                        !!sym(upper_or) <= OR_UPPER_LIMIT_PALB2 & !!sym(lower_or) >= OR_LOWER_LIMIT & !!sym(gene) == "PALB2" ~ 1,
                                        !!sym(upper_or) <= OR_UPPER_LIMIT_ATMCHEK2 & !!sym(lower_or) >= OR_LOWER_LIMIT & !!sym(gene) == "ATM" ~ 1,
                                        !!sym(upper_or) <= OR_UPPER_LIMIT_ATMCHEK2 & !!sym(lower_or) >= OR_LOWER_LIMIT & !!sym(gene) == "CHEK2" ~ 1,
                                        .default = 0)) -> df
  }
  
  df %>%
    mutate(include = if_else(zero_obs+single_obs+low_coverage+discrepant_evidence+moderate == 0, "Y", "N")) %>%
    select(zero_obs,single_obs,low_coverage,discrepant_evidence,moderate,include) %>%
    rename_with( ~ paste0(.x, "_", dataset)) -> output
  
  output_df <- cbind(df, output)
  return(output_df)
}

# plotting output evidence strengths
plot_prep <- function(df, col, ev) {
  df %>% 
    mutate(!!sym(ev) := as.numeric(!!sym(ev))) %>%
    mutate(!!sym(col) := case_when(!!sym(ev) <= -8 ~ "B_VSTR",
                                   !!sym(ev) <= -4 ~ "B_STR",
                                   !!sym(ev) <= -2 ~ "B_MOD",
                                   !!sym(ev) <= -1 ~ "B_SUP",
                                   !!sym(ev) >= 8 ~ "P_VSTR",
                                   !!sym(ev) >= 4 ~ "P_STR",
                                   !!sym(ev) >= 2 ~ "P_MOD",
                                   !!sym(ev) >= 1 ~ "P_SUP",
                                   is.na(!!sym(ev)) ~ "None"), .after=!!sym(ev)) -> df
  return(df)
}

# case:control max denominator imbalance lookup (see Sup Table 14)
read_csv("Analysis/P2_enriched/denom_limits.csv") -> denom_limits

# identify the dataset imbalance size at which false positive evidence starts to occur (used to create denom_limits lookup file)
# thresholds where there are 0 case carriers and 1 or more control carriers (number of controls set at 1 by default)
ratio_threshold_controls <- function(controls=1, ci=0, dataset_ratio=10, target) {
  ## all cases = 0, all controls = 1, want to find the ratio at which acmg points are >0 (and >1 for a permissive slant)
  control_denom = c(seq(500,25000000,by=500))
  control_carrier = c(rep(controls,length(control_denom)))
  case_denom = c(rep(10000,length(control_denom)))
  case_carrier = c(rep(0,length(control_denom)))
  ratio = c(control_denom/case_denom[1])
  gene = c(rep("BRCA1",length(control_denom))) #BRCA1 is just a placeholder, the actual gene does not matter as the 'target' variable sets the actual target odds ratio
  ratio_threshold <- data.frame(hgvs="TEST",case_carrier, case_denom, control_carrier, control_denom, ratio, gene)
  
  ratio_threshold_stats <- batch_stats(ratio_threshold, "hgvs", "case_carrier", "case_denom", "control_carrier", "control_denom", target_odds_BRCA1=target, ci=ci)
  ratio_zero <- ratio_threshold_stats[which.min(abs(ratio_threshold_stats$ACMG_points-0)), 6]
  return(ratio_zero)
}

# thresholds where there are 0 control carriers and 1 or more case carriers (number of case carriers set at 1 by default)
ratio_threshold_case <- function(cases=1, ci=0, dataset_ratio=10, target=4) {
  case_denom = c(seq(500,2500000,by=500))
  case_carrier = c(rep(cases,length(case_denom)))
  control_denom = c(rep(10000,length(case_denom)))
  control_carrier = c(rep(0,length(case_denom)))
  ratio = c(case_denom/control_denom[1])
  gene = c(rep("BRCA1",length(case_denom))) #BRCA1 is just a placeholder, the actual gene does not matter as the 'target' variable sets the actual target odds ratio
  ratio_threshold <- data.frame(hgvs="TEST",control_carrier, control_denom, case_carrier, case_denom, ratio, gene)
  ratio_threshold_stats <- batch_stats(ratio_threshold, "hgvs", "case_carrier", "case_denom", "control_carrier", "control_denom", target_odds_BRCA1=target)
  ratio_zero <- ratio_threshold_stats[which.min(abs(ratio_threshold_stats$ACMG_points-0)), 6]
  return(ratio_zero)
}


# Sup Figure 5: Plotting the dataset imbalance sizes
ratio_threshold_plot_controls <- function(controls=1, ci=0, dataset_ratio=10, target, title="Control:Case plot") {
  ## all cases = 0, all controls = 1, want to find the ratio at which acmg points are >0 (and >1 for a permissive slant)
  control_denom = c(seq(500,250000,by=500))
  control_carrier = c(rep(controls,length(control_denom)))
  case_denom = c(rep(10000,length(control_denom)))
  case_carrier = c(rep(0,length(control_denom)))
  ratio = c(control_denom/case_denom[1])
  gene = c(rep("BRCA1",length(control_denom))) #BRCA1 is just a placeholder, the actual gene does not matter as the 'target' variable sets the actual target odds ratio
  ratio_threshold <- data.frame(hgvs="TEST",case_carrier, case_denom, control_carrier, control_denom, ratio, gene)
  
  ratio_threshold_stats <- batch_stats(ratio_threshold, "hgvs", "case_carrier", "case_denom", "control_carrier", "control_denom", target_odds_BRCA1=target, ci=ci)
  ratio_zero <- ratio_threshold_stats[which.min(abs(ratio_threshold_stats$ACMG_points-0)), 6]
  print(ratio_zero)
  # # uncomment the below line to include a plot point of your dataset ratio against the two thresholds
  dataset_points <- ratio_threshold_stats[which.min(abs(ratio_threshold_stats$ratio-dataset_ratio)), 13]
  
  ratio_threshold_stats %>%
    ggplot(aes(ratio)) +
    # ratio:points curve
    geom_line(aes(y=ACMG_points), color="black", linewidth=1) +
    # threshold boundary points=0
    geom_segment(x=0,y=0,xend=ratio_zero,yend=0) +
    geom_segment(x=ratio_zero,y=0,xend=ratio_zero,yend=min(ratio_threshold_stats$ACMG_points)) +
    geom_rect(mapping = aes(xmin=0,xmax=ratio_zero,ymax=0,ymin=min(ratio_threshold_stats$ACMG_points)),fill="#3182BD",alpha=0.005) +
    # # plot the dataset ratio; uncomment to include a plot point of your dataset ratio against the two thresholds
    geom_segment(x=0,y=dataset_points,xend=dataset_ratio,yend=dataset_points, linetype="dashed", linewidth=1, colour="red") +
    geom_point(x=dataset_ratio, y=dataset_points, size=5, colour="red") +
    geom_segment(x=dataset_ratio,y=dataset_points,xend=dataset_ratio,yend=min(ratio_threshold_stats$ACMG_points), linetype="dashed", linewidth=1, colour="red") +
    
    theme_bw() +
    theme(text=element_text(size=16)) +
    scale_y_continuous(breaks=seq(round(min(ratio_threshold_stats$ACMG_points)),4,1), limits=c(-12,3), expand=c(0,0)) +
    scale_x_continuous(breaks=seq(0,25,1), limits=c(0,12), expand=c(0.01,0)) +
    ggtitle("Total controls exceed total cases",subtitle=title) +
    xlab("Control:Case Ratio")
}

ratio_threshold_plot_case <- function(cases=1, ci=0, dataset_ratio=10, target=4, title="Case:control plot") {
  case_denom = c(seq(500,250000,by=500))
  case_carrier = c(rep(cases,length(case_denom)))
  control_denom = c(rep(10000,length(case_denom)))
  control_carrier = c(rep(0,length(case_denom)))
  ratio = c(case_denom/control_denom[1])
  gene = c(rep("BRCA1",length(case_denom))) #BRCA1 is just a placeholder, the actual gene does not matter as the 'target' variable sets the actual target odds ratio
  ratio_threshold <- data.frame(hgvs="TEST",control_carrier, control_denom, case_carrier, case_denom, ratio, gene)
  
  ratio_threshold_stats <- batch_stats(ratio_threshold, "hgvs", "case_carrier", "case_denom", "control_carrier", "control_denom", target_odds_BRCA1=target)
  ratio_zero <- ratio_threshold_stats[which.min(abs(ratio_threshold_stats$ACMG_points-0)), 6]
  print(ratio_zero)
  # # uncomment the below line to include a plot point of your dataset ratio against the two thresholds
  dataset_points <- ratio_threshold_stats[which.min(abs(ratio_threshold_stats$ratio-dataset_ratio)), 13]
  
  ratio_threshold_stats %>%
    ggplot(aes(ratio)) +
    # ratio:points curve
    geom_line(aes(y=ACMG_points), color="black", linewidth=1) +
    # threshold boundary points=0
    geom_segment(x=0,y=0,xend=ratio_zero,yend=0) +
    geom_segment(x=ratio_zero,y=0,xend=ratio_zero,yend=3) +
    # # plot the dataset ratio; uncomment to include a plot point of your dataset ratio against the two thresholds
    geom_segment(x=0,y=dataset_points,xend=dataset_ratio,yend=dataset_points, linetype="dashed", linewidth=1, colour="red") +
    geom_point(x=dataset_ratio, y=dataset_points, size=5, colour="red") +
    geom_segment(x=dataset_ratio,y=dataset_points,xend=dataset_ratio,yend=3, linetype="dashed", linewidth=1, colour="red") +
    
    theme_bw() +
    theme(text=element_text(size=14)) +
    scale_y_continuous(breaks=seq(-12,4,1), limits=c(-12,3), expand=c(0,0)) +
    scale_x_continuous(breaks=seq(0,25,1), limits=c(0,12), expand=c(0.01,0)) +
    coord_trans(y="reverse") +
    ggtitle("Total cases exceed total controls",subtitle=title) +
    xlab("Case:Control Ratio")
}

# DATASET GENERATION ----
## BA1 Lookup File ----
## Creates a lookup file using gnomAD v4.1 for variants which meet stand-alone BA1 criteria for benignity

## read in gnomAD v4.1 exome+genome files
brca1 <- read_tsv("Data sources/gnomAD_v4.1/filtered_variants_BRCA1.vcf", comment="##")
brca2 <- read_tsv("Data sources/gnomAD_v4.1/filtered_variants_BRCA2.vcf", comment="##")
palb2 <- read_tsv("Data sources/gnomAD_v4.1/filtered_variants_PALB2.vcf", comment="##")
atm   <- read_tsv("Data sources/gnomAD_v4.1/filtered_variants_ATM.vcf", comment="##")
chek2 <- read_tsv("Data sources/gnomAD_v4.1/filtered_variants_CHEK2.vcf", comment="##")
gnomAD_all <- rbind(brca1, brca2, palb2, atm, chek2)

# Quality filter and data cleaning
## AC0: Genotype quality too low, allele count is 0 after removing low quality genotypes
## AS VSQR: Failed VSQR filter
## RF: Flag for failing the random forest model quality check
## Exomes filtered = the exome dataset failed one of the above tests (usually AC0). Same for Genomes filtered.
## If both filtered, remove variant, as failed both datasets. 
gnomAD_all %>%
  filter(FILTER!="BOTH_FILTERED") %>%
  mutate(gene = case_when(grepl("17",`#CHROM`) ~ "BRCA1",
                          grepl("13",`#CHROM`) ~ "BRCA2",
                          grepl("16",`#CHROM`) ~ "PALB2",
                          grepl("11",`#CHROM`) ~ "ATM",
                          grepl("22",`#CHROM`) ~ "CHEK2"),
         `#CHROM` = gsub("chr","", `#CHROM`),
         ID = paste("chr",`#CHROM`, POS, REF, ALT, sep="_"),
         `#CHROM` = as.double(`#CHROM`)) -> gnomAD_all_genes

write_csv(gnomAD_all_genes, "Data sources/gnomAD_v4.1/MTAF_VCEPs/gnomAD_all_genes.csv")

gnomAD_all_genes %>%
  mutate(
    # extract Grpmax FAF (which per gnomAD website does not include bottlenecked populations)
    grpmax_faf_v4 = case_when(FILTER=="PASS" ~ str_extract(INFO, "(?<=faf95_max_joint=).+?(?=;)"),
                              FILTER=="EXOMES_FILTERED" ~ str_extract(INFO, "(?<=faf95_max_genomes=).+?(?=;)"),
                              FILTER=="GENOMES_FILTERED" ~ str_extract(INFO, "(?<=faf95_max_exomes=).+?(?=;)")),
    across(ends_with("v4"), ~ as.numeric(.x))) -> gnomAD_all_genes_FAF

write_csv(gnomAD_all_genes_FAF, "Data sources/gnomAD_v4.1/MTAF_VCEPs/gnomAD_all_genes_prelift.csv")

gnomAD_all_genes_FAF %>%
  # Identify variants which are more common than the MTAF for that gene
  mutate(BA1_filter_gnomad = case_when(gene=="BRCA1" & grpmax_faf_v4>BRCA1_MTAF ~ "BA1", 
                                       gene=="BRCA2" & grpmax_faf_v4>BRCA2_MTAF ~ "BA1", 
                                       gene=="PALB2" & grpmax_faf_v4>PALB2_MTAF ~ "BA1", 
                                       gene=="ATM" & grpmax_faf_v4>ATM_MTAF ~ "BA1", 
                                       gene=="CHEK2" & grpmax_faf_v4>CHEK2_MTAF ~ "BA1", 
                                       .default="Rare")) %>%
  filter(BA1_filter_gnomad == "BA1") -> gnomAD_BA1

write_csv(gnomAD_BA1, "Data sources/gnomAD_v4.1/MTAF_VCEPs/gnomAD_BA1_filter_prelift.csv")

# liftover to build 37 using https://www.ensembl.org/Homo_sapiens/Tools/AssemblyConverter?db=core (Ensembl release v113)
# most of the datasets for this analysis use GrCh37, this stage allows merging with downstream datasets
read_tsv("Data sources/gnomAD_v4.1/MTAF_VCEPs/gnomAD_BA1_filter_postlift.vcf", comment="##") -> gnomAD_BA1

gnomAD_BA1 %>% 
  mutate(ID=paste(`#CHROM`,POS,REF,ALT,sep="-")) %>% 
  select(ID) %>% 
  write_csv("Data sources/gnomAD_v4.1/MTAF_VCEPs/vv_prep.csv")

# run liftover variants in the vv_prep csv file through variant validator https://variantvalidator.org/service/validate/batch/ 
read_tsv("Data sources/gnomAD_v4.1/MTAF_VCEPs/vv_output.txt",comment="#") -> hgvs_annotations #variant validator output file
hgvs_annotations %>%
  # use MANE select transcript and produce column using HGVS nomenclature
  filter(`Select transcript`=="MANE") %>%
  separate_wider_delim("HGVS_transcript",delim=":",names_sep="part") %>%
  select("Input","gene"="Gene_Symbol","hgvs_cdna"="HGVS_transcriptpart2") -> hgvs_annotations

gnomAD_BA1 %>%
  # annotate BA1 lookup file with HGVS nomenclature and chromosome position information in format to match other datasets
  mutate(Input=paste(`#CHROM`,POS,REF,ALT,sep="-")) %>% 
  left_join(hgvs_annotations, by="Input") %>%
  mutate(ID=paste0("chr",`#CHROM`,"_",POS,"_",REF,"_",ALT)) %>% 
  mutate(rare_flag="BA1_eligible") %>%
  select("ID","gene","hgvs_cdna","rare_flag") -> gnomAD_BA1_list

# output final BA1 lookup file
write_csv(gnomAD_BA1_list, "Data sources/gnomAD_v4.1/gnomAD_BA1_filter_final.csv")


## UK Biobank ----
## White female breast cancer cases
UKB_case_1 <- read_tsv("Data Sources/UKBiobank_Aug25/UKB_WES_full/BRCA1_BRCA2.white-female-counts.BC_stratified.Aug25.vcf", comment="##")
UKB_case_2 <- read_tsv("Data sources/UKBiobank_Aug25/UKB_WES_full/ATM-CHEK2-PALB2.white-female-counts.BC_stratified.Aug25.vcf", comment="##")

UKB_case_1 %>% select(-all_of(setdiff(colnames(UKB_case_1),colnames(UKB_case_2)))) -> UKB_case_1
UKB_case_2 %>% select(-all_of(setdiff(colnames(UKB_case_2),colnames(UKB_case_1)))) -> UKB_case_2

UKB_case <- rbind(UKB_case_1, UKB_case_2)

UKB_case %>% 
  # convert allele counts to counts of individuals
  mutate(White_BC_female_carriers = white_female_BC_AC - white_female_BC_hom,
         White_BC_female_total = white_female_BC_AN/2) %>%
  select(-starts_with("white_", ignore.case=FALSE)) -> UKB_case

# Total non-breast cancer controls (male and female)
UKB_control1 <- read_tsv("Data sources/UKBiobank_Aug25/UKB_WES_full/BRCA1_1.UKB_variants.counts.tallied.Aug25.vcf")
UKB_control2 <- read_tsv("Data sources/UKBiobank_Aug25/UKB_WES_full/BRCA1_2.UKB_variants.counts.tallied.Aug25.vcf")
UKB_control3 <- read_tsv("Data sources/UKBiobank_Aug25/UKB_WES_full/BRCA2_1.UKB_variants.counts.tallied.Aug25.vcf")
UKB_control4 <- read_tsv("Data sources/UKBiobank_Aug25/UKB_WES_full/BRCA2_2.UKB_variants.counts.tallied.Aug25.vcf")
UKB_control5 <- read_tsv("Data sources/UKBiobank_Aug25/UKB_WES_full/ATM_1.UKB_variants.counts.tallied.Aug25.vcf")
UKB_control6 <- read_tsv("Data sources/UKBiobank_Aug25/UKB_WES_full/CHEK2_1.UKB_variants.counts.tallied.Aug25.vcf")
UKB_control7 <- read_tsv("Data sources/UKBiobank_Aug25/UKB_WES_full/PALB2_1.UKB_variants.counts.tallied.Aug25.vcf")

UKB_control <- rbind(UKB_control1, UKB_control2, UKB_control3, UKB_control4, UKB_control5, UKB_control6, UKB_control7)
UKB_control %>%
  # create column of total tested (all non-carriers, all heterozygous carriers, and all homozygous carriers (het_hom = het carriers + hom carriers))
  mutate(White_nonBC_total = White_nonBC_het_hom+White_nonBC_non_carriers) %>%
  select(1:5, White_nonBC_total_carriers=White_nonBC_het_hom, White_nonBC_total) -> UKB_control

UKB_case %>% full_join(UKB_control, by=c("#CHROM","POS","ID","REF","ALT")) -> UKB

## variantvalidator check
# note: chromosomal coordiantes used in case HGVSc is different when using different transcripts
# VariantValidator run in batch with the 5 chosen transcipts as options
UKB %>% mutate(ID = gsub("_",":",ID)) -> UKB_annotated

write_csv(UKB_annotated, "Data sources/UKBiobank_Aug25/ukb_pervv.csv")
read_tsv("Data sources/UKBiobank_Aug25/ukb_postvv.txt",skip=2) -> ukb_vv #variant validator output

left_join(UKB_annotated, ukb_vv, by=c("ID"="Input")) -> UKB

# filter out variants which fail variantvalidator check
UKB %>% filter(!is.na(HGVS_transcript)) -> UKB_vv

# excluding PTVs which are excluded from the BRIDGES dataset (in preparation for later combination of datasets)
UKB_vv %>%
  separate_wider_delim(HGVS_transcript, names=c("gene","hgvs"),delim=":") %>%
  filter(White_BC_female_carriers > 0 | White_nonBC_total_carriers > 0) %>%
  ## amending variant validator HGVS for this variant due to alternate allele presence, which should be T>C as reported in UKB.
  mutate(hgvs = str_replace(hgvs,"c.7397=","c.7397T>C")) %>%
  ## exclude PTVs in the last exon
  mutate(final_exon = case_when(SYMBOL=="BRCA1" & EXON=="23/23" ~ 1,
                                SYMBOL=="BRCA2" & EXON=="27/27" ~ 1,
                                SYMBOL=="ATM" & EXON=="63/63" ~ 1,
                                SYMBOL=="CHEK2" & EXON=="15/15" ~ 1,
                                SYMBOL=="PALB2" & EXON=="13/13" ~ 1,
                                .default=0)) %>%
  filter((final_exon==0 & grepl(paste(consequence_filter,collapse="|"),Consequence)) | grepl("missense",Consequence)) %>%
  ## exclude splice variants affecting the penultimate exon for genes where exon skipping not confirmed pathogenic at those sites (per Dorling et al., 2021)
  filter((SYMBOL=="BRCA2" & !grepl(paste("c.9502-1G","c.9502-2A","c.9648+1G","c.9648+2T",sep="|"),HGVSc)) | SYMBOL %in% c("BRCA1","PALB2","ATM","CHEK2")) %>%
  filter((SYMBOL=="CHEK2" & !grepl(paste("c.1462-1G","c.1462-2A","c.1542+1G","c.1542+2T",sep="|"),HGVSc)) | SYMBOL %in% c("BRCA1","PALB2","ATM","BRCA2")) %>%
  ## exclude 7 canonical BRCA1 VUS variants per ENIGMA: c.594-2A>C, c.4096+1G>A, c.4096+2T>C, c.4096+1G>A, c.4186-2A>G, c.4358-1G>C, c.4358-2del
  filter((SYMBOL=="BRCA1" & !grepl(paste("c.594-2A>C","c.4096+1G>A","c.4096+2T>C","c.4096+1G>A","c.4186-2A>G","c.4358-1G>C","c.4358-2del",sep="|"),HGVSc)) | SYMBOL %in% c("CHEK2","PALB2","ATM","BRCA2")) %>%
  select(SYMBOL, hgvs, ID, Consequence, starts_with("White")) -> UKB_filtered

write_csv(UKB_filtered, "Data sources/UKBiobank_Aug25/UKB_femaleBC_totalnonBC_ptvmiss.csv")


# Split UKB controls across NDRS and UKB breast cancer cases
## separation of UKB controls into proportional segments (segment 1: for use with NDRS cases, segment 2: for use with UKB cases)
## aiming for best use of control data, division of controls should get the case/control ratio as close to 1 as possible for both case sets to reduce impact of imbalance (see Rowlands et al., 2024)
## Sup Figure 4: A case:control ratio of 6.62 is the lowest ratio obtainable across NDRS cases and UKB cases.
## Therefore, 122,232 controls are retained for UKB cases, and 297,141 controls are used for NDRS cases, and carriers are sampled where 29.1% are retained for UKB and 70.8% are 

UKB_filtered %>% mutate(White_nonBC_noncarrier = White_nonBC_total-White_nonBC_total_carriers) -> UKB_filtered

UKB_filtered %>%
  rowwise() %>%
  mutate(White_nonBC_total_noncarrier_UKBpart = round(White_nonBC_total*0.2914637),
         White_nonBC_total_noncarrier_NDRSpart = White_nonBC_total-White_nonBC_total_noncarrier_UKBpart,
         White_nonBC_total_carriers_UKBpart = sample_controls(White_nonBC_total_carriers),
         White_nonBC_total_carriers_NDRSpart = White_nonBC_total_carriers - White_nonBC_total_carriers_UKBpart) -> UKB_filtered_split

write_csv(UKB_filtered_split, "Data Sources/UKBiobank Data/WES/UKB_BRCA_Filtered_split_ptvmiss.csv")


## CARRIERS ----
carriers <- read_csv("Data sources/CARRIERS/CARRIERS_per_variant.csv")

# reformat nomeclature columns to align with other datasets
carriers %>%
  mutate(unique_id = paste0(CAVA_GENE, ":", c., ":", p.), .before = CAVA_CSN) %>%
  mutate(ID = paste0("chr",`#CHROM`,"_",GRCh37Location,"_",REF,"_",ALT), .before=unique_id) %>%
  select(-Sample_ID) %>%
  group_by(unique_id, ID, CAVA_SO, `case control status`, `race/ethnicity`) %>%
  summarise(count = n()) %>%
  pivot_wider(id_cols=c(unique_id,ID,CAVA_SO,`race/ethnicity`), names_from=`case control status`, values_from=`count`) %>%
  separate_wider_delim(unique_id, delim = ":", names=c("gene", "c.", "p.")) %>%
  mutate(c. = if_else((grepl("p.", p.) & !grepl("NA", p.)), c., paste(c., p., sep="_")),
         p. = if_else((grepl("p.", p.) & !grepl("NA", p.)), p., NA),
         c. = gsub("_NA","",c.),
         p. = gsub(" \\(p.Phe3065_Gln3066delinsLeuTer\\)","",p.),
         Case = replace_na(Case, 0),
         Control = replace_na(Control, 0)) %>%
  # this is needed because of the BRCA2 c.9196 variant with strange p. nomenclature, which we amended above
  group_by (gene, c., p., ID, CAVA_SO, `race/ethnicity`) %>%
  summarise (Case = sum(Case),
             Control = sum(Control)) -> carriers_per_var

# no VEP transcript annotations needed as carriers dataset already using build 37 versions of MANE select transcripts, and consequence already annotated

## variantvalidator check
carriers_per_var %>%
  mutate(vep_hgvs = case_when(gene == "BRCA1" ~ "NM_007294.3",
                              gene == "BRCA2" ~ "NM_000059.3",
                              gene == "PALB2" ~ "NM_024675.3",
                              gene == "ATM" ~ "NM_000051.3",
                              gene == "CHEK2" ~ "NM_007194.3"),
         vep_hgvs = paste0(vep_hgvs,":",c.)) -> carriers_per_var

write_csv(carriers_per_var, "Data sources/CARRIERS/carriers_prevv.csv")

read_tsv("Data sources/CARRIERS/carriers_postvv.txt",skip=2) %>% distinct() -> carriers_vv #variant validator output file

left_join(carriers_per_var, carriers_vv, by=c("vep_hgvs"="Input")) -> carriers_per_var

# filter out variants which fail variantvalidator check
carriers_per_var %>% filter(!is.na(HGVS_transcript)) -> carriers_per_var_vv

# excluding PTVs which are excluded from the BRIDGES dataset (in preparation for later combination of datasets)
carriers_per_var_vv %>%
  separate_wider_delim(HGVS_transcript, names=c("hgvs_transcript","hgvs"),delim=":") %>%
  mutate (Case_denom = 32247,
          Control_denom = 32544, .after=Control) %>%
  filter(`race/ethnicity` == "Non_Hispanic_White") %>%
  filter (gene %in% c("BRCA1","BRCA2","PALB2","ATM","CHEK2")) %>%
  ## exclude PTVs in the last exon
  mutate(final_exon = case_when(gene=="BRCA1" & GRCh37_POS <= 41197819 ~ 1,
                                gene=="BRCA2" & GRCh37_POS >= 32972299 ~ 1,
                                gene=="ATM" & GRCh37_POS >= 108236052 ~ 1,
                                gene=="CHEK2" & GRCh37_POS <= 29083974 ~ 1,
                                gene=="PALB2" & GRCh37_POS <= 23614990 ~ 1,
                                .default=0)) %>%
  filter((final_exon==0 & grepl(paste(consequence_filter,collapse="|"),CAVA_SO)) | grepl("missense",CAVA_SO)) %>%
  ## exclude splice variants affecting the penultimate exon for genes where exon skipping not confirmed pathogenic at those sites
  ### BRCA2: 9502-1,-2, 9648+1,+2. CHEK2: c.1462-1,-2, 1542+1,+2
  filter((gene=="BRCA2" & !grepl(paste("c.9502-1G","c.9502-2A","c.9648+1G","c.9648+2T",sep="|"),c.)) | gene %in% c("BRCA1","PALB2","ATM","CHEK2")) %>%
  filter((gene=="CHEK2" & !grepl(paste("c.1462-1G","c.1462-2A","c.1542+1G","c.1542+2T",sep="|"),c.)) | gene %in% c("BRCA1","PALB2","ATM","BRCA2")) %>%
  ## exclude 7 canonical BRCA1 VUS variants per ENIGMA: c.594-2A>C, c.4096+1G>A, c.4096+2T>C, c.4096+1G>A, c.4186-2A>G, c.4358-1G>C, c.4358-2del
  filter((gene=="BRCA1" & !grepl(paste("c.594-2A>C","c.4096+1G>A","c.4096+2T>C","c.4096+1G>A","c.4186-2A>G","c.4358-1G>C","c.4358-2del",sep="|"),c.)) | gene %in% c("CHEK2","PALB2","ATM","BRCA2")
  ) -> carriers_per_var_clean

write_csv (carriers_per_var_clean, "Data sources/CARRIERS/carriers_per_var_clean_ptvmiss.csv")


## BRIDGES ----
read_csv("Data sources/BRIDGES Data/ATM_CHEK2_BRCA1_BRCA2_PALB2_missense.csv") -> BRIDGES_missense
read_csv("Data sources/BRIDGES Data/ATM_CHEK2_BRCA1_BRCA2_PALB2_truncating.csv") -> BRIDGES_ptv

# VEP transcript annotations not needed as bridges dataset already using build 37 versions of MANE select transcripts, and consequence already annotated

# generate list of IDs for variant validator 
## (need hgvs nomenclature for merging with other datasets, genomic coords not suitable because the coords in UKB aren't always accurate/appropriate, like chr13_32911461_TGAACAAATGGGCA_TGGACAAATGGGCA for a missense variant)
# note: chromosomal coordiantes used in case HGVSc is different when using different transcripts
# VariantValidator run in batch with the 5 chosen transcipts as options
BRIDGES_missense %>%
  filter(Cases>0 | Controls >0) %>%
  mutate(VV_ID = ID,
         VV_ID = str_replace_all(VV_ID,"_","-"),
         VV_ID = str_replace_all(VV_ID,"chr","")) -> BRIDGES_missense
BRIDGES_missense %>% select(VV_ID) -> VV_IDs
write_csv(VV_IDs, "Data sources/BRIDGES Data/missense_vv_ids.csv")

BRIDGES_ptv %>%
  filter(Cases>0 | Controls >0) %>%
  mutate(VV_ID = ID,
         VV_ID = str_replace_all(VV_ID,"_","-"),
         VV_ID = str_replace_all(VV_ID,"chr","")) -> BRIDGES_ptv
BRIDGES_ptv %>% select(VV_ID) -> VV_IDs
write_csv(VV_IDs, "Data sources/BRIDGES Data/truncating_vv_ids.csv")

# read in variant validator outputs and clean up BRIDGES dataset for merging
read_tsv("Data sources/BRIDGES Data/missense_vv_output.txt",comment="#") -> BRIDGES_vv_missense
read_tsv("Data sources/BRIDGES Data/truncating_vv_output.txt",comment="#") -> BRIDGES_vv_ptv

BRIDGES_missense %>% 
  left_join(BRIDGES_vv_missense, by=c("VV_ID"="Input")) %>%
  # remove variants which fail variant validator check
  filter(!is.na(HGVS_transcript)) %>%
  separate_wider_delim(HGVS_transcript, names=c("gene","hgvs"),delim=":") %>%
  mutate(Case_denom = 42062, Control_denom = 44035, .after=Controls) %>%
  mutate(Gene_Symbol = case_when(grepl("chr11", ID) ~ "ATM",
                                 grepl("chr17", ID) ~ "BRCA1",
                                 grepl("chr13", ID) ~ "BRCA2",
                                 grepl("chr16", ID) ~ "PALB2",
                                 grepl("chr22", ID) ~ "CHEK2"), .before=ID) -> BRIDGES_processed_missense
write_csv(BRIDGES_processed_missense, "Data sources/BRIDGES Data/BRIDGES_missense_processed.csv")

BRIDGES_ptv %>% 
  left_join(BRIDGES_vv_ptv, by=c("VV_ID"="Input")) %>%
  # remove variants which fail variant validator check
  filter(!is.na(HGVS_transcript)) %>%
  separate_wider_delim(HGVS_transcript, names=c("gene","hgvs"),delim=":") %>%
  mutate(Case_denom = 42062, Control_denom = 44035, .after=Controls) %>%
  mutate(Gene_Symbol = case_when(grepl("chr11", ID) ~ "ATM",
                                 grepl("chr17", ID) ~ "BRCA1",
                                 grepl("chr13", ID) ~ "BRCA2",
                                 grepl("chr16", ID) ~ "PALB2",
                                 grepl("chr22", ID) ~ "CHEK2"), .before=ID) -> BRIDGES_processed_trunc
write_csv(BRIDGES_processed_trunc, "Data sources/BRIDGES Data/BRIDGES_truncating_processed.csv")

BRIDGES <- rbind(BRIDGES_processed_trunc,BRIDGES_processed_missense) 
write_csv(BRIDGES,"Data sources/BRIDGES Data/bridges_ptvmiss.csv")


## NDRS-UK ----
read_csv("Data Sources/NDRS_femalebreast_2025/output_1_for_non_snv.csv") -> ndrs_indels #annotated, but not pre-filtered for incorrect nomenclature
read_csv("Data Sources/NDRS_femalebreast_2025/output_to_classify.csv") -> ndrs_singles #pre-filtered for incorrect nomenclature

# nomenclature fixes for indels (removes extraneous lettering and regex errors)
ndrs_indels %>%
  mutate(hgvs_cdna = gsub("del]", "del", hgvs_cdna),
         hgvs_cdna = gsub("del_","del",hgvs_cdna),
         hgvs_cdna = gsub("bp","",hgvs_cdna),
         hgvs_cdna = gsub("c.3847_c.3848heT_dell","c.3847_3848del",hgvs_cdna),
         hgvs_cdna = gsub("c.n", "c.", hgvs_cdna),
         hgvs_cdna = gsub("heT_", "", hgvs_cdna),
         hgvs_cdna = gsub("heTdel", "del", hgvs_cdna),
         hgvs_cdna = gsub("del-p", "del", hgvs_cdna),
         hgvs_cdna = gsub("dup-p", "del", hgvs_cdna),
         hgvs_cdna = gsub("c./", "c.", hgvs_cdna),
         hgvs_cdna = gsub("1)2", "1_2", hgvs_cdna),
         hgvs_cdna = gsub("7-1", "7_1", hgvs_cdna),
         hgvs_cdna = gsub("\\(\\)", "", hgvs_cdna),
         hgvs_cdna = gsub("del9ins", "delins", hgvs_cdna),
         hgvs_cdna = gsub("c.3680_3681del\\(leu1227Glnfs*5\\)", "c.3680_3681del", hgvs_cdna),
         hgvs_cdna = gsub("c.10095delinsGAATTATATCT\\(p.ser3366fs\\)", "c.10095delinsGAATTATATCT", hgvs_cdna),
         hgvs_cdna = gsub("(r645fsx15)", "", hgvs_cdna),
         hgvs_cdna = gsub("c.2800_01del", "c.2800_2801del", hgvs_cdna),
         hgvs_cdna = gsub("c.3847_48del", "c.3847_3848del", hgvs_cdna),
         hgvs_cdna = gsub("c.1942-1945del","c.1942_1945del",hgvs_cdna),
         hgvs_cdna = gsub("c.4065_68del", "c.4065_4068del", hgvs_cdna),
         hgvs_cdna = gsub("c.4065_4068del(p.Asn1355fs)","c.4065_4068del",hgvs_cdna),
         hgvs_cdna = gsub("c.4065-4068del","c.4065_4068del",hgvs_cdna),
         hgvs_cdna = gsub("c.6280_81del", "c.6280_6281del", hgvs_cdna),
         hgvs_cdna = gsub("c.6763_64insA", "c.6763_6764insA", hgvs_cdna),
         hgvs_cdna = gsub("c.983_86del", "c.983_986del", hgvs_cdna),
         hgvs_cdna = gsub("c.5266dup-p", "c.5266dupC", hgvs_cdna),
         hgvs_cdna = gsub("c.5266dupC,c.288C>T", "c.5266dupC", hgvs_cdna),
         hgvs_cdna = gsub("c.2157del\\(\\)", "c.2157del", hgvs_cdna),
         hgvs_cdna = gsub("c.3847_3847del", "c.3847_3848del", hgvs_cdna),
         hgvs_cdna = gsub("c.6275_6275del", "c.6275_6276del", hgvs_cdna),
         hgvs_cdna = gsub("c.3481_3491del11\\*\\*","c.3481_3491del11",hgvs_cdna),
         hgvs_cdna = gsub("c.891_899delinsGATACTTAGc.891_899delinsGATACTTCAG", "c.891_899delinsGATACTTCAG", hgvs_cdna),
         hgvs_cdna = gsub("c.3036-3039del","c.3036_3039del",hgvs_cdna)
  ) -> ndrs_indels

# Initial formatting of NDRS data to group based on cDNA HGVS
ndrs_singles %>%
  select(gene, hgvs_cdna, WhiteCount) %>%
  mutate(WhiteCount = replace_na(WhiteCount, 0)) %>%
  group_by (gene, hgvs_cdna) %>%
  summarise (across (everything (), ~ sum(.x))) %>%
  mutate(vartype="missense") %>%
  filter(WhiteCount != 0) -> ndrs_singles_grouped

ndrs_indels %>%
  select(gene, hgvs_cdna, WhiteCount=`Full screen count`) %>%
  mutate(WhiteCount = replace_na(WhiteCount, 0)) %>%
  group_by (gene, hgvs_cdna) %>%
  summarise (across (everything (), ~ sum(.x))) %>%
  mutate(vartype="PTV") %>%
  filter(WhiteCount != 0) -> ndrs_indels_grouped

rbind(ndrs_singles_grouped, ndrs_indels_grouped) -> ndrs

# VEP annotation - converts HGVS to VEP format using same transcripts as other datasets (see Supplementary Methods)
ndrs %>%
  mutate(vep_hgvs = if_else(gene == "BRCA1",
                            # BRCA1 transcript
                            paste0("ENST00000357654.3", ":", hgvs_cdna),
                            # BRCA2 transcript
                            paste0("ENST00000544455.1", ":", hgvs_cdna))) -> ndrs

write_csv(ndrs, "Data Sources/NDRS_femalebreast_2025/ndrs_vep_file.csv")

# VEP call results
read_tsv("Data Sources/NDRS_femalebreast_2025/vep_annotation_output_preVEPcleaning.txt") -> vep_annotations

vep_annotations %>%
  filter(grepl("ENST00000357654.3|ENST00000544455.1", Feature)) %>%
  select("#Uploaded_variation", "Location", "Consequence", "SYMBOL", "EXON", "INTRON") %>%
  filter(grepl("BRCA", SYMBOL)) -> vep_annotations

left_join(ndrs, vep_annotations, join_by("vep_hgvs" == "#Uploaded_variation")) -> ndrs_annotated

# variants not identified by VEP:
ndrs_annotated %>% filter(is.na(Consequence)) -> ndrs_excluded

# variants which were identified by VEP:
filter(ndrs_annotated, !is.na(Consequence)) -> ndrs_annotated

# variant validator file (standardise the hgvs nomenclature to match the other datasets)
ndrs_annotated %>%
  mutate(vep_hgvs = case_when(gene == "BRCA1" ~ "NM_007294.3",
                              gene == "BRCA2" ~ "NM_000059.3",
                              gene == "PALB2" ~ "NM_024675.3",
                              gene == "ATM" ~ "NM_000051.3",
                              gene == "CHEK2" ~ "NM_007194.3"),
         vep_hgvs = paste0(vep_hgvs,":",hgvs_cdna)) -> ndrs_annotated

write_csv(ndrs_annotated, "Data Sources/NDRS_femalebreast_2025/ndrs_prevv.csv")
read_tsv("Data sources/NDRS_femalebreast_2025/ndrs_postvv.txt",skip=2) %>% distinct() -> ndrs_vv #variant validator output
left_join(ndrs_annotated, ndrs_vv, by=c("vep_hgvs"="Input")) -> ndrs_annotated

# filter out variants which fail variantvalidator check
ndrs_annotated %>% filter(!is.na(HGVS_transcript)) -> ndrs_annotated_vv

# excluding PTVs which are excluded from the BRIDGES dataset (in preparation for later combination of datasets)
ndrs_annotated_vv %>%
  separate_wider_delim(HGVS_transcript, names=c("hgvs_transcript","hgvs"), delim=":") %>%
  # sum counts across variants which are now identical ie have same nomenclature + consequence
  group_by (gene, hgvs, Consequence, EXON, INTRON) %>%
  summarise (WhiteCount = sum(WhiteCount)) %>%
  mutate(denominator = 44917) %>%
  ## exclude PTVs in the last exon and restrict to just PTVs and missense variants
  mutate(final_exon = case_when(gene=="BRCA1" & EXON=="23/23" ~ 1,
                                gene=="BRCA2" & EXON=="27/27" ~ 1,
                                .default=0)) %>%
  filter((final_exon==0 & grepl(paste(consequence_filter,collapse="|"),Consequence)) | grepl("missense",Consequence)) %>%
  ## exclude splice variants affecting the penultimate exon for genes where exon skipping not confirmed pathogenic at those sites
  ### BRCA2: 9502-1,-2, 9648+1,+2. CHEK2: c.1462-1,-2, 1542+1,+2
  filter((gene=="BRCA2" & !grepl(paste("c.9502-1G","c.9502-2A","c.9648+1G","c.9648+2T",sep="|"),hgvs)) | gene %in% c("BRCA1","PALB2","ATM","CHEK2")) %>%
  ## exclude 7 canonical BRCA1 VUS variants per ENIGMA: c.594-2A>C, c.4096+1G>A, c.4096+2T>C, c.4096+1G>A, c.4186-2A>G, c.4358-1G>C, c.4358-2del
  filter((gene=="BRCA1" & !grepl(paste("c.594-2A>C","c.4096+1G>A","c.4096+2T>C","c.4096+1G>A","c.4186-2A>G","c.4358-1G>C","c.4358-2del",sep="|"),hgvs)) | gene %in% c("CHEK2","PALB2","ATM","BRCA2")
  ) -> ndrs_annotated_cleaned

# merge with partition of UKB controls
read_csv("Data Sources/UKBiobank Data/WES/UKB_BRCA_Filtered_split_ptvmiss.csv") -> UKB_filtered_split

ndrs_annotated_cleaned %>%
  ungroup() %>%
  full_join(UKB_filtered_split, by=c("gene"="SYMBOL","hgvs")) %>%
  filter (gene %in% c("BRCA1", "BRCA2")) %>%
  select(gene,hgvs,WhiteCount,denominator,White_nonBC_total_carriers_NDRSpart,White_nonBC_total_noncarrier_NDRSpart) -> ndrs_ukb_combined

write_csv(ndrs_ukb_combined, "Data Sources/NDRS_femalebreast_2025/ndrs_ptvmiss_whiteonly.csv")


## AMBRY + gnomAD ----
### Ambry ----
read_csv("Data sources/Ambry/file_1.csv") -> file1
read_csv("Data sources/Ambry/file_2.csv") -> file2
read_csv("Data sources/Ambry/file_3.csv") -> file3
read_csv("Data sources/Ambry/file_4.csv") -> file4
read_csv("Data sources/Ambry/file_5.csv") -> file5
read_csv("Data sources/Ambry/file_6.csv") -> file6
bind_rows(file1,file2,file3,file4,file5,file6) -> Ambry

Ambry %>%
  distinct() %>%
  filter(symbol %in% c("ATM","CHEK2","PALB2","BRCA1","BRCA2")) %>%
  mutate(nucleotide_id = case_when(symbol=="ATM" ~ "NM_000051.4",
                                   symbol=="CHEK2" ~ "NM_007194.4",
                                   symbol=="PALB2" ~ "NM_024675.4",
                                   symbol=="BRCA1" ~ "NM_007294.4",
                                   symbol=="BRCA2" ~ "NM_000059.4")) %>%
  mutate(ID=paste0(nucleotide_id,":",c_variant),.before=gender) %>%
  mutate(gene_ID=paste0(symbol,":",c_variant),.before=gender) %>%
  mutate(ID=gsub("DEL","del",ID)) %>%
  mutate(ID=gsub("DUP","dup",ID)) %>%
  mutate(ID=gsub("INS","ins",ID)) %>%
  mutate(ID=gsub("Ex","EX",ID)) %>%
  mutate(gene_ID=gsub("DEL","del",gene_ID)) %>%
  mutate(gene_ID=gsub("DUP","dup",gene_ID)) %>%
  mutate(gene_ID=gsub("INS","ins",gene_ID)) %>%
  mutate(gene_ID=gsub("Ex","EX",gene_ID)) %>%
  mutate(genotype=gsub("N/A",NA,genotype)) %>%
  mutate(phenotype=case_when(gender=="female" & !is.na(`PHx Breast`) ~ "female_breast",
                             gender=="female" & !is.na(`PHx OvCa`) ~ "female_ovarian",
                             gender=="male" & !is.na(`PHx Breast`) ~ "male_breast",
                             gender=="male" & !is.na(`PHx OvCa`) ~ "male_ovarian",
                             (gender=="unknown" | is.na(gender)) & !is.na(`PHx Breast`) ~ "other_breast",
                             (gender=="unknown" | is.na(gender)) & !is.na(`PHx OvCa`) ~ "other_ovarian")) %>%
  mutate(fhx = case_when(!is.na(`FHx Breast`) &  !is.na(`FHX OvCa`) ~ "fhx_breastovarian",
                         !is.na(`FHx Breast`) ~ "fhx_breast",
                         !is.na(`FHX OvCa`) ~ "fhx_ovarian")) -> ambry_interim

# column of total number of variant observations
ambry_interim %>%
  select(ID,gene_ID,genotype,phenotype,fhx) %>%
  group_by(ID,gene_ID) %>%
  summarise(n=n()) %>%
  ungroup() -> ambry_totalno

# number of each gentoype observed
ambry_interim %>%
  select(ID,gene_ID,genotype,phenotype,fhx) %>%
  group_by(ID,gene_ID,genotype) %>%
  summarise(n=n()) %>%
  mutate(genotype=replace_na(genotype, "PHASE_UNK")) %>%
  ungroup() %>%
  pivot_wider(id_cols=c(ID,gene_ID),names_from=genotype,values_from=n) -> ambry_genotypes

# number of each cancer phenotype observed
ambry_interim %>%
  select(ID,gene_ID,genotype,phenotype,fhx) %>%
  group_by(ID,gene_ID,phenotype) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  pivot_wider(id_cols=c(ID,gene_ID),names_from=phenotype,values_from=n) -> ambry_phenotypes

# combine into one summary table
ambry_totalno %>%
  left_join(ambry_genotypes,by=c("ID","gene_ID")) %>%
  left_join(ambry_phenotypes,by=c("ID","gene_ID")) %>%
  mutate(breast_denom=187642,.after=female_breast) -> ambry_summary

write_csv(ambry_summary,"Data sources/Ambry/ambry_per_variant.csv")

# add VEP annotations. Use gene_id column to submit data through VEP.

read_tsv("Data sources/Ambry/VEP_build38.txt") -> build38_VEP #VEP output file
build38_VEP %>% 
  filter(Feature %in% c("ENST00000675843.1","ENST00000357654.9","ENST00000404276.6","ENST00000261584.9","ENST00000380152.8")) %>%
  select(`#Uploaded_variation`, Consequence, EXON) -> build38_VEP

ambry_summary %>% 
  full_join(build38_VEP,by=c("gene_ID"="#Uploaded_variation")) %>%
  # filter for variants seen at least once in female breast cancer patients
  filter(!is.na(female_breast)) -> ambry_summary_VEP

# run through variant validator (confirm nomenclature is standard/accurate before merging)
## NOTE: VV run using Build 38 transcripts here using the 'ID' column!!

write_csv(ambry_summary_VEP, "Data sources/Ambry/ambry_pre_vv.csv")

read_tsv("Data sources/Ambry/ambry_postvv.txt",skip=2) %>% distinct() -> ambry_vv #variant validator output file

left_join(ambry_summary_VEP, ambry_vv, by=c("ID"="Input")) -> ambry_summary_VEP

# filter out variants which fail variantvalidator check
ambry_summary_VEP %>% filter(!is.na(HGVS_transcript)) -> ambry_summary_VEP_VV

# note: liftover done as part of variantvalidator check in build 38 (provides coordinates for both build 37 and 38 when providing gene and hgvs)
# re-load build 37 coordinates with correct transcripts to re-generate hgvs in case any difference between builds
ambry_summary_VEP_VV %>%
  mutate(ID_build37 = paste0("chr",GRCh37_CHR,":",GRCh37_POS,":",GRCh37_REF,":",GRCh37_ALT),.before=ID) -> ambry_summary_VEP_VV
write_csv(ambry_summary_VEP_VV, "Data sources/Ambry/ambry_per_variant_VEP_VV_prelift.csv")

read_tsv("Data sources/Ambry/ambry_per_variant_VEP_VV_postlift.txt", skip=2) -> ambry_vep_build37 #variant validator output file using build 37 genomic coordinates

ambry_summary_VEP_VV %>%
  left_join(ambry_vep_build37, by=c("ID_build37"="Input")) %>%
  # some variants appear duplicated because they were listed in the Ambry file under both old and new HGVS nomenclature (so left join appears many-to-many)
  # nomenclature has been standardized using variantvalidator, therefore can be combined. 
  distinct() %>%
  select(-ID,-gene_ID) %>%
  separate_wider_delim(HGVS_transcript.y, names=c("gene","hgvs"), delim=":") %>%
  group_by (ID_build37, gene, hgvs) %>%
  summarise (female_breast = sum(female_breast),
             breast_denom=breast_denom) %>%
  distinct() -> ambry_counts_fixed

write_csv(ambry_counts_fixed, "Data sources/Ambry/ambry_counts_fixed.csv")


### gnomAD v4.1 ----
### Load in exome data files. NOTE: gnomAD v4.1 is in build 38
brca1 <- read_tsv("Data sources/gnomAD_v4.1/exomes_BRCA1.vcf", comment="##")
brca2 <- read_tsv("Data sources/gnomAD_v4.1/exomes_BRCA2.vcf", comment="##")
palb2 <- read_tsv("Data sources/gnomAD_v4.1/exomes_PALB2.vcf", comment="##")
atm   <- read_tsv("Data sources/gnomAD_v4.1/exomes_ATM.vcf", comment="##")
chek2 <- read_tsv("Data sources/gnomAD_v4.1/exomes_CHEK2.vcf", comment="##")
gnomAD_all <- rbind(brca1, brca2, palb2, atm, chek2)

## AC0: Genotype quality too low, allele count is 0 after removing low quality genotypes
## AS VSQR: Failed VSQR filter
## Filter to retain only variants which pass
gnomAD_all %>%
  filter(FILTER=="PASS") %>%
  mutate(gene = case_when(grepl("17",`#CHROM`) ~ "BRCA1",
                          grepl("13",`#CHROM`) ~ "BRCA2",
                          grepl("16",`#CHROM`) ~ "PALB2",
                          grepl("11",`#CHROM`) ~ "ATM",
                          grepl("22",`#CHROM`) ~ "CHEK2"),
         `#CHROM` = gsub("chr","", `#CHROM`),
         ID = paste(`#CHROM`, POS, REF, ALT, sep="-"),
         `#CHROM` = as.double(`#CHROM`)) -> gnomAD_exome_all_genes

write_csv(gnomAD_exome_all_genes, "Data sources/gnomAD_v4.1/gnomAD_exome_all_genes.csv")

gnomAD_exome_all_genes %>%
  mutate(
    # extract AC and AN for female NFE non-UKB partition - for Ambry data prep only
    AC_NFE_v4 = str_extract(INFO, "(?<=AC_non_ukb_nfe=)[:digit:]+"),
    AN_NFE_v4 = str_extract(INFO, "(?<=AN_non_ukb_nfe=)[:digit:]+"),
    HOM_NFE_v4 = str_extract(INFO, "(?<=nhomalt_non_ukb_nfe=)[:digit:]+"),
    grpmax_faf_NFE_v4 = str_extract(INFO, "(?<=faf95_non_ukb_nfe=).+?(?=;)"),
    across(ends_with("v4"), ~ as.numeric(.x))) %>%
  
  # remove variants if not seen (AC_NFE_v4 = 0)
  filter(AC_NFE_v4 != 0) %>%
  
  # convert counts to individuals
  mutate(AC_NFE_v4 = AC_NFE_v4-HOM_NFE_v4,
         AN_NFE_v4 = AN_NFE_v4/2) -> gnomAD_prevep

write_csv(gnomAD_prevep, "Data sources/gnomAD_v4.1/gnomAD_exomes_prevep.csv")

# annotate gnomAD data with VEP annotations including HGVS nomenclature
read_tsv("Data sources/gnomAD_v4.1/vep_annotation_exomes.txt") -> VEP #VEP output file

VEP %>%
  filter(Feature %in% c("ENST00000675843.1","ENST00000357654.9","ENST00000404276.6","ENST00000261584.9","ENST00000380152.8")) %>%
  separate_wider_delim(HGVSc, delim=":", names=c("Transcript","hgvs")) %>%
  select(`#Uploaded_variation`, hgvs, Consequence, EXON) -> VEP

gnomAD_prevep %>%
  left_join(VEP,by=c("ID"="#Uploaded_variation")) -> gnomAD_postvep

write_csv(gnomAD_postvep, "Data sources/gnomAD_v4.1/gnomAD_exomes_postvep.csv")

# liftover gnomAD to build 37
gnomAD_postvep %>% select(1:7) %>% mutate(INFO = ".") -> gnomAD_preliftover
write_csv(gnomAD_preliftover, "Data sources/gnomAD_v4.1/gnomAD_exomes_preliftover.csv")

read_tsv("Data sources/gnomAD_v4.1/gnomAD_exomes_preliftover.vcf", col_names=FALSE) -> gnomAD_liftover # liftover output file

gnomAD_postvep %>%
  left_join(gnomAD_liftover, by=c("ID"="X3")) %>%
  rename("CHROM_BUILD37"="X1", "POS_BUILD37"="X2","REF_BUILD37"="X4","ALT_BUILD37"="X5","CHROM_BUILD38"="#CHROM","POS_BUILD38"="POS","ID_BUILD38"="ID","REF_BUILD38"="REF","ALT_BUILD38"="ALT") %>%
  select(-X6,-X7,-X8) %>%
  mutate(gene_ID = paste0(gene,":",hgvs)) -> gnomAD_annotated

write_csv(gnomAD_annotated, "Data sources/gnomAD_v4.1/gnomAD_exomes_annotated.csv")

# final variant validator check using build 37 coordinates (in case hgvs format from build 37 VEP is different; ensures all datasets use the same origin for hgvs nomenclature)
gnomAD_annotated %>% 
  mutate(ID_build37 = paste0("chr",CHROM_BUILD37,":",POS_BUILD37,":",REF_BUILD37,":",ALT_BUILD37)) -> gnomAD_annotated_prevv

write_csv(gnomAD_annotated_prevv, "Data sources/gnomAD_v4.1/gnomAD_annotated_prevv.csv")
read_tsv("Data sources/gnomAD_v4.1/gnomAD_annotated_postvv.txt",skip=2) -> gnomAD_annotated_postvv #variant validator output file

left_join(gnomAD_annotated_prevv, gnomAD_annotated_postvv, by=c("ID_build37"="Input")) -> gnomAD_annotated_postvv

# filter out variants which fail variantvalidator check
gnomAD_annotated_postvv %>% filter(!is.na(HGVS_transcript))-> gnomAD_annotated_postvv
write_csv(gnomAD_annotated_postvv,"Data sources/gnomAD_v4.1/gnomAD_clean.csv")


### Combine ----
read_csv("Data sources/Ambry/ambry_per_variant_VEP_VV_prelift.csv") -> ambry_summary_VEP_VV
read_tsv("Data sources/Ambry/ambry_per_variant_VEP_VV_postlift.txt", skip=2) -> ambry_vep_build37
read_csv("Data sources/Ambry/ambry_counts_fixed.csv") -> ambry_counts_fixed
read_csv("Data sources/gnomAD_v4.1/gnomAD_clean.csv") -> gnomAD_annotated_postvv

# excluding PTVs from Ambry which are excluded from the BRIDGES dataset (in preparation for later combination of datasets)
ambry_summary_VEP_VV %>%
  ## exclude PTVs in the last exon
  mutate(final_exon = case_when(Gene_Symbol=="BRCA1" & EXON=="23/23" ~ 1,
                                Gene_Symbol=="BRCA2" & EXON=="27/27" ~ 1,
                                Gene_Symbol=="ATM" & EXON=="63/63" ~ 1,
                                Gene_Symbol=="CHEK2" & EXON=="15/15" ~ 1,
                                Gene_Symbol=="PALB2" & EXON=="13/13" ~ 1,
                                .default=0)) %>%
  filter((final_exon==0 & grepl(paste(consequence_filter,collapse="|"),Consequence)) | grepl("missense",Consequence)) %>%
  
  # join with counts and VV annotations
  left_join(ambry_vep_build37, by=c("ID_build37"="Input")) %>%
  select(-ID,-gene_ID,-n:-Warnings.x) %>%
  distinct() %>%
  left_join(ambry_counts_fixed,by="ID_build37") %>%
  filter(!is.na(HGVS_transcript.y)) %>% 
  
  ## exclude splice variants affecting the penultimate exon for genes where exon skipping not confirmed pathogenic at those sites
  ### BRCA2: 9502-1,-2, 9648+1,+2. CHEK2: c.1462-1,-2, 1542+1,+2
  filter((gene=="NM_000059.3" & !grepl(paste("c.9502-1G","c.9502-2A","c.9648+1G","c.9648+2T",sep="|"),hgvs)) | gene %in% c("NM_007294.3","NM_024675.3","NM_000051.3","NM_007194.3")) %>%
  filter((gene=="NM_007194.3" & !grepl(paste("c.1462-1G","c.1462-2A","c.1542+1G","c.1542+2T",sep="|"),hgvs)) | gene %in% c("NM_007294.3","NM_024675.3","NM_000051.3","NM_000059.3")) %>%
  ## exclude 7 canonical BRCA1 VUS variants per ENIGMA: c.594-2A>C, c.4096+1G>A, c.4096+2T>C, c.4096+1G>A, c.4186-2A>G, c.4358-1G>C, c.4358-2del
  filter((gene=="NM_007294.3" & !grepl(paste("c.594-2A>C","c.4096+1G>A","c.4096+2T>C","c.4096+1G>A","c.4186-2A>G","c.4358-1G>C","c.4358-2del",sep="|"),hgvs)) | gene %in% c("NM_007194.3","NM_024675.3","NM_000051.3","NM_000059.3")
  ) %>%
  mutate(gene = case_when(gene=="NM_007294.3" ~ "BRCA1",
                          gene=="NM_000059.3" ~ "BRCA2",
                          gene=="NM_024675.3" ~ "PALB2",
                          gene=="NM_000051.3" ~ "ATM",
                          gene=="NM_007194.3" ~ "CHEK2")) -> ambry_summary_VEP

# excluding PTVs from gnomAD which are excluded from the BRIDGES dataset (in preparation for later combination of datasets)
gnomAD_annotated_postvv %>%
  select(-gene,-hgvs) %>%
  separate_wider_delim(HGVS_transcript, names=c("gene","hgvs"), delim=":") %>%
  mutate(gene = case_when(gene=="NM_007294.3" ~ "BRCA1",
                          gene=="NM_000059.3" ~ "BRCA2",
                          gene=="NM_024675.3" ~ "PALB2",
                          gene=="NM_000051.3" ~ "ATM",
                          gene=="NM_007194.3" ~ "CHEK2")) %>%
  ## exclude PTVs in the last exon
  mutate(final_exon = case_when(gene=="BRCA1" & EXON=="23/23" ~ 1,
                                gene=="BRCA2" & EXON=="27/27" ~ 1,
                                gene=="ATM" & EXON=="63/63" ~ 1,
                                gene=="CHEK2" & EXON=="15/15" ~ 1,
                                gene=="PALB2" & EXON=="13/13" ~ 1,
                                .default=0)) %>%
  filter((final_exon==0 & grepl(paste(consequence_filter,collapse="|"),Consequence)) | grepl("missense",Consequence)) %>%
  ## exclude splice variants affecting the penultimate exon for genes where exon skipping not confirmed pathogenic at those sites
  ### BRCA2: 9502-1,-2, 9648+1,+2. CHEK2: c.1462-1,-2, 1542+1,+2
  filter((gene=="BRCA2" & !grepl(paste("c.9502-1G","c.9502-2A","c.9648+1G","c.9648+2T",sep="|"),hgvs)) | gene %in% c("BRCA1","PALB2","ATM","CHEK2")) %>%
  filter((gene=="CHEK2" & !grepl(paste("c.1462-1G","c.1462-2A","c.1542+1G","c.1542+2T",sep="|"),hgvs)) | gene %in% c("BRCA1","PALB2","ATM","BRCA2")) %>%
  ## exclude 7 canonical BRCA1 VUS variants per ENIGMA: c.594-2A>C, c.4096+1G>A, c.4096+2T>C, c.4096+1G>A, c.4186-2A>G, c.4358-1G>C, c.4358-2del
  filter((gene=="BRCA1" & !grepl(paste("c.594-2A>C","c.4096+1G>A","c.4096+2T>C","c.4096+1G>A","c.4186-2A>G","c.4358-1G>C","c.4358-2del",sep="|"),hgvs)) | gene %in% c("CHEK2","PALB2","ATM","BRCA2")
  ) %>%
  
  # merge gnomAD counts with Ambry counts
  full_join(ambry_summary_VEP, by=c("gene","hgvs")) %>%
  # c,7397 variant is the same variant as in the UKB dataset which appears incorrectly in the variant validator nomenclature due to alternative allele presence. Amended here to match with UKB dataset.
  mutate(hgvs = str_replace_all(hgvs,"c.7397=","c.7397T>C")) -> ambry_gnomAD_merged

write_csv(ambry_gnomAD_merged, "Data sources/Ambry/ambry_ptvmiss_merged.csv")


# ANALYSIS ----
## Collate variant data ----
UKB <- read_csv("Data Sources/UKBiobank Data/WES/UKB_BRCA_Filtered_split_ptvmiss.csv")
CARRIERS <- read_csv("Data sources/CARRIERS/carriers_per_var_clean_ptvmiss.csv")
BRIDGES <- read_csv("Data sources/BRIDGES Data/bridges_ptvmiss.csv")
NDRS <- read_csv("Data Sources/NDRS_femalebreast_2025/ndrs_ptvmiss_whiteonly.csv")
AMBRY <- read_csv("Data sources/Ambry/ambry_ptvmiss_merged.csv")

NDRS %>%
  full_join(AMBRY, join_by("hgvs", "gene")) %>%
  full_join(UKB, join_by("hgvs", "gene"=="SYMBOL")) %>%
  full_join(BRIDGES, join_by("hgvs", "gene"=="Gene_Symbol")) %>%
  full_join(CARRIERS, join_by("hgvs", "gene")) -> joined

# remove very common variants
read_csv("Data sources/gnomAD_v4.1/gnomAD_BA1_filter_final.csv") -> gnomAD_BA1_list

joined %>%
  left_join(gnomAD_BA1_list, by=c("gene","hgvs"="hgvs_cdna")) %>%
  filter(is.na(rare_flag)) -> joined_rare 

write_csv(joined_rare, "Analysis/P2_enriched/ptvmiss_rare.csv")

joined_rare %>%
  # set empty values, add or estimate denominators for all variants
  select(gene,hgvs,
         NDRS_case_carrier=WhiteCount,NDRS_case_total=denominator,NDRS_control_carrier=White_nonBC_total_carriers_NDRSpart.x, NDRS_control_total=White_nonBC_total_noncarrier_NDRSpart.x,
         AMBRY_case_carrier=female_breast,AMBRY_case_total=breast_denom,AMBRY_control_carrier=AC_NFE_v4,AMBRY_control_total=AN_NFE_v4,
         UKB_case_carrier=White_BC_female_carriers,UKB_case_total=White_BC_female_total,UKB_control_carrier=White_nonBC_total_carriers_UKBpart,UKB_control_total=White_nonBC_total_noncarrier_UKBpart, # for partitioned controls - UKB split
         BRIDGES_case_carrier=Cases,BRIDGES_case_total=Case_denom.x,BRIDGES_control_carrier=Controls,BRIDGES_control_total=Control_denom.x,
         CARRIERS_case_carrier=Case,CARRIERS_case_total=Case_denom.y,CARRIERS_control_carrier=Control,CARRIERS_control_total=Control_denom.y) %>%
  mutate(across(3:22, ~ replace_na(.x, 0))) %>%
  # Replace 0 values in the denominators with average allele denominator for that dataset
  mutate(BRIDGES_control_total=44035,
         BRIDGES_case_total=42062,
         UKB_control_total = if_else(UKB_control_carrier==0 & UKB_control_total==0,122232,UKB_control_total), #for partitioned controls (UKB partition)
         UKB_case_total = if_else(UKB_case_carrier==0 & UKB_case_total==0,18477,UKB_case_total),
         NDRS_control_total = if_else(NDRS_control_carrier==0 & NDRS_control_total==0,297141,NDRS_control_total), #for partitioned controls (NDRS partition)
         NDRS_case_total = 44917,
         AMBRY_control_total = if_else(AMBRY_control_carrier==0 & AMBRY_control_total==0,175054,AMBRY_control_total),
         AMBRY_case_total = 187642,
         CARRIERS_case_total=32247,
         CARRIERS_control_total=32544) -> joined_names

write_csv(joined_names, "Analysis/P2_enriched/ptvmissvar_rare_cleaned.csv")

## Odds ratios and LRs ----
### Set column names for calculated values
ukbiobank <- c("p-value_ukbiobank","OR_ukbiobank","Lower_CI_ukbiobank","Upper_CI_ukbiobank","LR_path_ukbiobank","ACMG_points_ukbiobank")
carriers <- c("p-value_carriers","OR_carriers","Lower_CI_carriers","Upper_CI_carriers","LR_path_carriers","ACMG_points_carriers")
bridges <- c("p-value_bridges","OR_bridges","Lower_CI_bridges","Upper_CI_bridges","LR_path_bridges","ACMG_points_bridges")
ndrs <- c("p-value_ndrs","OR_ndrs","Lower_CI_ndrs","Upper_CI_ndrs","LR_path_ndrs","ACMG_points_ndrs")
ambry <- c("p-value_ambry","OR_ambry","Lower_CI_ambry","Upper_CI_ambry","LR_path_ambry","ACMG_points_ambry")

# run LR and OR calculations (function = batch_stats)
read_csv("Analysis/P2_enriched/ptvmissvar_rare_cleaned.csv") -> joined_names_clean

joined_names_clean %>%
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

write_csv(stats_target_4, "Analysis/X_UpdatedCalcChey_Apr26/ptvmissvar_LR_oddsAF_ormeta.csv")

## Parameter filtering ----
### add discrepancy thresholds
stats_target_4 %>%
  mutate(across(starts_with("Target_Odds_"), ~ round(.x,2))) %>%
  left_join(denom_limits, by=c("UKB_control_carrier"="obs_control","UKB_case_carrier"="obs_case","Target_Odds_ukbiobank"="target_odds")) %>% rename(control_case_ratio_ukbiobank=control_case_ratio, case_control_ratio_ukbiobank=case_control_ratio) %>%
  left_join(denom_limits, by=c("BRIDGES_control_carrier"="obs_control","BRIDGES_case_carrier"="obs_case","Target_Odds_bridges"="target_odds")) %>% rename(control_case_ratio_bridges=control_case_ratio, case_control_ratio_bridges=case_control_ratio) %>%
  left_join(denom_limits, by=c("CARRIERS_control_carrier"="obs_control","CARRIERS_case_carrier"="obs_case","Target_Odds_carriers"="target_odds")) %>% rename(control_case_ratio_carriers=control_case_ratio, case_control_ratio_carriers=case_control_ratio) %>%
  left_join(denom_limits, by=c("NDRS_control_carrier"="obs_control","NDRS_case_carrier"="obs_case","Target_Odds_ndrs"="target_odds")) %>% rename(control_case_ratio_ndrs=control_case_ratio, case_control_ratio_ndrs=case_control_ratio) %>%
  left_join(denom_limits, by=c("AMBRY_control_carrier"="obs_control","AMBRY_case_carrier"="obs_case","Target_Odds_ambry"="target_odds")) %>% rename(control_case_ratio_ambry=control_case_ratio, case_control_ratio_ambry=case_control_ratio) -> stats_target_4b

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

# single obs inclusion rule: include single obs if there is any data in another dataset (and that other dataset is NOT benign evidence from NDRS as we're discounting this)
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
  mutate(ACMG_combo=log(LR_combo,2.08)) %>%
  # remove variants which are not present in any dataset
  mutate(UKB_zero_from_part = if_else(zero_obs_ukbiobank==1 & zero_obs_bridges==1 & zero_obs_carriers==1 & zero_obs_ndrs==1 & zero_obs_ambry==1, 1, 0)) %>%
  filter(UKB_zero_from_part==0) %>%
  mutate(constituent_data = gsub("^$|^ $",NA,constituent_data)) -> case_counts_with_combo

plot_prep(case_counts_with_combo, "UK Biobank", "ACMG_points_ukbiobank")-> case_counts_with_combo
plot_prep(case_counts_with_combo, "BRIDGES", "ACMG_points_bridges")-> case_counts_with_combo
plot_prep(case_counts_with_combo, "CARRIERS", "ACMG_points_carriers")-> case_counts_with_combo
plot_prep(case_counts_with_combo, "NDRS", "ACMG_points_ndrs")-> case_counts_with_combo
plot_prep(case_counts_with_combo, "AMBRY", "ACMG_points_ambry")-> case_counts_with_combo
plot_prep(case_counts_with_combo, "All Datasets", "ACMG_combo")-> case_counts_with_combo

# output final dataset
write_csv(case_counts_with_combo, "Analysis/X_UpdatedCalcChey_Apr26//ptvmissvar_LR_oddsAF_ormeta_filtered.csv")


# Tables and Additional Analysis ----
## Table 3 ----

# set for either missense or PTVs
read_csv("Analysis/X_UpdatedCalcChey_Apr26/ptvmissvar_LR_oddsAF_orbase_filtered.csv") -> case_counts_with_combo
case_counts_with_combo %>% filter(Consequence=="Missense") -> case_counts_with_combo
case_counts_with_combo %>% filter(Consequence=="PTV") -> case_counts_with_combo

# total number of variants in each dataset (pre-filtering)
case_counts_with_combo %>%
  # filter(gene=="BRCA2") %>%
  mutate(`UK Biobank` = if_else(UKB_case_carrier == 0 & UKB_control_carrier == 0,"zero",`UK Biobank`)) %>%
  mutate(BRIDGES = if_else(BRIDGES_case_carrier==0 & BRIDGES_control_carrier==0,"zero",BRIDGES)) %>%
  mutate(CARRIERS = if_else(CARRIERS_case_carrier==0 & CARRIERS_control_carrier==0,"zero",CARRIERS)) %>%
  mutate(NDRS = if_else(NDRS_case_carrier==0 & NDRS_control_carrier==0,"zero",NDRS)) %>%
  mutate(AMBRY = if_else(AMBRY_case_carrier==0 & AMBRY_control_carrier==0,"zero",AMBRY)) %>%
  
  select(`UK Biobank`, BRIDGES, CARRIERS, NDRS, AMBRY, `All Datasets`) %>%
  mutate(across(everything(), ~ replace_na(.x, "None"))) %>%
  pivot_longer(c(`UK Biobank`, BRIDGES, CARRIERS, NDRS, AMBRY, `All Datasets`)) %>%
  group_by(name) %>%
  filter(value!="zero") %>%
  summarise(n=n())

order <- c("P_VSTR","P_STR","P_MOD","P_SUP","None","B_SUP","B_MOD","B_STR","B_VSTR")
data <- c("UK Biobank","BRIDGES","CARRIERS","NDRS","AMBRY","All Datasets")

case_counts_with_combo %>%
  filter(gene=="BRCA1") %>%
  # remove variants which met criteria for likely reduced penetrance (remove data in all datasets, not just datasets which meet this criteria)
  mutate(`All Datasets` = if_else((moderate_ukbiobank==1|moderate_bridges==1|moderate_carriers==1|moderate_ndrs==1|moderate_ambry==1), "pen", `All Datasets`)) %>%
  filter(`All Datasets` == "pen") -> case_counts_mod_penetrance

case_counts_with_combo %>%
  #filter(gene=="CHEK2") %>%
  # remove variants which met criteria for likely reduced penetrance (remove data in all datasets, not just datasets which meet this criteria)
  mutate(`All Datasets` = if_else((moderate_ukbiobank==1|moderate_bridges==1|moderate_carriers==1|moderate_ndrs==1|moderate_ambry==1), "pen", `All Datasets`)) %>%
  filter(`All Datasets` != "pen" | is.na(`All Datasets`)) %>%
  # remove data which failed parameter checks for specific datasets
  ## note: use include_DATA rather than include_DATA_combo field, as these cols consider only that dataset in a vaccuum (and the _combo fields account for obs in the other datasets)
  mutate(`UK Biobank` = if_else(include_ukbiobank == "N","RM",`UK Biobank`)) %>%
  mutate(BRIDGES = if_else(include_bridges == "N","RM",BRIDGES)) %>%
  mutate(CARRIERS = if_else(include_carriers == "N","RM",CARRIERS)) %>%
  mutate(NDRS = if_else((zero_obs_ndrs+single_obs_ndrs+discrepant_evidence_ndrs+moderate_ndrs>0),"RM",NDRS)) %>% #to make sure we count benign evidence from NDRS set, even if that isn't ultimately used
  mutate(AMBRY = if_else(include_ambry == "N","RM",AMBRY)) %>%
  # this column uses constituent data field as this accounts for all data including where there is a single obs in one dataset and another obs in another dataset.
  mutate(`All Datasets` = if_else(is.na(constituent_data), "RM",`All Datasets`)) %>%
  # remove data which is present for variants that have not been seen in specific datasets
  mutate(`UK Biobank` = if_else(UKB_case_carrier == 0 & UKB_control_carrier == 0,"zero",`UK Biobank`)) %>%
  mutate(BRIDGES = if_else(BRIDGES_case_carrier==0 & BRIDGES_control_carrier==0,"zero",BRIDGES)) %>%
  mutate(CARRIERS = if_else(CARRIERS_case_carrier==0 & CARRIERS_control_carrier==0,"zero",CARRIERS)) %>%
  mutate(NDRS = if_else(NDRS_case_carrier==0 & NDRS_control_carrier==0,"zero",NDRS)) %>%
  mutate(AMBRY = if_else(AMBRY_case_carrier==0 & AMBRY_control_carrier==0,"zero",AMBRY)) -> case_summary_table

case_summary_table %>%
  select(gene, hgvs,`UK Biobank`, BRIDGES, CARRIERS, NDRS, AMBRY, `All Datasets`) %>%
  mutate(across(everything(), ~ replace_na(.x, "None"))) %>%
  pivot_longer(c(`UK Biobank`, BRIDGES, CARRIERS, NDRS, AMBRY, `All Datasets`)) %>%
  group_by(gene, name, value) %>%
  filter(value != "RM" & value != "zero") %>% 
  summarise(n=n()) -> preplot

# get breakdown of combined evidence given
preplot %>% filter(name=='All Datasets')

# sums variants with evidence towards/against pathogenicity, can be filtered by gene in above code block 
preplot %>%
  mutate(value = case_when(grepl("P_",value) ~ "P_",
                           grepl("B_",value) ~ "B_",
                           grepl("None",value) ~ "None")) %>%
  group_by(gene,name,value) %>%
  summarise(n=sum(n)) %>%
  print(n=100)

## Table S3 and S4 ----
# sums variants removed due to singlet observation
case_counts_with_combo %>%
  #filter(gene=="PALB2") %>%
  mutate(`UK Biobank` = if_else(single_obs_ukbiobank==1,"single",`UK Biobank`)) %>%
  mutate(BRIDGES = if_else(single_obs_bridges==1,"single",BRIDGES)) %>%
  mutate(CARRIERS = if_else(single_obs_carriers==1,"single",CARRIERS)) %>%
  mutate(NDRS = if_else(single_obs_ndrs==1 & ACMG_points_ndrs>0,"single",NDRS)) %>%
  mutate(AMBRY = if_else(single_obs_ambry==1,"single",AMBRY)) %>%
  mutate(`All Datasets` = if_else((grepl("single",`UK Biobank`)|grepl("single",BRIDGES)|grepl("single",CARRIERS)|(grepl("single",NDRS)&ACMG_points_ndrs>0)|grepl("single",AMBRY)) & is.na(constituent_data), "single", `All Datasets`)) %>%
  select(gene, `UK Biobank`, BRIDGES, CARRIERS, NDRS, AMBRY, `All Datasets`) %>%
  mutate(across(everything(), ~ replace_na(.x, "None"))) %>%
  pivot_longer(c(`UK Biobank`, BRIDGES, CARRIERS, NDRS, AMBRY, `All Datasets`)) %>%
  group_by(gene, name) %>%
  filter(value=="single") %>%
  summarise(n=n()) %>% print(n=100)

# sums variants getting benign evidence in NDRS (ACMG points < 0)
case_counts_with_combo %>%
  mutate(NDRS = if_else(ACMG_points_ndrs<0,"benignexcl",NDRS)) %>%
  mutate(`All Datasets` = if_else(ACMG_points_ndrs<0 & (zero_obs_ukbiobank==1&zero_obs_bridges==1&zero_obs_carriers==1&zero_obs_ambry==1), "benignexcl", `All Datasets`)) %>%
  select(gene, NDRS, `All Datasets`) %>%
  mutate(across(everything(), ~ replace_na(.x, "None"))) %>%
  pivot_longer(c(NDRS, `All Datasets`)) %>%
  group_by(gene, name) %>%
  filter(value=="benignexcl") %>%
  summarise(n=n()) %>% print(n=100)

# sums variants removed due to discrepancy
case_counts_with_combo %>%
  # filter(gene=="BRCA2") %>%
  mutate(`UK Biobank` = if_else(discrepant_evidence_ukbiobank==1,"dis",`UK Biobank`)) %>%
  mutate(BRIDGES = if_else(discrepant_evidence_bridges==1 & single_obs_bridges==0,"dis",BRIDGES)) %>%
  mutate(CARRIERS = if_else(discrepant_evidence_carriers==1 & single_obs_carriers==0,"dis",CARRIERS)) %>%
  mutate(NDRS = if_else(discrepant_evidence_ndrs==1 & single_obs_ndrs==0,"dis",NDRS)) %>%
  mutate(AMBRY = if_else(discrepant_evidence_ambry==1 & single_obs_ambry==0,"dis",AMBRY)) %>%
  mutate(`All Datasets` = if_else((grepl("dis",`UK Biobank`)|grepl("dis",BRIDGES)|grepl("dis",CARRIERS)|grepl("dis",NDRS)|grepl("dis",AMBRY))
                                  & single_obs_ukbiobank==0 & single_obs_bridges==0 & single_obs_carriers==0 & single_obs_ndrs==0 & single_obs_ambry==0
                                  & is.na(constituent_data), "dis", `All Datasets`)) %>%
  select(gene,UKB_case_carrier,UKB_control_carrier ,  `UK Biobank`, BRIDGES, CARRIERS, NDRS, AMBRY, `All Datasets`) %>%
  mutate(across(everything(), ~ replace_na(.x, "None"))) %>%
  pivot_longer(c(`UK Biobank`, BRIDGES, CARRIERS, NDRS, AMBRY, `All Datasets`)) %>%
  group_by(gene, UKB_case_carrier,UKB_control_carrier,name) %>%
  filter(value=="dis") %>%
  summarise(n=n())

# sums variants removed due to suspected moderate penetrance
case_counts_with_combo %>%
  #filter(gene=="BRCA2") %>%
  mutate(`UK Biobank` = if_else(moderate_ukbiobank==1,"pen",`UK Biobank`)) %>%
  mutate(BRIDGES = if_else(moderate_bridges==1,"pen",BRIDGES)) %>%
  mutate(CARRIERS = if_else(moderate_carriers==1,"pen",CARRIERS)) %>%
  mutate(NDRS = if_else(moderate_ndrs==1,"pen",NDRS)) %>%
  mutate(AMBRY = if_else(moderate_ambry==1,"pen",AMBRY)) %>%
  mutate(`All Datasets` = if_else((grepl("pen",`UK Biobank`)|grepl("pen",BRIDGES)|grepl("pen",CARRIERS)|grepl("pen",NDRS)|grepl("pen",AMBRY)), "pen", `All Datasets`)) %>%
  select(gene, `UK Biobank`, BRIDGES, CARRIERS, NDRS, AMBRY, `All Datasets`) %>%
  mutate(across(everything(), ~ replace_na(.x, "None"))) %>%
  pivot_longer(c(`UK Biobank`, BRIDGES, CARRIERS, NDRS, AMBRY, `All Datasets`)) %>%
  group_by(gene,name) %>%
  filter(value=="pen") %>%
  summarise(n=n())

## Figure 2 ----
# Boxplot breakdown of only variants in the combined dataset
case_summary_table %>%
  select(gene, hgvs, ACMG_combo,`UK Biobank`, BRIDGES, CARRIERS, NDRS, AMBRY, `All Datasets`) %>%
  mutate(across(everything(), ~ replace_na(.x, "None"))) %>%
  pivot_longer(c(`UK Biobank`, BRIDGES, CARRIERS, NDRS, AMBRY, `All Datasets`)) %>%
  filter(name=="All Datasets" & value != "RM" & value != "zero") %>%
  group_by(gene) %>%
  mutate(avg=mean(ACMG_combo)) %>%
  mutate(med=median(ACMG_combo)) %>%
  ungroup() -> boxplot_casecounts

# statistical tests for significantly different groups (per-gene comparison)
kruskal_test(ACMG_combo ~ gene, data=boxplot_casecounts)
dunn_test(ACMG_combo ~ gene, data=boxplot_casecounts, p.adjust.method="bonferroni") -> dunn
dunn %>% add_xy_position(x="gene") -> dunndunn
# added for prettier plotting
dunndunn %>% mutate(y.position=c(176,195,291,214,310,233,329,253,348,272)) -> dunndunn #if using missense variants
dunndunn %>% mutate(y.position=c(176,195,253,272,291,214,233,310,329,348)) -> dunndunn #uf using ptvs

# boxplot with statistical comparisons plotted
boxplot_casecounts %>%
  mutate(ACMG_combo = if_else(ACMG_combo <= -160, -160, ACMG_combo)) %>%
  mutate(ACMG_combo = if_else(ACMG_combo >= 160, 160, ACMG_combo)) %>%
  ggplot(aes(x=gene, y=ACMG_combo)) +
  geom_boxplot(width=0.5) +
  coord_flip() +
  stat_pvalue_manual(dunndunn, label="p {scales::pvalue(p.adj)}", bracket.size=0.05, bracket.nudge.y=0, tip.length=0.01, hide.ns=TRUE, coord.flip=TRUE) +
  geom_hline(yintercept=0, linetype="dashed", color="red", size=1) +
  geom_text(aes(y=avg, label=paste0("median:", round(med,2)),vjust=4,size=20),check_overlap=TRUE) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        text = element_text(size=16),
        axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "none") +
  scale_y_continuous(breaks=seq(-1060,160,10), limits=c(-160,320), expand = c(0.01, 0)) +
  ggtitle("PS4-LLR applied for variants in each gene") +
  ylab("PS4-LLR") +
  xlab("Gene")

## Figure S3 ----
## BRCA1 Sankey
read_csv("Analysis/Reduced Penetrance/ReducedPenetrance_forSankey_missense_new.csv") -> vcep_redpen
vcep_redpen$x <- factor(vcep_redpen$x, levels=c(2,4))
vcep_redpen[order(vcep_redpen$x),] -> vcep_redpen
vcep_redpen$next_x <- factor(vcep_redpen$next_x, levels=c(2,4))
vcep_redpen[order(vcep_redpen$next_x),] -> vcep_redpen
vcep_redpen$node <- factor(vcep_redpen$node, levels=c("B_VSTR","B_STR","B_MOD","B_SUP","None","P_SUP","P_MOD","P_STR","P_VSTR"))
vcep_redpen[order(vcep_redpen$node),] -> vcep_redpen
vcep_redpen$next_node <- factor(vcep_redpen$next_node, levels=c("B_VSTR","B_STR","B_MOD","B_SUP","None","P_SUP","P_MOD","P_STR","P_VSTR"))
vcep_redpen[order(vcep_redpen$next_node),] -> vcep_redpen
vcep_redpen %>% filter(gene=="BRCA1") %>% select(-1,-2) -> df

ggplot(df, aes(x=x, next_x=next_x, node=node, next_node=next_node, fill=node, label=node)) +
  geom_sankey(flow.alpha = 0.75, node.color = 1) +
  geom_sankey_label() +
  theme_sankey(base_size=22) +
  theme(text=element_text(size=20),
        plot.title = element_text(hjust = 0.5)) +
  geom_sankey_label(size = 4, color=c("white","white","black","black","black","black","black","white","white",
                                      "white","white","black","black","black","black","black","white","white")) +
  xlab("Target odds of association") +
  scale_fill_manual(values = colour_palette) +
  ggtitle("BRCA1") +
  theme(legend.position = "none") -> BRCA1

## BRCA2 sankey
read_csv("Analysis/Reduced Penetrance/ReducedPenetrance_forSankey_missense_new.csv") -> vcep_redpen
vcep_redpen$x <- factor(vcep_redpen$x, levels=c(2,4))
vcep_redpen[order(vcep_redpen$x),] -> vcep_redpen
vcep_redpen$next_x <- factor(vcep_redpen$next_x, levels=c(2,4))
vcep_redpen[order(vcep_redpen$next_x),] -> vcep_redpen
vcep_redpen$node <- factor(vcep_redpen$node, levels=c("B_VSTR","B_STR","B_MOD","B_SUP","None","P_SUP","P_MOD","P_STR","P_VSTR"))
vcep_redpen[order(vcep_redpen$node),] -> vcep_redpen
vcep_redpen$next_node <- factor(vcep_redpen$next_node, levels=c("B_VSTR","B_STR","B_MOD","B_SUP","None","P_SUP","P_MOD","P_STR","P_VSTR"))
vcep_redpen[order(vcep_redpen$next_node),] -> vcep_redpen
vcep_redpen %>% filter(gene=="BRCA2") %>% select(-1,-2) -> df

ggplot(df, aes(x=x, next_x=next_x, node=node, next_node=next_node, fill=node, label=node)) +
  geom_sankey(flow.alpha = 0.75, node.color = 1) +
  geom_sankey_label() +
  theme_sankey(base_size=22) +
  theme(text=element_text(size=20),
        plot.title = element_text(hjust = 0.5)) +
  geom_sankey_label(size = 4, color=c("white","white","black","black","black","black","black","white","white",
                                      "white","white","black","black","black","black","black","white","white")) +
  xlab("Target odds of association") +
  scale_fill_manual(values = colour_palette) +
  ggtitle("BRCA2") +
  theme(legend.position = "none") -> BRCA2

## PALB2 sankey
read_csv("Analysis/Reduced Penetrance/ReducedPenetrance_forSankey_missense_new.csv") -> vcep_redpen
vcep_redpen$x <- factor(vcep_redpen$x, levels=c(2,4))
vcep_redpen[order(vcep_redpen$x),] -> vcep_redpen
vcep_redpen$next_x <- factor(vcep_redpen$next_x, levels=c(2,4))
vcep_redpen[order(vcep_redpen$next_x),] -> vcep_redpen
vcep_redpen$node <- factor(vcep_redpen$node, levels=c("B_VSTR","B_STR","B_MOD","B_SUP","None","P_SUP","P_MOD","P_STR","P_VSTR"))
vcep_redpen[order(vcep_redpen$node),] -> vcep_redpen
vcep_redpen$next_node <- factor(vcep_redpen$next_node, levels=c("B_VSTR","B_STR","B_MOD","B_SUP","None","P_SUP","P_MOD","P_STR","P_VSTR"))
vcep_redpen[order(vcep_redpen$next_node),] -> vcep_redpen
vcep_redpen %>% filter(gene=="PALB2") %>% select(-1,-2) -> df

ggplot(df, aes(x=x, next_x=next_x, node=node, next_node=next_node, fill=node, label=node)) +
  geom_sankey(flow.alpha = 0.75, node.color = 1) +
  geom_sankey_label() +
  theme_sankey(base_size=22) +
  theme(text=element_text(size=20),
        plot.title = element_text(hjust = 0.5)) +
  geom_sankey_label(size = 4, color=c("white","white","black","black","black","black","black","white","white",
                                      "white","white","black","black","black","black","black","white","white")) +
  xlab("Target odds of association") +
  scale_fill_manual(values = colour_palette) +
  ggtitle("PALB2") +
  theme(legend.position = "none") -> PALB2

## ATM sankey
read_csv("Analysis/Reduced Penetrance/ReducedPenetrance_forSankey_missense_new.csv") -> vcep_redpen
vcep_redpen$x <- factor(vcep_redpen$x, levels=c(2,4))
vcep_redpen[order(vcep_redpen$x),] -> vcep_redpen
vcep_redpen$next_x <- factor(vcep_redpen$next_x, levels=c(2,4))
vcep_redpen[order(vcep_redpen$next_x),] -> vcep_redpen
vcep_redpen$node <- factor(vcep_redpen$node, levels=c("B_VSTR","B_STR","B_MOD","B_SUP","None","P_SUP","P_MOD","P_STR","P_VSTR"))
vcep_redpen[order(vcep_redpen$node),] -> vcep_redpen
vcep_redpen$next_node <- factor(vcep_redpen$next_node, levels=c("B_VSTR","B_STR","B_MOD","B_SUP","None","P_SUP","P_MOD","P_STR","P_VSTR"))
vcep_redpen[order(vcep_redpen$next_node),] -> vcep_redpen
vcep_redpen %>% filter(gene=="ATM") %>% filter(!is.na(x)) %>% select(-1,-2) -> df

ggplot(df, aes(x=x, next_x=next_x, node=node, next_node=next_node, fill=node, label=node)) +
  geom_sankey(flow.alpha = 0.75, node.color = 1) +
  geom_sankey_label() +
  theme_sankey(base_size=22) +
  theme(text=element_text(size=20),
        plot.title = element_text(hjust = 0.5)) +
  geom_sankey_label(size = 4, color=c("white","white","black","black","black","black","black","white","white",
                                      "white","white","black","black","black","black","black","white","white")) +
  xlab("Target odds of association") +
  scale_fill_manual(values = colour_palette) +
  theme(legend.position = "none") +
  ggtitle("ATM") -> ATM

## CHEK2 sankey
read_csv("Analysis/Reduced Penetrance/ReducedPenetrance_forSankey_missense_new.csv") -> vcep_redpen
vcep_redpen$x <- factor(vcep_redpen$x, levels=c(2,4))
vcep_redpen[order(vcep_redpen$x),] -> vcep_redpen
vcep_redpen$next_x <- factor(vcep_redpen$next_x, levels=c(2,4))
vcep_redpen[order(vcep_redpen$next_x),] -> vcep_redpen
vcep_redpen$node <- factor(vcep_redpen$node, levels=c("B_VSTR","B_STR","B_MOD","B_SUP","None","P_SUP","P_MOD","P_STR","P_VSTR"))
vcep_redpen[order(vcep_redpen$node),] -> vcep_redpen
vcep_redpen$next_node <- factor(vcep_redpen$next_node, levels=c("B_VSTR","B_STR","B_MOD","B_SUP","None","P_SUP","P_MOD","P_STR","P_VSTR"))
vcep_redpen[order(vcep_redpen$next_node),] -> vcep_redpen
vcep_redpen %>% filter(gene=="CHEK2") %>% filter(!is.na(x)) %>% select(-1,-2) -> df

ggplot(df, aes(x=x, next_x=next_x, node=node, next_node=next_node, fill=node, label=node)) +
  geom_sankey(flow.alpha = 0.75, node.color = 1) +
  geom_sankey_label() +
  theme_sankey(base_size=22) +
  theme(text=element_text(size=20),
        plot.title = element_text(hjust = 0.5)) +
  geom_sankey_label(size = 4, color=c("white","white","black","black","black","black","black","white","white",
                                      "white","white","black","black","black","black","black","white","white")) + 
  xlab("Target odds of association") +
  scale_fill_manual(values = colour_palette) +
  theme(legend.position = "none") +
  ggtitle("CHEK2") -> CHEK2

ggarrange(BRCA1,BRCA2,PALB2,ATM,CHEK2,labels=c("A","B","C","D","E"),ncol=3,nrow=2)


## Figure 3 ----
### ClinVar validation plotting 
# Bar plot using ClinVar B/LB and P/LP variants to show concordance, all ClinVar classified 
read_csv("Analysis/Validation_datasets/clinvar_geneset_total.csv") -> clinvar_allstar
## VEP consequence for ClinVar variants
read_tsv("Analysis/Validation_datasets/clinvar_vep.txt") -> clinvar_vep
clinvar_vep %>% filter(Feature %in% c("ENST00000278616.4","ENST00000357654.3","ENST00000380152.3","ENST00000328354.6","ENST00000261584.4")) -> clinvar_vep

clinvar_allstar %>%
  mutate(ID = paste(transcript,hgvs_cdna,sep=":")) %>%
  left_join(clinvar_vep,by=c("ID"="#Uploaded_variation")) %>%
  select(gene,hgvs_cdna,hgvs_prot,Type,ClinicalSignificance,ReviewStatus,ClinVar_Stars,Consequence) %>%
  mutate(ClinicalSignificance = gsub("Pathogenic/Likely pathogenic","Likely pathogenic/Pathogenic",ClinicalSignificance)) -> clinvar_allstar
clinvar_allstar %>% 
  filter(grepl("missense",Consequence)) %>%
  filter(ClinVar_Stars > 0) %>%
  filter(!grepl("#N/A",ClinicalSignificance) & !grepl("Uncertain",ClinicalSignificance) & !grepl("not provided",ClinicalSignificance) & !grepl("Conflict",ClinicalSignificance)) %>% 
  write_csv("Analysis/P2_enriched/clinvar_missense_classified.csv")

clinvar_allstar %>% 
  filter(grepl("missense",Consequence)) %>%
  filter(ClinVar_Stars>0) %>%
  filter(grepl("#N/A",ClinicalSignificance) | grepl("Uncertain",ClinicalSignificance) | grepl("Conflict",ClinicalSignificance))

clinvar_allstar %>%
  filter(grepl("missense",Consequence)) %>%
  inner_join(case_counts_with_combo,by=c("gene", "hgvs_cdna"="hgvs"))

# all clinvar missense variants merged with PS4-LLR, even if clinvar variant is not in the PS4_LLR combined dataset
clinvar_allstar %>% 
  filter(grepl("missense",Consequence)) %>%
  filter(ClinVar_Stars > 0) %>%
  filter(!grepl("#N/A",ClinicalSignificance) & !grepl("Uncertain",ClinicalSignificance) & !grepl("not provided",ClinicalSignificance) & !grepl("Conflict",ClinicalSignificance)) %>%
  left_join(case_counts_with_combo, by=c("gene","hgvs_cdna"="hgvs")) %>%
  group_by(gene, ClinicalSignificance,`All Datasets`) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  pivot_wider(names_from=`All Datasets`,values_from=n) %>%
  select(gene,ClinicalSignificance,P_VSTR,P_STR,P_MOD,P_SUP,`NA`,B_SUP,B_MOD,B_STR,B_VSTR) %>%
  print(n=100)
clinvar_allstar %>% 
  filter(grepl("missense",Consequence)) %>%
  filter(ClinVar_Stars>0) %>%
  filter(grepl("#N/A",ClinicalSignificance) | grepl("Uncertain",ClinicalSignificance) | grepl("Conflict",ClinicalSignificance)) %>%
  left_join(case_counts_with_combo, by=c("gene","hgvs_cdna"="hgvs")) %>%
  group_by(gene, ClinicalSignificance,`All Datasets`) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  pivot_wider(names_from=`All Datasets`,values_from=n) %>%
  select(gene,ClinicalSignificance,P_VSTR,P_STR,P_MOD,P_SUP,`NA`,B_SUP,B_MOD,B_STR,B_VSTR) %>%
  print(n=100)

inner_join(case_counts_with_combo, clinvar_allstar, by=c("gene", "hgvs"="hgvs_cdna")) -> vcep_with_clinvar

write_csv(vcep_with_clinvar, "Analysis/P2_enriched/missense_LR_oddsAF_palb2hi_filtered_withClinVar.csv")
read_csv("Analysis/P2_enriched/missense_LR_oddsAF_palb2hi_filtered_withClinVar.csv") -> vcep_with_clinvar

# number of clinvar variants, any classification, which have data in the combined dataset
vcep_with_clinvar %>%
  # select(1,2,120:127) %>%
  filter(ClinVar_Stars > 0) %>%
  filter(ACMG_combo != 0)

# number of VUS/conflicting clinvar variants which have data in combined dataset (ie are not filtered out)
vcep_with_clinvar %>%
  # select(1,2,120:127) %>%
  filter(grepl("#N/A",ClinicalSignificance) | grepl("Uncertain",ClinicalSignificance) | grepl("Conflict",ClinicalSignificance)) %>%
  filter(ClinVar_Stars > 0) %>%
  filter(ACMG_combo != 0)

# number of 1* or above classified clinvar variants (LP/P/LB/B) which have data in the combined dataset
vcep_with_clinvar %>%
  # select(1,2,120:127) %>%
  #filter(grepl("missense",Consequence)) %>%
  filter(ClinVar_Stars > 0) %>%
  filter(!grepl("#N/A",ClinicalSignificance) & !grepl("Uncertain",ClinicalSignificance) & !grepl("not provided",ClinicalSignificance) & !grepl("Conflict",ClinicalSignificance)) %>%
  filter(Include == "Yes") -> vcep_with_clinvar
vcep_with_clinvar %>% group_by(ClinicalSignificance) %>% summarise(n=n()) %>% print(n=100)
vcep_with_clinvar %>% group_by(ClinicalSignificance,`All Datasets`) %>% summarise(n=n()) %>% print(n=100)

vcep_with_clinvar %>% filter(grepl("B_",`All Datasets`) & grepl("athogeni",ClinicalSignificance)) -> vcep_with_clinvar_PLP_benPS4LLR
vcep_with_clinvar %>% filter(grepl("P_",`All Datasets`) & grepl("enign",ClinicalSignificance)) -> vcep_with_clinvar_BLB_pathPS4LLR

vcep_with_clinvar %>%
  group_by(gene, ClinicalSignificance,`All Datasets`) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  pivot_wider(names_from=`All Datasets`,values_from=n) %>%
  select(gene,ClinicalSignificance,P_VSTR,P_STR,P_MOD,P_SUP,`NA`,B_SUP,B_MOD,B_STR,B_VSTR) %>%
  print(n=100)

vcep_with_clinvar %>%
  select(ACMG_combo, ClinicalSignificance) %>%
  ggplot(aes(x=ACMG_combo)) +
  geom_histogram(aes(fill=ClinicalSignificance),binwidth=5, color="black", center=0.5, position="stack") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        text = element_text(size=14)) +
  scale_fill_manual(values=c("#130857","#444cc7","#2f9fd4","#cc6600","#B2182B","#5e0101")) + 
  scale_x_continuous(breaks=seq(-200,200,10), limits=c(-200, 200), expand = c(0.01, 0)) +
  scale_y_continuous(breaks=seq(0,50,5), limits=c(0,50), expand=c(0.01,0)) +
  xlab("PS4-LLR Evidence") +
  ylab("No. variants") +
  labs(fill="Clinical Significance") +
  ggtitle("PS4-LLR (Log Likelihood Ratio) for variants classified in ClinVar")

## Figure S2 ----
# plot of variants with at least 2* classifications
vcep_with_clinvar %>% 
  filter(ClinVar_Stars >= 2) %>% 
  filter(!grepl("#N/A",ClinicalSignificance) & !grepl("Uncertain",ClinicalSignificance) & !grepl("not provided",ClinicalSignificance) & !grepl("Conflict",ClinicalSignificance)) %>%
  filter(ACMG_combo != 0) -> vcep_with_clinvar

vcep_with_clinvar %>% group_by(gene) %>% summarise(n=n()) %>% print(n=100)
vcep_with_clinvar %>% group_by(ClinicalSignificance,`All Datasets`) %>% summarise(n=n()) %>% print(n=100)
vcep_with_clinvar %>%
  filter(gene %in% c("BRCA1","BRCA2")) %>%
  group_by(ClinicalSignificance,`All Datasets`) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  pivot_wider(names_from=`All Datasets`,values_from=n) %>%
  select(ClinicalSignificance,P_VSTR,P_STR,P_MOD,P_SUP,`NA`,B_SUP,B_MOD,B_STR,B_VSTR) %>%
  print(n=100)

vcep_with_clinvar %>%
  filter(gene=="BRCA1") %>%
  select(ACMG_combo, ClinicalSignificance) %>%
  ggplot(aes(x=ACMG_combo)) +
  geom_histogram(aes(fill=ClinicalSignificance),binwidth=15, color="black", center=0.5, position="stack") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        legend.position = "none",
        text = element_text(size=12)) +
  scale_fill_manual(values=c("#130857","#444cc7","#2f9fd4","#cc6600","#B2182B","#5e0101")) + 
  scale_x_continuous(breaks=seq(-1100,200,50), limits=c(-1100, 200), expand = c(0.01, 0)) +
  scale_y_continuous(breaks=seq(0,50,5), limits=c(0,50), expand=c(0.01,0)) +
  xlab("Log Likelihood Ratio (LLR)") +
  ylab("No. variants") +
  ggtitle("Case-control LLR for variants classified in ClinVar \n(BRCA1)") -> brca1
vcep_with_clinvar %>%
  filter(gene=="BRCA2") %>%
  select(ACMG_combo, ClinicalSignificance) %>%
  ggplot(aes(x=ACMG_combo)) +
  geom_histogram(aes(fill=ClinicalSignificance),binwidth=15, color="black", center=0.5, position="stack") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        legend.position = "none",
        text = element_text(size=12)) +
  scale_fill_manual(values=c("#130857","#444cc7","#2f9fd4","#cc6600","#B2182B","#5e0101")) + 
  scale_x_continuous(breaks=seq(-1100,200,50), limits=c(-1100, 200), expand = c(0.01, 0)) +
  scale_y_continuous(breaks=seq(0,50,5), limits=c(0,50), expand=c(0.01,0)) +
  xlab("Log Likelihood Ratio (LLR)") +
  ylab("No. variants") +
  ggtitle("Case-control LLR for variants classified in ClinVar \n(BRCA2)") -> brca2
vcep_with_clinvar %>%
  filter(gene=="PALB2") %>%
  select(ACMG_combo, ClinicalSignificance) %>%
  ggplot(aes(x=ACMG_combo)) +
  geom_histogram(aes(fill=ClinicalSignificance),binwidth=15, color="black", center=0.5, position="stack") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        legend.position = "none",
        text = element_text(size=12)) +
  scale_fill_manual(values=c("#130857","#5e0101")) + 
  scale_x_continuous(breaks=seq(-1100,200,50), limits=c(-1100, 200), expand = c(0.01, 0)) +
  scale_y_continuous(breaks=seq(0,50,5), limits=c(0,50), expand=c(0.01,0)) +
  xlab("Log Likelihood Ratio (LLR)") +
  ylab("No. variants") +
  ggtitle("Case-control LLR for variants classified in ClinVar \n(PALB2)") -> palb2
vcep_with_clinvar %>%
  filter(gene=="ATM") %>%
  select(ACMG_combo, ClinicalSignificance) %>%
  ggplot(aes(x=ACMG_combo)) +
  geom_histogram(aes(fill=ClinicalSignificance),binwidth=15, color="black", center=0.5, position="stack") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        legend.position = "none",
        text = element_text(size=12)) +
  scale_fill_manual(values=c("#130857","#444cc7","#cc6600","#B2182B","#5e0101")) + 
  scale_x_continuous(breaks=seq(-1100,200,50), limits=c(-1100, 200), expand = c(0.01, 0)) +
  scale_y_continuous(breaks=seq(0,50,5), limits=c(0,50), expand=c(0.01,0)) +
  xlab("Log Likelihood Ratio (LLR)") +
  ylab("No. variants") +
  ggtitle("Case-control LLR for variants classified in ClinVar \n(ATM)") -> atm
vcep_with_clinvar %>%
  filter(gene=="CHEK2") %>%
  select(ACMG_combo, ClinicalSignificance) %>%
  ggplot(aes(x=ACMG_combo)) +
  geom_histogram(aes(fill=ClinicalSignificance),binwidth=15, color="black", center=0.5, position="stack") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        legend.position = "none",
        text = element_text(size=12)) +
  scale_fill_manual(values=c("#cc6600","#B2182B")) + 
  scale_x_continuous(breaks=seq(-1100,200,50), limits=c(-1100, 200), expand = c(0.01, 0)) +
  scale_y_continuous(breaks=seq(0,50,5), limits=c(0,50), expand=c(0.01,0)) +
  xlab("Log Likelihood Ratio (LLR)") +
  ylab("No. variants") +
  ggtitle("Case-control LLR for variants classified in ClinVar \n(CHEK2)") -> chek2

ggarrange(brca1,brca2,palb2,atm,chek2,labels=c("A","B","C","D","E"),ncol=2,nrow=3)


# Sensitivity/Specificity
vcep_with_clinvar %>% filter(ACMG_combo>=1) %>% filter(grepl("athogenic",ClinicalSignificance)) #81
vcep_with_clinvar %>% filter(ACMG_combo<=-1) %>% filter(grepl("athogenic",ClinicalSignificance)) 
vcep_with_clinvar %>% filter(ACMG_combo<=-8) %>% filter(grepl("athogenic",ClinicalSignificance))
vcep_with_clinvar %>% filter(ACMG_combo<=1 & ACMG_combo>=-1) %>% filter(grepl("athogenic",ClinicalSignificance)) 
vcep_with_clinvar %>% filter(grepl("athogenic",ClinicalSignificance)) -> path_mean
mean(path_mean$ACMG_combo)

vcep_with_clinvar %>% filter(ACMG_combo<=-1) %>% filter(grepl("enign",ClinicalSignificance))
vcep_with_clinvar %>% filter(ACMG_combo<=1 & ACMG_combo>=-1) %>% filter(grepl("enign",ClinicalSignificance))
vcep_with_clinvar %>% filter(ACMG_combo>=1) %>% filter(grepl("enign",ClinicalSignificance))
vcep_with_clinvar %>% filter(ACMG_combo>=8) %>% filter(grepl("enign",ClinicalSignificance))
vcep_with_clinvar %>% filter(grepl("enign",ClinicalSignificance)) -> ben_mean
mean(ben_mean$ACMG_combo)

sensitivity = 81/(81+33) #71.3%
specificity = 208/(208+13) #91.7%
accuracy = (81+208) / (335) #86.2%

library(ROCR)
pred=prediction()


# Bar plot using ClinVar B/LB and P/LP variants to show concordance, 2* classifications in BRCA1 and BRCA2
read_csv("Analysis/P2_enriched/missense_LR_oddsAF_palb2hi_filtered_withClinVar.csv") -> vcep_with_clinvar

vcep_with_clinvar %>%
  filter(!grepl("#N/A",ClinicalSignificance) & !grepl("Uncertain",ClinicalSignificance) & !grepl("not provided",ClinicalSignificance) & !grepl("Conflict",ClinicalSignificance)) %>%
  filter(ClinVar_Stars >= 2) -> vcep_with_clinvar

vcep_with_clinvar %>% filter(ACMG_combo>1) %>% filter(gene=="BRCA1") %>% filter(grepl("athogenic",ClinicalSignificance)) #29
vcep_with_clinvar %>% filter(`All Datasets`=="P_SUP") %>% filter(gene=="BRCA1"|gene=="BRCA2") %>% filter(grepl("enign",ClinicalSignificance)) #7
vcep_with_clinvar %>% filter(`All Datasets`=="P_MOD") %>% filter(gene=="BRCA1"|gene=="BRCA2") %>% filter(grepl("enign",ClinicalSignificance)) #2
vcep_with_clinvar %>% filter(`All Datasets`=="P_STR") %>% filter(gene=="BRCA1"|gene=="BRCA2") %>% filter(grepl("enign",ClinicalSignificance)) #7
vcep_with_clinvar %>% filter(`All Datasets`=="P_VSTR") %>% filter(gene=="BRCA1"|gene=="BRCA2") %>% filter(grepl("enign",ClinicalSignificance)) #13
vcep_with_clinvar %>% filter(ACMG_combo != 0) %>% filter(gene=="BRCA1"|gene=="BRCA2") %>% filter(grepl("enign",ClinicalSignificance)) #37

vcep_with_clinvar %>% filter(ACMG_combo< -1) %>% filter(gene=="BRCA2"|gene=="BRCA1") %>% filter(grepl("enign",ClinicalSignificance)) #100
vcep_with_clinvar %>% filter(`All Datasets`=="P_SUP") %>% filter(gene=="BRCA1") %>% filter(grepl("enign",ClinicalSignificance)) #1
vcep_with_clinvar %>% filter(`All Datasets`=="P_MOD") %>% filter(gene=="BRCA1") %>% filter(grepl("enign",ClinicalSignificance)) #6
vcep_with_clinvar %>% filter(`All Datasets`=="P_STR") %>% filter(gene=="BRCA1") %>% filter(grepl("enign",ClinicalSignificance)) #6
vcep_with_clinvar %>% filter(`All Datasets`=="P_VSTR") %>% filter(gene=="BRCA1") %>% filter(grepl("enign",ClinicalSignificance)) #87
vcep_with_clinvar %>% filter(ACMG_combo != 0) %>% filter(gene=="BRCA2") %>% filter(grepl("enign",ClinicalSignificance)) #106


vcep_with_clinvar %>%
  filter(gene=="BRCA1") %>%
  select(ACMG_combo, ClinicalSignificance) %>%
  ggplot(aes(x=ACMG_combo)) +
  geom_histogram(aes(fill=ClinicalSignificance),binwidth=15, color="black", center=0.5, position="stack") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        text = element_text(size=12)) +
  scale_fill_manual(values=c("#130857","#444cc7","#2f9fd4","#cc6600","#B2182B","#5e0101")) + 
  scale_x_continuous(breaks=seq(-1100,200,50), limits=c(-1100, 200), expand = c(0.01, 0)) +
  scale_y_continuous(breaks=seq(0,30,2), limits=c(0,30), expand=c(0.01,0)) +
  xlab("Log Likelihood Ratio (LLR)") +
  ylab("No. variants") +
  ggtitle("Case-control LLR for variants classified in ClinVar \n(BRCA1)") -> brca1
vcep_with_clinvar %>%
  filter(gene=="BRCA2") %>%
  select(ACMG_combo, ClinicalSignificance) %>%
  ggplot(aes(x=ACMG_combo)) +
  geom_histogram(aes(fill=ClinicalSignificance),binwidth=15, color="black", center=0.5, position="stack") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        text = element_text(size=12)) +
  scale_fill_manual(values=c("#130857","#444cc7","#2f9fd4","#cc6600","#B2182B","#5e0101")) + 
  scale_x_continuous(breaks=seq(-1100,200,50), limits=c(-1100, 200), expand = c(0.01, 0)) +
  scale_y_continuous(breaks=seq(0,30,2), limits=c(0,30), expand=c(0.01,0)) +
  xlab("Log Likelihood Ratio (LLR)") +
  ylab("No. variants") +
  ggtitle("Case-control LLR for variants classified in ClinVar \n(BRCA2)") -> brca2

ggarrange(brca1,brca2,labels=c("A","B"),ncol=2,nrow=1)


### ClinVar VUS and Conflicting variants ----
# total no. missense VUS and conflicting variants in ClinVar
# clinvar_allstar %>%
#   mutate(ID = paste(gene, hgvs_cdna, sep=":")) %>%
#   # inner_join(case_counts_with_combo, clinvar_allstar, by=c("gene", "hgvs_cdna")) %>%
#   # left_join(clinvar_vep, by=c("ID"="#Uploaded_variation")) %>%
#   filter(grepl("missense",Consequence)) %>%
#   # filter(gene %in% c("BRCA1","BRCA2")) %>%
#   filter(!grepl("#N/A",ClinicalSignificance) & ((grepl("Uncertain",ClinicalSignificance) 
#                                                  | grepl("Conflicting",ClinicalSignificance)
#   ))) -> vcep_with_clinvar

read_csv("Analysis/P2_enriched/missense_LR_oddsAF_palb2hi_filtered_withClinVar.csv") -> vcep_with_clinvar

vcep_with_clinvar %>%
  # select(1,2,105,106,110:112) %>%
  # filter(gene %in% c("BRCA1","BRCA2")) %>%
  filter(ACMG_combo!=0) %>%
  filter(ClinVar_Stars > 0) %>%
  filter(!grepl("#N/A",ClinicalSignificance) & ((grepl("Uncertain",ClinicalSignificance) 
                                                 | grepl("Conflicting",ClinicalSignificance)
  ))) %>%
  mutate(ClinicalSignificance = case_when(ClinicalSignificance=="Conflicting classifications of pathogenicity; risk factor" ~ "Conflicting",
                                          ClinicalSignificance=="Conflicting classifications of pathogenicity" ~ "Conflicting",
                                          ClinicalSignificance=="Uncertain significance" ~ "Uncertain/VUS")) %>%
  filter(!is.na(constituent_data)) -> vcep_with_clinvar_2

write_csv(vcep_with_clinvar_2, "Analysis/P2_enriched/clinvar_missense_vus_with_data.csv")

vcep_with_clinvar_2 %>%
  filter(gene=="BRCA1"|gene=="BRCA2") %>%
  select(ACMG_combo, ClinicalSignificance) %>%
  ggplot(aes(x=ACMG_combo)) +
  geom_histogram(aes(fill=ClinicalSignificance),binwidth=15, color="black", center=0.5, position="stack") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        text = element_text(size=20)) +
  scale_fill_manual(values=c("green","purple")) + 
  scale_x_continuous(breaks=seq(-1100,200,50), limits=c(-1100, 200), expand = c(0.01, 0)) +
  scale_y_continuous(breaks=seq(0,1700,100), limits=c(0,1700), expand=c(0.01,0)) +
  xlab("EPs / Log Likelihood Ratio (LLR)") +
  ggtitle("Case-control LLR for ClinVar VUS and \nConflicting classifications")

