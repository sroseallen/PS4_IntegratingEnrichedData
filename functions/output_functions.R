### Enrichment Factor calculations (+ confidence intervals)
EF_calculator <- function(PTV_carriers_lab, tot_tested_lab, PTV_carriers_unselected, tot_tested_unselected) {
  EF <- (PTV_carriers_lab/tot_tested_lab) * (tot_tested_unselected/PTV_carriers_unselected)
  LCI <- exp(log(EF)-(1.96*sqrt(1/PTV_carriers_lab - 1/tot_tested_lab + 1/PTV_carriers_unselected - 1/tot_tested_unselected)))
  UCI <- exp(log(EF)+(1.96*sqrt(1/PTV_carriers_lab - 1/tot_tested_lab + 1/PTV_carriers_unselected - 1/tot_tested_unselected)))
  return(list(EF,LCI,UCI))
}

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
# set seed to ensure results can be replicated; seed used for main paper analysis=123
set.seed(123)

sample_controls <- function(carriers, prob=c(0.2914637,0.7085363)) {
  assign <- sample(c("A", "B"), size=carriers, replace=TRUE, prob=prob) #A = ukb, B = ndrs
  A <- sum(str_count(assign, "A"))
  B <- sum(str_count(assign, "B"))
  return (A)
}

# batch generation of OR and LR stats (run per-variant)
batch_stats <- function(df, hgvs, col1, col2, col3, col4, new_colnames=c("p-value","OR","Lower_CI","Upper_CI","LR_path","ACMG_points","Target_Odds"), target_odds_BRCA1=4, target_odds_BRCA2=4, target_odds_PALB2=4, target_odds_ATM=2, target_odds_CHEK2=2, ci=0) {
  
  newdata_df <- data.frame(matrix(ncol=7,nrow=0))
  colnames(newdata_df) <- new_colnames
  
  # pull out case carrier and control carrier counts, and total tested counts for cases and controls
  for (i in 1:nrow(df)) {
    df[i,col1] -> a #case carrier
    df[i,col2] -> b #control carrier
    df[i,col3] -> c #total cases
    df[i,col4] -> d #total controls
    df[i, "gene"] -> gene_name
    
    a <- as.numeric(a)
    b <- as.numeric(b)
    c <- as.numeric(c)
    d <- as.numeric(d)
    
    # calculate OR, p-value, and confidence intervals (and decide if haldane correction is necessary for OR calculation)
    if (a != 0 & c == 0) {
      new_OR <- all_ORstat(a,(b-a),c,(d-c),haldane=TRUE,critical_value=1.96)
    } else {
      new_OR <- all_ORstat(a,(b-a),c,(d-c),haldane=FALSE,critical_value=1.96)
    }
    
    ## variants which have a count that prevents proper calculation of LR with PS4-LR-Calc + will need to be entered as "NC - not calculable" in the function output, see Sup Methods
    brca1_highcount <- c("c.1487G>A")
    brca2_highcount <- c("c.7544C>T")
    palb2_highcount <- c()
    atm_highcount <- c("c.6067G>A","c.998C>T","c.5558A>T","c.5071A>C","c.1229T>C","c.3925G>A")
    chek2_highcount <- c("c.1100del","c.470T>C")
    
    # calculate LR and ACMG points
    if (gene_name == "BRCA1") {
      ## if variant count too high to calculate, output that here (see Supp methods)
      if (df[i,hgvs] %in% brca1_highcount & "LR_path_ambry" %in% new_colnames) {
        new_LR <- data.frame(LR_path = c("NC - not calculable"), ACMG_points = c("NC - not calculable"))
        ## if there are no variant carriers, output NA
      } else if (a==0 & c==0) {
        new_LR <- data.frame(LR_path = c("NA"), ACMG_points = c("NA"))
        ## otherwise, calculate the LR per updated PS4-LR-Calc methodology (Rowlands et al., 2025 + see Supp methods)
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
    
    # add all stats to main df
    newdata_df <- rbind(newdata_df, new_stats)
  }
  cbind(df, newdata_df)
}

# filtering datasets based on determined parameters (single obs, discrepancies between datasets, potential moderate penetrance)
parameter_filtering <- function(df, dataset, case, control, case_denom, control_denom, points, upper_or, lower_or, control_case_thresh, case_control_thresh, gene="gene", OR_UPPER_LIMIT_BRCA=4, OR_UPPER_LIMIT_PALB2=3, OR_UPPER_LIMIT_ATMCHEK2=2, OR_LOWER_LIMIT=1, EXCLUDE_ZERO=TRUE, EXCLUDE_SINGLE=TRUE, EXCLUDE_DISCREPANT=TRUE, EXCLUDE_MODERATE=TRUE) {
  df %>% mutate(zero_obs = 0, single_obs=0, discrepant_evidence=0, low_coverage=0, moderate=0) -> df
  
  ## Parameter 1.1: Zero observation across cases and controls
  if (EXCLUDE_ZERO==TRUE) {
    df %>% mutate(zero_obs := if_else((!!sym(case) == 0 & !!sym(control) == 0),1,0)) -> df
  }
  ## Parameter 1.2: Single observation across cases and controls
  if (EXCLUDE_SINGLE==TRUE) {
    df %>% mutate(single_obs := if_else((!!sym(case) == 1 & !!sym(control) == 0) | (!!sym(case) == 0 & !!sym(control) == 1),1,0)) -> df
  }
  
  # Parameter 2: Case:control ratio exceeds limit for appropriate use (i.e. would give false positive or false negative evidence)
  if (EXCLUDE_DISCREPANT==TRUE) {
    df %>% mutate(discrepant_evidence := case_when((!!sym(control_denom) / !!sym(case_denom)) > !!sym(control_case_thresh) ~ 1,
                                                   (!!sym(case_denom) / !!sym(control_denom)) > !!sym(case_control_thresh) ~ 1,
                                                   .default = 0)) -> df
  }
  
  ## Parameter 3: Variants of likely low penetrance (ie target OR is higher than actual OR, set to target OR of 4) (see Supp methods, Supplementary Table 14)
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

# plotting categorical evidence strengths per ACMG v3.0 boundaries
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


# identify the dataset imbalance size at which false positive evidence starts to occur (used to create denom_limits lookup file)
# thresholds where there are 0 case carriers and 1 or more control carriers (number of controls set at 1 by default)
ratio_threshold_controls <- function(controls=1, ci=0, dataset_ratio=10, target) {
  ## all cases = 0, all controls = 1, want to find the ratio at which acmg points are >0 (and >1 for a permissive slant)
  control_denom = c(seq(500,2500000,by=500))
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

