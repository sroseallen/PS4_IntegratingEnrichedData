library(tidyverse)
library(stringr)

### FISHERS TEST (2-tailed) (non-Haldane)
p_value_fishers <- function(a, b, c, d) {
  fisher_matrix <- matrix (c(a,b,c,d), nrow=2, byrow=TRUE)
  return (fisher.test(fisher_matrix)$p.value)
}

### CHI SQUARED TEST
p_value_chi <- function (a,b,c,d) {
  chisq <- chisq.test(x=matrix(c(a,b,c,d), nrow=2, ncol=2), correct=TRUE)[[3]]
  return(chisq)
}

### OR CALCULATION (non-Haldane)
odds <- function(a,b,c,d) {
  ((a/b)/(c/d))
}

lower_ci <- function(a,b,c,d,critical_value=1.96) {
  OR <- odds(a,b,c,d)
  exp(log(OR) - (critical_value*sqrt(1/a + 1/b + 1/c + 1/d)))
}
upper_ci <- function(a,b,c,d,critical_value=1.96) {
  OR <- odds(a,b,c,d)
  exp(log(OR) + (critical_value*sqrt(1/a + 1/b + 1/c + 1/d)))
}


all_ORstat <- function(a,b,c,d,haldane=FALSE,critical_value=1.96) {
  if (c==0 & haldane==FALSE) {
    warning("No observation of this variant in controls (c=0), no odds ratio generated.\nHighly suggest to re-run this function with the 'haldane=TRUE' argument.")
  }
  if (a==0 & haldane==FALSE) {
    warning("No observation of this variant in cases (a=0), no odds ratio generated.")
  }
  
  p_calc <- p_value_fishers(a,b,c,d)
  
  if(haldane == TRUE) {
    a=a+0.5
    b=b+0.5
    c=c+0.5
    d=d+0.5
    odds_calc <- odds(a,b,c,d)
    lower_calc <- lower_ci(a,b,c,d,critical_value)
    upper_calc <- upper_ci(a,b,c,d,critical_value)
    
  } else if (haldane == FALSE) {
    if (a != 0) {
      odds_calc <- odds(a,b,c,d)
      lower_calc <- lower_ci(a,b,c,d,critical_value)
      upper_calc <- upper_ci(a,b,c,d,critical_value)
    } else if (a == 0) {
      odds_calc <- 0
      lower_calc <- 0
      upper_calc <- 0
    }
  }
  
  df <- data.frame(p_value=p_calc, odds_ratio=odds_calc, lower_ci_95=lower_calc, upper_ci_95=upper_calc)
  return(df)
}
 
# # Example cases: Generating for a single variant
# all_ORstat(516,58,(48826-516),(50703-58))
# all_ORstat(516,0,(48826-516),(50703-58), haldane=TRUE)
# 
# # if using a spreadsheet
# f <- read_csv("BRIDGES_UKB_Combinations.csv") 
# 
# for (i in 1:length(f$hgvs_cdna)) {
#   pull(f[i,27]) -> a
#   pull(f[i,29]) -> b
#   pull(f[i,28]) - a -> c
#   pull(f[i,30]) - b -> d
#   
#   if(is.na(b)) {b <- 0}
#   if(is.na(d)) {d <- 0}
#   
#   f[i,37] <- odds(a,b,c,d)
#   f[i,38] <- lower_ci(a,b,c,d)
#   f[i,39] <- upper_ci(a,b,c,d)
#   f[i,40] <- p_value_fishers(a,b,c,d)
#   
#   if(b==0 & a!=0){
#     a+0.5 -> a_haldane
#     b+0.5 -> b_haldane
#     c+0.5 -> c_haldane
#     d+0.5 -> d_haldane
#     
#     f[i,37] <- odds(a_haldane,b_haldane,c_haldane,d_haldane)
#     f[i,38] <- lower_ci(a_haldane,b_haldane,c_haldane,d_haldane)
#     f[i,39] <- upper_ci(a_haldane,b_haldane,c_haldane,d_haldane)
#   }
# }
# 
# write.csv(f, "output.csv")
