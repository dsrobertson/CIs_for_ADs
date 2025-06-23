### Code for Section 4 (Figures 3): values of coverage and CI width for the
### standard (Wald) CI as the value of p_CE varies.

################################################################################

source('MUSEC_trial.R')

N = 10^6  # Number of bootstrap replicates

### Run simulation

phat_CE_vec = phat_CE + seq(from = -0.07, to = 0.14, by = 0.01)

stop1_vec_Wald = rep(NA, length(phat_CE_vec))

CI_coverage_Wald = CI_coverage_stop_early_Wald = 
  CI_width_Wald = CI_width_stop_early_Wald = CI_width_se_Wald = 
  CI_width_se_stop_early_Wald =   CI_coverage_continue_Wald = 
  CI_width_continue_Wald = CI_width_se_continue_Wald =
  CI_consistency_Wald = CI_consistency_stop_early_Wald = CI_consistency_continue_Wald = 
  rep(NA, length(phat_CE_vec))

standard_CI = matrix(nrow = N, ncol = 2)

set.seed(7)

Sys.time()

for(i in 1:length(phat_CE_vec)){
  
  print(i)
  
  phat_CE = phat_CE_vec[i]
  
  s1_1 = rbinom(N, size = n1_1, prob = phat_CE)
  s0_1 = rbinom(N, size = n0_1, prob = phat_P)
  
  ptilde_1 = (s0_1 + s1_1)/(n0_1 + n1_1)
  Ihat_1 = 1/(ptilde_1*(1-ptilde_1)*(1/n0_1 + 1/n1_1))
  Z1 = (s1_1/n1_1 - s0_1/n0_1)*sqrt(Ihat_1)

  s1_2 = s1_1 + rbinom(N, size = n1_2-n1_1, prob = phat_CE)
  s0_2 = s0_1 + rbinom(N, size = n0_2-n0_1, prob = phat_P)
  
  ptilde_2 = (s0_2 + s1_2)/(n0_2 + n1_2)
  Ihat_2 = 1/(ptilde_2*(1-ptilde_2)*(1/n0_2 + 1/n1_2))
  Z2 = (s1_2/n1_2 - s0_2/n0_2)*sqrt(Ihat_2) 
  
  
  # Stage 1
  phat_CE_boot = s1_1/n1_1
  phat_P_boot = s0_1/n0_1
  
  MLE_boot = phat_CE_boot - phat_P_boot
  
  MLE_se = sqrt(phat_CE_boot*(1-phat_CE_boot)/n1_1 + 
                  phat_P_boot*(1-phat_P_boot)/n0_1)

  standard_CI1 = cbind(MLE_boot - qnorm(1-alpha/2)*MLE_se,
                   MLE_boot + qnorm(1-alpha/2)*MLE_se)
  
  # Stage 2
  phat_CE_boot = s1_2/n1_2
  phat_P_boot = s0_2/n0_2
  
  MLE_boot = phat_CE_boot - phat_P_boot
  
  MLE_se = sqrt(phat_CE_boot*(1-phat_CE_boot)/n1_2 + 
                  phat_P_boot*(1-phat_P_boot)/n0_2)
  
  standard_CI2 = cbind(MLE_boot - qnorm(1-alpha/2)*MLE_se,
                   MLE_boot + qnorm(1-alpha/2)*MLE_se)

  
  standard_CI[Z1 > e1,] = standard_CI1[Z1 > e1,]
  standard_CI[Z1 <= e1,] = standard_CI2[Z1 <= e1,]
  
  stop1 = rej2 = rep(0, N)
  
  stop1[Z1 > e1] = 1
  
  rej2[Z2 > e2] = 1
  
  rej = cbind(stop1, rej2)
  
  ################
  
  delta = phat_CE - phat_P
  
  stop1_vec_Wald[i] = mean(stop1)
  
  
  ### Overall (unconditional)
  
  CI_coverage_Wald[i] = c(mean(standard_CI[,1] <= delta & standard_CI[,2] >= delta))
  
  CI_width_Wald[i] = c(mean(standard_CI[,2] - standard_CI[,1]))
  
  CI_width_se_Wald[i] = c(sd(standard_CI[,2] - standard_CI[,1]))
  
  CI_consistency_Wald[i] = c(mean((standard_CI[,1] > 0 & rowSums(rej)) |
                                    (standard_CI[,1] < 0 & !(rowSums(rej)))))
  
  
  ### Means for trials that stop early
  
  CI_coverage_stop_early_Wald[i] = c(mean(standard_CI[stop1==1,1] <= delta &
                                     standard_CI[stop1==1,2] >= delta))
  
  CI_width_stop_early_Wald[i] = c(mean(standard_CI[stop1==1,2] - standard_CI[stop1==1,1]))
  
  CI_width_se_stop_early_Wald[i] = c(sd(standard_CI[stop1==1,2] - standard_CI[stop1==1,1]))
  
  CI_consistency_stop_early_Wald[i] = c(mean(standard_CI[stop1==1,1] > 0))
  
  
  ### Means for trials that continue
  
  CI_coverage_continue_Wald[i] = c(mean(standard_CI[stop1==0,1] <= delta &
                                   standard_CI[stop1==0,2] >= delta))
  
  CI_width_continue_Wald[i] = c(mean(standard_CI[stop1==0,2] - standard_CI[stop1==0,1]))
  
  CI_width_se_continue_Wald[i] = c(sd(standard_CI[stop1==0,2] - standard_CI[stop1==0,1]))
  
  CI_consistency_continue_Wald[i] = c(mean((standard_CI[stop1==0,1] > 0 & rej2[stop1==0]) |
                                           (standard_CI[stop1==0,1] < 0 & !(rej2[stop1==0]))))
  
}

Sys.time()


save(stop1_vec_Wald, CI_coverage_Wald, CI_coverage_stop_early_Wald,
     CI_coverage_continue_Wald, CI_width_Wald, CI_width_se_Wald, 
     CI_width_stop_early_Wald, CI_width_se_stop_early_Wald,
     CI_width_continue_Wald, CI_width_se_continue_Wald,
     file = 'sim_results_Wald.Rdata')
