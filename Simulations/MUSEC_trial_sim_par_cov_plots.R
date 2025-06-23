### Code for Section 4 (Figures 2 and 3): Plots of probability of early 
### stopping, coverage and CI width as the value of p_CE varies.

################################################################################

source('MUSEC_trial.R')

N = 10^5  # Number of bootstrap replicates
M = 10^4   # Number of inner bootstrap replicates

### Set up parallel computing

library(doParallel)
library(doRNG)
library(data.table)

mc = 32; cl = makeCluster(mc)

registerDoParallel(cl)
registerDoRNG(7)

### Run simulation

phat_CE_vec = phat_CE + seq(from = -0.07, to = 0.14, by = 0.01)


CI_names = c('Wald test (standard)', 'Unconditional exact', 
             'RCI', 'Adjusted asymptotic', 'Unconditional bootstrap',
             'Conditional exact', 'Restricted conditional',
             'Conditional likelihood', 'Penalised likelihood')


CI_coverage = CI_coverage_stop_early = 
  CI_width = CI_width_stop_early = CI_width_se = CI_width_se_stop_early = 
  matrix(nrow = 9, ncol = length(phat_CE_vec))

rownames(CI_coverage) = rownames(CI_coverage_stop_early) = CI_names


CI_names = c('Wald test (standard)', 'Unconditional exact', 
             'RCI', 'Adjusted asymptotic', 'Unconditional bootstrap',
             'Conditional exact', 'Restricted conditional',
             'Conditional likelihood')

CI_coverage_continue = CI_width_continue = CI_width_se_continue = 
  matrix(nrow = 8, ncol = length(phat_CE_vec))

rownames(CI_coverage_continue) = CI_names

stop1_vec = rep(NA, length(phat_CE_vec))

Sys.time()

for(i in 1:length(phat_CE_vec)){
 
  print(i)
  
  phat_CE = phat_CE_vec[i]
   
  CI_results = foreach (j = 1:N,
                        .combine = function(x,y)rbindlist(list(x,y)),
                        .packages = c('mvtnorm','DPQ','cubature'),
                        .errorhandling='pass') %dopar% {
                          
                          s1_1 = rbinom(1, size = n1_1, prob = phat_CE)
                          s0_1 = rbinom(1, size = n0_1, prob = phat_P)
                          
                          ptilde_1 = (s0_1 + s1_1)/(n0_1 + n1_1)
                          Ihat_1 = 1/(ptilde_1*(1-ptilde_1)*(1/n0_1 + 1/n1_1))
                          Z1 = (s1_1/n1_1 - s0_1/n0_1)*sqrt(Ihat_1)
                          
                          if(Z1 > e1 | Z1 < f1){
                            
                            stop1 = 1
                            rej2 = 0
                            
                            s1_2 = s0_2 = NA
                            
                            
                            # Standard  (Wald)
                            phat_CE_boot = s1_1/n1_1
                            phat_P_boot = s0_1/n0_1
                            
                            MLE_boot = phat_CE_boot - phat_P_boot
                            
                            MLE_se = sqrt(phat_CE_boot*(1-phat_CE_boot)/n1_1 + 
                                            phat_P_boot*(1-phat_P_boot)/n0_1)
                            
                            standard_CI = c(MLE_boot - qnorm(1-alpha/2)*MLE_se,
                                            MLE_boot + qnorm(1-alpha/2)*MLE_se)
                        
                            # Exact
                            lb = (qnorm(alpha/2)+Z1)/sqrt(Ihat_1)
                            ub = (qnorm(1-alpha/2)+Z1)/sqrt(Ihat_1)
                            
                            exact_CI = c(lb, ub)
                            
                            
                            # RCI
                            RCI = c(Z1/sqrt(Ihat_1) - ck_alpha/sqrt(Ihat_1)*sqrt(2),
                                    Z1/sqrt(Ihat_1) + ck_alpha/sqrt(Ihat_1)*sqrt(2))
                            
                            # Adjusted asymptotic
                            V1 = Ihat_1
                            
                            theta = Z1/sqrt(V1)
                            
                            fn1 = function(z) {
                              (1/sqrt(V1))*dnorm((z-theta*V1)/sqrt(V1))
                            }
                            
                            E_Z_int1 = function(x) {
                              x*fn1(x)
                            }
                            
                            E_Z = integrate(E_Z_int1, lower = e1, upper = Inf)$value
                            
                            E_Z_sqrtV = (1/sqrt(V1))*integrate(E_Z_int1, lower = e1, upper = Inf)$value
                            
                            E_Z2_int1 = function(x) {
                              (x^2)*fn1(x)
                            }
                            
                            E_Z2_V = (1/V1)*integrate(E_Z2_int1, lower = e1, upper = Inf)$value
                            
                            E_V_int1 = function(x) {
                              fn1(x)
                            }
                            
                            E_sqrtV = sqrt(V1)*cubintegrate(E_V_int1, lower = e1, upper = Inf)$integral
                            
                            E_V = V1*cubintegrate(E_V_int1, lower = e1, upper = Inf)$integral
                            
                            mu_adj = E_Z_sqrtV - theta*E_sqrtV
                            
                            if(E_Z2_V - 2*theta*E_Z + theta^2*E_V - mu_adj^2 > 0){
                              
                              sigma_adj = sqrt(E_Z2_V - 2*theta*E_Z + theta^2*E_V - mu_adj^2)
                              
                              adjusted_asym_est = theta - mu_adj/sqrt(V1)
                              
                              adjusted_asym_CI = c(adjusted_asym_est - qnorm(1-alpha/2)*sigma_adj/sqrt(V1),
                                                   adjusted_asym_est + qnorm(1-alpha/2)*sigma_adj/sqrt(V1))
                              
                            } else {
                              
                              adjusted_asym_CI = c(NA, NA)
                            }
                            
                            # Unconditional bootstrap
                            
                            MLEsim = rep(NA, M)
                            
                            s1_1_boot = rbinom(M, size = n1_1, prob = phat_CE_boot)
                            s0_1_boot = rbinom(M, size = n0_1, prob = phat_P_boot)
                            
                            ptilde_1_boot = (s0_1_boot + s1_1_boot)/(n0_1 + n1_1)
                            Ihat_1_boot = 1/(ptilde_1_boot*(1-ptilde_1_boot)*(1/n0_1 + 1/n1_1))
                            Z1_boot = (s1_1_boot/n1_1 - s0_1_boot/n0_1)*sqrt(Ihat_1_boot) 
                            
                            s1_2_boot = s1_1_boot + rbinom(M, size = n1_2-n1_1, prob = phat_CE_boot)
                            s0_2_boot = s0_1_boot + rbinom(M, size = n0_2-n0_1, prob = phat_P_boot)
                            
                            ptilde_2_boot = (s0_2_boot + s1_2_boot)/(n0_2 + n1_2)
                            Ihat_2_boot = 1/(ptilde_2_boot*(1-ptilde_2_boot)*(1/n0_2 + 1/n1_2))
                            Z2_boot = (s1_2_boot/n1_2 - s0_2_boot/n0_2)*sqrt(Ihat_2_boot) 
                            
                            MLEsim[Z1_boot > e1] = Z1_boot[Z1_boot > e1]/sqrt(Ihat_1_boot[Z1_boot > e1])
                            MLEsim[Z1_boot <= e1] = Z2_boot[Z1_boot <= e1]/sqrt(Ihat_2_boot[Z1_boot <= e1])
                            
                            bootstrap_uncond_CI = c(quantile(MLEsim, probs = 0.025),
                                                    quantile(MLEsim, probs = 0.975))
                            
                            # Conditional exact
                            
                            cb_optim  = function(theta, z, I1, u, q) {
                              
                              log_S1     = pnorm(u, theta*sqrt(I1), 1,
                                                 lower.tail = FALSE, log.p = TRUE)
                              
                              num    = pnorm(z, theta*sqrt(I1), lower.tail = FALSE, log.p = TRUE)
                              
                              int    = num - log_S1
                              
                              (int - log(q))^2
                            }
                            
                            clb = optimise(cb_optim, c(standard_CI[1]-100, standard_CI[1]+1),
                                           z = Z1, I1 = Ihat_1,
                                           u = e1, q = 0.025)$minimum
                            
                            cub = optimise(cb_optim, c(standard_CI[2]-1, standard_CI[2]+1),
                                           z = Z1, I1 = Ihat_1,
                                           u = e1, q = 0.975)$minimum
                            
                            exact_cond_CI = c(max(clb,-1), min(cub,1))
                            
                            rlb = (e1 - qnorm(0.975))/sqrt(Ihat_1)
                            
                            if(rlb < cub){
                              
                              r_cond_CI = c(max(clb, rlb), cub)
                              
                            } else {
                              
                              r_cond_CI = c(NA, NA)
                            }
                            
                            
                            # Conditional (penalised) likelihood
                            
                            plog_lik = function(theta, lambda, z, u, I1){
                              
                              -0.5*(z - theta*sqrt(I1))^2 -
                                lambda*pnorm(u, theta*sqrt(I1), 1,
                                             lower.tail = FALSE, log.p = TRUE)
                            }
                            
                            pLE = function(lambda, z, u, I1){
                              
                              optimise(plog_lik,
                                       c(standard_CI[1]-1, standard_CI[2]+1),
                                       lambda = lambda, z = z,
                                       u = u, I1 = I1, maximum = TRUE)$maximum
                              
                            }
                            
                            
                            cMLEsim = pMLEsim = rep(NA, M)
                            boot_count = 0
                            
                            while(boot_count < M){
                              
                              s1_1_boot = rbinom(1, size = n1_1, prob = phat_CE_boot)
                              s0_1_boot = rbinom(1, size = n0_1, prob = phat_P_boot)
                              
                              ptilde_1_boot = (s0_1_boot + s1_1_boot)/(n0_1 + n1_1)
                              Ihat_1_boot = 1/(ptilde_1_boot*(1-ptilde_1_boot)*(1/n0_1 + 1/n1_1))
                              Z1_boot = (s1_1_boot/n1_1 - s0_1_boot/n0_1)*sqrt(Ihat_1_boot)
                              
                              if(Z1_boot > e1){
                                
                                boot_count = boot_count + 1
                                
                                cMLEsim[boot_count] = optimise(plog_lik,
                                                               c(MLE_boot-100,MLE_boot+1),
                                                               lambda = 1,
                                                               z = Z1_boot, I1 = Ihat_1_boot,
                                                               u = e1, maximum = TRUE)$maximum
                                
                                lambda_star = uniroot(pLE, c(0,1), z = e1, u = e1,
                                                      I1 = Ihat_1_boot)$root
                                
                                
                                pMLEsim[boot_count] = optimise(plog_lik,
                                                               c(MLE_boot-1,MLE_boot+1),
                                                               lambda = lambda_star, 
                                                               z = Z1_boot, I1 = Ihat_1_boot,
                                                               u = e1, maximum = TRUE)$maximum
                                
                              }
                            }
                            
                            clik_CI = c(max(quantile(cMLEsim, probs = 0.025),-1),
                                        min(quantile(cMLEsim, probs = 0.975),1))
                            
                            plik_CI = c(quantile(pMLEsim, probs = 0.025),
                                        quantile(pMLEsim, probs = 0.975))
                            
                          } else {
                            
                            stop1 = 0
                            
                            s1_2 = s1_1 + rbinom(1, size = n1_2-n1_1, prob = phat_CE)
                            s0_2 = s0_1 + rbinom(1, size = n0_2-n0_1, prob = phat_P)
                            
                            ptilde_2 = (s0_2 + s1_2)/(n0_2 + n1_2)
                            Ihat_2 = 1/(ptilde_2*(1-ptilde_2)*(1/n0_2 + 1/n1_2))
                            Z2 = (s1_2/n1_2 - s0_2/n0_2)*sqrt(Ihat_2) 
                            
                            if(Z2 > e2){
                              rej2 = 1
                            } else {
                              rej2 = 0
                            }
                            
                            # Standard (Wald)
                            phat_CE_boot = s1_2/n1_2
                            phat_P_boot = s0_2/n0_2
                            
                            MLE_boot = phat_CE_boot - phat_P_boot
                            
                            MLE_se = sqrt(phat_CE_boot*(1-phat_CE_boot)/n1_2 + 
                                            phat_P_boot*(1-phat_P_boot)/n0_2)
                            
                            standard_CI = c(MLE_boot - qnorm(1-alpha/2)*MLE_se,
                                            MLE_boot + qnorm(1-alpha/2)*MLE_se)

                            # Exact
                            cov_matrix = matrix(c(1, sqrt(Ihat_1/Ihat_2),
                                                  sqrt(Ihat_1/Ihat_2), 1),
                                                nrow = 2, byrow = TRUE)
                            
                            pval_fn = function(delta) {
                              
                              mu = delta*sqrt(c(Ihat_1,Ihat_2))
                              
                              pnorm(delta*sqrt(Ihat_1) - e1) + 
                                pmvnorm(lower=c(f1, Z2), upper=c(e1, Inf),
                                        mean = mu, sigma = cov_matrix)[1]
                            }
                            
                            obj_fn = function(p, q){
                              L = (pval_fn(delta=p) - q)^2
                            }
                            
                            lb = optimize(obj_fn, q = alpha/2,
                                          interval = c(standard_CI[1]-0.5,
                                                       standard_CI[1]+0.5))$minimum
                            
                            ub = optimize(obj_fn, q = 1 - alpha/2,
                                          interval = c(standard_CI[2]-0.5,
                                                       standard_CI[2]+0.5))$minimum
                            
                            exact_CI = c(lb, ub)
                            
                            
                            # RCI
                            RCI = c(Z2/sqrt(Ihat_2) - ck_alpha/sqrt(Ihat_2),
                                    Z2/sqrt(Ihat_2) + ck_alpha/sqrt(Ihat_2))
                            
                            
                            # Adjusted asymptotic
                            
                            V1 = Ihat_1
                            V2 = Ihat_2
                            
                            I2 = V2 - V1
                            
                            if(I2 > 0){
                              
                              theta = Z2/sqrt(V2)
                              
                              fn1 = function(z) {
                                (1/sqrt(V1))*dnorm((z-theta*V1)/sqrt(V1))
                              }
                              
                              fn2_int = function(x, z) {
                                (1/sqrt(I2))*dnorm((z-x-theta*I2)/sqrt(I2))*fn1(x)
                              }
                              
                              fn2 = function(z) {
                                sapply(z, function (z_i) cubintegrate(fn2_int, lower = f1, upper = e1,
                                                                      z = z_i)$integral)
                              }
                              
                              E_Z_int1 = function(x) {
                                x*fn1(x)
                              }
                              
                              E_Z_int2 = function(x) {
                                x*fn2(x)
                              }
                              
                              E_Z = integrate(E_Z_int1, lower = e1, upper = Inf)$value +
                                cubintegrate(E_Z_int2, lower = -Inf, upper = Inf)$integral
                              
                              E_Z_sqrtV = (1/sqrt(V1))*integrate(E_Z_int1, lower = e1, upper = Inf)$value +
                                (1/sqrt(V2))*cubintegrate(E_Z_int2, lower = -Inf, upper = Inf)$integral
                              
                              
                              E_Z2_int1 = function(x) {
                                (x^2)*fn1(x)
                              }
                              
                              E_Z2_int2 = function(x) {
                                (x^2)*fn2(x)
                              }
                              
                              E_Z2_V = (1/V1)*integrate(E_Z2_int1, lower = e1, upper = Inf)$value +
                                (1/V2)*cubintegrate(E_Z2_int2, lower = -Inf, upper = Inf)$integral
                              
                              E_V_int1 = function(x) {
                                fn1(x)
                              }
                              
                              E_V_int2 = function(x) {
                                fn2(x)
                              }
                              
                              E_sqrtV = sqrt(V1)*integrate(E_V_int1, lower = e1, upper = Inf)$value +
                                sqrt(V2)*cubintegrate(E_V_int2, lower = -Inf, upper = Inf)$integral
                              
                              E_V = V1*integrate(E_V_int1, lower = e1, upper = Inf)$value +
                                V2*cubintegrate(E_V_int2, lower = -Inf, upper = Inf)$integral
                              
                              mu_adj = E_Z_sqrtV - theta*E_sqrtV
                              
                              if(E_Z2_V - 2*theta*E_Z + theta^2*E_V - mu_adj^2 > 0){
                                
                                sigma_adj = sqrt(E_Z2_V - 2*theta*E_Z + theta^2*E_V - mu_adj^2)
                                
                                adjusted_asym_est = theta - mu_adj/sqrt(V2)
                                
                                adjusted_asym_CI = c(adjusted_asym_est - qnorm(1-alpha/2)*sigma_adj/sqrt(V2),
                                                     adjusted_asym_est + qnorm(1-alpha/2)*sigma_adj/sqrt(V2))
                                
                              } else {
                                
                                adjusted_asym_CI = c(NA, NA)
                              }
                              
                            } else {
                              
                              adjusted_asym_CI = c(NA, NA)
                            }
                            
                            # Unconditional bootstrap
                            MLEsim = rep(NA, M)
                            
                            s1_1_boot = rbinom(M, size = n1_1, prob = phat_CE_boot)
                            s0_1_boot = rbinom(M, size = n0_1, prob = phat_P_boot)
                            
                            ptilde_1_boot = (s0_1_boot + s1_1_boot)/(n0_1 + n1_1)
                            Ihat_1_boot = 1/(ptilde_1_boot*(1-ptilde_1_boot)*(1/n0_1 + 1/n1_1))
                            Z1_boot = (s1_1_boot/n1_1 - s0_1_boot/n0_1)*sqrt(Ihat_1_boot) 
                            
                            s1_2_boot = s1_1_boot + rbinom(M, size = n1_2-n1_1, prob = phat_CE_boot)
                            s0_2_boot = s0_1_boot + rbinom(M, size = n0_2-n0_1, prob = phat_P_boot)
                            
                            ptilde_2_boot = (s0_2_boot + s1_2_boot)/(n0_2 + n1_2)
                            Ihat_2_boot = 1/(ptilde_2_boot*(1-ptilde_2_boot)*(1/n0_2 + 1/n1_2))
                            Z2_boot = (s1_2_boot/n1_2 - s0_2_boot/n0_2)*sqrt(Ihat_2_boot) 
                            
                            MLEsim[Z1_boot > e1] = Z1_boot[Z1_boot > e1]/sqrt(Ihat_1_boot[Z1_boot > e1])
                            MLEsim[Z1_boot <= e1] = Z2_boot[Z1_boot <= e1]/sqrt(Ihat_2_boot[Z1_boot <= e1])
                            
                            bootstrap_uncond_CI = c(quantile(MLEsim, probs = 0.025),
                                                    quantile(MLEsim, probs = 0.975))
                            
                            # Conditional exact CIs
                            
                            if(I2 > 0) {
                              
                              f2  = function(z, theta, I1, I2, u) {
                                
                                -0.5*(z - theta*sqrt(I2))^2 +
                                  pnorm(u, z*sqrt(I1/I2), sqrt((I2 - I1)/I2), log.p = TRUE) - 
                                  log(sqrt(2*pi))
                                
                              }
                              
                              fcond = function(x, theta, I1, I2, u){
                                
                                log_S2 = pnorm(u, theta*sqrt(I1), log.p = TRUE)
                                exp(f2(x, theta, I1, I2, u) - log_S2)
                              }
                              
                              
                              cb_optim = function(theta, z, I1, I2, u, q) {
                                
                                int = integrate(fcond, lower = -Inf, upper = z,
                                                theta = theta, I1 = I1, I2 = I2, u = u)$value
                                
                                (int - q)^2
                              }
                              
                              clb = optimise(cb_optim,
                                             interval = c(standard_CI[1]-0.5,
                                                          standard_CI[1]+0.5),
                                             z = Z2, I1 = Ihat_1,
                                             I2 = Ihat_2, u = e1, q = 0.975)$minimum
                              
                              cub = optimise(cb_optim,
                                             interval = c(standard_CI[2]-0.5,
                                                          standard_CI[2]+0.5),
                                             z = Z2, I1 = Ihat_1,
                                             I2 = Ihat_2, u = e1, q = 0.025)$minimum
                              
                              exact_cond_CI = c(clb, cub)
                              
                              rub = (e1 - qnorm(0.025))/sqrt(Ihat_1)
                              
                              r_cond_CI = c(clb, min(cub, rub))
                              
                            } else {
                              
                              exact_cond_CI = r_cond_CI = c(NA, NA)
                            }
                            
                            
                            # Conditional likelihood
                            
                            cMLE_obj = function(theta, z, I1, I2, u){
                              
                              (z/sqrt(I2) - theta + (sqrt(I1)*dnorm(u - theta*sqrt(I1)))/
                                 (I2*pnorm(u - theta*sqrt(I1))))^2
                              
                            }
                            
                            cMLEsim = rep(NA, M)
                            boot_count = 0
                            
                            while(boot_count < M){
                              
                              s1_1_boot = rbinom(1, size = n1_1, prob = phat_CE_boot)
                              s0_1_boot = rbinom(1, size = n0_1, prob = phat_P_boot)
                              
                              ptilde_1_boot = (s0_1_boot + s1_1_boot)/(n0_1 + n1_1)
                              Ihat_1_boot = 1/(ptilde_1_boot*(1-ptilde_1_boot)*(1/n0_1 + 1/n1_1))
                              Z1_boot = (s1_1_boot/n1_1 - s0_1_boot/n0_1)*sqrt(Ihat_1_boot)
                              
                              if(Z1_boot <= e1){
                                
                                boot_count = boot_count + 1
                                
                                s1_2_boot = s1_1_boot + rbinom(1, size = n1_2-n1_1, prob = phat_CE_boot)
                                s0_2_boot = s0_1_boot + rbinom(1, size = n0_2-n0_1, prob = phat_P_boot)
                                
                                ptilde_2_boot = (s0_2_boot + s1_2_boot)/(n0_2 + n1_2)
                                Ihat_2_boot = 1/(ptilde_2_boot*(1-ptilde_2_boot)*(1/n0_2 + 1/n1_2))
                                Z2_boot = (s1_2_boot/n1_2 - s0_2_boot/n0_2)*sqrt(Ihat_2_boot)
                                
                                cMLEsim[boot_count] = optimise(cMLE_obj,
                                                               c(Z2_boot/sqrt(Ihat_2_boot)-1,
                                                                 Z2_boot/sqrt(Ihat_2_boot)+1),
                                                               z = Z2_boot, I1 = Ihat_1_boot,
                                                               I2 = Ihat_2_boot,
                                                               u = e1)$minimum
                                
                              }
                            }
                            
                            clik_CI = c(quantile(cMLEsim, probs = 0.025),
                                        quantile(cMLEsim, probs = 0.975))
                            
                            plik_CI = clik_CI
                            
                          }
                          
                          return(list(rej = c(stop1, rej2),
                                      s0 = c(s0_1, s0_2),
                                      s1 = c(s1_1, s1_2),
                                      standard_CI = standard_CI,
                                      exact_CI = exact_CI,
                                      RCI = RCI,
                                      adjusted_asym_CI = adjusted_asym_CI,
                                      bootstrap_uncond_CI = bootstrap_uncond_CI,
                                      exact_cond_CI = exact_cond_CI,
                                      r_cond_CI = r_cond_CI,
                                      clik_CI = clik_CI,
                                      plik_CI = plik_CI))
                        }

  rej = matrix(CI_results$rej, ncol = 2, byrow = TRUE)
  stop1 = rej[,1]
  rej2 = rej[,2]
  
  s0 = matrix(CI_results$s0, ncol = 2, byrow = TRUE)
  s1 = matrix(CI_results$s1, ncol = 2, byrow = TRUE)
  
  standard_CI = matrix(CI_results$standard_CI, ncol = 2, byrow = TRUE)
  exact_CI = matrix(CI_results$exact_CI, ncol = 2, byrow = TRUE)
  RCI = matrix(CI_results$RCI, ncol = 2, byrow = TRUE)
  adjusted_asym_CI = matrix(CI_results$adjusted_asym_CI, ncol = 2, byrow = TRUE)
  bootstrap_uncond_CI = matrix(CI_results$bootstrap_uncond_CI, ncol = 2, byrow = TRUE)
  exact_cond_CI = matrix(CI_results$exact_cond_CI, ncol = 2, byrow = TRUE)
  r_cond_CI = matrix(CI_results$r_cond_CI, ncol = 2, byrow = TRUE)
  clik_CI = matrix(CI_results$clik_CI, ncol = 2, byrow = TRUE)
  plik_CI = matrix(CI_results$plik_CI, ncol = 2, byrow = TRUE)
  
  ################
  
  delta = phat_CE - phat_P
  
  stop1_vec[i] = mean(stop1)
  
  
  ### Overall (unconditional)
  
  CI_coverage[,i] = c(mean(standard_CI[,1] <= delta & standard_CI[,2] >= delta),
                      mean(exact_CI[,1] <= delta & exact_CI[,2] >= delta),
                      mean(RCI[,1] <= delta & RCI[,2] >= delta),
                      sum(adjusted_asym_CI[,1] <= delta & adjusted_asym_CI[,2] >= delta, na.rm = T)/N,
                      mean(bootstrap_uncond_CI[,1] <= delta & bootstrap_uncond_CI[,2] >= delta),
                      mean(exact_cond_CI[,1] <= delta & exact_cond_CI[,2] >= delta),
                      sum(r_cond_CI[,1] <= delta & r_cond_CI[,2] >= delta, na.rm = TRUE)/N,
                      mean(clik_CI[,1] <= delta & clik_CI[,2] >= delta),
                      mean(plik_CI[,1] <= delta & plik_CI[,2] >= delta))
  
  CI_width[,i] = c(mean(standard_CI[,2] - standard_CI[,1]),
                   mean(exact_CI[,2] - exact_CI[,1]),
                   mean(RCI[,2] - RCI[,1]),
                   mean(adjusted_asym_CI[,2] - adjusted_asym_CI[,1], na.rm = TRUE),
                   mean(bootstrap_uncond_CI[,2] - bootstrap_uncond_CI[,1]),
                   mean(exact_cond_CI[,2] - exact_cond_CI[,1]),
                   mean(r_cond_CI[,2] - r_cond_CI[,1], na.rm = TRUE),
                   mean(clik_CI[,2] - clik_CI[,1]),
                   mean(plik_CI[,2] - plik_CI[,1]))

  CI_width_se[,i] = c(sd(standard_CI[,2] - standard_CI[,1]),
                      sd(exact_CI[,2] - exact_CI[,1]),
                      sd(RCI[,2] - RCI[,1]),
                      sd(adjusted_asym_CI[,2] - adjusted_asym_CI[,1], na.rm = TRUE),
                      sd(bootstrap_uncond_CI[,2] - bootstrap_uncond_CI[,1]),
                      sd(exact_cond_CI[,2] - exact_cond_CI[,1]),
                      sd(r_cond_CI[,2] - r_cond_CI[,1], na.rm = TRUE),
                      sd(clik_CI[,2] - clik_CI[,1]),
                      sd(plik_CI[,2] - plik_CI[,1]))
  
  
  ### Means for trials that stop early
  
  CI_coverage_stop_early[,i] = c(mean(standard_CI[stop1==1,1] <= delta & standard_CI[stop1==1,2] >= delta),
                                 mean(exact_CI[stop1==1,1] <= delta & exact_CI[stop1==1,2] >= delta),
                                 mean(RCI[stop1==1,1] <= delta & RCI[stop1==1,2] >= delta),
                                 mean(adjusted_asym_CI[stop1==1,1] <= delta & adjusted_asym_CI[stop1==1,2] >= delta),
                                 mean(bootstrap_uncond_CI[stop1==1,1] <= delta & bootstrap_uncond_CI[stop1==1,2] >= delta),
                                 mean(exact_cond_CI[stop1==1,1] <= delta & exact_cond_CI[stop1==1,2] >= delta),
                                 sum(r_cond_CI[stop1==1,1] <= delta & r_cond_CI[stop1==1,2] >= delta, na.rm  = T)/sum(stop1),
                                 mean(clik_CI[stop1==1,1] <= delta & clik_CI[stop1==1,2] >= delta),
                                 mean(plik_CI[stop1==1,1] <= delta & plik_CI[stop1==1,2] >= delta))
  
  CI_width_stop_early[,i] = c(mean(standard_CI[stop1==1,2] - standard_CI[stop1==1,1]),
                              mean(exact_CI[stop1==1,2] - exact_CI[stop1==1,1]),
                              mean(RCI[stop1==1,2] - RCI[stop1==1,1]),
                              mean(adjusted_asym_CI[stop1==1,2] - adjusted_asym_CI[stop1==1,1]),
                              mean(bootstrap_uncond_CI[stop1==1,2] - bootstrap_uncond_CI[stop1==1,1]),
                              mean(exact_cond_CI[stop1==1,2] - exact_cond_CI[stop1==1,1]),
                              mean(r_cond_CI[stop1==1,2] - r_cond_CI[stop1==1,1], na.rm = TRUE),
                              mean(clik_CI[stop1==1,2] - clik_CI[stop1==1,1]),
                              mean(plik_CI[stop1==1,2] - plik_CI[stop1==1,1]))
  
  CI_width_se_stop_early[,i] = c(sd(standard_CI[stop1==1,2] - standard_CI[stop1==1,1]),
                                 sd(exact_CI[stop1==1,2] - exact_CI[stop1==1,1]),
                                 sd(RCI[stop1==1,2] - RCI[stop1==1,1]),
                                 sd(adjusted_asym_CI[stop1==1,2] - adjusted_asym_CI[stop1==1,1]),
                                 sd(bootstrap_uncond_CI[stop1==1,2] - bootstrap_uncond_CI[stop1==1,1]),
                                 sd(exact_cond_CI[stop1==1,2] - exact_cond_CI[stop1==1,1]),
                                 sd(r_cond_CI[stop1==1,2] - r_cond_CI[stop1==1,1], na.rm = TRUE),
                                 sd(clik_CI[stop1==1,2] - clik_CI[stop1==1,1]),
                                 sd(plik_CI[stop1==1,2] - plik_CI[stop1==1,1]))
  
  ### Means for trials that continue
  
  CI_coverage_continue[,i] = c(mean(standard_CI[stop1==0,1] <= delta &
                                      standard_CI[stop1==0,2] >= delta),
                               mean(exact_CI[stop1==0,1] <= delta &
                                      exact_CI[stop1==0,2] >= delta),
                               mean(RCI[stop1==0,1] <= delta &
                                      RCI[stop1==0,2] >= delta),
                               sum(adjusted_asym_CI[stop1==0,1] <= delta &
                                      adjusted_asym_CI[stop1==0,2] >= delta, na.rm = T)/sum(stop1==0),
                               mean(bootstrap_uncond_CI[stop1==0,1] <= delta &
                                      bootstrap_uncond_CI[stop1==0,2] >= delta),
                               sum(exact_cond_CI[stop1==0,1] <= delta &
                                      exact_cond_CI[stop1==0,2] >= delta, na.rm = TRUE)/sum(stop1==0),
                               sum(r_cond_CI[stop1==0,1] <= delta &
                                      r_cond_CI[stop1==0,2] >= delta, na.rm = TRUE)/sum(stop1==0),
                               mean(clik_CI[stop1==0,1] <= delta &
                                      clik_CI[stop1==0,2] >= delta))
  
  
  CI_width_continue[,i] =  c(mean(standard_CI[stop1==0,2] - standard_CI[stop1==0,1]),
                             mean(exact_CI[stop1==0,2] - exact_CI[stop1==0,1]),
                             mean(RCI[stop1==0,2] - RCI[stop1==0,1]),
                             mean(adjusted_asym_CI[stop1==0,2] - adjusted_asym_CI[stop1==0,1], na.rm = T),
                             mean(bootstrap_uncond_CI[stop1==0,2] - bootstrap_uncond_CI[stop1==0,1]),
                             mean(exact_cond_CI[stop1==0,2] - exact_cond_CI[stop1==0,1], na.rm = TRUE),
                             mean(r_cond_CI[stop1==0,2] - r_cond_CI[stop1==0,1], na.rm = TRUE),
                             mean(clik_CI[stop1==0,2] - clik_CI[stop1==0,1]))
  
  CI_width_se_continue[,i] = c(sd(standard_CI[stop1==0,2] - standard_CI[stop1==0,1]),
                               sd(exact_CI[stop1==0,2] - exact_CI[stop1==0,1]),
                               sd(RCI[stop1==0,2] - RCI[stop1==0,1]),
                               sd(adjusted_asym_CI[stop1==0,2] - adjusted_asym_CI[stop1==0,1], na.rm = T),
                               sd(bootstrap_uncond_CI[stop1==0,2] - bootstrap_uncond_CI[stop1==0,1]),
                               sd(exact_cond_CI[stop1==0,2] - exact_cond_CI[stop1==0,1], na.rm = TRUE),
                               sd(r_cond_CI[stop1==0,2] - r_cond_CI[stop1==0,1], na.rm = TRUE),
                               sd(clik_CI[stop1==0,2] - clik_CI[stop1==0,1]))
}

Sys.time()
stopCluster(cl)


save(phat_CE_vec, stop1_vec, CI_coverage, CI_coverage_stop_early,
     CI_coverage_continue, CI_width, CI_width_se, CI_width_stop_early,
     CI_width_se_stop_early, CI_width_continue, CI_width_se_continue,
     file = 'sim_results_cov_plot.Rdata')


#####################
library(ggplot2)
library(tidyr)
library(dplyr)


source('MUSEC_trial.R')
load('sim_results_cov_plot.Rdata')
load('sim_results_Wald.Rdata')

alpha = 0.05
N = 10^5

### Probability of early stopping

df <- data.frame(phat_CE_vec, stop1_vec)

ggplot(df, aes(x = phat_CE_vec, y = stop1_vec)) +
  geom_line(color = 'blue') +
  scale_x_continuous(name = expression(p[CE])) +
  scale_y_continuous(name = 'Probability of early stopping', limits = c(0, 1)) + 
  theme_minimal()

cbbPalette = c("#000000", "#D55E00", "#009E73",
               "#56B4E9", "#A020F0", "#CC79A7", "#CC79A7",
               "#F0E442", "#999999")

leg_label = c('Wald test (standard)', 'Unconditional final', 
              'RCI', 'Adjusted asymptotic', 'Bootstrap',
              'Conditional final', 'Restricted conditional',
              'Conditional likelihood', 'Penalised likelihood')


### Unconditional coverage

CI_coverage[1,] = CI_coverage_Wald

# Monte Carlo error
MC_error = matrix(nrow = 9, ncol = length(phat_CE_vec))

for(i in 1:9){
  MC_error[i,] = qnorm(1-alpha/2)*sqrt(CI_coverage[i,]*(1 - CI_coverage[i,])/N)
}

MC_error[1,] = MC_error[1,]/sqrt(10)

# Create plots
CI_coverage_df <- data.frame(phat_CE_vec, 
                        lower_bound1 = CI_coverage[1,] - MC_error[1,],
                        upper_bound1 = CI_coverage[1,] + MC_error[1,],
                        lower_bound2 = CI_coverage[2,] - MC_error[2,],
                        upper_bound2 = CI_coverage[2,] + MC_error[2,],
                        lower_bound3 = CI_coverage[3,] - MC_error[3,],
                        upper_bound3 = CI_coverage[3,] + MC_error[3,],
                        lower_bound4 = CI_coverage[4,] - MC_error[4,],
                        upper_bound4 = CI_coverage[4,] + MC_error[4,],
                        lower_bound5 = CI_coverage[5,] - MC_error[5,],
                        upper_bound5 = CI_coverage[5,] + MC_error[5,],
                        lower_bound6 = CI_coverage[6,] - MC_error[6,],
                        upper_bound6 = CI_coverage[6,] + MC_error[6,],
                        lower_bound7 = CI_coverage[7,] - MC_error[7,],
                        upper_bound7 = CI_coverage[7,] + MC_error[7,],
                        lower_bound8 = CI_coverage[8,] - MC_error[8,],
                        upper_bound8 = CI_coverage[8,] + MC_error[8,],
                        lower_bound9 = CI_coverage[9,] - MC_error[9,],
                        upper_bound9 = CI_coverage[9,] + MC_error[9,])


CI_coverage_long <- CI_coverage_df %>%
  pivot_longer(
    cols = -phat_CE_vec,
    names_to = c(".value", "group"),
    names_pattern = "(lower|upper)_bound(\\d+)"
  ) %>%
  mutate(mean_value = (lower + upper) / 2) %>%
  mutate(group = factor(group, levels = seq_along(leg_label), labels = leg_label))

p1 = ggplot(CI_coverage_long, aes(x = phat_CE_vec, y = mean_value, color = group, linetype = group, group = group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.75) +
  scale_color_manual(values = cbbPalette, name = "CI method") +
  scale_fill_manual(values = cbbPalette, name = "CI method") +
  scale_linetype_manual(values = c(1,rep(1,5),2,1,1), name = "CI method") +
  geom_hline(yintercept = 0.95, color = "red", linetype = 2, linewidth = 1) +
  scale_y_continuous(limits = c(0.9, 1), breaks = seq(0.90, 1.00, by = 0.02)) + 
  labs(x = expression(p[CE]), y = "Coverage") +
  theme_minimal() +
  theme(
    legend.position = "bottom", # Specify position
    legend.box = "horizontal"  # Needed for multi-row layout sometimes
  ) +
  guides(color = guide_legend(nrow = 2, title.position = "top", 
                              title.hjust = 0.5))

### Unconditional CI width

CI_width[1,] = CI_width_Wald

CI_width_df <- as.data.frame(t(CI_width)) 
colnames(CI_width_df) <- leg_label 
CI_width_df$phat_CE <- phat_CE_vec 

CI_width_df <- CI_width_df %>%
  pivot_longer(cols = -phat_CE, 
               names_to = "Row", 
               values_to = "CI_width") %>%
  mutate(Row = factor(Row, levels = leg_label))  # Set factor levels to match leg_label order

p2 = ggplot(CI_width_df, aes(x = phat_CE, y = CI_width, color = Row, linetype = Row, group = Row)) +
  geom_line(linewidth = 0.75) +
  scale_color_manual(values = cbbPalette, name = "CI method") +
  scale_linetype_manual(values = c(1,rep(1,5),2,1,1), name = "CI method") +
  labs(x = expression(p[CE]), y = "CI width") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, by = 0.1))


### Coverage conditional on stopping at stage 1

CI_coverage_stop_early[1,] = CI_coverage_stop_early_Wald


# Monte Carlo error
MC_error = matrix(nrow = 9, ncol = length(phat_CE_vec))

for(i in 1:9){
  MC_error[i,] = qnorm(1-alpha/2)*sqrt(CI_coverage_stop_early[i,]*(1 - CI_coverage_stop_early[i,])/(N*stop1_vec))
}

MC_error[1,] = MC_error[1,]/sqrt(10)

# Create plots
CI_coverage_stop_early_df <- data.frame(phat_CE_vec, 
                             lower_bound1 = CI_coverage_stop_early[1,] - MC_error[1,],
                             upper_bound1 = CI_coverage_stop_early[1,] + MC_error[1,],
                             lower_bound2 = CI_coverage_stop_early[2,] - MC_error[2,],
                             upper_bound2 = CI_coverage_stop_early[2,] + MC_error[2,],
                             lower_bound3 = CI_coverage_stop_early[3,] - MC_error[3,],
                             upper_bound3 = CI_coverage_stop_early[3,] + MC_error[3,],
                             lower_bound4 = CI_coverage_stop_early[4,] - MC_error[4,],
                             upper_bound4 = CI_coverage_stop_early[4,] + MC_error[4,],
                             lower_bound5 = CI_coverage_stop_early[5,] - MC_error[5,],
                             upper_bound5 = CI_coverage_stop_early[5,] + MC_error[5,],
                             lower_bound6 = CI_coverage_stop_early[6,] - MC_error[6,],
                             upper_bound6 = CI_coverage_stop_early[6,] + MC_error[6,],
                             lower_bound7 = CI_coverage_stop_early[7,] - MC_error[7,],
                             upper_bound7 = CI_coverage_stop_early[7,] + MC_error[7,],
                             lower_bound8 = CI_coverage_stop_early[8,] - MC_error[8,],
                             upper_bound8 = CI_coverage_stop_early[8,] + MC_error[8,],
                             lower_bound9 = CI_coverage_stop_early[9,] - MC_error[9,],
                             upper_bound9 = CI_coverage_stop_early[9,] + MC_error[9,])


CI_coverage_stop_early_long <- CI_coverage_stop_early_df %>%
  pivot_longer(
    cols = -phat_CE_vec,
    names_to = c(".value", "group"),
    names_pattern = "(lower|upper)_bound(\\d+)"
  ) %>%
  mutate(mean_value = (lower + upper) / 2) %>%
  mutate(group = factor(group, levels = seq_along(leg_label), labels = leg_label))

p3 = ggplot(CI_coverage_stop_early_long, aes(x = phat_CE_vec, y = mean_value, color = group, linetype = group, group = group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.75) +
  scale_color_manual(values = cbbPalette, name = "CI method") +
  scale_fill_manual(values = cbbPalette, name = "CI method") +
  scale_linetype_manual(values = c(1,rep(1,5),2,1,1), name = "CI method") +
  geom_hline(yintercept = 0.95, color = "red", linetype = 2, linewidth = 1) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) + 
  labs(x = expression(p[CE]), y = "Coverage") +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(fill = NA)),
         fill = "none",
         linetype = guide_legend(override.aes = list(fill = NA)))



### CI width conditional on stopping at stage 1

CI_width_stop_early[1,] = CI_width_stop_early_Wald

CI_width_stop_early_df <- as.data.frame(t(CI_width_stop_early)) 
colnames(CI_width_stop_early_df) <- leg_label 
CI_width_stop_early_df$phat_CE <- phat_CE_vec 

CI_width_stop_early_df <- CI_width_stop_early_df %>%
  pivot_longer(cols = -phat_CE, 
               names_to = "Row", 
               values_to = "CI_width") %>%
  mutate(Row = factor(Row, levels = leg_label))  # Set factor levels to match leg_label order

p4 = ggplot(CI_width_stop_early_df, aes(x = phat_CE, y = CI_width, color = Row, linetype = Row, group = Row)) +
  geom_line(linewidth = 0.75) +
  scale_color_manual(values = cbbPalette, name = "CI method") +
  scale_linetype_manual(values = c(1,rep(1,5),2,1,1), name = "CI method") +
  labs(x = expression(p[CE]), y = "CI width") +
  theme_minimal() +
  scale_y_continuous(limits = c(0.1, 1.2), breaks = seq(0.1, 1.2, by = 0.1))


### Coverage conditional on continuing to stage 2
leg_label = c('Wald test (standard)', 'Unconditional exact', 
              'RCI', 'Adjusted asymptotic', 'Bootstrap',
              'Conditional exact', 'Restricted conditional',
              'Conditional likelihood')


CI_coverage_continue[1,] = CI_coverage_continue_Wald


# Monte Carlo error
MC_error = matrix(nrow = 8, ncol = length(phat_CE_vec))

for(i in 1:8){
  MC_error[i,] = qnorm(1-alpha/2)*sqrt(CI_coverage_continue[i,]*
                                         (1 - CI_coverage_continue[i,])/(N*(1-stop1_vec)))
}

MC_error[1,] = MC_error[1,]/sqrt(10)

# Create plots
CI_coverage_continue_df <- data.frame(phat_CE_vec, 
                             lower_bound1 = CI_coverage_continue[1,] - MC_error[1,],
                             upper_bound1 = CI_coverage_continue[1,] + MC_error[1,],
                             lower_bound2 = CI_coverage_continue[2,] - MC_error[2,],
                             upper_bound2 = CI_coverage_continue[2,] + MC_error[2,],
                             lower_bound3 = CI_coverage_continue[3,] - MC_error[3,],
                             upper_bound3 = CI_coverage_continue[3,] + MC_error[3,],
                             lower_bound4 = CI_coverage_continue[4,] - MC_error[4,],
                             upper_bound4 = CI_coverage_continue[4,] + MC_error[4,],
                             lower_bound5 = CI_coverage_continue[5,] - MC_error[5,],
                             upper_bound5 = CI_coverage_continue[5,] + MC_error[5,],
                             lower_bound6 = CI_coverage_continue[6,] - MC_error[6,],
                             upper_bound6 = CI_coverage_continue[6,] + MC_error[6,],
                             lower_bound7 = CI_coverage_continue[7,] - MC_error[7,],
                             upper_bound7 = CI_coverage_continue[7,] + MC_error[7,],
                             lower_bound8 = CI_coverage_continue[8,] - MC_error[8,],
                             upper_bound8 = CI_coverage_continue[8,] + MC_error[8,])


CI_coverage_continue_long <- CI_coverage_continue_df %>%
  pivot_longer(
    cols = -phat_CE_vec,
    names_to = c(".value", "group"),
    names_pattern = "(lower|upper)_bound(\\d+)"
  ) %>%
  mutate(mean_value = (lower + upper) / 2) %>%
  mutate(group = factor(group, levels = seq_along(leg_label), labels = leg_label))

p5 = ggplot(CI_coverage_continue_long, aes(x = phat_CE_vec, y = mean_value, color = group, linetype = group, group = group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.75) +
  scale_color_manual(values = cbbPalette, name = "CI method") +
  scale_fill_manual(values = cbbPalette, name = "CI method") +
  scale_linetype_manual(values = c(1,rep(1,5),2,1,1), name = "CI method") +
  geom_hline(yintercept = 0.95, color = "red", linetype = 2, linewidth = 1) +
  scale_y_continuous(limits = c(0.55, 1), breaks = seq(0.5, 1, by = 0.1)) + 
  labs(x = expression(p[CE]), y = "Coverage") +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(fill = NA)),
         fill = "none",
         linetype = guide_legend(override.aes = list(fill = NA)))



### CI width conditional on continuing to stage 2

CI_width_continue[1,] = CI_width_continue_Wald

CI_width_continue_df <- as.data.frame(t(CI_width_continue)) 
colnames(CI_width_continue_df) <- leg_label 
CI_width_continue_df$phat_CE <- phat_CE_vec 

CI_width_continue_df <- CI_width_continue_df %>%
  pivot_longer(cols = -phat_CE, 
               names_to = "Row", 
               values_to = "CI_width") %>%
  mutate(Row = factor(Row, levels = leg_label))  # Set factor levels to match leg_label order

p6 = ggplot(CI_width_continue_df, aes(x = phat_CE, y = CI_width, color = Row, linetype = Row, group = Row)) +
  geom_line(linewidth = 0.75) +
  scale_color_manual(values = cbbPalette, name = "CI method") +
  scale_linetype_manual(values = c(1,rep(1,5),2,1,1), name = "CI method") +
  labs(x = expression(p[CE]), y = "CI width") +
  theme_minimal() +
  scale_y_continuous(limits = c(0.1, 0.35), breaks = seq(0, 0.3, by = 0.1))




##################################################################
# Single plot

library(cowplot)

shared_legend <- get_legend(p1)
shared_legend <- cowplot::get_plot_component(p1, 'guide-box-bottom', return_all = FALSE)

row1 <- plot_grid(p1 + theme(legend.position = "none"), # Remove legend from p1 now
                  p2  + theme(legend.position = "none"),
                  ncol = 2,
                  align = 'vh',
                  axis = 'tblr') # Adjust alignment as needed

row2 <- plot_grid(p3  + theme(legend.position = "none"),
                  p4  + theme(legend.position = "none"),
                  ncol = 2,
                  align = 'vh',
                  axis = 'tblr')

row3 <- plot_grid(p5  + theme(legend.position = "none"),
                  p6 + theme(legend.position = "none"),
                  ncol = 2,
                  align = 'vh',
                  axis = 'tblr')


title_row1 <- ggdraw() +
  draw_label("Unconditional", fontface = 'bold', hjust = 0.5, size=12)

title_row2 <- ggdraw() +
  draw_label("Conditional on stopping at stage 1", fontface = 'bold', hjust = 0.5, size=12)

title_row3 <- ggdraw() +
  draw_label("Conditional on continuing to stage 2", fontface = 'bold', hjust = 0.5, size=12)

plots_with_titles <- plot_grid(
  title_row1, row1,
  title_row2, row2,
  title_row3, row3,
  ncol = 1,
  rel_heights = c(0.1, 1, 0.1, 1, 0.1, 1)
)

final_plot <- plot_grid(
  plots_with_titles,
  shared_legend,
  ncol = 1,
  rel_heights = c(1, 0.1) 
)

print(final_plot)
