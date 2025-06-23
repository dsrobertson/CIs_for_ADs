### Code for Section 4 (Tables 3-6 and A.3.1 - A.3.2): Performance of 
### different CI methods for fixed values of p_CE.

################################################################################

source('MUSEC_trial.R')

N = 10^5  # Number of bootstrap replicates
M = 10^4  # Number of inner bootstrap replicates

### Set up parallel computing

library(doParallel)
library(doRNG)
library(data.table)

mc = 32; cl = makeCluster(mc)

registerDoParallel(cl)
registerDoRNG(7)

### Run simulation

#phat_CE = phat_CE + 0.08 # Uncomment out for Tables 6, A.3.1 and A.3.2

Sys.time()

CI_results = foreach (i = 1:N,
                      .combine = function(x,y)rbindlist(list(x,y)),
                      .packages = c('mvtnorm','cubature'),
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
                          
                          
                          # Standard
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
                          
                          exact_cond_CI = c(clb,cub)
                          
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
                          
                          clik_CI = c(quantile(cMLEsim, probs = 0.025),
                                      quantile(cMLEsim, probs = 0.975))
                          
                          # clik_CI = c(max(quantile(cMLEsim, probs = 0.025),-1),
                          #             min(quantile(cMLEsim, probs = 0.975),1))

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
                          
                          # Standard
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

Sys.time()

stopCluster(cl)

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

save(rej, s0, s1, standard_CI, exact_CI, RCI, adjusted_asym_CI, bootstrap_uncond_CI,
     exact_cond_CI, r_cond_CI, clik_CI, plik_CI, file = 'sim_results.Rdata')


################

stop1 = rej[,1]
rej2 = rej[,2]

delta = phat_CE - phat_P

mean(stop1)


### Overall (unconditional)

c(sum(is.na(adjusted_asym_CI[,1])), sum(is.na(adjusted_asym_CI[,2])))
c(sum(is.na(exact_cond_CI[,1])), sum(is.na(exact_cond_CI[,2])))
c(sum(is.na(r_cond_CI[,1])), sum(is.na(r_cond_CI[,2])))

CI_names = c('Wald test (standard)', 'Unconditional exact', 
             'RCI', 'Adjusted asymptotic', 'Bootstrap', 'Conditional exact',
             'Restricted conditional','Conditional likelihood',
             'Penalised likelihood')

CI_coverage = c(mean(standard_CI[,1] <= delta & standard_CI[,2] >= delta),
                mean(exact_CI[,1] <= delta & exact_CI[,2] >= delta),
                mean(RCI[,1] <= delta & RCI[,2] >= delta),
                sum(adjusted_asym_CI[,1] <= delta & adjusted_asym_CI[,2] >= delta, na.rm = TRUE)/10^5,
                mean(bootstrap_uncond_CI[,1] <= delta & bootstrap_uncond_CI[,2] >= delta),
                mean(exact_cond_CI[,1] <= delta & exact_cond_CI[,2] >= delta),
                sum(r_cond_CI[,1] <= delta & r_cond_CI[,2] >= delta, na.rm = TRUE)/10^5,
                mean(clik_CI[,1] <= delta & clik_CI[,2] >= delta),
                mean(plik_CI[,1] <= delta & plik_CI[,2] >= delta))

names(CI_coverage) = CI_names
print(round(CI_coverage,3))

CI_consistency = c(mean((standard_CI[,1] > 0 & rowSums(rej)) | (standard_CI[,1] < 0 & !(rowSums(rej)))),
                   mean((exact_CI[,1] > 0 & rowSums(rej)) | (exact_CI[,1] < 0 & !(rowSums(rej)))),
                   mean((RCI[,1] > 0 & rowSums(rej)) | (RCI[,1] < 0 & !(rowSums(rej)))),
                   sum((adjusted_asym_CI[,1] > 0 & rowSums(rej)) | (adjusted_asym_CI[,1] < 0 & !(rowSums(rej))), na.rm = TRUE)/10^5,
                   mean((bootstrap_uncond_CI[,1] > 0 & rowSums(rej)) | (bootstrap_uncond_CI[,1] < 0 & !(rowSums(rej)))),
                   mean((exact_cond_CI[,1] > 0 & rowSums(rej)) | (exact_cond_CI[,1] < 0 & !(rowSums(rej)))),
                   sum((r_cond_CI[,1] > 0 & rowSums(rej)) | (r_cond_CI[,1] < 0 & !(rowSums(rej))), na.rm = TRUE)/10^5,
                   mean((clik_CI[,1] > 0 & rowSums(rej)) | (clik_CI[,1] < 0 & !(rowSums(rej)))),
                   mean((plik_CI[,1] > 0 & rowSums(rej)) | (plik_CI[,1] < 0 & !(rowSums(rej)))))

names(CI_consistency) = CI_names
print(round(CI_consistency,3))

CI_coverage_lower = c(mean(standard_CI[,1] > delta),
                      mean(exact_CI[,1] > delta),
                      mean(RCI[,1] > delta ),
                      sum(adjusted_asym_CI[,1] > delta, na.rm = TRUE)/10^5,
                      mean(bootstrap_uncond_CI[,1] > delta),
                      mean(exact_cond_CI[,1] > delta),
                      sum(r_cond_CI[,1] > delta, na.rm = TRUE)/10^5,
                      mean(clik_CI[,1] > delta),
                      mean(plik_CI[,1] > delta))

CI_coverage_upper = c(mean(standard_CI[,2] < delta),
                      mean(exact_CI[,2] < delta),
                      mean(RCI[,2] < delta),
                      sum(adjusted_asym_CI[,2] < delta, na.rm = TRUE)/10^5,
                      mean(bootstrap_uncond_CI[,2] < delta),
                      mean(exact_cond_CI[,2] < delta),
                      sum(r_cond_CI[,2] < delta, na.rm = TRUE)/10^5,
                      mean(clik_CI[,2] < delta),
                      mean(plik_CI[,2] < delta))

CI_coverage = cbind(CI_coverage_lower, CI_coverage_upper)
rownames(CI_coverage) = CI_names
print(round(CI_coverage,3))

CI_width = c(mean(standard_CI[,2] - standard_CI[,1]),
             mean(exact_CI[,2] - exact_CI[,1]),
             mean(RCI[,2] - RCI[,1]),
             mean(adjusted_asym_CI[,2] - adjusted_asym_CI[,1], na.rm = TRUE),
             mean(bootstrap_uncond_CI[,2] - bootstrap_uncond_CI[,1]),
             mean(exact_cond_CI[,2] - exact_cond_CI[,1]),
             mean(r_cond_CI[,2] - r_cond_CI[,1], na.rm = TRUE),
             mean(clik_CI[,2] - clik_CI[,1]),
             mean(plik_CI[,2] - plik_CI[,1]))

CI_width_se = c(sd(standard_CI[,2] - standard_CI[,1]),
                sd(exact_CI[,2] - exact_CI[,1]),
                sd(RCI[,2] - RCI[,1]),
                sd(adjusted_asym_CI[,2] - adjusted_asym_CI[,1], na.rm = TRUE),
                sd(bootstrap_uncond_CI[,2] - bootstrap_uncond_CI[,1]),
                sd(exact_cond_CI[,2] - exact_cond_CI[,1]),
                sd(r_cond_CI[,2] - r_cond_CI[,1], na.rm = TRUE),
                sd(clik_CI[,2] - clik_CI[,1]),
                sd(plik_CI[,2] - plik_CI[,1]))

CI_width = cbind(CI_width, CI_width_se)
rownames(CI_width) = CI_names
print(round(CI_width,3))


### Means for trials that stop early

CI_coverage_stop_early = c(mean(standard_CI[stop1==1,1] <= delta & standard_CI[stop1==1,2] >= delta),
                           mean(exact_CI[stop1==1,1] <= delta & exact_CI[stop1==1,2] >= delta),
                           mean(RCI[stop1==1,1] <= delta & RCI[stop1==1,2] >= delta),
                           mean(adjusted_asym_CI[stop1==1,1] <= delta & adjusted_asym_CI[stop1==1,2] >= delta),
                           mean(bootstrap_uncond_CI[stop1==1,1] <= delta & bootstrap_uncond_CI[stop1==1,2] >= delta),
                           mean(exact_cond_CI[stop1==1,1] <= delta & exact_cond_CI[stop1==1,2] >= delta),
                           sum(r_cond_CI[stop1==1,1] <= delta & r_cond_CI[stop1==1,2] >= delta, na.rm = TRUE)/sum(stop1),
                           mean(clik_CI[stop1==1,1] <= delta & clik_CI[stop1==1,2] >= delta),
                           mean(plik_CI[stop1==1,1] <= delta & plik_CI[stop1==1,2] >= delta))

names(CI_coverage_stop_early) = CI_names
print(round(CI_coverage_stop_early,3))

CI_consistency_stop_early = c(mean(standard_CI[stop1==1,1] > 0),
                              mean(exact_CI[stop1==1,1] > 0),
                              mean(RCI[stop1==1,1] > 0),
                              mean(adjusted_asym_CI[stop1==1,1] > 0),
                              mean(bootstrap_uncond_CI[stop1==1,1] > 0),
                              mean(exact_cond_CI[stop1==1,1] > 0),
                              sum(r_cond_CI[stop1==1,1] > 0, na.rm = T)/sum(stop1),
                              mean(clik_CI[stop1==1,1] > 0),
                              mean(plik_CI[stop1==1,1] > 0))

names(CI_consistency_stop_early) = CI_names
print(round(CI_consistency_stop_early,3))


CI_coverage_lower_stop_early = c(mean(standard_CI[stop1 == 1,1] > delta ),
                                 mean(exact_CI[stop1 == 1,1] > delta),
                                 mean(RCI[stop1 == 1,1] > delta ),
                                 mean(adjusted_asym_CI[stop1 == 1,1] > delta),
                                 mean(bootstrap_uncond_CI[stop1 == 1,1] > delta),
                                 mean(exact_cond_CI[stop1 == 1,1] > delta),
                                 sum(r_cond_CI[stop1 == 1,1] > delta, na.rm = T)/sum(stop1),
                                 mean(clik_CI[stop1 == 1,1] > delta),
                                 mean(plik_CI[stop1 == 1,1] > delta))

CI_coverage_upper_stop_early = c(mean(standard_CI[stop1 == 1,2] < delta),
                                 mean(exact_CI[stop1 == 1,2] < delta),
                                 mean(RCI[stop1 == 1,2] < delta),
                                 mean(adjusted_asym_CI[stop1 == 1,2] < delta),
                                 mean(bootstrap_uncond_CI[stop1 == 1,2] < delta),
                                 mean(exact_cond_CI[stop1 == 1,2] < delta),
                                 sum(r_cond_CI[stop1 == 1,2] < delta, na.rm = T)/sum(stop1),
                                 mean(clik_CI[stop1 == 1,2] < delta),
                                 mean(plik_CI[stop1 == 1,2] < delta))

CI_coverage_stop_early = cbind(CI_coverage_lower_stop_early, CI_coverage_upper_stop_early)
rownames(CI_coverage_stop_early) = CI_names
print(round(CI_coverage_stop_early,3))


CI_width_stop_early = c(mean(standard_CI[stop1==1,2] - standard_CI[stop1==1,1]),
                        mean(exact_CI[stop1==1,2] - exact_CI[stop1==1,1]),
                        mean(RCI[stop1==1,2] - RCI[stop1==1,1]),
                        mean(adjusted_asym_CI[stop1==1,2] - adjusted_asym_CI[stop1==1,1]),
                        mean(bootstrap_uncond_CI[stop1==1,2] - bootstrap_uncond_CI[stop1==1,1]),
                        mean(exact_cond_CI[stop1==1,2] - exact_cond_CI[stop1==1,1]),
                        mean(r_cond_CI[stop1==1,2] - r_cond_CI[stop1==1,1], na.rm = T),
                        mean(clik_CI[stop1==1,2] - clik_CI[stop1==1,1]),
                        mean(plik_CI[stop1==1,2] - plik_CI[stop1==1,1]))

CI_width_se_stop_early = c(sd(standard_CI[stop1==1,2] - standard_CI[stop1==1,1]),
                           sd(exact_CI[stop1==1,2] - exact_CI[stop1==1,1]),
                           sd(RCI[stop1==1,2] - RCI[stop1==1,1]),
                           sd(adjusted_asym_CI[stop1==1,2] - adjusted_asym_CI[stop1==1,1]),
                           sd(bootstrap_uncond_CI[stop1==1,2] - bootstrap_uncond_CI[stop1==1,1]),
                           sd(exact_cond_CI[stop1==1,2] - exact_cond_CI[stop1==1,1]),
                           sd(r_cond_CI[stop1==1,2] - r_cond_CI[stop1==1,1], na.rm = T),
                           sd(clik_CI[stop1==1,2] - clik_CI[stop1==1,1]),
                           sd(plik_CI[stop1==1,2] - plik_CI[stop1==1,1]))

CI_width_stop_early = cbind(CI_width_stop_early, CI_width_se_stop_early)
rownames(CI_width_stop_early) = CI_names
print(round(CI_width_stop_early,3))


### Means for trials that continue

c(sum(is.na(adjusted_asym_CI[stop1==0,1])), sum(is.na(adjusted_asym_CI[stop1==0,2])))
c(sum(is.na(r_cond_CI[stop1==0,1])), sum(is.na(r_cond_CI[stop1==0,2])))

CI_names = c('Wald test (standard)', 'Unconditional exact', 
             'RCI', 'Adjusted asymptotic', 'Bootstrap', 'Conditional exact',
             'Restricted conditional','Conditional likelihood')


CI_coverage_continue = c(mean(standard_CI[stop1==0,1] <= delta &
                                standard_CI[stop1==0,2] >= delta),
                         mean(exact_CI[stop1==0,1] <= delta &
                                exact_CI[stop1==0,2] >= delta),
                         mean(RCI[stop1==0,1] <= delta &
                                RCI[stop1==0,2] >= delta),
                         sum(adjusted_asym_CI[stop1==0,1] <= delta &
                                adjusted_asym_CI[stop1==0,2] >= delta, na.rm = T)/sum(stop1==0),
                         mean(bootstrap_uncond_CI[stop1==0,1] <= delta &
                                bootstrap_uncond_CI[stop1==0,2] >= delta),
                         mean(exact_cond_CI[stop1==0,1] <= delta &
                                exact_cond_CI[stop1==0,2] >= delta),
                         sum(r_cond_CI[stop1==0,1] <= delta &
                                r_cond_CI[stop1==0,2] >= delta, na.rm = T)/sum(stop1 == 0),
                         mean(clik_CI[stop1==0,1] <= delta &
                                clik_CI[stop1==0,2] >= delta))

names(CI_coverage_continue) = CI_names
print(round(CI_coverage_continue,3))

CI_consistency_continue = c(mean((standard_CI[stop1==0,1] > 0 & rej2[stop1==0]) | (standard_CI[stop1==0,1] < 0 & !(rej2[stop1==0]))),
                            mean((exact_CI[stop1==0,1] > 0 & rej2[stop1==0]) | (exact_CI[stop1==0,1] < 0 & !(rej2[stop1==0]))),
                            mean((RCI[stop1==0,1] > 0 & rej2[stop1==0] | (RCI[stop1==0,1] < 0 & !(rej2[stop1==0])))),
                            sum((adjusted_asym_CI[stop1==0,1] > 0 & rej2[stop1==0]) | (adjusted_asym_CI[stop1==0,1] < 0 & !(rej2[stop1==0])), na.rm = T)/sum(stop1==0),
                            mean((bootstrap_uncond_CI[stop1==0,1] > 0 & rej2[stop1==0]) | (bootstrap_uncond_CI[stop1==0,1] < 0 & !(rej2[stop1==0]))),
                            mean((exact_cond_CI[stop1==0,1] > 0 & rej2[stop1==0]) | (exact_cond_CI[stop1==0,1] < 0 & !(rej2[stop1==0]))),
                            sum((r_cond_CI[stop1==0,1] > 0 & rej2[stop1==0]) | (r_cond_CI[stop1==0,1] < 0 & !(rej2[stop1==0])), na.rm=T)/sum(stop1==0),
                            mean((clik_CI[stop1==0,1] > 0 & rej2[stop1==0]) | (clik_CI[stop1==0,1] < 0 & !(rej2[stop1==0]))))

names(CI_consistency_continue) = CI_names
print(round(CI_consistency_continue,3))

CI_coverage_lower_continue = c(mean(standard_CI[stop1==0,1] > delta ),
                               mean(exact_CI[stop1==0,1] > delta),
                               mean(RCI[stop1==0,1] > delta ),
                               sum(adjusted_asym_CI[stop1==0,1] > delta, na.rm = T)/sum(stop1==0),
                               mean(bootstrap_uncond_CI[stop1==0,1] > delta),
                               mean(exact_cond_CI[stop1==0,1] > delta),
                               sum(r_cond_CI[stop1==0,1] > delta, na.rm=T)/sum(stop1==0),
                               mean(clik_CI[stop1==0,1] > delta))

CI_coverage_upper_continue = c(mean(standard_CI[stop1==0,2] < delta),
                               mean(exact_CI[stop1==0,2] < delta),
                               mean(RCI[stop1==0,2] < delta),
                               sum(adjusted_asym_CI[stop1==0,2] < delta, na.rm = T)/sum(stop1==0),
                               mean(bootstrap_uncond_CI[stop1==0,2] < delta),
                               mean(exact_cond_CI[stop1==0,2] < delta),
                               sum(r_cond_CI[stop1==0,2] < delta, na.rm=T)/sum(stop1==0),
                               mean(clik_CI[stop1==0,2] < delta))

CI_coverage_continue = cbind(CI_coverage_lower_continue, CI_coverage_upper_continue )
rownames(CI_coverage_continue) = CI_names
print(round(CI_coverage_continue,3))


CI_width_continue = c(mean(standard_CI[stop1==0,2] - standard_CI[stop1==0,1]),
                      mean(exact_CI[stop1==0,2] - exact_CI[stop1==0,1]),
                      mean(RCI[stop1==0,2] - RCI[stop1==0,1]),
                      mean(adjusted_asym_CI[stop1==0,2] - adjusted_asym_CI[stop1==0,1], na.rm = T),
                      mean(bootstrap_uncond_CI[stop1==0,2] - bootstrap_uncond_CI[stop1==0,1]),
                      mean(exact_cond_CI[stop1==0,2] - exact_cond_CI[stop1==0,1]),
                      mean(r_cond_CI[stop1==0,2] - r_cond_CI[stop1==0,1], na.rm=T),
                      mean(clik_CI[stop1==0,2] - clik_CI[stop1==0,1]))

CI_width_se_continue = c(sd(standard_CI[stop1==0,2] - standard_CI[stop1==0,1]),
                         sd(exact_CI[stop1==0,2] - exact_CI[stop1==0,1]),
                         sd(RCI[stop1==0,2] - RCI[stop1==0,1]),
                         sd(adjusted_asym_CI[stop1==0,2] - adjusted_asym_CI[stop1==0,1], na.rm = T),
                         sd(bootstrap_uncond_CI[stop1==0,2] - bootstrap_uncond_CI[stop1==0,1]),
                         sd(exact_cond_CI[stop1==0,2] - exact_cond_CI[stop1==0,1]),
                         sd(r_cond_CI[stop1==0,2] - r_cond_CI[stop1==0,1], na.rm=T),
                         sd(clik_CI[stop1==0,2] - clik_CI[stop1==0,1]))

CI_width_continue = cbind(CI_width_continue, CI_width_se_continue)
rownames(CI_width_continue) = CI_names
print(round(CI_width_continue,3))
