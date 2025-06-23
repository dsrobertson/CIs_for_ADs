### Code for Tables A.2.2 and A.2.4: calculation of confidence intervals
### for the MUSEC trial at the (approximate) stopping boundaries

################################################################################


original = gsDesign(k=2, test.type = 1, beta = 0.2, sfu='OF', n.fix = 400)
e1 = original$upper$bound[1]
e2 = original$upper$bound[2]

alpha = 0.05

################################################################################
############ Stop at (approximate) boundary at stage 1

### Interim data

n1_1 = 101   # Number of subjects randomized to CE arm
n0_1 = 97    # Number of subjects randomized to placebo

s1_1 = seq(0,101, by = 1)
s0_1 = seq(0,97, by = 1)

combine_s1_s0 = expand.grid(s1 = s1_1, s0 = s0_1)

Z1 = rep(NA, nrow(combine_s1_s0))


for(i in 1:nrow(combine_s1_s0)){
  
  s1_1 = combine_s1_s0[i,1]
  s0_1 = combine_s1_s0[i,2]
  
  ptilde_1 = (s0_1 + s1_1)/(n0_1 + n1_1)
  Ihat_1 = 1/(ptilde_1*(1-ptilde_1)*(1/n0_1 + 1/n1_1))
  Z1[i] = (s1_1/n1_1 - s0_1/n0_1)*sqrt(Ihat_1)   # Wald test statistic
  
}

j = which.min((Z1 - e1)^2)

# combine_s1_s0[j,]
# print(c(Z1[j], e1))

s1_1 = combine_s1_s0[j,1]
s0_1 = combine_s1_s0[j,2]

ptilde_1 = (s0_1 + s1_1)/(n0_1 + n1_1)
Ihat_1 = 1/(ptilde_1*(1-ptilde_1)*(1/n0_1 + 1/n1_1))
Z1 = (s1_1/n1_1 - s0_1/n0_1)*sqrt(Ihat_1) 

##### Calculate point estimates and CIs

# Standard
phat_CE = s1_1/n1_1
phat_P = s0_1/n0_1

MLE = phat_CE - phat_P

MLE_se = sqrt(phat_CE*(1-phat_CE)/n1_1 + 
                phat_P*(1-phat_P)/n0_1)

standard_CI = c(MLE - qnorm(1-alpha/2)*MLE_se,
                MLE + qnorm(1-alpha/2)*MLE_se)


# Exact
lb = (qnorm(alpha/2)+Z1)/sqrt(Ihat_1)
ub = (qnorm(1-alpha/2)+Z1)/sqrt(Ihat_1)

MUE = (qnorm(0.5)+Z1)/sqrt(Ihat_1)

exact_CI = c(lb, ub)


# RCI

ck_alpha = e2

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

M = 10^6 

MLEsim = rep(NA, M)

s1_1_boot = rbinom(M, size = n1_1, prob = phat_CE)
s0_1_boot = rbinom(M, size = n0_1, prob = phat_P)

ptilde_1_boot = (s0_1_boot + s1_1_boot)/(n0_1 + n1_1)
Ihat_1_boot = 1/(ptilde_1_boot*(1-ptilde_1_boot)*(1/n0_1 + 1/n1_1))
Z1_boot = (s1_1_boot/n1_1 - s0_1_boot/n0_1)*sqrt(Ihat_1_boot) 

n1_2 = 143   
n0_2 = 134  

s1_2_boot = s1_1_boot + rbinom(M, size = n1_2-n1_1, prob = phat_CE)
s0_2_boot = s0_1_boot + rbinom(M, size = n0_2-n0_1, prob = phat_P)

ptilde_2_boot = (s0_2_boot + s1_2_boot)/(n0_2 + n1_2)
Ihat_2_boot = 1/(ptilde_2_boot*(1-ptilde_2_boot)*(1/n0_2 + 1/n1_2))
Z2_boot = (s1_2_boot/n1_2 - s0_2_boot/n0_2)*sqrt(Ihat_2_boot) 

MLEsim[Z1_boot > e1] = Z1_boot[Z1_boot > e1]/sqrt(Ihat_1_boot[Z1_boot > e1])
MLEsim[Z1_boot <= e1] = Z2_boot[Z1_boot <= e1]/sqrt(Ihat_2_boot[Z1_boot <= e1])

bootstrap_uncond = mean(MLEsim)

bootstrap_uncond_CI = c(quantile(MLEsim, probs = 0.025),
                        quantile(MLEsim, probs = 0.975))




### Randomisation-based CI

set.seed(7)

stage1_outcomes = sample(c(rep(1, s1_1 + s0_1), rep(0, n0_1 + n1_1 - s1_1 - s0_1)))

stage1_treatment = c(rep(0, n0_1), rep(1, n1_1))

n1 = n0_1 + n1_1


N = 10^6   # Bootstrap replicates

more_extreme_vec = s0_1_boot = s1_1_boot = rep(0, N)

dqset.seed(7)

for(i in 1:N){
  
  stage1_assignment = dqsample(stage1_treatment, size = n1)
  
  s0_1_boot[i] = sum(stage1_outcomes[stage1_assignment == 0])
  s1_1_boot[i] = sum(stage1_outcomes[stage1_assignment == 1])
  
}


ptilde_1_boot = (s0_1_boot + s1_1_boot)/(n0_1 + n1_1)
Ihat_1_boot = 1/(ptilde_1_boot*(1-ptilde_1_boot)*(1/n0_1 + 1/n1_1))
Z1_boot = (s1_1_boot/n1_1 - s0_1_boot/n0_1)*sqrt(Ihat_1_boot) 


more_extreme_vec[Z1_boot > Z1] = 1

Sys.time()

p_adjusted = sum(more_extreme_vec)/N

E_adjusted = MLE_se*qnorm(1 - p_adjusted)

random_CI = c(E_adjusted - MLE_se*qnorm(1-alpha/2),
              E_adjusted + MLE_se*qnorm(1-alpha/2))


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

cMUE = optimise(cb_optim, c(standard_CI[2]-100, standard_CI[2]+1),
                z = Z1, I1 = Ihat_1,
                u = e1, q = 0.5)$minimum

exact_cond_CI = c(clb, cub)

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


cMLE = optimise(plog_lik, c(MLE-100,MLE+1),
                lambda = 1,
                z = Z1, I1 = Ihat_1,
                u = e1, maximum = TRUE)$maximum

pLE = function(lambda, z, u, I1){
  
  optimise(plog_lik,
           c(standard_CI[1]-1, standard_CI[2]+1),
           lambda = lambda, z = z,
           u = u, I1 = I1, maximum = TRUE)$maximum
  
}

lambda_star = uniroot(pLE, c(0,1), z = e1, u = e1,
                      I1 = Ihat_1)$root


pMLE = optimise(plog_lik, c(MLE-1,MLE+1),
                lambda = lambda_star, 
                z = Z1, I1 = Ihat_1,
                u = e1, maximum = TRUE)$maximum



cMLEsim = pMLEsim = rep(NA, M)
boot_count = 0

while(boot_count < M){
  
  s1_1_boot = rbinom(1, size = n1_1, prob = phat_CE)
  s0_1_boot = rbinom(1, size = n0_1, prob = phat_P)
  
  ptilde_1_boot = (s0_1_boot + s1_1_boot)/(n0_1 + n1_1)
  Ihat_1_boot = 1/(ptilde_1_boot*(1-ptilde_1_boot)*(1/n0_1 + 1/n1_1))
  Z1_boot = (s1_1_boot/n1_1 - s0_1_boot/n0_1)*sqrt(Ihat_1_boot)
  
  if(Z1_boot > e1){
    
    boot_count = boot_count + 1
    
    cMLEsim[boot_count] = optimise(plog_lik, c(MLE-100,MLE+1),
                                   lambda = 1,
                                   z = Z1_boot, I1 = Ihat_1_boot,
                                   u = e1, maximum = TRUE)$maximum
    
    lambda_star = uniroot(pLE, c(0,1), z = e1, u = e1,
                          I1 = Ihat_1_boot)$root
    
    
    pMLEsim[boot_count] = optimise(plog_lik, c(MLE-1,MLE+1),
                                   lambda = lambda_star, 
                                   z = Z1_boot, I1 = Ihat_1_boot,
                                   u = e1, maximum = TRUE)$maximum
    
  }
}

clik_CI = c(quantile(cMLEsim, probs = 0.025),
            quantile(cMLEsim, probs = 0.975))

plik_CI = c(quantile(pMLEsim, probs = 0.025),
            quantile(pMLEsim, probs = 0.975))



###


point.estimates = rbind(MLE, MUE, NA, adjusted_asym_est, 
                        bootstrap_uncond, E_adjusted, cMUE, cMUE, cMLE, pMLE)

point.estimates = round(point.estimates, digits = 3)


CIs_raw = rbind(standard_CI, exact_CI, RCI, adjusted_asym_CI, 
                bootstrap_uncond_CI, random_CI, exact_cond_CI, r_cond_CI,
                clik_CI, plik_CI)

CIs = round(CIs_raw, digits = 3)

CI_names = rownames(CIs) = rownames(point.estimates) = 
  c('Wald test (standard)', 'Unconditional exact', 
    'RCI', 'Adjusted asymptotic', 'Unconditional bootstrap',
    'Randomisation-based', 'Conditional exact', 'Restricted conditional',
    'Conditional likelihood', 'Penalised likelihood')


print(CIs)

print(point.estimates)

print(round(CIs_raw[,2] - CIs_raw[,1], digits = 3))


################################################################################
############ Stop at (approximate) boundary at stage 2

### Final data

n1_2 = 143   # Number of subjects randomized to CE arm
n0_2 = 134   # Number of subjects randomized to placebo

s1_2 = seq(0,143, by = 1)
s0_2 = seq(0,134, by = 1)

combine_s1_s0 = expand.grid(s1 = s1_2, s0 = s0_2)

Z2 = rep(NA, nrow(combine_s1_s0))


for(i in 1:nrow(combine_s1_s0)){
  
  s1_2 = combine_s1_s0[i,1]
  s0_2 = combine_s1_s0[i,2]
  
  ptilde_2 = (s0_2 + s1_2)/(n0_2 + n1_2)
  Ihat_2 = 1/(ptilde_2*(1-ptilde_2)*(1/n0_2 + 1/n1_2))
  Z2[i] = (s1_2/n1_2 - s0_2/n0_2)*sqrt(Ihat_2)   # Wald test statistic
  
}

j = which.min((Z2 - e2)^2)

# combine_s1_s0[j,]
# print(c(Z2[j], e2))


### Trial data

s1_2 = unname(combine_s1_s0[j,1])
s0_2 = unname(combine_s1_s0[j,2])

ptilde_2 = (s0_2 + s1_2)/(n0_2 + n1_2)
Ihat_2 = 1/(ptilde_2*(1-ptilde_2)*(1/n0_2 + 1/n1_2))
Z2 = (s1_2/n1_2 - s0_2/n0_2)*sqrt(Ihat_2) 


s1_1 = 45 
s0_1 = 30

# c(s1_2 - s1_1, n1_2 - n1_1)
# c(s0_2 - s0_1, n0_2 - n0_1)

ptilde_1 = (s0_1 + s1_1)/(n0_1 + n1_1)
Ihat_1 = 1/(ptilde_1*(1-ptilde_1)*(1/n0_1 + 1/n1_1))
Z1 = (s1_1/n1_1 - s0_1/n0_1)*sqrt(Ihat_1)   # Wald test statistic

# c(Z1, e1)


##### Calculate point estimates and CIs

alpha = 0.05

### Standard/naive (Wald)

phat_CE = s1_2/n1_2
phat_P = s0_2/n0_2

MLE = phat_CE - phat_P
MLE_se = sqrt(phat_CE*(1-phat_CE)/n1_2 + phat_P*(1-phat_P)/n0_2)

standard_CI = c(MLE - qnorm(1-alpha/2)*MLE_se,
                MLE + qnorm(1-alpha/2)*MLE_se)


### Exact (unconditional)

library(mvtnorm)

cov_matrix = matrix(c(1, sqrt(Ihat_1/Ihat_2), sqrt(Ihat_1/Ihat_2), 1),
                    nrow = 2, byrow = TRUE)

e1 = original$upper$bound[1]
f1 = -Inf

pval_fn = function(delta) {
  
  mu = delta*sqrt(c(Ihat_1,Ihat_2))
  
  pnorm(delta*sqrt(Ihat_1) - e1) + 
    pmvnorm(lower=c(f1, Z2), upper=c(e1, Inf),
            mean = mu, sigma = cov_matrix)[1]
}

obj_fn = function(p, q){
  L = (pval_fn(delta=p) - q)^2
}

lb = optimize(obj_fn, q = alpha/2, interval = c(0,1), tol = 1e-10)$minimum
ub = optimize(obj_fn, q = 1 - alpha/2, interval = c(0,1), tol = 1e-10)$minimum

exact_CI = c(lb, ub)   # Exact (final) confidence interval
MUE = optimize(obj_fn, q = 0.5, interval = c(0,1), tol = 1e-10)$minimum   # Median unbiased estimate



### RCI following Jennison and Turnbull book

ck_alpha = e2 = gsDesign(k=2, test.type = 2, alpha = 0.025, sfu='OF')$upper$bound[2]

## Stage 1
c(Z1/sqrt(Ihat_1) - ck_alpha/sqrt(Ihat_1)*sqrt(2),
  Z1/sqrt(Ihat_1) + ck_alpha/sqrt(Ihat_1)*sqrt(2))

# Stage 2
RCI = c(Z2/sqrt(Ihat_2) - ck_alpha/sqrt(Ihat_2),
        Z2/sqrt(Ihat_2) + ck_alpha/sqrt(Ihat_2))



### Parametric bootstrap (unconditional)

N = 10^6   # Number of bootstrap replicates

MLEsim = rep(NA, N)

set.seed(7)

s1_1_boot = rbinom(N, size = n1_1, prob = phat_CE)
s0_1_boot = rbinom(N, size = n0_1, prob = phat_P)

ptilde_1_boot = (s0_1_boot + s1_1_boot)/(n0_1 + n1_1)
Ihat_1_boot = 1/(ptilde_1_boot*(1-ptilde_1_boot)*(1/n0_1 + 1/n1_1))
Z1_boot = (s1_1_boot/n1_1 - s0_1_boot/n0_1)*sqrt(Ihat_1_boot) 

s1_2_boot = s1_1_boot + rbinom(N, size = n1_2-n1_1, prob = phat_CE)
s0_2_boot = s0_1_boot + rbinom(N, size = n0_2-n0_1, prob = phat_P)

ptilde_2_boot = (s0_2_boot + s1_2_boot)/(n0_2 + n1_2)
Ihat_2_boot = 1/(ptilde_2_boot*(1-ptilde_2_boot)*(1/n0_2 + 1/n1_2))
Z2_boot = (s1_2_boot/n1_2 - s0_2_boot/n0_2)*sqrt(Ihat_2_boot) 

MLEsim[Z1_boot > e1] = Z1_boot[Z1_boot > e1]/sqrt(Ihat_1_boot[Z1_boot > e1])
MLEsim[Z1_boot <= e1] = Z2_boot[Z1_boot <= e1]/sqrt(Ihat_2_boot[Z1_boot <= e1])

bootstrap_uncond_CI = c(quantile(MLEsim, probs = 0.025),
                        quantile(MLEsim, probs = 0.975))

bootstrap_uncond = mean(MLEsim)



### Adjusted asymptotic (unconditional)

V1 = Ihat_1
V2 = Ihat_2

I2 = V2 - V1

theta = Z2/sqrt(Ihat_2)

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
  integrate(E_Z_int2, lower = -Inf, upper = Inf)$value

E_Z_sqrtV = (1/sqrt(V1))*integrate(E_Z_int1, lower = e1, upper = Inf)$value +
  (1/sqrt(V2))*integrate(E_Z_int2, lower = -Inf, upper = Inf)$value


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
  sqrt(V2)*integrate(E_V_int2, lower = -Inf, upper = Inf)$value

E_V = V1*integrate(E_V_int1, lower = e1, upper = Inf)$value +
  V2*integrate(E_V_int2, lower = -Inf, upper = Inf)$value

mu_adj = E_Z_sqrtV - theta*E_sqrtV

sigma_adj = sqrt(E_Z2_V - 2*theta*E_Z + theta^2*E_V - mu_adj^2)

adjusted_asym_est = theta - mu_adj/sqrt(V2)

adjusted_asym_CI = c(adjusted_asym_est - qnorm(1-alpha/2)*sigma_adj/sqrt(V2),
                     adjusted_asym_est + qnorm(1-alpha/2)*sigma_adj/sqrt(V2))



### Randomisation-based CI

set.seed(7)

stage1_outcomes = sample(c(rep(1, s1_1 + s0_1), rep(0, n0_1 + n1_1 - s1_1 - s0_1)))

stage2_outcomes = sample(c(rep(1, s1_2 - s1_1 + s0_2 - s0_1),
                           rep(0, (n0_2 - n0_1) - (s0_2 - s0_1) + (n1_2 - n1_1) - (s1_2 - s1_1))))

stage1_treatment = c(rep(0, n0_1), rep(1, n1_1))
stage2_treatment = c(rep(0, n0_2 - n0_1), rep(1, n1_2 - n1_1))

n1 = n0_1 + n1_1
n2 = n0_2 - n0_1 + n1_2 - n1_1

N = 10^6   # Bootstrap replicates

more_extreme_vec = s0_1_boot = s1_1_boot = s1_2_boot = s0_2_boot = rep(0, N)

dqset.seed(7)

for(i in 1:N){
  
  stage1_assignment = dqsample(stage1_treatment, size = n1)
  
  s0_1_boot[i] = sum(stage1_outcomes[stage1_assignment == 0])
  s1_1_boot[i] = sum(stage1_outcomes[stage1_assignment == 1])
  
  stage2_assignment = dqsample(stage2_treatment, size = n2)
  
  s1_2_boot[i] = s1_1_boot[i] + sum(stage2_outcomes[stage2_assignment == 1])
  s0_2_boot[i] = s0_1_boot[i] + sum(stage2_outcomes[stage2_assignment == 0])
  
}


ptilde_1_boot = (s0_1_boot + s1_1_boot)/(n0_1 + n1_1)
Ihat_1_boot = 1/(ptilde_1_boot*(1-ptilde_1_boot)*(1/n0_1 + 1/n1_1))
Z1_boot = (s1_1_boot/n1_1 - s0_1_boot/n0_1)*sqrt(Ihat_1_boot) 

ptilde_2_boot = (s0_2_boot + s1_2_boot)/(n0_2 + n1_2)
Ihat_2_boot = 1/(ptilde_2_boot*(1-ptilde_2_boot)*(1/n0_2 + 1/n1_2))
Z2_boot = (s1_2_boot/n1_2 - s0_2_boot/n0_2)*sqrt(Ihat_2_boot) 


more_extreme_vec[Z1_boot > e1] = 1
more_extreme_vec[Z1_boot <= e1 & Z2_boot > Z2] = 1

Sys.time()

p_adjusted = sum(more_extreme_vec)/N

E_adjusted = MLE_se*qnorm(1 - p_adjusted)

random_CI = c(E_adjusted - MLE_se*qnorm(1-alpha/2),
              E_adjusted + MLE_se*qnorm(1-alpha/2))


############################################################################
### Exact (conditional)

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


clb = optimise(cb_optim, c(-0.5,1), z = Z2, I1 = Ihat_1,
               I2 = Ihat_2, u = e1, q = 0.975)$minimum

cub = optimise(cb_optim, c(-0.5,1), z = Z2, I1 = Ihat_1,
               I2 = Ihat_2, u = e1, q = 0.025)$minimum


exact_cond_CI = c(clb, cub)

cMUE = optimise(cb_optim, c(-0.5,1), z = Z2, I1 = Ihat_1,
                I2 = Ihat_2, u = e1, q = 0.5)$minimum


## Restricted conditional

rub = (e1 - qnorm(0.025))/sqrt(Ihat_1)


r_cond_CI = c(clb, min(cub, rub))


## Conditional MLE

cMLE_obj = function(theta, z, I1, I2, u){
  
  (z/sqrt(I2) - theta + (sqrt(I1)*dnorm(u - theta*sqrt(I1)))/
     (I2*pnorm(u - theta*sqrt(I1))))^2
  
}


cMLE = optimise(cMLE_obj, c(MLE-1,MLE+1),
                z = Z2, I1 = Ihat_1, I2 = Ihat_2,
                u = e1)$minimum


N = 10^6   # Number of bootstrap replicates

cMLEsim = rep(NA, N)

boot_count = 0

set.seed(7)

while(boot_count < N){
  
  s1_1_boot = rbinom(1, size = n1_1, prob = phat_CE)
  s0_1_boot = rbinom(1, size = n0_1, prob = phat_P)
  
  ptilde_1_boot = (s0_1_boot + s1_1_boot)/(n0_1 + n1_1)
  Ihat_1_boot = 1/(ptilde_1_boot*(1-ptilde_1_boot)*(1/n0_1 + 1/n1_1))
  Z1_boot = (s1_1_boot/n1_1 - s0_1_boot/n0_1)*sqrt(Ihat_1_boot)
  
  if(Z1_boot <= e1){
    
    boot_count = boot_count + 1
    
    s1_2_boot = s1_1_boot + rbinom(1, size = n1_2-n1_1, prob = phat_CE)
    s0_2_boot = s0_1_boot + rbinom(1, size = n0_2-n0_1, prob = phat_P)
    
    ptilde_2_boot = (s0_2_boot + s1_2_boot)/(n0_2 + n1_2)
    Ihat_2_boot = 1/(ptilde_2_boot*(1-ptilde_2_boot)*(1/n0_2 + 1/n1_2))
    Z2_boot = (s1_2_boot/n1_2 - s0_2_boot/n0_2)*sqrt(Ihat_2_boot)
    
    cMLEsim[boot_count] = optimise(cMLE_obj,
                                   c(Z2_boot/sqrt(Ihat_2_boot)-1,Z2_boot/sqrt(Ihat_2_boot)+1),
                                   z = Z2_boot, I1 = Ihat_1_boot, I2 = Ihat_2_boot,
                                   u = e1)$minimum
  }
}

clik_CI = c(quantile(cMLEsim, probs = 0.025),
            quantile(cMLEsim, probs = 0.975))


############

point.estimates = rbind(MLE, MUE, NA, adjusted_asym_est, 
                        bootstrap_uncond, E_adjusted, cMUE, cMUE, cMLE)

point.estimates = round(point.estimates, digits = 3)

CIs_raw = rbind(standard_CI, exact_CI, RCI, adjusted_asym_CI, 
                bootstrap_uncond_CI, random_CI, exact_cond_CI, r_cond_CI, clik_CI)

CIs = round(CIs_raw, digits = 3)

CI_names = rownames(point.estimates) = rownames(CIs) = 
  c('Wald test (standard)', 'Unconditional exact', 
    'RCI', 'Adjusted asymptotic', 'Unconditional bootstrap',
    'Randomisation-based', 'Conditional exact', 'Restricted conditional',
    'Conditional likelihood')

print(point.estimates)

print(CIs)

print(round(CIs_raw[,2] - CIs_raw[,1], 3))

