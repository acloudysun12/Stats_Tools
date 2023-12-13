library(mvtnorm)
library(MASS)
library(dplyr)
source("C:/tfs_projects/Research/Exam IC Research/avs/create_ib_data.R")

?optim

data_sim = rnorm(200, mean = 1, sd = 1)

neg_log_likeli_data_calc = function(vec_param){
  vec_likeli = log(dnorm(data_sim, mean = vec_param[1], sd = exp(vec_param[2])))
  return(-sum(vec_likeli))
}

a = optim(par = c(1.5, 1), fn = neg_log_likeli_data_calc, method = "BFGS", lower = -Inf, upper = Inf) # try different starts
a$par

neg_log_likeli_data_calc(c(a$par[1], exp(a$par[2]))) 
neg_log_likeli_data_calc(c(1, exp(0)))


A = matrix(round(runif(4, -3, 3), 1), nrow = 2)
A[2,1] = 0
cov = t(A)%*%A
eigen(cov)
data_sim2 = mvrnorm(n = 200, mu = c(0, 3), Sigma = cov)

dmvnorm(x = c(-.1, 3.2), mean = c(0, 3), sigma = as.matrix(cov))


neg_log_likeli_mult_dim_calc = function(vec_param){
  U_sim = matrix(c(vec_param[3], 0, vec_param[4:5]), nrow = 2, byrow = FALSE)
  cov_sim = t(U_sim) %*% U_sim
  vec_likeli = dmvnorm(x = data_sim2, mean = vec_param[1:2], sigma = cov_sim, log = TRUE)
  return(-sum(vec_likeli))
}

b = optim(par = runif(5, -5, 5), fn = neg_log_likeli_mult_dim_calc, method = "BFGS",  control = list(trace = 1, ndeps = rep(0.0001, 5)))
b$par
b$value





### ANALYSIS with SIMULATED DATA

mix1 = rnorm(500, mean = 0, sd = 0.5)
mix2 = rnorm(500, mean = 4, sd = 0.5)
binom_rv = rbinom(n = 500, 1, prob = 0.5)
data_sim = ifelse(binom_rv == 1, mix1, mix2)
hist(data_sim)

mu = c(0, 1)
log_sd = c(0, 0)
lik = cbind(dnorm(data_sim, mean = mu[1], sd = exp(log_sd[1])), dnorm(data_sim, mean = mu[2], sd = exp(log_sd[2])))
log_lik = log(lik)
pz_x_t = lik/rowSums(lik)
# ln_pz_x_t = log(pz_x_t)
# ln_pz_t = log(cbind(rep(0.5, 500), rep(0.5, 500)))
# neg_weighted_ll = -pz_x_t*(log_lik + ln_pz_t - ln_pz_x_t)
# sum(neg_weighted_ll)

vec_param = c(0, 0.5, 4, 0.5, 0.5) # mean1, log_sd1, mean2, log_sd2, prop


# E STEP
calc_E = function(vec_param, dat){
  wlik1 = dnorm(dat, mean = vec_param[1], sd = exp(vec_param[2]))*vec_param[5]
  wlik2 = dnorm(dat, mean = vec_param[3], sd = exp(vec_param[4]))*(1-vec_param[5])
  wlik = cbind(wlik1, wlik2)
  new_pz_x_t = wlik/rowSums(wlik) # this soasqlves for E step, feeds into M step
  return(new_pz_x_t)
}
pz_x_t = calc_E(vec_param)



pz_x_t_iter = pz_x_t
dat = data_sim
# M Step analytically
calc_M = function(pz_x_t_iter, dat){
  mu1 = sum(dat*pz_x_t_iter[,1])/sum(pz_x_t_iter[,1])
  mu2 = sum(dat*pz_x_t_iter[,2])/sum(pz_x_t_iter[,2])
  sig_sq_1 = sum((dat - mu1)^2*pz_x_t_iter[,1])/sum(pz_x_t_iter[,1])
  sig_sq_2 = sum((dat - mu2)^2*pz_x_t_iter[,2])/sum(pz_x_t_iter[,2])
  prop_mix = sum(pz_x_t_iter[,1])/length(dat)
  new_param_vec = c(mu1, log(sqrt(sig_sq_1)), mu2, log(sqrt(sig_sq_2)), prop_mix)
  return(new_param_vec)
}



# With analytical M step first
vec_params_mtx = matrix(nrow = 0, ncol = 5)
weighted_lls = numeric(0)
vec_param = c(0, 0.1, 1, 0.5, 0.5) # mean1, log_sd1, mean2, log_sd2, prop
vec_params_mtx = rbind(vec_params_mtx, vec_param)
new_ll_weighted = compute_ll(data_sim, vec_param)
weighted_lls = c(weighted_lls, new_ll_weighted)
for (i in 1:5){
  pz_x_t = calc_E(vec_param, data_sim) # first calculate posteriors (E step)
  vec_param = calc_M(pz_x_t_iter = pz_x_t, dat = data_sim)
  new_ll_weighted = compute_ll(data_sim, vec_param)
  vec_params_mtx = rbind(vec_params_mtx, vec_param)
  weighted_lls = c(weighted_lls, new_ll_weighted)
}
weighted_lls
vec_params_mtx


# M STEP numerical
# E STEP
calc_E2 = function(vec_param, dat){
  prob = exp(vec_param[5])/(1+exp(vec_param[5]))
  wlik1 = dnorm(dat, mean = vec_param[1], sd = exp(vec_param[2]))*prob
  wlik2 = dnorm(dat, mean = vec_param[3], sd = exp(vec_param[4]))*(1-prob)
  wlik = cbind(wlik1, wlik2)
  new_pz_x_t = wlik/rowSums(wlik) # this soasqlves for E step, feeds into M step
  return(new_pz_x_t)
}

# Now with non-analytical M 
vec_params_mtx = matrix(nrow = 0, ncol = 5)
weighted_lls = numeric(0)
vec_param = c(0, 0.1, 1, 0.5, 0) # mean1, log_sd1, mean2, log_sd2, log-odds(prop)
vec_params_mtx = rbind(vec_params_mtx, vec_param)
weighted_lls = c(weighted_lls, compute_ll2(data_sim, vec_param))
for (i in 1:5){
  pz_x_t = calc_E2(vec_param, data_sim) # first calculate posteriors (E step)
  ln_pz_x_t = log(pz_x_t)
  mle_results = optim(par = vec_param, fn = neg_weighted_ll_calc, method = "BFGS", control = list(ndeps = rep(0.0001, 5))) # now find MLEs (M step)
  vec_param = mle_results$par
  new_ll_weighted = compute_ll2(data_sim, vec_param)
  vec_params_mtx = rbind(vec_params_mtx, vec_param)
  weighted_lls = c(weighted_lls, new_ll_weighted)
}
weighted_lls
vec_params_mtx

###### RESULTS TIE OUT!!!!! YES!


compute_ll = function(sample_data, vec_param){
  ll1_weighted = dnorm(sample_data, mean = vec_param[1], sd = exp(vec_param[2]))*vec_param[5]
  ll2_weighted = dnorm(sample_data, mean = vec_param[3], sd = exp(vec_param[4]))*(1-vec_param[5])
  return(sum(log(ll1_weighted + ll2_weighted)))
}

compute_ll2 = function(sample_data, vec_param){
  prob = exp(vec_param[5])/(1+exp(vec_param[5]))
  ll1_weighted = dnorm(sample_data, mean = vec_param[1], sd = exp(vec_param[2]))*prob
  ll2_weighted = dnorm(sample_data, mean = vec_param[3], sd = exp(vec_param[4]))*(1-prob)
  return(sum(log(ll1_weighted + ll2_weighted)))
}

neg_weighted_ll_calc = function(vec_param){
  prob = exp(vec_param[5])/(1+exp(vec_param[5])) # log-odds to prob
  lik1 = dnorm(data_sim, mean = vec_param[1], sd = exp(vec_param[2]))
  lik2 = dnorm(data_sim, mean = vec_param[3], sd = exp(vec_param[4]))
  lik = cbind(lik1, lik2)
  log_lik = log(lik)
  pz_t = cbind(rep(prob, length(data_sim)), rep(1 - prob, length(data_sim)))
  ln_pz_t = log(pz_t)
  neg_weighted_ll = -pz_x_t*(log_lik + ln_pz_t)
  return(sum(neg_weighted_ll))
}



#### REAL DATA, ENC_RATIO
length(df_mdl_all$enc_ratio)
hist(df_mdl_all$enc_ratio, main = "Hist ENC Ratio (Entire pop)")

data_sim = df_mdl_all$enc_ratio
# Now with numerical calc M 
vec_params_mtx = matrix(nrow = 0, ncol = 5)
weighted_lls = numeric(0)
vec_param = c(0, 0.1, 1, 0.5, 0) # mean1, log_sd1, mean2, log_sd2, log-odds(prop)
vec_params_mtx = rbind(vec_params_mtx, vec_param)
weighted_lls = c(weighted_lls, compute_ll2(data_sim, vec_param))
for (i in 1:25){
  pz_x_t = calc_E2(vec_param, data_sim) # first calculate posteriors (E step)
  ln_pz_x_t = log(pz_x_t)
  mle_results = optim(par = vec_param, fn = neg_weighted_ll_calc, method = "BFGS", control = list(ndeps = rep(0.0001, 5))) # now find MLEs (M step)
  vec_param = mle_results$par
  new_ll_weighted = compute_ll2(data_sim, vec_param)
  vec_params_mtx = rbind(vec_params_mtx, vec_param)
  weighted_lls = c(weighted_lls, new_ll_weighted)
}
weighted_lls
vec_params_mtx[25,]

cbind(mix = round(pz_x_t[,1], 0), tgt = ib_pop_all$TARGET, ct = 1) %>% as.data.frame() %>% group_by(mix, tgt) %>% dplyr::summarize(sum_ct = sum(ct))
pred_obj0 = prediction(predictions = pz_x_t[,2], labels = ib_pop_all$TARGET)
perf_obj0 = performance(prediction.obj = pred_obj1, measure = "auc")
auc_val0 = as.numeric(perf_obj1@y.values)
auc_val0


#### ENC_RATIO but with 3 mixtures!!!
# data_sim = df_mdl_all$enc_ratio[ib_pop_all$window_idx <= 2017.25]
data_sim = df_mdl_all$enc_ratio
data_sample = data_sim
calc_E3 = function(vec_param, dat){
  prob1 = exp(vec_param[7])/(1+exp(vec_param[7]))
  prob2 = min(1 - prob1, exp(vec_param[8])/(1+exp(vec_param[8]))) # bound by the highest it can be = 1 - prob1
  prob3 = max(1 - prob1 - prob2, 0)
  wlik1 = dnorm(dat, mean = vec_param[1], sd = exp(vec_param[2]))*prob1
  wlik2 = dnorm(dat, mean = vec_param[3], sd = exp(vec_param[4]))*prob2
  wlik3 = dnorm(dat, mean = vec_param[5], sd = exp(vec_param[6]))*prob3
  wlik = cbind(wlik1, wlik2, wlik3)
  new_pz_x_t = wlik/rowSums(wlik) # this soasqlves for E step, feeds into M step
  return(new_pz_x_t)
}
compute_ll3 = function(dat, vec_param){
  prob1 = exp(vec_param[7])/(1+exp(vec_param[7]))
  prob2 = min(1 - prob1, exp(vec_param[8])/(1+exp(vec_param[8]))) # bound by the highest it can be = 1 - prob1
  prob3 = max(1 - prob1 - prob2, 0)
  wlik1 = dnorm(dat, mean = vec_param[1], sd = exp(vec_param[2]))*prob1
  wlik2 = dnorm(dat, mean = vec_param[3], sd = exp(vec_param[4]))*prob2
  wlik3 = dnorm(dat, mean = vec_param[5], sd = exp(vec_param[6]))*prob3
  return(sum(log(wlik1 + wlik2 + wlik3)))
}
neg_weighted_ll_calc = function(vec_param){
  prob1 = exp(vec_param[7])/(1+exp(vec_param[7]))
  prob2 = min(1 - prob1, exp(vec_param[8])/(1+exp(vec_param[8]))) # bound by the highest it can be = 1 - prob1
  prob3 = max(1 - prob1 - prob2, 0)
  lik1 = dnorm(data_sim, mean = vec_param[1], sd = exp(vec_param[2]))
  lik2 = dnorm(data_sim, mean = vec_param[3], sd = exp(vec_param[4]))
  lik3 = dnorm(data_sim, mean = vec_param[5], sd = exp(vec_param[6]))
  lik = cbind(lik1, lik2, lik3)
  log_lik = log(lik)
  pz_t = cbind(rep(prob1, length(data_sim)), rep(prob2, length(data_sim)), rep(prob3, length(data_sim)))
  ln_pz_t = log(pz_t)
  neg_weighted_ll = -pz_x_t*(log_lik + ln_pz_t)
  return(sum(neg_weighted_ll))
}
# Now with numerical calc M 
vec_params_mtx = matrix(nrow = 0, ncol = 8)
weighted_lls = numeric(0)
vec_param = c(0, 0.1, 0.5, 0.3, 1, 0.5, -1, 0) # mean1, log_sd1, mean2, log_sd2, mean3, log_sd3, log-odds(prop1), log-odds(prop2)
vec_params_mtx = rbind(vec_params_mtx, vec_param)
weighted_lls = c(weighted_lls, compute_ll3(data_sim, vec_param))
for (i in 1:25){
  pz_x_t = calc_E3(vec_param, data_sim) # first calculate posteriors (E step)
  ln_pz_x_t = log(pz_x_t)
  mle_results = optim(par = vec_param, fn = neg_weighted_ll_calc, method = "BFGS", control = list(ndeps = rep(0.0001, 8))) # now find MLEs (M step)
  vec_param = mle_results$par
  new_ll_weighted = compute_ll3(data_sim, vec_param)
  vec_params_mtx = rbind(vec_params_mtx, vec_param)
  weighted_lls = c(weighted_lls, new_ll_weighted)
}
weighted_lls
vec_params_mtx[25,]

argmax = function(vec){return(as.numeric(which(vec == max(vec))))}
cbind(mix = apply(X = pz_x_t, 1, argmax), tgt = ib_pop_all$TARGET[ib_pop_all$window_idx <= 2017.25], ct = 1) %>% as.data.frame() %>% group_by(mix, tgt) %>% dplyr::summarize(sum_ct = sum(ct))

pred_obj1 = prediction(predictions = pz_x_t[,1], labels = ib_pop_all$TARGET)
perf_obj1 = performance(prediction.obj = pred_obj1, measure = "auc")
auc_val1 = as.numeric(perf_obj1@y.values)
auc_val1

pz_x_t_val = calc_E3(vec_params_mtx[25,], ib_pop_all$enc_ratio[ib_pop_all$window_idx >= 2017.25]) # first calculate posteriors (E step)
pred_obj2 = prediction(predictions = pz_x_t_val[,1], labels = ib_pop_all$TARGET[ib_pop_all$window_idx >= 2017.25])
perf_obj2 = performance(prediction.obj = pred_obj2, measure = "auc")
auc_val2 = as.numeric(perf_obj2@y.values)
auc_val2





#### WHAT IF SOME BINARY VARIABLES INSTEAD?

# Consider data x as 3-d vector of random variables both ~Bern(p something)
xvec = cbind(df_mdl_all$bsa_other_flag, df_mdl_all$facts_inv_other_1yr_cnt_is_0, df_mdl_all$ib_rt_com_is_0)
sample_data = xvec
data_sim = xvec

colMeans(xvec)

compute_ll_binom = function(wliks){
  return(sum(log(wliks[,1] + wliks[,2])))
}

# function for optimization
neg_w_ll_calc_binom = function(vec_param){
  p11 = exp(vec_param[1])/(1+exp(vec_param[1]))
  p12 = exp(vec_param[2])/(1+exp(vec_param[2]))
  p13 = exp(vec_param[3])/(1+exp(vec_param[3]))
  p21 = exp(vec_param[4])/(1+exp(vec_param[4]))
  p22 = exp(vec_param[5])/(1+exp(vec_param[5]))
  p23 = exp(vec_param[6])/(1+exp(vec_param[6]))
  prob = exp(vec_param[7])/(1+exp(vec_param[7]))
  lik1 = p11^(sample_data[,1])*(1-p11)^(1-sample_data[,1])*p12^(sample_data[,2])*(1-p12)^(1-sample_data[,2])*p13^(sample_data[,3])*(1-p13)^(1-sample_data[,3])
  lik2 = p21^(sample_data[,1])*(1-p21)^(1-sample_data[,1])*p22^(sample_data[,2])*(1-p22)^(1-sample_data[,2])*p23^(sample_data[,3])*(1-p23)^(1-sample_data[,3])
  lik = cbind(lik1, lik2)
  log_lik = log(lik)
  pz_t = cbind(rep(prob, nrow(data_sim)), rep(1 - prob, nrow(data_sim)))
  ln_pz_t = log(pz_t)
  neg_w_ll_binom = -pz_x_t*(log_lik + ln_pz_t)
  return(sum(neg_w_ll_binom))
}

# E STEP
calc_E_binom = function(vec_param, dat){
  p11 = exp(vec_param[1])/(1+exp(vec_param[1]))
  p12 = exp(vec_param[2])/(1+exp(vec_param[2]))
  p13 = exp(vec_param[3])/(1+exp(vec_param[3]))
  p21 = exp(vec_param[4])/(1+exp(vec_param[4]))
  p22 = exp(vec_param[5])/(1+exp(vec_param[5]))
  p23 = exp(vec_param[6])/(1+exp(vec_param[6]))
  prob = exp(vec_param[7])/(1+exp(vec_param[7]))
  wlik1 = p11^(dat[,1])*(1-p11)^(1-dat[,1])*p12^(dat[,2])*(1-p12)^(1-dat[,2])*p13^(dat[,3])*(1-p13)^(1-dat[,3])*prob
  wlik2 = p21^(dat[,1])*(1-p21)^(1-dat[,1])*p22^(dat[,2])*(1-p22)^(1-dat[,2])*p23^(dat[,3])*(1-p23)^(1-dat[,3])*(1-prob)
  wlik = cbind(wlik1, wlik2)
  new_pz_x_t = wlik/rowSums(wlik) # this soasqlves for E step, feeds into M step
  return(list(new_pz_x_t = new_pz_x_t, new_wlik = wlik))
}

# Now with numerical calc M 
vec_params_mtx = matrix(nrow = 0, ncol = 7)
weighted_lls = numeric(0)
# vec_param = c(-1, -1, -1, 1, 1, 1, 0) # logit(p11), logit(p12), logit(p21), logit(p22), logit(prop)
vec_param = runif(7, -3, 3)
for (i in 1:10){
  pz_x_t = calc_E_binom(vec_param, data_sim)$new_pz_x_t # first calculate posteriors (E step)
  ln_pz_x_t = log(pz_x_t)
  mle_results = optim(par = vec_param, fn = neg_w_ll_calc_binom, method = "BFGS", control = list(ndeps = rep(0.0001, 7))) # now find MLEs (M step)
  vec_param = mle_results$par
  wlik = calc_E_binom(vec_param, data_sim)$new_wlik # with new EM step, record results
  new_ll_weighted = compute_ll_binom(wlik)
  vec_params_mtx = rbind(vec_params_mtx, vec_param)
  weighted_lls = c(weighted_lls, compute_ll_binom(calc_E_binom(vec_param, sample_data)$new_wlik))
}
weighted_lls
vec_params_mtx

cbind(mix = round(pz_x_t[,2], 0)+1, tgt = ib_pop_all$TARGET, ct = 1) %>% as.data.frame() %>% group_by(mix, tgt) %>% dplyr::summarize(cts = sum(ct))






library(plotly)
df_mdl_all$tgt = ib_pop_all$TARGET
df_mdl_all$tgt <- as.factor(df_mdl_all$tgt)
plot(df_mdl_all$enc_ratio, df_mdl_all$facts_inv_fin_cnt_log, pch = 16, cex = 0.5)

p = plot_ly(df_mdl_all, x = ~enc_ratio, y = ~current_prn_relation_cnt_sum, z = ~facts_inv_fin_cnt_log, color = ~tgt,  
            marker = list(symbol = 'circle', sizemode = 'diameter'), colors = c('#BF382A', '#0C4B8E')) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'enc'),
                      yaxis = list(title = 'prn'),
                      zaxis = list(title = 'fin')))




# Entropy
0.05*log(0.05) + 0.95*log(0.95)
0.5*log(0.5) + 0.5*log(0.5)
0.49*log(0.49) + 0.51*log(0.51)
0.01*log(0.01)
0.00001*log(0.00001) + 0.9999*log(0.9999)


