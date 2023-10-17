###### sparse linear regression ######
library(MASS)
library(Matrix)
library(picasso)
library(glmnet)
library(doSNOW)
library(foreach)
library(abind)
library(lars)
library(parallel)
n_cores = detectCores()
cluster = makeCluster(25, outfile="")
registerDoSNOW(cluster)
source("/home/roysapta/IHT/IHT.R")
source("/home/roysapta/IHT/IHT_screening.R")

## parameters
p = 2000 
k = 0.9 # sample size parameter
#r = 10  # signal streangth paramter r = 1 + delta

q0 = 1 # number of spiky signlas = 0, 1, 2, 3
snr = 6

gamma = 10
nrep = 200

eta0 = 0.5 #gradient descent step size


s = 13
n = floor(p^k)



# signal and noise independent design
r_list = seq(6, 9.5, 0.5)
for (k1 in 1:length(r_list)){
  r = r_list[k1]
  tp = sqrt(r)
  betamin = sqrt(2* r *log(p)/n)
  a = betamin
  S_x = diag(p)
  
  
  sigma_e = 1
  
  ####### selecting lambda_lasso list ########
  #  data_i <- data_lm_newscheme(n = n, mu_x = rep(0,p), S_x = S_x,
  #                              s = s, betamin = betamin, sigma_e = sigma_e, tp = tp, q = q0)
  
  
  #  X = data_i$X
  #  Y = data_i$Y 
  #  mod = glmnet(x = X, y = Y, alpha = 1, standardize = F, intercept = F)
  #  lambda_list = mod$lambda
  
  ############
  
  df = foreach(i = 1:nrep, .packages = c("MASS", "Matrix", "picasso", "glmnet","foreach", "lars"), .combine = '+')%dopar%{
    set.seed(i)
    data_i <- data_lm_q_spiky_fixed_snr(n = n, mu_x = rep(0,p), S_x = S_x,
                                        s = s, betamin = betamin, sigma_e = sigma_e, snr = snr, q = q0)
    
    
    X = data_i$X
    Y = data_i$Y  
    beta = data_i$beta
    epsilon = data_i$epsilon
    supp_true = data_i$supp_true
    #s_budget_max = 2*s
    
    # lasso
    lasso_mod = glmnet(x=X, y = Y, standardize = F, intercept = F, nlambda = 500)
    O = 1*(lasso_mod$beta != 0); O1 = colSums(O)
    id = min(which(O1 == s))
    if(id < Inf){
      prop_recovery_lasso = prod(sort(supp_true) == sort(which(lasso_mod$beta[,id]!=0)))
    }
    else{
      prop_recovery_lasso = 0
    }
    cat("nrep = ", i , "lasso complete of inner iter=",k1,"\n")
    
    ## IHT + screening 
    iht_screening_result = iht_one_cv_screening(X = X, Y= Y, beta0 = rep(0,p), gamma = gamma, a = a, v = 1, maxiter = 1e3, prec = 1e-7, s_list = 1:(2*s), eta_list = eta0)
    supp_est  = sort(iht_screening_result$ alpha_vec_sort_idx[1:s])
    prop_recovery_iht_screening = prod(sort(supp_true) == supp_est)
    
    cat("nrep = ", i , "iht+screening complete of main iter", k1,"\n")
    
    ## marginal screening 
    ms_result = ms_ordered_select(X = X, Y = Y, supp_true = supp_true, scale_X = F, s_level = s)
    prop_recovery_ms = 1 * I(ms_result$TPR ==1 && ms_result$FDR == 0)
    
    cat("nrep = ", i , "MS complete of main iter",  k1,"\n")
    
    M = cbind(prop_recovery_iht_screening = prop_recovery_iht_screening, prop_recovery_iht_screening_sd = prop_recovery_iht_screening^2,
              prop_recovery_ms = prop_recovery_ms, prop_recovery_ms_sd = prop_recovery_ms^2,
              prop_recovery_lasso = prop_recovery_lasso, prop_recovery_lasso_sd = prop_recovery_lasso^2)
    M = M/nrep
  }
  
  df = data.frame(df)
  df$prop_recovery_iht_screening_sd = sqrt( df$prop_recovery_iht_screening_sd - df$prop_recovery_iht_screening^2)
  #df$f_score_path_iht_screening_sd = sqrt( df$f_score_path_iht_screening_sd - df$f_score_path_iht_screening^2)
  
  df$prop_recovery_ms_sd = sqrt( df$prop_recovery_ms_sd -  df$prop_recovery_ms^2)
  #df$f_score_path_ms_sd = sqrt( df$f_score_path_ms_sd -  df$f_score_path_ms^2)
  
  df$prop_recovery_lasso_sd = sqrt( df$prop_recovery_lasso_sd -  df$prop_recovery_lasso^2)
  #df$f_score_path_lasso_sd = sqrt( df$f_score_path_lasso_sd -  df$f_score_path_lasso^2)
  #df$s_budget = 0:(2*s)
  
  w = 1.96/sqrt(nrep)
  
  df$prop_recovery_lasso_low = df$prop_recovery_lasso - w*df$prop_recovery_lasso_sd
  df$prop_recovery_lasso_up = df$prop_recovery_lasso + w*df$prop_recovery_lasso_sd
  
  
  df$prop_recovery_iht_screening_low = df$prop_recovery_iht_screening - w*df$prop_recovery_iht_screening_sd
  df$prop_recovery_iht_screening_up = df$prop_recovery_iht_screening + w*df$prop_recovery_iht_screening_sd
  
  df$prop_recovery_ms_low = df$prop_recovery_ms - w*df$prop_recovery_ms_sd
  df$prop_recovery_ms_up = df$prop_recovery_ms + w*df$prop_recovery_ms_sd
  
  
  save_df <- function (df, q0){
    if (q0 == 0){
      my_path = paste("/home/roysapta/IHT/cv_glmnet/everything_cv/varying_heterogeneity/error_path/iht_one/q_spiky/hetero/hetero_lars_error_path_r", r, "s", s,"p", p, "spiky", q0, "grad-step", eta0,"nrep",nrep,".csv", sep = "_")
      write.csv(df, file = my_path, row.names = F ) }
    else {
      my_path = paste("/home/roysapta/IHT/cv_glmnet/everything_cv/varying_heterogeneity/error_path/iht_one/q_spiky/hetero/hetero_lars_error_path_r", r, "s", s,"p", p, "spiky", q0, "grad-step", eta0,"nrep",nrep,".csv", sep = "_")
      write.csv(df, file = my_path, row.names = F )}
  }
  save_df(df, q0)
}
stopCluster(cluster)
