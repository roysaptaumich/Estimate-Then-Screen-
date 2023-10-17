###### sparse linear regression ######
library(MASS)
library(Matrix)
library(picasso)
library(glmnet)
library(doSNOW)
library(foreach)
library(abind)
library(parallel)
n_cores = detectCores()
cluster = makeCluster(25, outfile="")
registerDoSNOW(cluster)
source("/home/roysapta/IHT/IHT.R")
source("/home/roysapta/IHT/IHT_screening.R")

## parameters
p = 1000 
k = 0.9 # sample size parameter
delta_list = seq(0.5, 10, 0.5) # signal streangth paramter r = 1 + delta
l = length(delta_list)

# mu = 0.8 # equi correlation 
# phi = 0.8 # AR(1) correlation
gamma = 0.7


nrep = 200

# TPR_iht_one_screening_list = c()
# FDR_iht_one_screening_list = c()


TPR_iht_screening_list = c()
FDR_iht_screening_list = c()
TPR_iht_screening_list_sd = c()
FDR_iht_screening_list_sd = c()
prop_recovery_iht_screening_list = c()

TPR_iht_screening_no_list = c()
FDR_iht_screening_no_list = c()
TPR_iht_screening_no_list_sd = c()
FDR_iht_screening_no_list_sd = c()
prop_recovery_iht_screening_no_list = c()

TPR_ms_list = c()
FDR_ms_list = c()
TPR_ms_list_sd = c()
FDR_ms_list_sd = c()
prop_recovery_ms_list = c()

# TPR_iht = c()
# FDR_iht =c()
# 
TPR_lasso = c()
FDR_lasso= c()
TPR_lasso_sd = c()
FDR_lasso_sd = c()
prop_recovery_lasso_list = c()
# TPR_scad = c()
# FDR_scad = c()




## compute TPR and FPR

####################
####################

for (delta in delta_list){
  
  
  s = 4*floor(2* log (p))
  n = floor(p^k)
  iter  = which(delta_list == delta)
  
  cat(iter, "th main iteration has started","\n")
  
  # signal and noise independent design
  r = (1 + delta)
  snr = r
  betamin = sqrt(2* r *log(p)/n)
  a = betamin
  S_x = diag(p)
  
  # # signa and noise AR(1) design
  # # m_phi = 0.5 * (1 - phi)/(1 + phi)
  # # A0 = 5 * (phi/ (m_phi * sqrt(1 - phi^2))) * (v/gamma)^0.5 
  # # betamin  = A0 * sqrt(s * log(p) / n)
  # r = (2)*(1 + sqrt(1-v))^2
  # betamin = sqrt(2* r *log(p)/n)
  # a = betamin
  # S_x = phi^abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - 
  #                 (1:p - 1))
  
  sigma_e = 1
  
  q0 = 0.8 # probability of the weak signal
  
  
  
  # ## decide lambda values in picasso (scad)
  # lambda_list_scad <- lambda_select_scad(n = n, mu_x = rep(0,p), S_x = S_x, 
  #                                        s = s, betamin = betamin, sigma_e = sigma_e, ntest = ceiling(nrep/10))
  
  ## decide lambda values for glmnet (lasso)
  #lambda_list_lasso <- lambda_select_lasso(n = n, mu_x = rep(0,p), S_x = S_x, 
  #                                         s = s, betamin = betamin, sigma_e = sigma_e, ntest = 1, 
  #                                         data_gen_method = 2  , pkg = "glmnet", tp = tp, q = q0)
  
  
  
  error_list = foreach (i = 1:nrep, .packages = c("MASS", "Matrix", "picasso", "glmnet","foreach")) %dopar% {
    
    
    ## generate data
    data_i <- data_lm2_newscheme_SNR(n = n, mu_x = rep(0,p), S_x = S_x,
                                s = s, betamin = betamin, sigma_e = sigma_e, snr = snr)
    
    #data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
    #                  s = s, betamin = betamin, sigma_e = sigma_e, tp = tp, q = q0)
    
    X = data_i$X
    Y = data_i$Y  
    beta = data_i$beta
    epsilon = data_i$epsilon
    supp_true = data_i$supp_true
    
    
    
    # ## optimization with nonconvex penalty 
    # scad_result = scad_select(X = X, Y = Y, supp_true = supp_true, lambda_list = lambda_list_scad)
    # TPR_scad_rep = scad_result$TPR_cv
    # FDR_scad_rep = scad_result$FDR_cv
    # cat("nrep = ", i , "scad complete", "main iter =", iter, "of", l ,"\n")
    # 
    ## optimization with convex penalty
    lasso_result = lasso_cv_select(X = X, Y = Y, supp_true = supp_true, lambda_list = NULL)
    TPR_lasso_rep = lasso_result$TPR_cv
    FDR_lasso_rep = lasso_result$FDR_cv
    prop_recovery_lasso_rep = 1* I(TPR_lasso_rep ==1 && FDR_lasso_rep == 0)
    
    cat("nrep = ", i , "lasso complete", "main iter =", iter, "of", l ,"\n")
    
    
    ## IHT + screening 
    iht_screening_result = iht_cv_grid_screening_select(X = X, Y= Y,  supp_true = supp_true, beta0 = rep(0,p), gamma = gamma, a = a, v = 1, maxiter = 1e3, prec = 1e-7, lim_max = s+2)
    TPR_iht_screening_list_rep = iht_screening_result$TPR
    FDR_iht_screening_list_rep = iht_screening_result$FDR
    prop_recovery_iht_screening_rep = 1* I(TPR_iht_screening_list_rep ==1 && FDR_iht_screening_list_rep == 0)
    s_hat = iht_screening_result$s_hat
    
    cat("nrep = ", i , "iht+screening complete",  "main iter =", iter , "of", l ,"\n")
    
    
    ## IHT + screening no oracle
    #iht_screening_result_no = iht_cv_grid_screening_select_no_oracle(X = X, Y= Y,  supp_true = supp_true, beta0 = rep(0,p), gamma = gamma, a = a, v = 1, maxiter = 1e3, prec = 1e-7, lim_max = s+2)
    TPR_iht_screening_list_rep_no = iht_screening_result$TPR_no
    FDR_iht_screening_list_rep_no = iht_screening_result$FDR_no
    prop_recovery_iht_screening_rep_no = 1* I(TPR_iht_screening_list_rep_no ==1 && FDR_iht_screening_list_rep_no == 0)
    #s_hat = iht_screening_result$s_hat
    
    cat("nrep = ", i , "iht+screening no oracle complete",  "main iter =", iter , "of", l ,"\n")
    
    
    ## marginal screening 
    ms_result = ms_ordered_select_cv(X = X, Y = Y, supp_true = supp_true, scale_X = F, s_level_list = 1:(s+2))
    TPR_ms_list_rep = ms_result$TPR.1se
    FDR_ms_list_rep = ms_result$FDR.1se
    prop_recovery_ms_rep = 1* I(TPR_ms_list_rep ==1 && FDR_ms_list_rep == 0)
    
    cat("nrep = ", i , "MS complete",  "main iter =", iter , "of", l ,"\n")
    
    ## IHT one stage
    # iht_one_screening_result = iht_one_screening_select(X = X, Y= Y, s_list = s , eta_list = 0.1 , supp_true = supp_true, beta0 = rep(0,p), gamma = 0.7, a = a, v = v, maxiter = 1e3, prec = 1e-7)
    # TPR_iht_one_screening_list_rep = iht_one_screening_result$TPR
    # FDR_iht_one_screening_list_rep = iht_one_screening_result$FDR
    # 
    # cat("nrep = ", i , "iht one stage complete",  "main iter =", iter , "of", l ,"\n")
    
    ## summarize results
    
    return(list(FDR_lasso_rep = FDR_lasso_rep, TPR_lasso_rep = TPR_lasso_rep,
                FDR_ms_list_rep = FDR_ms_list_rep, TPR_ms_list_rep = TPR_ms_list_rep,
                FDR_iht_screening_list_rep = FDR_iht_screening_list_rep, TPR_iht_screening_list_rep = TPR_iht_screening_list_rep,
                FDR_iht_screening_list_rep_no = FDR_iht_screening_list_rep_no, TPR_iht_screening_list_rep_no = TPR_iht_screening_list_rep_no,
                prop_recovery_lasso_rep = prop_recovery_lasso_rep, prop_recovery_ms_rep = prop_recovery_ms_rep,
                prop_recovery_iht_screening_rep = prop_recovery_iht_screening_rep, prop_recovery_iht_screening_rep_no = prop_recovery_iht_screening_rep_no))
    #FDR_iht_one_screening_list_rep = FDR_iht_one_screening_list_rep, TPR_iht_one_screening_list_rep = TPR_iht_one_screening_list_rep))
    
  }
  
  
  error_list2 = do.call(rbind, error_list)
  
  # TPR_iht_one_screening_list = c(apply(do.call(rbind, error_list2[,8]), 2, mean), TPR_iht_one_screening_list)
  # FDR_iht_one_screening_list = c(apply(do.call(rbind, error_list2[,7]), 2, mean), FDR_iht_one_screening_list)
  
  
  TPR_lasso = c(TPR_lasso, apply(do.call(rbind, error_list2[,2]), 2, mean)  )
  FDR_lasso= c(FDR_lasso, apply(do.call(rbind, error_list2[,1]), 2, mean) )
  TPR_lasso_sd = c(TPR_lasso_sd, apply(do.call(rbind, error_list2[,2]), 2, sd)  )
  FDR_lasso_sd = c(FDR_lasso_sd , apply(do.call(rbind, error_list2[,1]), 2, sd) )
  
  TPR_ms_list = c( TPR_ms_list, apply(do.call(rbind, error_list2[,4]), 2, mean))
  FDR_ms_list = c(FDR_ms_list, apply(do.call(rbind, error_list2[,3]), 2, mean))
  TPR_ms_list_sd = c( TPR_ms_list_sd, apply(do.call(rbind, error_list2[,4]), 2, sd))
  FDR_ms_list_sd = c(FDR_ms_list_sd, apply(do.call(rbind, error_list2[,3]), 2, sd))
  
  TPR_iht_screening_list = c(TPR_iht_screening_list, apply(do.call(rbind, error_list2[,6]), 2, mean))
  FDR_iht_screening_list = c(FDR_iht_screening_list, apply(do.call(rbind, error_list2[,5]), 2, mean))
  TPR_iht_screening_list_sd = c(TPR_iht_screening_list_sd, apply(do.call(rbind, error_list2[,6]), 2, sd))
  FDR_iht_screening_list_sd = c(FDR_iht_screening_list_sd, apply(do.call(rbind, error_list2[,5]), 2, sd))
  
  TPR_iht_screening_no_list = c(TPR_iht_screening_no_list, apply(do.call(rbind, error_list2[,8]), 2, mean))
  FDR_iht_screening_no_list = c(FDR_iht_screening_no_list, apply(do.call(rbind, error_list2[,7]), 2, mean))
  TPR_iht_screening_no_list_sd = c(TPR_iht_screening_no_list_sd, apply(do.call(rbind, error_list2[,8]), 2, sd))
  FDR_iht_screening_no_list_sd = c(FDR_iht_screening_no_list_sd, apply(do.call(rbind, error_list2[,7]), 2, sd))
  
  
  
  
  prop_recovery_lasso_list = c(prop_recovery_lasso_list, apply(do.call(rbind, error_list2[,9]), 2, mean))
  prop_recovery_ms_list = c(prop_recovery_ms_list, apply(do.call(rbind, error_list2[,10]), 2, mean))
  prop_recovery_iht_screening_list = c(prop_recovery_iht_screening_list, apply(do.call(rbind, error_list2[,11]), 2, mean))
  prop_recovery_iht_screening_no_list = c(prop_recovery_iht_screening_no_list, apply(do.call(rbind, error_list2[,12]), 2, mean))
  
}


library(plotly)
df <- data.frame(r = (1 + delta_list) , FDR_lasso = FDR_lasso, FDR_ms_list = FDR_ms_list, FDR_iht = FDR_iht_screening_list, FDR_iht_no = FDR_iht_screening_no_list,  
                 FDR_lasso_sd = FDR_lasso_sd, FDR_ms_list_sd = FDR_ms_list_sd, FDR_iht_sd = FDR_iht_screening_list_sd, FDR_iht_no_sd = FDR_iht_screening_no_list_sd,
                 TPR_lasso = TPR_lasso, TPR_ms_list = TPR_ms_list, TPR_iht = TPR_iht_screening_list, TPR_iht_no = TPR_iht_screening_no_list,
                 TPR_lasso_sd = TPR_lasso_sd, TPR_ms_list_sd = TPR_ms_list_sd, TPR_iht_sd = TPR_iht_screening_list_sd, TPR_iht_no_sd = TPR_iht_screening_no_list_sd,
                 prop_recovery_lasso_list = prop_recovery_lasso_list, prop_recovery_ms_list = prop_recovery_ms_list,
                 prop_recovery_iht_screening_list = prop_recovery_iht_screening_list, prop_recovery_iht_screening_no_list = prop_recovery_iht_screening_no_list)
w = 1.96/sqrt(nrep)
df$FDR_lasso_low = df$FDR_lasso - w*df$FDR_lasso_sd
df$FDR_lasso_up = df$FDR_lasso + w*df$FDR_lasso_sd
df$TPR_lasso_low = df$TPR_lasso - w*df$TPR_lasso_sd
df$TPR_lasso_up = df$TPR_lasso + w*df$TPR_lasso_sd

df$FDR_iht_low = df$FDR_iht - w*df$FDR_iht_sd
df$FDR_iht_up = df$FDR_iht + w*df$FDR_iht_sd
df$TPR_iht_low = df$TPR_iht - w*df$TPR_iht_sd
df$TPR_iht_up = df$TPR_iht + w*df$TPR_iht_sd

df$FDR_iht_no_low = df$FDR_iht_no - w*df$FDR_iht_no_sd
df$FDR_iht_no_up = df$FDR_iht_no + w*df$FDR_iht_no_sd
df$TPR_iht_no_low = df$TPR_iht_no - w*df$TPR_iht_no_sd
df$TPR_iht_no_up = df$TPR_iht_no + w*df$TPR_iht_no_sd

df$FDR_ms_list_low = df$FDR_ms_list - w*df$FDR_ms_list_sd
df$FDR_ms_list_up = df$FDR_ms_list + w*df$FDR_ms_list_sd
df$TPR_ms_list_low = df$TPR_ms_list - w*df$TPR_ms_list_sd
df$TPR_ms_list_up = df$TPR_ms_list + w*df$TPR_ms_list_sd



# VARYING SNR
my_path = paste("/home/roysapta/IHT/cv_glmnet/everything_cv/varying_SNR_p_1000/df_fdr_tpr_iht_ms_lasso_CV_hetero-type2_newscheme_varying_SNR_p",p,"s",s , nrep,"rep.csv", sep = "_")

write.csv(df, file = my_path, row.names = F)
