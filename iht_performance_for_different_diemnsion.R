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

#df = read.csv("/home/roysapta/IHT/cv_glmnet/everything_cv/df_fdr_tpr_iht_ms_lasso_CV_hetero-type2_p1000_200_rep.csv")
#data.frame(df$r-1, df$prop_recovery_iht_screening_list)

delta_list = 9-1
r_list = 1+delta_list
p_list = 500*(1:6)
nrep = 100
O = matrix(0, nrow = length(r_list), ncol = length(p_list))
NO = matrix(0, nrow = length(r_list), ncol = length(p_list))
 for(i in 1:length(delta_list)) {
  r = r_list[i]
  for(j in 1:length(p_list)){
  p = p_list[j]
  k = 0.9
  
  s =   floor(2 * log(p)) #floor(p^(1-v))
  n = floor(p^k)
  
  gamma = 0.7
  
  # Identity covariance matrix
  S_x = diag(p)
  
  
  
  betamin = sqrt(2 *r * log(p)/n) 
  a = betamin
  #q0 = 0.8 # proportion of homogeneity
  sigma_e = 1
  snr = r
  
  error_list = foreach (k = 1:nrep, .packages = c("MASS", "Matrix", "picasso", "glmnet","foreach")) %dopar% {
    
  
  #data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
  #                 s = s, betamin = betamin, sigma_e = sigma_e, tp, q = q0)
  
  
  data_i <- data_lm2_newscheme_SNR(n = n, mu_x = rep(0,p), S_x = S_x,
                     s = s, betamin = betamin, sigma_e = sigma_e, snr = snr)
  X = data_i$X
  Y = data_i$Y  
  beta = data_i$beta
  epsilon = data_i$epsilon
  supp_true = data_i$supp_true
  
  iht_screening_result = iht_cv_grid_screening_select(X = X, Y= Y,  supp_true = supp_true, beta0 = rep(0,p), gamma = gamma, a = a, v = 1, maxiter = 1e3, prec = 1e-7, lim_max = s + 2)
  TPR_iht_screening_list_rep = iht_screening_result$TPR
  FDR_iht_screening_list_rep = iht_screening_result$FDR
  prop_recovery_iht_screening_rep = 1* I(TPR_iht_screening_list_rep ==1 && FDR_iht_screening_list_rep == 0)
  
  iht_screening_no_result = iht_cv_grid_screening_select_no_oracle(X = X, Y= Y,  supp_true = supp_true, beta0 = rep(0,p), gamma = gamma, a = a, v = 1, maxiter = 1e3, prec = 1e-7, lim_max = s + 2)
  TPR_iht_screening_no_list_rep = iht_screening_no_result$TPR
  FDR_iht_screening_no_list_rep = iht_screening_no_result$FDR
  prop_recovery_iht_no_screening_rep = 1* I(TPR_iht_screening_no_list_rep ==1 && FDR_iht_screening_no_list_rep == 0)
  
  
  
  cat(k,"th nrep for", j,"th p", i,"th r", "\n")
  return(list(prop_recovery_iht_screening_rep = prop_recovery_iht_screening_rep, prop_recovery_iht_no_screening_rep = prop_recovery_iht_no_screening_rep))
  }
  error_list2 = do.call(rbind, error_list)
  
  O[i, j] = apply(do.call(rbind, error_list2[,1]), 2, mean)
  NO[i, j] = apply(do.call(rbind, error_list2[,2]), 2, mean)
  #val = error_list/nrep
}
#val2 = df1
}

rnames = as.character(r_list)
cnames = as.character(p_list)

rownames(O) = rnames; colnames(O) = cnames
rownames(NO) = rnames; colnames(NO) = cnames


write.table(O, file = paste("/home/roysapta/IHT/cv_glmnet/everything_cv/varying_dimension_and_SNR/iht_performance_for_different_dimension_oracle_mat_snr", r_list, sep = "_" ))
write.table(NO, file = paste("/home/roysapta/IHT/cv_glmnet/everything_cv/varying_dimension_and_SNR/iht_performance_for_different_dimension_nonoracle_mat_snr10", r_list, sep = "_"))


#m<-read.table("/home/roysapta/IHT/cv_glmnet/everything_cv/varying_dimension_and_SNR/iht_performance_for_different_dimension_mat_small") 
#m <- as.matrix(t(m)[ncol(m):1,])

#require(plotly)
#c1 = as.character(6:11)
#fig <- plot_ly(z =abs(m[nrow(m):1,]), 
#               x = c1, y = p_list,
#               type = "heatmap") %>% #colors = colorRamp(c("yellow", "green"))
#  layout(xaxis = list(title = "r"), yaxis = list(title = "p"))
#fig

