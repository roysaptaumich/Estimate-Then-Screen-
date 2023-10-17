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
delta_list = seq(0.5, 16, 0.5) # signal streangth paramter r = 1 + delta
l = length(delta_list)
q0 = 0.8 # probability of the weak signal
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
  
  
  s = floor(2* log (p))
  n = floor(p^k)
  iter  = which(delta_list == delta)
  
  cat(iter, "th main iteration has started","\n")
  
  # signal and noise independent design
  r = (1 + delta)
  tp = sqrt(r)
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
  
  
  
  
  
  # ## decide lambda values in picasso (scad)
  # lambda_list_scad <- lambda_select_scad(n = n, mu_x = rep(0,p), S_x = S_x, 
  #                                        s = s, betamin = betamin, sigma_e = sigma_e, ntest = ceiling(nrep/10))
  
  ## decide lambda values for glmnet (lasso)
  #lambda_list_lasso <- lambda_select_lasso(n = n, mu_x = rep(0,p), S_x = S_x, 
  #                                         s = s, betamin = betamin, sigma_e = sigma_e, ntest = 1, 
  #                                         data_gen_method = 2  , pkg = "glmnet", tp = tp, q = q0)
  
  
  
  error_list = foreach (i = 1:nrep, .packages = c("MASS", "Matrix", "picasso", "glmnet","foreach")) %dopar% {
    
    
    ## generate data
    data_i <- data_lm_newscheme(n = n, mu_x = rep(0,p), S_x = S_x,
                       s = s, betamin = betamin, sigma_e = sigma_e, tp = tp)
    
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
    iht_screening_result = iht_cv_grid_screening_select(X = X, Y= Y,  supp_true = supp_true, beta0 = rep(0,p), gamma = gamma, a = a, v = 1, maxiter = 1e3, prec = 1e-7, lim_max = 15)
    TPR_iht_screening_list_rep = iht_screening_result$TPR
    FDR_iht_screening_list_rep = iht_screening_result$FDR
    prop_recovery_iht_screening_rep = 1* I(TPR_iht_screening_list_rep ==1 && FDR_iht_screening_list_rep == 0)
    s_hat = iht_screening_result$s_hat
    
    cat("nrep = ", i , "iht+screening complete",  "main iter =", iter , "of", l ,"\n")
    
    
    ## IHT + screening no oracle
    iht_screening_result_no = iht_cv_grid_screening_select_no_oracle(X = X, Y= Y,  supp_true = supp_true, beta0 = rep(0,p), gamma = gamma, a = a, v = 1, maxiter = 1e3, prec = 1e-7, lim_max = 15)
    TPR_iht_screening_list_rep_no = iht_screening_result_no$TPR
    FDR_iht_screening_list_rep_no = iht_screening_result_no$FDR
    prop_recovery_iht_screening_rep_no = 1* I(TPR_iht_screening_list_rep_no ==1 && FDR_iht_screening_list_rep_no == 0)
    #s_hat = iht_screening_result$s_hat
    
    cat("nrep = ", i , "iht+screening no oracle complete",  "main iter =", iter , "of", l ,"\n")
    
    
    ## marginal screening 
    ms_result = ms_ordered_select_cv(X = X, Y = Y, supp_true = supp_true, scale_X = F, s_level_list = 1:15)
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
  prop_recovery_iht_screening_no_list = c(prop_recovery_iht_screening_no_list, apply(do.call(rbind, error_list2[,11]), 2, mean))
  
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

# VARYING HETEROGENEITY
save_df <- function (df, q0){
  if (q0 == 1){
    my_path = paste("/home/roysapta/IHT/cv_glmnet/everything_cv/varying_heterogeneity/df_fdr_tpr_iht_ms_lasso_CV_homo_newscheme_p", p,"s",s,nrep,"rep.csv", sep = "_")
    write.csv(df, file = my_path, row.names = F ) }
  else {
    my_path = paste("/home/roysapta/IHT/cv_glmnet/everything_cv/varying_heterogeneity/df_fdr_tpr_iht_ms_lasso_CV_hetero_newscheme_p", p,"s",s,nrep,"rep.csv", sep = "_")
    write.csv(df, file = my_path, row.names = F )}
}
save_df(df, q0)

# VARYING SNR
#write.csv(df, file = "/home/roysapta/IHT/cv_glmnet/everything_cv/varying_SNR_p_1000/df_fdr_tpr_iht_ms_lasso_CV_hetero-type2_newscheme_tp=sqrt(r)_no_oracle_p1000_200_rep.csv", row.names = F)


# VARYING dimension and SNR

#write.csv(df, file = "/home/roysapta/IHT/cv_glmnet/everything_cv/df_fdr_tpr_iht_ms_lasso_CV_hetero-type2_newscheme_tp=sqrt(r)_no_oracle_p1000_200_rep.csv", row.names = F)


fig2 <- plot_ly(data = df, x = ~r) %>%
  add_trace(y = ~FDR_lasso, mode = 'lines', name = 'LASSO (FDR)', type = 'scatter', line = list(color = "red"))%>%
  add_trace(y = ~FDR_ms_list, mode = 'lines', name = 'MS (FDR)', type = 'scatter', line = list(color = "blue"))%>%
  add_trace(y = ~FDR_iht, mode = 'lines', name = 'ETS (FDR)', type = 'scatter', line = list(color = "green"))%>%
  #add_trace(y = ~FDR_iht_one_screening_list, mode = 'lines', name = 'ETS 1 (FDR)', type = 'scatter', line = list(color = "pink"))%>%
  
  #add_trace(y = ~TPR_iht_one_screening_list, mode = 'lines', name = 'ETS 1 (TPR)', linetype = I("dash"), line = list(color = "pink"))
  layout(legend=list(title=list(text='<b> </b>')) , yaxis= list(title = "FDR"), xaxis = list(title='r ')) %>%
  add_ribbons(ymin = ~FDR_lasso_low,
              ymax = ~FDR_lasso_up,
              showlegend = FALSE,
              opacity = 0.3) %>%
  add_ribbons(ymin = ~FDR_ms_list_low,
              ymax = ~FDR_ms_list_up,
              showlegend = FALSE,
              opacity = 0.3) %>%
  add_ribbons(ymin = ~FDR_iht_low,
              ymax = ~FDR_iht_up,
              showlegend = FALSE,
              fillcolor = 'rgba(7, 100, 181, 0.2)',
              opacity = 0.3)




fig2

fig3 <- plot_ly(data = df, x=~r) %>%
  add_trace(y = ~TPR_lasso, mode = 'lines', name = 'LASSO (TPR)', type = 'scatter', line = list(color = "red"))%>%
  add_trace(y = ~TPR_ms_list, mode = 'lines', name = 'MS (TPR)', type = 'scatter', line = list(color = "blue"))%>%
  add_trace(y = ~TPR_iht, mode = 'lines', name = 'ETS (TPR)', type = 'scatter', line = list(color = "green"))%>%
  layout(legend=list(title=list(text='<b> </b>')) , yaxis= list(title = "TPR"), xaxis = list(title='r ')) %>%
  add_ribbons(ymin = ~TPR_lasso_low,
              ymax = ~TPR_lasso_up,
              showlegend = FALSE,
              opacity = 0.3) %>%
  add_ribbons(ymin = ~TPR_ms_list_low,
              ymax = ~TPR_ms_list_up,
              showlegend = FALSE,
              opacity = 0.3) %>%
  add_ribbons(ymin = ~TPR_iht_low,
              ymax = ~TPR_iht_up,
              showlegend = FALSE,
              fillcolor = 'rgba(7, 100, 181, 0.2)',
              opacity = 0.3)

fig3

# plot(1-v_list, FDR_lasso, type='l', 
#      ylim = c(0.0,1),
#      xlab = "1-v", ylab = "FDR")
# #points(FDR_scad_cv, TPR_scad_cv)
# lines(1-v_list, FDR_iht_screening_list, col= "blue")
# lines(1-v_list, FDR_scad, col="red")
# #lines (1-v_list, FDR_iht, col = "blue")
# 
# 
# #points(FDR_lasso_cv, TPR_lasso_cv, col="red")
# #lines(FDR_iht2[, 1] , TPR_iht2[, 1], col="green")
# #lines(FDR_iht2[, 2] , TPR_iht2[, 2], col="blue")
# #lines(FDR_ms, TPR_ms, col="blueviolet")
# legend(x=0.51, y=0.4, legend=c("Lasso","ETS","SCAD"),
#        lty=1, col=c("black","blue", "red"), cex=0.9, text.width = 5)
# #legend(x=0.5, y=0, legend=c("SCAD-cv", "Lasso-cv"), pch=1, 
# #       col=c("black", "red"), bty="n")
# #legend(x = 0.7, y=0, legend = c("SIS"), lty=1, col = c("blueviolet"), bty = "n", ncol = 1)