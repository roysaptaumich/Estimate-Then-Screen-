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
#r = 6  # signal streangth paramter r = 1 + delta

q0 = 0.8 # probability of the weak signal

gamma = 0.7
nrep = 100




s = 1*floor(2* log (p))
n = floor(p^k)



# signal and noise independent design
r_list = 1 + seq(0.5,9,0.5)
for (k1 in 1:length(r_list)){
r = r_list[k1]
tp = sqrt(r)
betamin = sqrt(2* r *log(p)/n)
a = betamin
S_x = diag(p)


sigma_e = 1

####### selecting lambda_lasso list ########
data_i <- data_lm_newscheme(n = n, mu_x = rep(0,p), S_x = S_x,
                            s = s, betamin = betamin, sigma_e = sigma_e, tp = tp, q = q0)


X = data_i$X
Y = data_i$Y 
mod = glmnet(x = X, y = Y, alpha = 1, standardize = F, intercept = F, nlambda = 500)
lambda_list = mod$lambda

############

df = foreach(i = 1:nrep, .packages = c("MASS", "Matrix", "picasso", "glmnet","foreach"), .combine = '+')%dopar%{
  data_i <- data_lm_newscheme(n = n, mu_x = rep(0,p), S_x = S_x,
                              s = s, betamin = betamin, sigma_e = sigma_e, tp = tp, q = q0)

  
  X = data_i$X
  Y = data_i$Y  
  beta = data_i$beta
  epsilon = data_i$epsilon
  supp_true = data_i$supp_true
  
  
  # lasso
  lasso_result = glmnet(x= X, y = Y, alpha = 1, standardize = F, intercept = F, lambda = lambda_list)
  beta_mat = lasso_result$beta
  s_budget = colSums(beta_mat!=0); s_budget_max = 2*s; s_unique_budget = unique(s_budget)[which(unique(s_budget)<= s_budget_max)]
  FDR_path_lasso = rep(0, (s_budget_max+1))
  TPR_path_lasso = rep(0, (s_budget_max+1))
  
  for (k in 1:length(s_unique_budget)) {
    id  = which(s_budget == s_unique_budget[k])[1]
    set = which(beta_mat[ , id] != 0)
    FDR_path_lasso[s_unique_budget[k]] = length(setdiff(set, supp_true))/ max(1, length(set))
    TPR_path_lasso[s_unique_budget[k]] = length(intersect(set, supp_true))/ length(supp_true)
    
  }
  prop_recovery_path_lasso =  I(TPR_path_lasso == 1)* I(FDR_path_lasso == 0)
  cat("nrep = ", i , "lasso complete main iter = ", k1, "\n")
  
  ## IHT + screening 
  iht_screening_result = iht_cv_grid_screening_select(X = X, Y= Y,  supp_true = supp_true, beta0 = rep(0,p), gamma = gamma, a = a, v = 1, maxiter = 1e3, prec = 1e-7, lim_max = s+2)
  TPR_path_iht_screening = iht_screening_result$TPR_path[1:(s_budget_max+1)]
  FDR_path_iht_screening = iht_screening_result$FDR_path[1:(s_budget_max+1)]
  prop_recovery_path_iht_screening = I(TPR_path_iht_screening ==1) * I(FDR_path_iht_screening == 0)
  #s_hat = iht_screening_result$s_hat
  
  cat("nrep = ", i , "iht+screening complete main iter = ", k1, "\n")
  
  ## marginal screening 
  ms_result = ms_ordered_select(X = X, Y = Y, supp_true = supp_true, scale_X = F, s_level = s)
  TPR_path_ms = ms_result$TPR_path[1:(s_budget_max+1)]
  FDR_path_ms = ms_result$FDR_path[1:(s_budget_max+1)]
  prop_recovery_path_ms = I(TPR_path_ms ==1) * I(FDR_path_ms == 0)
  
  cat("nrep = ", i , "MS complete mai iter =",  k1, "\n")
  
  M = cbind(TPR_path_iht_screening = TPR_path_iht_screening, FDR_path_iht_screening = FDR_path_iht_screening, FDR_path_iht_screening_sd = FDR_path_iht_screening^2, TPR_path_iht_screening_sd = TPR_path_iht_screening^2,prop_recovery_path_iht_screening = prop_recovery_path_iht_screening,
                 TPR_path_ms = TPR_path_ms, FDR_path_ms = FDR_path_ms, FDR_path_ms_sd = FDR_path_ms^2 , TPR_path_ms_sd = TPR_path_ms^2 ,prop_recovery_path_ms = prop_recovery_path_ms,
                 TPR_path_lasso = TPR_path_lasso, FDR_path_lasso = FDR_path_lasso, FDR_path_lasso_sd = FDR_path_lasso^2 , TPR_path_lasso_sd = TPR_path_lasso^2,prop_recovery_path_lasso = prop_recovery_path_lasso)
  M = M/nrep
}
df = data.frame(df)
df$FDR_path_iht_screening_sd = sqrt( df$FDR_path_iht_screening_sd - df$FDR_path_iht_screening^2)
df$FDR_path_ms_sd = sqrt( df$FDR_path_ms_sd -  df$FDR_path_ms^2)
df$FDR_path_lasso_sd = sqrt( df$FDR_path_lasso_sd -  df$FDR_path_lasso^2)
df$s_budget = 0:(2*s)

w = 1.96/sqrt(nrep)

df$FDR_lasso_low = df$FDR_path_lasso - w*df$FDR_path_lasso_sd
df$FDR_lasso_up = df$FDR_path_lasso + w*df$FDR_path_lasso_sd
df$TPR_lasso_low = df$TPR_path_lasso - w*df$TPR_path_lasso_sd
df$TPR_lasso_up = df$TPR_path_lasso + w*df$TPR_path_lasso_sd

df$FDR_iht_low = df$FDR_path_iht_screening - w*df$FDR_path_iht_screening_sd
df$FDR_iht_up = df$FDR_path_iht_screening + w*df$FDR_path_iht_screening_sd
df$TPR_iht_low = df$TPR_path_iht_screening - w*df$TPR_path_iht_screening_sd
df$TPR_iht_up = df$TPR_path_iht_screening + w*df$TPR_path_iht_screening_sd


df$FDR_ms_list_low = df$FDR_path_ms - w*df$FDR_path_ms_sd
df$FDR_ms_list_up = df$FDR_path_ms + w*df$FDR_path_ms_sd
df$TPR_ms_list_low = df$TPR_path_ms - w*df$TPR_path_ms_sd
df$TPR_ms_list_up = df$TPR_path_ms + w*df$TPR_path_ms_sd


save_df <- function (df, q0){
  if (q0 == 1){
    my_path = paste("/home/roysapta/IHT/cv_glmnet/everything_cv/varying_heterogeneity/error_path/with_lasso/homo/homo_lars_error_path_r", r, "s", s,"p", p, "nrep",nrep,".csv", sep = "_")
    write.csv(df, file = my_path, row.names = F ) }
  else {
    my_path = paste("/home/roysapta/IHT/cv_glmnet/everything_cv/varying_heterogeneity/error_path/with_lasso/hetero/hetero_lars_error_path_r", r, "s", s,"p", p, "nrep",nrep,".csv", sep = "_")
    write.csv(df, file = my_path, row.names = F )}
}
save_df(df, q0)
}
stopCluster(cluster)

require(plotly)

fig2 <- plot_ly(data = df, x = ~s_budget) %>%
  add_trace(y = ~FDR_path_lasso, mode = 'lines', name = 'LASSO (FDR)', type = 'scatter', line = list(color = "red"))%>%
  add_trace(y = ~FDR_path_ms, mode = 'lines', name = 'MS (FDR)', type = 'scatter', line = list(color = "blue"))%>%
  add_trace(y = ~FDR_path_iht_screening, mode = 'lines', name = 'ETS (FDR)', type = 'scatter', line = list(color = "green"))%>%
  #add_trace(y = ~FDR_iht_one_screening_list, mode = 'lines', name = 'ETS 1 (FDR)', type = 'scatter', line = list(color = "pink"))%>%
  
  #add_trace(y = ~TPR_iht_one_screening_list, mode = 'lines', name = 'ETS 1 (TPR)', linetype = I("dash"), line = list(color = "pink"))
  layout(legend=list(title=list(text='<b> </b>')) , yaxis= list(title = "FDR"), xaxis = list(title='sparsity budget ')) %>%
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

fig3 <- plot_ly(data = df, x=~s_budget) %>%
  add_trace(y = ~TPR_path_lasso, mode = 'lines', name = 'LASSO (TPR)', type = 'scatter', line = list(color = "red"))%>%
  add_trace(y = ~TPR_path_ms, mode = 'lines', name = 'MS (TPR)', type = 'scatter', line = list(color = "blue"))%>%
  add_trace(y = ~TPR_path_iht_screening, mode = 'lines', name = 'ETS (TPR)', type = 'scatter', line = list(color = "green"))%>%
  layout(legend=list(title=list(text='<b> </b>')) , yaxis= list(title = "TPR"), xaxis = list(title='sparsity budget ')) %>%
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
