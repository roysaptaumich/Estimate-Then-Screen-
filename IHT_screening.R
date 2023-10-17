###### sparse linear regression ######
library(MASS)
library(Matrix)
library(picasso)
library(glmnet)
library(doSNOW)
library(foreach)
library(abind)
source("/home/roysapta/IHT/IHT.R")

iht_screening <- function(X, Y, s_list, g_list, beta0 = rep(0,p), gamma = 0.7,a, v, maxiter = 1e3, prec = 1e-7){
  ## the iterative hard thresholding algorithm
  ## Input: X, Y: data
  ##        scale_X (boolean): whether to scale X before training
  ##        s_list: list of possible values of parameter s
  ##        g_list (same length as s_list): list of possible values of parameter g
  ##        beta0: initial value of beta
  ##        maxiter: maximum number of iterations 
  ##        prec: precision where the iterations stop
  ##        gamma: partition parameter. Determines the size of 1st subsample and 2nd subsample
  ##        a: minimum signal strength
  ##        v: true sparsity parameter
  ## Output: iht_screening_result: list containing beta (coefficients), grad (gradient), steps_list (number of iterations performed)
  
  n = dim(X)[1]
  p = dim(X)[2]
  #X = scale(X)
  
  support_estimate = c()
  # partitoning the data
  n1 = floor(gamma * n)
  subset = sample(1:n, n1)
  X1 = X[subset,]; Y1 = as.matrix(Y[subset, ], ncol=1)
  X2 = X[-subset,]; Y2 = as.matrix(Y[-subset, ], ncol=1)
  
  # first setp IHT for 1st subsample
  model_iht2 = iht_cor2(X= X1, Y = Y1, scale_X = FALSE, s_list=s_list, g_list = g_list, beta0=beta0, maxiter = maxiter, prec = prec)
  beta_iht2 = model_iht2$beta
  grad_iht2 = model_iht2$grad
  steps_iht2 = model_iht2$steps_list
  
  
  
  # screening stage
  for (j in 1:p){
    feature_vec = as.matrix(X2[,j], ncol=1)
    alpha_selector = t(feature_vec) %*% (Y2 - X2[,-j] %*% as.matrix(beta_iht2[-j,], ncol=1))/norm(feature_vec,"2")
    threshold = a*norm(feature_vec, "2")/2 + v*log(p)/(a* norm(feature_vec,"2"))
    
    
    #print(c(abs(alpha_selector), threshold, norm(feature_vec,"2")))
    
    if(abs(alpha_selector)> threshold){
      support_estimate = c(j,support_estimate)
      #print("support entry success")
    }
    else{
      support_estimate = support_estimate
      #print("support entry failure")
    }
  }
  
  iht_screening_result = list(beta = beta_iht2, grad = grad_iht2, steps_list = steps_iht2, support_estimate = support_estimate)
  return(iht_screening_result)
  
  
}


iht_screening_select <- function(X, Y, s_list,  g_list, supp_true, beta0 = rep(0,p), gamma = 0.1,a, v, maxiter = 1e3, prec = 1e-7){
  ## the iterative hard thresholding algorithm
  ## Input: X, Y: data
  ##        scale_X (boolean): whether to scale X before training
  ##        s_list: list of possible values of parameter s
  ##        g_list (same length as s_list): list of possible values of parameter g
  ##        supp_true: true support of beta
  ##        beta0: initial value of beta
  ##        maxiter: maximum number of iterations 
  ##        prec: precision where the iterations stop
  ##        gamma: partition parameter. Determines the size of 1st subsample and 2nd subsample
  ##        a: minimum signal strength
  ##        v: true sparsity parameter
  ## Output: iht_screening_select: returns TPR and FPR
  
  model_iht_screening <- iht_screening(X, Y, s_list  = s_list, g_list = g_list, beta0 = rep(0,p), gamma = 0.7, a = a, v = v, maxiter = 1e3, prec = 1e-7)
  supp_beta_iht = model_iht_screening$support_estimate
  FDR = length(setdiff(supp_beta_iht, supp_true))/length(supp_beta_iht); FDR
  
  TPR = length(intersect(supp_beta_iht, supp_true))/length(supp_true); TPR
  
  iht_screening_select_result = list(TPR=TPR, FDR=FDR)
  return(iht_screening_select_result)
}


###### IHT CV + screen

iht_cv_screening <- function(X, Y, beta0 = rep(0,p), gamma = 0.7,a, v, maxiter = 1e3, prec = 1e-7){
  ## the iterative hard thresholding algorithm
  ## Input: X, Y: data
  ##        scale_X (boolean): whether to scale X before training
  ##        s_list: list of possible values of parameter s
  ##        g_list (same length as s_list): list of possible values of parameter g
  ##        beta0: initial value of beta
  ##        maxiter: maximum number of iterations 
  ##        prec: precision where the iterations stop
  ##        gamma: partition parameter. Determines the size of 1st subsample and 2nd subsample
  ##        a: minimum signal strength
  ##        v: true sparsity parameter
  ## Output: iht_screening_result: list containing beta (coefficients), grad (gradient), steps_list (number of iterations performed)
  
  n = dim(X)[1]
  p = dim(X)[2]
  #X = scale(X)
  
  support_estimate = c()
  # partitoning the data
  n1 = floor(gamma * n)
  subset = sample(1:n, n1)
  X1 = X[subset,]; Y1 = as.matrix(Y[subset, ], ncol=1)
  X2 = X[-subset,]; Y2 = as.matrix(Y[-subset, ], ncol=1)
  
  # first setp IHT for 1st subsample with CV
  model_iht_cv = iht_cv(X = X1, Y = Y1, s_list = 1: 50, g_list= 1:50, beta0 = rep(0,p), nfold = 10, n_cv = 1, maxiter = 1e3, prec = 1e-7)
  s_list = model_iht_cv$s.min; g_list = model_iht_cv$g.min
  s_hat = sum(model_iht_cv$coef_min!=0)
  
  model_iht2 = iht_cor2(X = X1, Y = Y1, scale_X =  FALSE, s_list =  s_list, g_list = g_list, beta0 = rep(0,p), maxiter = 1e3, prec = 1e-7)
  beta_iht2 = model_iht2$beta
  #support_estimate = model_iht2$support_estimate
  # grad_iht2 = model_iht2$grad
  # steps_iht2 = model_iht2$steps_list
  
  #iht_cor2()
  # screening stage
  for (j in 1:p){
    feature_vec = as.matrix(X2[,j], ncol=1)
    alpha_selector = t(feature_vec) %*% (Y2 - X2[,-j] %*% as.matrix(beta_iht2[-j,], ncol=1))/norm(feature_vec,"2")
    threshold = a*norm(feature_vec, "2")/2 + v*log(p)/(a* norm(feature_vec,"2"))
    
    
    #print(c(abs(alpha_selector), threshold, norm(feature_vec,"2")))
    
    if(abs(alpha_selector)> threshold){
      support_estimate = c(j,support_estimate)
      #print("support entry success")
    }
    else{
      support_estimate = support_estimate
      #print("support entry failure")
    }
  }
  iht_screening_result = list(beta = beta_iht2, s_hat = s_hat, support_estimate = support_estimate, s_list = s_list , g_list = g_list)
  return(iht_screening_result)
  
  
}


iht_cv_screening_select <- function(X, Y,  supp_true, beta0 = rep(0,p), gamma = 0.1,a, v, maxiter = 1e3, prec = 1e-7){
  ## the iterative hard thresholding algorithm
  ## Input: X, Y: data
  ##        scale_X (boolean): whether to scale X before training
  ##        s_list: list of possible values of parameter s
  ##        g_list (same length as s_list): list of possible values of parameter g
  ##        supp_true: true support of beta
  ##        beta0: initial value of beta
  ##        maxiter: maximum number of iterations 
  ##        prec: precision where the iterations stop
  ##        gamma: partition parameter. Determines the size of 1st subsample and 2nd subsample
  ##        a: minimum signal strength
  ##        v: true sparsity parameter
  ## Output: iht_screening_select: returns TPR and FPR
  
  model_iht_screening <- iht_cv_screening(X = X, Y = Y, beta0 = rep(0,p), gamma = 0.7, a = a, v = v, maxiter = 1e3, prec = 1e-7)
  supp_beta_iht = model_iht_screening$support_estimate
  s_hat = model_iht_screening$s_hat
  FDR = length(setdiff(supp_beta_iht, supp_true))/length(supp_beta_iht); FDR
  
  TPR = length(intersect(supp_beta_iht, supp_true))/length(supp_true); TPR
  
  iht_screening_select_result = list(TPR=TPR, FDR=FDR, s_hat = s_hat)
  return(iht_screening_select_result)
}

##### IHT CV grid + screen

iht_cv_grid_screening <- function(X, Y, beta0 = rep(0,p), gamma = 0.7,a, v, maxiter = 1e3, prec = 1e-7, lim_max = 15){
  ## the iterative hard thresholding algorithm
  ## Input: X, Y: data
  ##        scale_X (boolean): whether to scale X before training
  ##        s_list: list of possible values of parameter s
  ##        g_list (same length as s_list): list of possible values of parameter g
  ##        beta0: initial value of beta
  ##        maxiter: maximum number of iterations 
  ##        prec: precision where the iterations stop
  ##        gamma: partition parameter. Determines the size of 1st subsample and 2nd subsample
  ##        a: minimum signal strength
  ##        v: true sparsity parameter
  ## Output: iht_screening_result: list containing beta (coefficients), grad (gradient), steps_list (number of iterations performed)
  
  n = dim(X)[1]
  p = dim(X)[2]
  #X = scale(X)
  
  support_estimate = c()
  # partitoning the data
  n1 = floor(gamma * n)
  subset = sample(1:n, n1)
  X1 = X[subset,]; Y1 = as.matrix(Y[subset, ], ncol=1)
  X2 = X[-subset,]; Y2 = as.matrix(Y[-subset, ], ncol=1)
  
  # first setp IHT for 1st subsample with CV
  model_iht_cv = iht_cv_grid(X = X1, Y = Y1, s_list = 1: lim_max, g_list= 1:lim_max, beta0 = rep(0,p), nfold = 10, n_cv = 1, maxiter = 1e3, prec = 1e-7)
  s_list = model_iht_cv$s.min; g_list = model_iht_cv$g.min
  s_hat = sum(model_iht_cv$coef_min!=0)
  
  model_iht2 = iht_cor2(X = X1, Y = Y1, scale_X =  FALSE, s_list =  s_list, g_list = g_list, beta0 = rep(0,p), maxiter = 1e3, prec = 1e-7)
  beta_iht2 = model_iht2$beta
  #support_estimate = model_iht2$support_estimate
  # grad_iht2 = model_iht2$grad
  # steps_iht2 = model_iht2$steps_list
  
  #iht_cor2()
  # screening stage
  for (j in 1:p){
    feature_vec = as.matrix(X2[,j], ncol=1)
    alpha_selector = t(feature_vec) %*% (Y2 - X2[,-j] %*% as.matrix(beta_iht2[-j,], ncol=1))/norm(feature_vec,"2")
    threshold = a*norm(feature_vec, "2")/2 + v*log(p)/(a* norm(feature_vec,"2"))
    
    
    #print(c(abs(alpha_selector), threshold, norm(feature_vec,"2")))
    
    if(abs(alpha_selector)> threshold){
      support_estimate = c(j,support_estimate)
      #print("support entry success")
    }
    else{
      support_estimate = support_estimate
      #print("support entry failure")
    }
  }
  
  alpha_vec = c()
  for (j in 1:p){
    feature_vec = as.matrix(X2[,j], ncol=1)
    alpha_selector = t(feature_vec) %*% (Y2 - X2[,-j] %*% as.matrix(beta_iht2[-j,], ncol=1))/norm(feature_vec,"2")
    #threshold = a*norm(feature_vec, "2")/2 + v*log(p)/(a* norm(feature_vec,"2"))
    alpha_vec  = c(alpha_vec, alpha_selector)
    
    
  }
  alpha_vec_sort_idx = order(abs(alpha_vec), decreasing = T)
  support_estimate_no = alpha_vec_sort_idx[1:s_list]
  
  iht_screening_result = list(beta = beta_iht2, s_hat = s_hat, support_estimate = support_estimate, support_estimate_no = support_estimate_no, s_list = s_list , g_list = g_list, alpha_vec_sort_idx = alpha_vec_sort_idx)
  return(iht_screening_result)
  
  
}


iht_cv_grid_screening_select <- function(X, Y,  supp_true, beta0 = rep(0,p), gamma = 0.1,a, v, maxiter = 1e3, prec = 1e-7, lim_max = 20){
  ## the iterative hard thresholding algorithm
  ## Input: X, Y: data
  ##        scale_X (boolean): whether to scale X before training
  ##        s_list: list of possible values of parameter s
  ##        g_list (same length as s_list): list of possible values of parameter g
  ##        supp_true: true support of beta
  ##        beta0: initial value of beta
  ##        maxiter: maximum number of iterations 
  ##        prec: precision where the iterations stop
  ##        gamma: partition parameter. Determines the size of 1st subsample and 2nd subsample
  ##        a: minimum signal strength
  ##        v: true sparsity parameter
  ## Output: iht_screening_select: returns TPR and FPR
  n = dim(X)[1]
  p = dim(X)[2]
  model_iht_screening <- iht_cv_grid_screening(X = X, Y = Y, beta0 = rep(0,p), gamma = 0.7, a = a, v = v, maxiter = 1e3, prec = 1e-7, lim_max = lim_max)
  supp_beta_iht = model_iht_screening$support_estimate
  s_hat = model_iht_screening$s_hat
  alpha_vec_sort_idx = model_iht_screening$alpha_vec_sort_idx
  FDR = length(setdiff(supp_beta_iht, supp_true))/max(1,length(supp_beta_iht)); FDR
  
  TPR = length(intersect(supp_beta_iht, supp_true))/length(supp_true); TPR
  
  supp_beta_iht = model_iht_screening$support_estimate_no
  #s_hat = model_iht_screening$s_hat
  FDR_no = length(setdiff(supp_beta_iht, supp_true))/length(supp_beta_iht); FDR_no
  
  TPR_no = length(intersect(supp_beta_iht, supp_true))/length(supp_true); TPR_no
  
  FDR_path = c(); TPR_path = c()
  for (j in 1:p) {
   set = alpha_vec_sort_idx[1:j]
   
   FDR_path[j] = length(setdiff(set, supp_true))/max(1,length(set))
   TPR_path[j] = length(intersect(set, supp_true))/length(supp_true)
  }
  
  iht_screening_select_result = list(TPR=TPR, FDR=FDR, FDR_no = FDR_no, TPR_no = TPR_no, FDR_path = c(0,FDR_path), TPR_path = c(0,TPR_path) ,s_hat = s_hat)
  return(iht_screening_select_result)
}
##### IHT CV no oracle grid + screen

iht_cv_grid_screening_no_oracle <- function(X, Y, beta0 = rep(0,p), gamma = 0.7,a, v, maxiter = 1e3, prec = 1e-7, lim_max = 15){
  ## the iterative hard thresholding algorithm
  ## Input: X, Y: data
  ##        scale_X (boolean): whether to scale X before training
  ##        s_list: list of possible values of parameter s
  ##        g_list (same length as s_list): list of possible values of parameter g
  ##        beta0: initial value of beta
  ##        maxiter: maximum number of iterations 
  ##        prec: precision where the iterations stop
  ##        gamma: partition parameter. Determines the size of 1st subsample and 2nd subsample
  ##        a: minimum signal strength
  ##        v: true sparsity parameter
  ## Output: iht_screening_result: list containing beta (coefficients), grad (gradient), steps_list (number of iterations performed)
  
  n = dim(X)[1]
  p = dim(X)[2]
  #X = scale(X)
  
  alpha_vec = c()
  # partitoning the data
  n1 = floor(gamma * n)
  subset = sample(1:n, n1)
  X1 = X[subset,]; Y1 = as.matrix(Y[subset, ], ncol=1)
  X2 = X[-subset,]; Y2 = as.matrix(Y[-subset, ], ncol=1)
  
  # first setp IHT for 1st subsample with CV
  model_iht_cv = iht_cv_grid(X = X1, Y = Y1, s_list = 1: lim_max, g_list= 1:lim_max, beta0 = rep(0,p), nfold = 10, n_cv = 1, maxiter = 1e3, prec = 1e-7)
  s_list = model_iht_cv$s.min; g_list = model_iht_cv$g.min
  s_hat = sum(model_iht_cv$coef_min!=0)
  
  model_iht2 = iht_cor2(X = X1, Y = Y1, scale_X =  FALSE, s_list =  s_list, g_list = g_list, beta0 = rep(0,p), maxiter = 1e3, prec = 1e-7)
  beta_iht2 = model_iht2$beta
  #support_estimate = model_iht2$support_estimate
  # grad_iht2 = model_iht2$grad
  # steps_iht2 = model_iht2$steps_list
  
  #iht_cor2()
  # screening stage
  for (j in 1:p){
    feature_vec = as.matrix(X2[,j], ncol=1)
    alpha_selector = t(feature_vec) %*% (Y2 - X2[,-j] %*% as.matrix(beta_iht2[-j,], ncol=1))/norm(feature_vec,"2")
    #threshold = a*norm(feature_vec, "2")/2 + v*log(p)/(a* norm(feature_vec,"2"))
    alpha_vec  = c(alpha_vec, alpha_selector)
    
    
  }
  alpha_vec_sort_idx = order(abs(alpha_vec), decreasing = T)
  support_estimate = alpha_vec_sort_idx[1:s_list]
  iht_screening_result = list(beta = beta_iht2, s_hat = s_hat, support_estimate = support_estimate, s_list = s_list , g_list = g_list)
  return(iht_screening_result)
  
  
}


iht_cv_grid_screening_select_no_oracle <- function(X, Y,  supp_true, beta0 = rep(0,p), gamma = 0.1,a, v, maxiter = 1e3, prec = 1e-7, lim_max = 15){
  ## the iterative hard thresholding algorithm
  ## Input: X, Y: data
  ##        scale_X (boolean): whether to scale X before training
  ##        s_list: list of possible values of parameter s
  ##        g_list (same length as s_list): list of possible values of parameter g
  ##        supp_true: true support of beta
  ##        beta0: initial value of beta
  ##        maxiter: maximum number of iterations 
  ##        prec: precision where the iterations stop
  ##        gamma: partition parameter. Determines the size of 1st subsample and 2nd subsample
  ##        a: minimum signal strength
  ##        v: true sparsity parameter
  ## Output: iht_screening_select: returns TPR and FPR
  
  model_iht_screening <- iht_cv_grid_screening_no_oracle(X = X, Y = Y, beta0 = rep(0,p), gamma = 0.7, a = a, v = v, maxiter = 1e3, prec = 1e-7, lim_max = lim_max)
  supp_beta_iht = model_iht_screening$support_estimate
  s_hat = model_iht_screening$s_hat
  FDR = length(setdiff(supp_beta_iht, supp_true))/length(supp_beta_iht); FDR
  
  TPR = length(intersect(supp_beta_iht, supp_true))/length(supp_true); TPR
  
  iht_screening_select_result = list(TPR=TPR, FDR=FDR, s_hat = s_hat)
  return(iht_screening_select_result)
}






#########

iht_one_screening <- function(X, Y, s_list, eta_list, beta0 = rep(0,p), gamma = 0.7,a, v, maxiter = 1e3, prec = 1e-7){
  ## the iterative hard thresholding algorithm
  ## Input: X, Y: data
  ##        scale_X (boolean): whether to scale X before training
  ##        s_list: list of possible values of parameter s
  ##        eta_list (same length as s_list): list of possible step size eta
  ##        beta0: initial value of beta
  ##        maxiter: maximum number of iterations 
  ##        prec: precision where the iterations stop
  ##        gamma: partition parameter. Determines the size of 1st subsample and 2nd subsample
  ##        a: minimum signal strength
  ##        v: true sparsity parameter
  ## Output: iht_one_screening_result: list containing beta (coefficients),  steps_list (number of iterations performed)
  
  n = dim(X)[1]
  p = dim(X)[2]
  #X = scale(X)
  
  support_estimate = c()
  # partitoning the data
  n1 = floor(gamma * n)
  subset = sample(1:n, n1)
  X1 = X[subset,]; Y1 = as.matrix(Y[subset, ], ncol=1)
  X2 = X[-subset,]; Y2 = as.matrix(Y[-subset, ], ncol=1)
  
  # first setp IHT for 1st subsample
  
  model_iht1 = iht_one(X1, Y1, scale_X = FALSE, s_list=s_list, eta_list = eta_list, beta0=beta0, maxiter = maxiter, prec = prec)
  beta_iht1 = model_iht1$beta
  steps_iht1 = model_iht1$steps_list
  
  
  
  # screening stage
  for (j in 1:p){
    feature_vec = as.matrix(X2[,j], ncol=1)
    alpha_selector = t(feature_vec) %*% (Y2 - X2[,-j] %*% as.matrix(beta_iht1[-j,], ncol=1))/norm(feature_vec,"2")
    threshold = a*norm(feature_vec, "2")/2 + v*log(p)/(a* norm(feature_vec,"2"))
    
    
    #print(c(abs(alpha_selector), threshold, norm(feature_vec,"2")))
    
    if(abs(alpha_selector)> threshold){
      support_estimate = c(j,support_estimate)
      #print("support entry success")
    }
    else{
      support_estimate = support_estimate
      #print("support entry failure")
    }
  }
  
  iht_one_screening_result = list(beta = beta_iht1, steps_list = steps_iht1, support_estimate = support_estimate)
  return(iht_one_screening_result)
  
  
}



iht_one_screening_select <- function(X, Y, s_list,  eta_list, supp_true, beta0 = rep(0,p), gamma = 0.1,a, v, maxiter = 1e3, prec = 1e-7){
  ## the iterative hard thresholding algorithm
  ## Input: X, Y: data
  ##        scale_X (boolean): whether to scale X before training
  ##        s_list: list of possible values of parameter s
  ##        g_list (same length as s_list): list of possible values of parameter g
  ##        supp_true: true support of beta
  ##        beta0: initial value of beta
  ##        maxiter: maximum number of iterations 
  ##        prec: precision where the iterations stop
  ##        gamma: partition parameter. Determines the size of 1st subsample and 2nd subsample
  ##        a: minimum signal strength
  ##        v: true sparsity parameter
  ## Output: iht_screening_select: returns TPR and FPR
  
  model_iht_screening <- iht_one_screening(X, Y, s_list = s_list, eta_list = eta_list, beta0 = rep(0,p), gamma = 0.7,a, v, maxiter = 1e3, prec = 1e-7)
  supp_beta_iht = which(model_iht_screening$beta != 0)
  FDR = length(setdiff(supp_beta_iht, supp_true))/length(supp_beta_iht)
  
  TPR = length(intersect(supp_beta_iht, supp_true))/length(supp_true)
  
  iht_screening_select_result = list(TPR=TPR, FDR=FDR)
  return(iht_screening_select_result)
}


#########

iht_one_cv_screening <- function(X, Y, s_list, eta_list, beta0 = rep(0,p), gamma = 0.7,a, v, maxiter = 1e3, prec = 1e-7){
  ## the iterative hard thresholding algorithm
  ## Input: X, Y: data
  ##        scale_X (boolean): whether to scale X before training
  ##        s_list: list of possible values of parameter s
  ##        eta_list (same length as s_list): list of possible step size eta
  ##        beta0: initial value of beta
  ##        maxiter: maximum number of iterations 
  ##        prec: precision where the iterations stop
  ##        gamma: partition parameter. Determines the size of 1st subsample and 2nd subsample
  ##        a: minimum signal strength
  ##        v: true sparsity parameter
  ## Output: iht_one_screening_result: list containing beta (coefficients),  steps_list (number of iterations performed)
  
  n = dim(X)[1]
  p = dim(X)[2]
  #X = scale(X)
  
  support_estimate = c()
  # partitoning the data
  if (gamma <1 && gamma >0){
  n1 = floor(gamma * n)
  subset = sample(1:n, n1)
  X1 = X[subset,]; Y1 = as.matrix(Y[subset, ], ncol=1)
  X2 = X[-subset,]; Y2 = as.matrix(Y[-subset, ], ncol=1)
  }
  else{
    X1 = X; X2 = X
    Y1 = Y; Y2 = Y
  }
  # first setp IHT for 1st subsample
  
  model_iht1 = iht_one_cv(X = X1, Y = Y1, s_list=s_list, eta_list = eta_list, beta0=beta0, maxiter = maxiter, prec = prec)
  beta_iht1 = model_iht1$coef_min
  steps_iht1 = model_iht1$steps_list
  s_hat = model_iht1$s.min
  
  
  # screening stage
#  for (j in 1:p){
#    feature_vec = as.matrix(X2[,j], ncol=1)
#    alpha_selector = t(feature_vec) %*% (Y2 - X2[,-j] %*% as.matrix(beta_iht1[-j,], ncol=1))/norm(feature_vec,"2")
#    threshold = a*norm(feature_vec, "2")/2 + v*log(p)/(a* norm(feature_vec,"2"))
    
    
    #print(c(abs(alpha_selector), threshold, norm(feature_vec,"2")))
    
#    if(abs(alpha_selector)> threshold){
#      support_estimate = c(j,support_estimate)
      #print("support entry success")
#    }
#    else{
#      support_estimate = support_estimate
      #print("support entry failure")
#    }
#  }
  
 # alpha_vec = c()
  fun = function(j){
    feature_vec = as.matrix(X2[,j], ncol=1)
    alpha_selector = t(feature_vec) %*% (Y2 - X2[,-j] %*% as.matrix(beta_iht1[-j,], ncol=1))/norm(feature_vec,"2")
    #threshold = a*norm(feature_vec, "2")/2 + v*log(p)/(a* norm(feature_vec,"2"))
    #alpha_vec  = c(alpha_vec, alpha_selector)
    return(alpha_selector)
    
  }
  alpha_vec = unlist(mclapply(1:p , FUN = fun))
  alpha_vec_sort_idx = order(abs(alpha_vec), decreasing = T)
  #support_estimate_no = alpha_vec_sort_idx[1:s_hat]
  
  iht_screening_result = list(beta = beta_iht1, s_hat = s_hat, support_estimate = support_estimate, s_list = s_list , eta_list = eta_list, alpha_vec_sort_idx = alpha_vec_sort_idx)
  return(iht_screening_result)
  
  # iht_one_screening_result = list(beta = beta_iht1, steps_list = steps_iht1, support_estimate = support_estimate)
  #return(iht_one_screening_result)
  
  
}



iht_one_screening_select <- function(X, Y, s_list,  eta_list, supp_true, beta0 = rep(0,p), gamma = 0.7,a, v, maxiter = 1e3, prec = 1e-7){
  ## the iterative hard thresholding algorithm
  ## Input: X, Y: data
  ##        scale_X (boolean): whether to scale X before training
  ##        s_list: list of possible values of parameter s
  ##        g_list (same length as s_list): list of possible values of parameter g
  ##        supp_true: true support of beta
  ##        beta0: initial value of beta
  ##        maxiter: maximum number of iterations 
  ##        prec: precision where the iterations stop
  ##        gamma: partition parameter. Determines the size of 1st subsample and 2nd subsample
  ##        a: minimum signal strength
  ##        v: true sparsity parameter
  ## Output: iht_screening_select: returns TPR and FPR
  
  model_iht_screening <- iht_one_cv_screening(X = X, Y = Y, s_list = s_list, eta_list = eta_list, beta0 = rep(0,p), gamma = 0.7,a = a, v = 1, maxiter = 1e3, prec = 1e-7)
  supp_beta_iht = model_iht_screening$support_estimate
  FDR = length(setdiff(supp_beta_iht, supp_true))/length(supp_beta_iht)
  
  TPR = length(intersect(supp_beta_iht, supp_true))/length(supp_true)
  
  #supp_beta_iht = model_iht_screening$support_estimate_no
  s_hat = model_iht_screening$s_hat
  #FDR_no = length(setdiff(supp_beta_iht, supp_true))/length(supp_beta_iht); FDR_no
  
  #TPR_no = length(intersect(supp_beta_iht, supp_true))/length(supp_true); TPR_no
  
  alpha_vec_sort_idx = model_iht_screening$alpha_vec_sort_idx
  FDR_path = c(); TPR_path = c()
  for (j in 1:p) {
    set = alpha_vec_sort_idx[1:j]
    
    FDR_path[j] = length(setdiff(set, supp_true))/max(1,length(set))
    TPR_path[j] = length(intersect(set, supp_true))/length(supp_true)
  }
  
  iht_screening_select_result = list(TPR=TPR, FDR=FDR, FDR_path = c(0,FDR_path), TPR_path = c(0,TPR_path) ,s_hat = s_hat)
  return(iht_screening_select_result)
  
  #iht_screening_select_result = list(TPR=TPR, FDR=FDR)
  #return(iht_screening_select_result)
}


lasso_screening <- function(X, Y, lambda, gamma = 0.7, lambda_val,a, v, do.CV = F){
  ## the iterative hard thresholding algorithm
  ## Input: X, Y: data
  ##        scale_X (boolean): whether to scale X before training
  ##        s_list: list of possible values of parameter s
  ##        eta_list (same length as s_list): list of possible step size eta
  ##        beta0: initial value of beta
  ##        maxiter: maximum number of iterations 
  ##        prec: precision where the iterations stop
  ##        gamma: partition parameter. Determines the size of 1st subsample and 2nd subsample
  ##        a: minimum signal strength
  ##        v: true sparsity parameter
  ## Output: iht_one_screening_result: list containing beta (coefficients),  steps_list (number of iterations performed)
  
  n = dim(X)[1]
  p = dim(X)[2]
  #X = scale(X)
  
  support_estimate = c()
  # partitoning the data
  if (gamma <1 && gamma >0){
    n1 = floor(gamma * n)
    subset = sample(1:n, n1)
    X1 = X[subset,]; Y1 = as.matrix(Y[subset, ], ncol=1)
    X2 = X[-subset,]; Y2 = as.matrix(Y[-subset, ], ncol=1)
  }
  else{
    X1 = X; X2 = X
    Y1 = Y; Y2 = Y
  }
  
  # first setp IHT for 1st subsample
  if(do.CV == T){
  model_lasso_cv = cv.glmnet(x = X1, y = Y1, intercept = F, standardize = F, alpha =1 , nfolds = 10)
  lambda_fit = model_lasso_cv$lambda.min
  beta_lasso = glmnet(x= X1, y = Y1, lambda = lambda_fit, standardize = F, intercept = F)$beta
  print(lambda_fit)
  }
  else{
    lambda_fit = lambda_val
    beta_lasso = glmnet(x= X1, y = Y1, lambda = lambda_fit, standardize = F, intercept = F)$beta
    
  }
  
  # screening stage
  #  for (j in 1:p){
  #    feature_vec = as.matrix(X2[,j], ncol=1)
  #    alpha_selector = t(feature_vec) %*% (Y2 - X2[,-j] %*% as.matrix(beta_iht1[-j,], ncol=1))/norm(feature_vec,"2")
  #    threshold = a*norm(feature_vec, "2")/2 + v*log(p)/(a* norm(feature_vec,"2"))
  
  
  #print(c(abs(alpha_selector), threshold, norm(feature_vec,"2")))
  
  #    if(abs(alpha_selector)> threshold){
  #      support_estimate = c(j,support_estimate)
  #print("support entry success")
  #    }
  #    else{
  #      support_estimate = support_estimate
  #print("support entry failure")
  #    }
  #  }
  
  #alpha_vec = c()
  #for (j in 1:p){
  #  feature_vec = as.matrix(X2[,j], ncol=1)
  #  alpha_selector = t(feature_vec) %*% (Y2 - X2[,-j] %*% as.matrix(beta_lasso[-j,], ncol=1))/norm(feature_vec,"2")
  #  #threshold = a*norm(feature_vec, "2")/2 + v*log(p)/(a* norm(feature_vec,"2"))
  #  alpha_vec  = c(alpha_vec, alpha_selector)
    
    
  #}
  fun = function(j){
    feature_vec = as.matrix(X2[,j], ncol=1)
    alpha_selector = t(feature_vec) %*% (Y2 - X2[,-j] %*% as.matrix(beta_lasso[-j,], ncol=1))/norm(feature_vec,"2")
    return(alpha_selector)
  }
  alpha_vec = unlist(mclapply(1:p, FUN = fun))
  alpha_vec_sort_idx = order(abs(alpha_vec), decreasing = T)
  #support_estimate_no = alpha_vec_sort_idx[1:s_hat]
  
  lasso_screening_result = list(beta = beta_lasso, support_estimate = support_estimate, alpha_vec_sort_idx = alpha_vec_sort_idx)
  return(lasso_screening_result)
  
  # iht_one_screening_result = list(beta = beta_iht1, steps_list = steps_iht1, support_estimate = support_estimate)
  #return(iht_one_screening_result)
  
  
}

