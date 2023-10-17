library(parallel)
fm_cov <- function(p, eig_list){
    ## generate a covariance matrix from the factor model
    ## Input: p: dimension
    ##        eig_list: vector of eigenvalues of the low rank component Sigma_b
    ## Output: Sigma_f: covariance matrix (Sigma_f = Sigma_b + diag(p))
    set.seed(2020)
    # construct Sigma_b
    K = length(eig_list)
    V_0 = matrix(rnorm(n = p * K, mean = 0, sd = 1), nrow = p, ncol = K)
    V = qr.Q(qr(V_0))
    Sigma_b = V %*% diag(eig_list) %*% t(V)
    # construct Sigma_u
    Sigma_u = diag(p)
    # obtain covariance
    Sigma_f = Sigma_b + Sigma_u
    return(Sigma_f)
}


data_lm_q_spiky_fixed_snr <- function(n, mu_x, S_x, s, betamin, sigma_e, snr = 10, q=1){
  ## generate linear model data
  ## Input: n: sample size
  ##        mu_x: mean of of each observation x_i
  ##        S_x: Covariance matrix of x_i
  ##        beta_min: minimum possible absolute value of entries of beta
  ##        sigma_e: standard deviation of the noise
  ##        q = probablity of signal being at the same level. q=1 means all signal are coming from minimum level.
  ##        tp = maximum signal magnitude
  ## Output: data_lm: a list containing X, Y, beta, supp_true (true support), epsilon (noise)
  p = length(mu_x)
  X = matrix(rnorm(n*p), ncol = p)
  
  beta0 = matrix(betamin, nrow = s, ncol = 1)
  if(q>0){
  spiky = (snr - (s-q)*a^2)/q
  beta0[1:q] = sqrt(spiky)
  }
  beta = matrix(0, nrow = p, ncol = 1)
  supp_true = sample(p, s)
  beta[supp_true] = beta0
  epsilon = as.matrix(rnorm(n, mean = 0, sd = sigma_e))
  Y = X %*% beta + epsilon
  data= list(X = X, Y = Y, beta = beta, supp_true = supp_true, epsilon = epsilon)
  return(data)
}

data_lm <- function(n, mu_x, S_x, s, betamin, sigma_e, tp, q=1){
  ## generate linear model data
  ## Input: n: sample size
  ##        mu_x: mean of of each observation x_i
  ##        S_x: Covariance matrix of x_i
  ##        beta_min: minimum possible absolute value of entries of beta
  ##        sigma_e: standard deviation of the noise
  ##        q = probablity of signal being at the same level. q=1 means all signal are coming from minimum level.
  ##        tp = maximum signal magnitude
  ## Output: data_lm: a list containing X, Y, beta, supp_true (true support), epsilon (noise)
  p = length(mu_x)
  X = mvrnorm(n, mu = mu_x, Sigma = S_x)
  choice_mat = as.matrix(rbinom(s, 1, q))
  beta0 = as.matrix((1 + rnorm(s)^2/n)^1) * choice_mat * betamin +  (1- choice_mat) * tp
  beta = matrix(0, nrow = p, ncol = 1)
  supp_true = sample(p, s)
  beta[supp_true] = beta0
  epsilon = as.matrix(rnorm(n, mean = 0, sd = sigma_e))
  Y = X %*% beta + epsilon
  data= list(X = X, Y = Y, beta = beta, supp_true = supp_true, epsilon = epsilon)
  return(data)
  }

data_lm2 <- function(n, mu_x, S_x, s, betamin, sigma_e, tp){
  ## generate linear model data
  ## Input: n: sample size
  ##        mu_x: mean of of each observation x_i
  ##        S_x: Covariance matrix of x_i
  ##        beta_min: minimum possible absolute value of entries of beta
  ##        sigma_e: standard deviation of the noise
  ##        q = probablity of signal being at the same level. q=1 means all signal are coming from minimum level.
  ##        tp = maximum signal magnitude
  ## Output: data_lm: a list containing X, Y, beta, supp_true (true support), epsilon (noise)
  p = length(mu_x)
  X = mvrnorm(n, mu = mu_x, Sigma = S_x)
  #choice_mat = as.matrix(rbinom(s, 1, q))
  beta0 = as.matrix((1 + rnorm(s)^2/n)^1)  * betamin 
  beta0[sample(s,1)] = tp
  beta = matrix(0, nrow = p, ncol = 1)
  supp_true = sample(p, s)
  beta[supp_true] = beta0
  epsilon = as.matrix(rnorm(n, mean = 0, sd = sigma_e))
  Y = X %*% beta + epsilon
  data = list(X = X, Y = Y, beta = beta, supp_true = supp_true, epsilon = epsilon)
  return(data)
}
data_lm_newscheme <- function(n, mu_x, S_x, s, betamin, sigma_e, tp, q=1){
  ## generate linear model data
  ## Input: n: sample size
  ##        mu_x: mean of of each observation x_i
  ##        S_x: Covariance matrix of x_i
  ##        beta_min: minimum possible absolute value of entries of beta
  ##        sigma_e: standard deviation of the noise
  ##        q = probablity of signal being at the same level. q=1 means all signal are coming from minimum level.
  ##        tp = maximum signal magnitude
  ## Output: data_lm: a list containing X, Y, beta, supp_true (true support), epsilon (noise)
  p = length(mu_x)
  X = mvrnorm(n, mu = mu_x, Sigma = S_x)
  choice_mat = as.matrix(rbinom(s, 1, q))
  beta0 = as.matrix((1 + rnorm(s)^2/n)^0.5) * choice_mat * betamin +  (1- choice_mat) * tp
  beta = matrix(0, nrow = p, ncol = 1)
  supp_true = sample(p, s)
  beta[supp_true] = beta0
  epsilon = as.matrix(rnorm(n, mean = 0, sd = sigma_e))
  Y = X %*% beta + epsilon
  data= list(X = X, Y = Y, beta = beta, supp_true = supp_true, epsilon = epsilon)
  return(data)
}

data_lm2_newscheme <- function(n, mu_x, S_x, s, betamin, sigma_e, tp){
  ## generate linear model data
  ## Input: n: sample size
  ##        mu_x: mean of of each observation x_i
  ##        S_x: Covariance matrix of x_i
  ##        beta_min: minimum possible absolute value of entries of beta
  ##        sigma_e: standard deviation of the noise
  ##        q = probablity of signal being at the same level. q=1 means all signal are coming from minimum level.
  ##        tp = maximum signal magnitude
  ## Output: data_lm: a list containing X, Y, beta, supp_true (true support), epsilon (noise)
  p = length(mu_x)
  X = mvrnorm(n, mu = mu_x, Sigma = S_x)
  #choice_mat = as.matrix(rbinom(s, 1, q))
  beta0 = as.matrix((1 + rnorm(s)^2/n)^0.5)  * betamin 
  beta0[sample(s,1)] = tp
  beta = matrix(0, nrow = p, ncol = 1)
  supp_true = sample(p, s)
  beta[supp_true] = beta0
  epsilon = as.matrix(rnorm(n, mean = 0, sd = sigma_e))
  Y = X %*% beta + epsilon
  data = list(X = X, Y = Y, beta = beta, supp_true = supp_true, epsilon = epsilon)
  return(data)
}

data_lm2_newscheme_SNR <- function(n, mu_x, S_x, s, betamin, sigma_e, snr){
  ## generate linear model data
  ## Input: n: sample size
  ##        mu_x: mean of of each observation x_i
  ##        S_x: Covariance matrix of x_i
  ##        beta_min: minimum possible absolute value of entries of beta
  ##        sigma_e: standard deviation of the noise
  ##        q = probablity of signal being at the same level. q=1 means all signal are coming from minimum level.
  ##        tp = maximum signal magnitude
  ## Output: data_lm: a list containing X, Y, beta, supp_true (true support), epsilon (noise)
  p = length(mu_x)
  X = mvrnorm(n, mu = mu_x, Sigma = S_x)
  #choice_mat = as.matrix(rbinom(s, 1, q))
  beta0 = as.matrix((1 + rnorm(s)^2/n)^0.5)  * betamin
  T1 = snr - sum(beta0[1:(s-1) ]^2)  
  beta0[s] = sqrt(T1)
  beta = matrix(0, nrow = p, ncol = 1)
  supp_true = sample(p, s)
  beta[supp_true] = beta0
  epsilon = as.matrix(rnorm(n, mean = 0, sd = sigma_e))
  Y = X %*% beta + epsilon
  data = list(X = X, Y = Y, beta = beta, supp_true = supp_true, epsilon = epsilon)
  return(data)
}



find_ind <- function(y,x, ind_x){
  # y = array whose indices to be found
  # x = mother array inside which y lies
  # ind_x = indices of y should be picked from ind_x
  
  id = unlist(apply(as.matrix(y), 1, function(c) which(as.matrix(x) == c)))
  return( intersect(id, ind_x))
}


iht_one <- function(X, Y, s_list, eta_list, beta0, maxiter = 1e3, prec = 1e-7){
  ## the iterative hard thresholding algorithm
  ## Input: X, Y: data
  ##        scale_X (boolean): whether to scale X before training
  ##        s_list: list of possible values of parameter s
  ##        eta_list (same length as s_list): list of possible step-sizes eta
  ##        beta0: initial value of beta
  ##        maxiter: maximum number of iterations 
  ##        prec: precision where the iterations stop
  ## Output: iht_cor_result: list containing beta (coefficients), steps_list (number of iterations performed)
  n = dim(X)[1]
  p = dim(X)[2]
#  if (scale_X) 
#    X = scale(X)
  #Sxx = t(X) %*% X / n
  Sxy = t(X) %*% Y / n
  s_list_len = length(s_list)
  steps_list = rep(0, s_list_len)
  beta_iht_one = matrix(0, nrow=p, ncol=s_list_len)
  
  
  for (k in 1:s_list_len){
    s_iht = s_list[k]
    eta_iht = eta_list
    t = 0
    betat = beta0
    while (t < maxiter){
      beta_old = betat
      grad = - Sxy + t(X)%*% (X %*% beta_old)/n
      beta_mid = beta_old - eta_iht * grad 
      betat = rep(0,p)
      beta_mid_sort = sort(abs(beta_mid), decreasing = TRUE)
      indt2 = which(abs(beta_mid)>= beta_mid_sort[s_iht])
      betat[indt2] = beta_mid[indt2]
      
    
      # identify convergence condition
      if (sum((betat - beta_old)^2) < prec * (sum((beta_old)^2) + 1)) 
        break
      t = t + 1
    }
    beta_iht_one[, k] = betat
    steps_list[k] = t
  }
  
  iht_one_result = list(beta=beta_iht_one,  steps_list=steps_list)
  return(iht_one_result)
}

iht_one_par <- function(X, Y, s_list, eta_list, beta0, maxiter = 1e3, prec = 1e-7){
  ## the iterative hard thresholding algorithm
  ## Input: X, Y: data
  ##        scale_X (boolean): whether to scale X before training
  ##        s_list: list of possible values of parameter s
  ##        eta_list (same length as s_list): list of possible step-sizes eta
  ##        beta0: initial value of beta
  ##        maxiter: maximum number of iterations 
  ##        prec: precision where the iterations stop
  ## Output: iht_cor_result: list containing beta (coefficients), steps_list (number of iterations performed)
  loc.env = new.env()
  n = dim(X)[1]
  p = dim(X)[2]
  #  if (scale_X) 
  #    X = scale(X)
  #Sxx = t(X) %*% X / n
  Sxy = t(X) %*% Y / n
  s_list_len = length(s_list)
  loc.env$steps_list = rep(0, s_list_len)
  loc.env$beta_iht_one = matrix(0, nrow=p, ncol=s_list_len)
  
  
  fun = function(k){
    s_iht = s_list[k]
    eta_iht = eta_list
    t = 0
    betat = beta0
    while (t < maxiter){
      beta_old = betat
      grad = - Sxy + t(X)%*% (X %*% beta_old)/n
      beta_mid = beta_old - eta_iht * grad 
      betat = rep(0,p)
      beta_mid_sort = sort(abs(beta_mid), decreasing = TRUE)
      indt2 = which(abs(beta_mid)>= beta_mid_sort[s_iht])
      betat[indt2] = beta_mid[indt2]
      
      
      # identify convergence condition
      if (sum((betat - beta_old)^2) < prec * (sum((beta_old)^2) + 1)) 
        break
      t = t + 1
    }
    loc.env$beta_iht_one[, k] = betat
    loc.env$steps_list[k] = t
  }
  mclapply(1:s_list_len, fun(k))
  iht_one_result = list(beta= loc.env$beta_iht_one,  steps_list= loc.env$steps_list)
  return(iht_one_result)
}


iht_cor2 <- function(X, Y, scale_X = TRUE, s_list, g_list, beta0, maxiter = 1e3, prec = 1e-7){
    ## the iterative hard thresholding algorithm
    ## Input: X, Y: data
    ##        scale_X (boolean): whether to scale X before training
    ##        s_list: list of possible values of parameter s
    ##        g_list (same length as s_list): list of possible values of parameter g
    ##        beta0: initial value of beta
    ##        maxiter: maximum number of iterations 
    ##        prec: precision where the iterations stop
    ## Output: iht_cor_result: list containing beta (coefficients), grad (gradient), steps_list (number of iterations performed)
    n = dim(X)[1]
    p = dim(X)[2]
    if (scale_X) 
        X = scale(X)
    Sxx = t(X) %*% X / n
    Sxy = t(X) %*% Y / n
    s_list_len = length(s_list)
    steps_list = rep(0, s_list_len)
    beta_iht_cor = matrix(0, nrow=p, ncol=s_list_len)
    grad_iht_cor = matrix(0, nrow=p, ncol=s_list_len)

    for (k in 1:s_list_len){
        s_iht = s_list[k]
        g_iht = g_list[k]
        t = 0
        betat = beta0
        while (t < maxiter){
            beta_old = betat
            grad = - Sxy + Sxx %*% beta_old
            indt1 = which(beta_old != 0)
            indt_null = which(beta_old == 0)
            #grad_1 = grad
            grad_1 = grad[indt_null] # this is edited line; 
            indt2 = order(-abs(grad_1))[1:g_iht]
            #grad_sort = sort(abs(grad_1), decreasing=TRUE)
            #indt2 = which(abs(grad_1) >= grad_sort[g_iht]) ; #cat("projected prev grad ind", indt2 ,"\n")
            indt2 = find_ind(grad_1[indt2], grad, indt_null);  # This is edited line
            indt = union(indt1, indt2)
       

            # refit 1
            Xt = X[, indt]
            betat = rep(0, p)
            betat[indt] = solve(t(Xt) %*% Xt) %*% (t(Xt) %*% Y)
            
            # truncation 
            #betat_sort = sort(abs(betat), decreasing=TRUE)
            indt0 = order(-abs(betat))[1:s_iht]
            #indt0 = which(abs(betat) >= betat_sort[s_iht])
            
            # refit 2
            Xt0 = X[, indt0]
            betat = rep(0, p)
            betat[indt0] = solve(t(Xt0) %*% Xt0) %*% (t(Xt0) %*% Y)
            
            # identify convergence condition
            if (sum((betat - beta_old)^2) < prec * (sum((beta_old)^2) + 1)) 
                break
            t = t + 1
        }
        beta_iht_cor[, k] = betat
        grad_iht_cor[, k] = - Sxy + Sxx %*% betat
        steps_list[k] = t
    }

    iht_cor_result = list(beta=beta_iht_cor, grad= grad_iht_cor, steps_list=steps_list)
    return(iht_cor_result)
}

marginal_screening <- function(X, Y, s_list, scale_X = TRUE){
    ## the iterative hard thresholding algorithm
    ## Input: X, Y: data
    ##        s_list: list of values of parameter s
    ##        scale_X (boolean): whether to scale X before training
    ## Output: ms_result: list of beta (coef), grad (gradient)
    n = dim(X)[1]
    p = dim(X)[2]
    beta0 = rep(0, p)
    if (scale_X) 
        X = scale(X)
    Sxx = t(X) %*% X / n
    Sxy = t(X) %*% Y / n
    s_list_len = length(s_list)
    beta_ms = matrix(0, nrow=p, ncol=s_list_len)
    grad_ms = matrix(0, nrow=p, ncol=s_list_len)
  
    for (k in 1:s_list_len){
        s_ms = s_list[k]
        t = 0
        #betat = beta0
        beta_old = beta0
        grad_1 =  - Sxy
        grad_sort = sort(abs(grad_1), decreasing=TRUE)
        indt = which(abs(grad_1) >= grad_sort[s_ms])

        # refit 1
        Xt = X[, indt]
        betat = rep(0, p)
        betat[indt] = solve(t(Xt) %*% Xt) %*% (t(Xt) %*% Y)
      
        beta_ms[, k] = betat
        grad_ms[, k] = - Sxy
    }
  
    ms_result = list(beta=beta_ms, grad=grad_ms)
    return(ms_result)
}

lambda_select_scad <- function(n, mu_x, S_x, s, betamin, sigma_e, nlambda = 100, lambda_ratio = 1.2, ntest, q=1){
    ## select a list of lambda for SCAD
    ## Input: mu_x: mean of x_i
    ##        S_x: variance of x_i
    ##        s: sparsity
    ##        betamin: minimum possible absolute value of beta
    ##        sigma_e: standard deviation of noise
    ##        nlambda: number of lambda's to be chosen
    ##        lambda_ratio: the ratio between adjacent lambda's
    ##        ntest: number of tests to perform to evaluate the range of lambda
    ## Output: lambda_list: list of lambda
    lambda_max = 0
    for (ntest in 1:ceiling(nrep/10)){
        data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                          s = s, betamin = betamin, sigma_e = sigma_e, q = q)
        X = data_i$X
        Y = data_i$Y
        model_scad = picasso(X, Y, nlambda=nlambda, 
                             method="scad", intercept=FALSE)
        lambda_scad = model_scad$lambda
        lambda_max = max(max(lambda_scad),lambda_max)
    }
    lambda_max = lambda_max * 2
    ratios = lambda_ratio ^ seq(0, nlambda-1)
    lambda_list = lambda_max / ratios
    return(lambda_list)
}

lambda_select_lasso <- function(n, mu_x, S_x, s, betamin, sigma_e, nlambda = 100, lambda_ratio = 1.2, ntest, tp, q = 1, pkg = "glmnet", data_gen_method = 2){
  ## select a list of lambda for LASSO
  ## Input: mu_x: mean of x_i
  ##        S_x: variance of x_i
  ##        s: sparsity
  ##        betamin: minimum possible absolute value of beta
  ##        sigma_e: standard deviation of noise
  ##        nlambda: number of lambda's to be chosen
  ##        lambda_ratio: the ratio between adjacent lambda's
  ##        ntest: number of tests to perform to evaluate the range of lambda
  ## Output: lambda_list: list of lambda
  lambda_max = 0
  for (ntest in 1:ceiling(nrep/10)){
    if (data_gen_method == 1){
      data_i <- data_lm(n = n, mu_x = rep(0,p), S_x = S_x,
                        s = s, betamin = betamin, sigma_e = sigma_e, tp , q)
    }
    else{
      data_i <- data_lm2(n = n, mu_x = rep(0,p), S_x = S_x,
                         s = s, betamin = betamin, sigma_e = sigma_e, tp )
    }
    
    X = data_i$X
    Y = data_i$Y
    if (pkg == "glmnet"){
      model_scad = glmnet(X, Y, nlambda=nlambda, 
                          family = "gaussian", standardize = F, intercept = F)
    }
    else{
      model_scad = picasso(X, Y, nlambda=nlambda, 
                           method="l1", intercept=FALSE)
    }
    lambda_scad = model_scad$lambda
    lambda_max = max(max(lambda_scad),lambda_max)
  }
  lambda_max = lambda_max * 2
  ratios = lambda_ratio ^ seq(0, nlambda-1)
  lambda_list = lambda_max / ratios
  return(lambda_list)
}

scad_select <- function(X, Y, supp_true, lambda_list, nfold = 10){
    ## variable selection using SCAD
    ## Input: X, Y: data
    ##        supp_true: the true support (not used during training)
    ##        lambda_list: list of lambda values for training
    ##        nfold: number of folds for cross validation
    ## Output: scad_result: list of TPR, FDR and their values at the best lambda value
    s = length(supp_true)
    n = dim(X)[1]
    p = dim(X)[2]
    
    ## TPR and FDR
    model_scad = picasso(X, Y, lambda=lambda_list, 
                         method="scad", intercept=FALSE)
    beta_scad = model_scad$beta
    beta_supp = apply(beta_scad, 2, function(c) sum(c != 0))
    TP = apply(beta_scad, 2, function(c) length(intersect(which(c != 0), supp_true)))
    FP = beta_supp - TP
    TPR = TP / s
    FN = s - TP
    Ham = FN + FP
    
    # deal with cases with zero supp
    beta_supp1 = beta_supp
    beta_supp1[which(beta_supp1 == 0)] = 1
    FDR = FP / beta_supp1
    
    ## cross validation 
    cv_score = rep(0, length(lambda_list))
    n0 = floor(n / nfold) * nfold
    ind_sample = matrix(sample(n0), nrow = nfold, ncol = n0 / nfold)
    for (k in 1:length(lambda_list)){
        lambda = lambda_list[k]
        for (j in 1:nfold){
            ind_test = ind_sample[j, ]
            ind_train = setdiff(c(1:n0), ind_test)
            X_test = X[ind_test, ]
            Y_test = Y[ind_test, ]
            X_train = X[ind_train, ]
            Y_train = Y[ind_train, ]
            model_train = picasso(X_train, Y_train, lambda=lambda, 
                                  method="scad", intercept=FALSE)
            beta_train = model_train$beta
            residual = as.matrix(Y_test - X_test %*%beta_train)
            cv_score[k] = cv_score[k] + sum((residual)^2)/nfold 
        }
    }
    ind_cv = which(cv_score == min(cv_score))
    TPR_cv = TPR[ind_cv]
    FDR_cv = FDR[ind_cv]

    scad_result = list(TPR = TPR, FDR = FDR, TPR_cv = TPR_cv, FDR_cv = FDR_cv, Ham_loss = Ham)
    return(scad_result)
}

lasso_select <- function(X, Y, supp_true, lambda_list, nfold = 10){
    ## variable selection using LASSO
    ## Input: X, Y: data
    ##        supp_true: the true support (not used during training)
    ##        lambda_list: list of lambda values for training
    ##        nfold: number of folds for cross validation
    ## Output: scad_result: list of TPR, FDR and their values at the best lambda value
  
    s = length(supp_true)
    n = dim(X)[1]
    p = dim(X)[2]
  
    ## TPR and FDR
    model_lasso = picasso(X, Y, lambda=lambda_list, 
                          method="l1", intercept=FALSE)
    beta_lasso = model_lasso$beta
    beta_supp = apply(beta_lasso, 2, function(c) sum(c != 0))
    TP = apply(beta_lasso, 2, function(c) length(intersect(which(c != 0), supp_true)))
    FP = beta_supp - TP
    TPR = TP / s
    FN = s - TP
    Ham = FN + FP
    # deal with cases with zero supp
    beta_supp1 = beta_supp
    beta_supp1[which(beta_supp1 == 0)] = 1
    FDR = FP / beta_supp1
    
    ## cross validation 
    cv_score = rep(0, length(lambda_list))
    n0 = floor(n / nfold) * nfold
    ind_sample = matrix(sample(n0), nrow = nfold, ncol = n0 / nfold)
    for (k in 1:length(lambda_list)){
        lambda = lambda_list[k]
        for (j in 1:nfold){
            ind_test = ind_sample[j, ]
            ind_train = setdiff(c(1:n0), ind_test)
            X_test = X[ind_test, ]
            Y_test = Y[ind_test, ]
            X_train = X[ind_train, ]
            Y_train = Y[ind_train, ]
            model_train = picasso(X_train, Y_train, lambda=lambda, 
                                  method="l1", intercept=FALSE)
            beta_train = model_train$beta
            residual = as.matrix(Y_test - X_test %*% beta_train)
            cv_score[k] = cv_score[k] + sum((residual)^2)/nfold 
        }
    }
    ind_cv = which(cv_score == min(cv_score))
    TPR_cv = TPR[ind_cv]
    FDR_cv = FDR[ind_cv]
    lambda_cv = lambda_list[ind_cv]
    
    beta_hat_cv = picasso(X, Y, lambda = lambda_cv, method = "l1")$beta 
  
    lasso_result = list(TPR = TPR, FDR = FDR, TPR_cv = TPR_cv, FDR_cv = FDR_cv, lambda_cv = lambda_cv, beta = beta_hat_cv  ,Ham_loss = Ham)
    return(lasso_result)
}

lasso_cv_select <- function(X, Y, supp_true, lambda_list, nfold = 10){
  ## variable selection using LASSO
  ## Input: X, Y: data
  ##        supp_true: the true support (not used during training)
  ##        lambda_list: list of lambda values for training
  ##        nfold: number of folds for cross validation
  ## Output: scad_result: list of TPR, FDR and their values at the best lambda value
  
  s = length(supp_true)
  n = dim(X)[1]
  p = dim(X)[2]
  
  ## TPR and FDR
  lasso_cv <- cv.glmnet(X, Y, alpha = 1, lambda = lambda_list,
                        standardize = FALSE, intercept = F, nfolds = nfold)
  #lasso_cv_picasso = picasso(X, Y, lambda = lasso_cv$lambda.1se, method = "l1", intercept = F, standardize = T)
  beta_lasso_cv = as.vector(glmnet(X,Y, family = "gaussian", standardize = F, intercept = F, lambda = lasso_cv$lambda.1se)$beta)  #as.vector(coef(lasso_cv, s = "lambda.min"))[-1]
  #beta_lasso_cv  = lasso_cv_picasso$beta
  supp_beta_lasso_cv = which(beta_lasso_cv != 0)
  FDR_lasso_cv_screen = length(setdiff(supp_beta_lasso_cv, supp_true))/max(1,length(supp_beta_lasso_cv)); FDR_lasso_cv_screen
  TPR_lasso_cv_screen = length(intersect(supp_beta_lasso_cv, supp_true))/length(supp_true); TPR_lasso_cv_screen
  
  lasso_result = list( TPR_cv = TPR_lasso_cv_screen, FDR_cv = FDR_lasso_cv_screen, lambda_cv = lasso_cv$lambda.1se, beta = beta_lasso_cv  )
  return(lasso_result)
}



iht2_select <- function(X, Y, supp_true, s_list, g_list, beta0, maxiter = 1e3, prec = 1e-7){
    ## variable selection using iterative hard thresholding
    ## Input: X, Y: data
    ##        supp_true: the true support (not used during training)
    ##        s_list: list of s values for training
    ##        g_list (same length as s_list): list of g values for training
    ##        beta0: initialization
    ##        maxiter: maximum number of iterations
    ##        prec: precision to determine when to stop the iterations
    ## Output: scad_result: list of TPR, FDR and number of iterations
  
    s = length(supp_true)
    p = dim(X)[2]
    model_iht2 = iht_cor2(X, Y, s_list=s_list, g_list=g_list, beta0=beta0, maxiter = maxiter, prec = prec)
    beta_iht2 = model_iht2$beta
    grad_iht2 = model_iht2$grad
    steps_iht2 = model_iht2$steps_list
    
    beta_order = matrix(0, nrow=p, ncol=length(s_list))
    grad_order = matrix(0, nrow=p, ncol=length(s_list))
    for (j in 1:length(s_list)){
        beta_order[, j] = rev(order(abs(beta_iht2[, j])))
        grad_iht2[beta_order[1:s_list[j], j], j] = 0
        grad_order[, j] = rev(order(abs(grad_iht2[, j])))
        beta_order[(s_list[j] + 1):p, j] = grad_order[1:(p - s_list[j]), j]
    }

    TP = matrix(0, nrow=p, ncol=length(s_list))
    FP = matrix(0, nrow=p, ncol=length(s_list))
    FN = matrix(0, nrow=p, ncol=length(s_list))
    Ham = matrix(0, nrow=p, ncol=length(s_list))
    TPR = matrix(0, nrow=p, ncol=length(s_list))
    FDR = matrix(0, nrow=p, ncol=length(s_list))
    for (j in 1:(length(s_list))){
        TP_ind = as.numeric(is.element(beta_order[, j], supp_true))
        TP[, j] = cumsum(TP_ind)
        FN[,j] = s - TP[,j]
        TPR[, j] = TP[, j] / s
        FP[, j] = (1:p) - TP[, j]
        FDR[, j] = FP[, j] / (1:p)
        Ham[, j] = FP[, j] + FN[, j] 
    }
    iht2_result = list(TPR=TPR, FDR=FDR, steps_iht2=steps_iht2, Ham_loss = Ham)
    return(iht2_result)
}


ms_ordered_select <- function(X, Y, supp_true, scale_X = TRUE, path = F, s_level){
  ## variable selection using SIS (marginal screening)
  ## Input: X, Y: data
  ##        supp_true: the true support (not used during training)
  ##        scale_X: whether to scale X before training
  ## Output: ms_result: list of supp_estimate, TPR and FDR 
  
  s = length(supp_true)
  n = dim(X)[1]
  p = dim(X)[2]
  if (scale_X) 
    X = scale(X)
  #Sxx = t(X) %*% X / n
  Sxy = t(X) %*% Y / n
  grad = - Sxy
  grad_sort_idx = order(abs(grad), decreasing = TRUE)
  supp_estimate = grad_sort_idx[1:s_level]
  
  FDR = length(setdiff(supp_estimate, supp_true))/length(supp_estimate)
  
  TPR = length(intersect(supp_estimate, supp_true))/length(supp_true)
  
  
  FDR_path = c(); TPR_path = c()
  if(path == T){
  for (j in 1:p){
    set = grad_sort_idx[1:j]
    
    FDR_path[j] = length(setdiff(set, supp_true))/max(1,length(set))
    TPR_path[j] = length(intersect(set, supp_true))/length(supp_true)
  }
  
  }
  ms_result = list(supp_estimate = supp_estimate, TPR=TPR, FDR=FDR, grad_sort_idx = grad_sort_idx, FDR_path = c(0,FDR_path), TPR_path = c(0,TPR_path))
  return(ms_result)
}

ms_ordered_select_cv <- function(X, Y, supp_true, scale_X = TRUE, nfold = 10, s_level_list){
  ## variable selection using SIS (marginal screening)
  ## Input: X, Y: data
  ##        supp_true: the true support (not used during training)
  ##        scale_X: whether to scale X before training
  ## Output: ms_result: list of supp_estimate, TPR and FDR 
  
  #cv_score = matrix(0, nrow = length(s_list), ncol = length(g_list))
  n0 = floor(n / nfold) * nfold
  ind_sample = matrix(sample(n0), nrow = nfold, ncol = n0 / nfold)
  
  s = length(supp_true)
  n = dim(X)[1]
  p = dim(X)[2]
  if (scale_X) 
    X = scale(X)
  #Sxx = t(X) %*% X / n
  #Sxy = t(X) %*% Y / n
  #grad = - Sxy
  #grad_sort_idx = order(abs(grad), decreasing = TRUE)
  
  cv_score_mat =  foreach(k = 1:length(s_level_list), .combine = 'rbind')%do%{
    #cat(k, "th s_level","\n")
    s_level = s_level_list[k]
    #M = grad_sort_idx[1:s_level]
    #M = sort(M)
    
    val =  foreach (j = 1:nfold, .combine = '+')%do%{
      #cat(j, "th fold", k, "th s_level","\n")
      ind_test = ind_sample[j, ]
      ind_train = setdiff(c(1:n0), ind_test)
      X1 = X[ind_train, ]
      Y1 = Y[ind_train, ]
      #Sxx = t(X1)%*%X1/n
      Sxy = t(X1)%*%Y1/n
      grad = -Sxy
      grad_sort_idx = order(abs(grad), decreasing = TRUE)
      M = grad_sort_idx[1:s_level]
      M = sort(M) 
      X_train  = X1[, M]
      Y_train = Y1
      X_test = X[ind_test , M]
      Y_test = Y[ind_test, ]
      X_train = X[ind_train, M]
      Y_train = Y[ind_train, ]
      beta_OLS = solve(t(X_train)%*%X_train, t(X_train)%*%Y_train)
      cv_val = sum((Y_test - X_test %*% beta_OLS)^2)/nfold
      vec = c(cv_val, (nfold^2) * cv_val^2/(nfold -1))
      #print(vec)
    }
    val  = val
    # cat(k,"th inner foreach complete")
  }
  #print("foreach complete")
  id_min = which.min(cv_score_mat[,1]); s_level.min = s_level_list[id_min]; cv.min = cv_score_mat[id_min, 1]
  support.min = sort(grad_sort_idx[1:s_level.min]) 
  TPR.min = length(intersect(support.min, supp_true))/length(supp_true)
  FDR.min = length(setdiff(support.min, supp_true))/ length(support.min)
  
  se.min = sqrt(cv_score_mat[id_min, 2] - cv_score_mat[id_min, 1]^2)
  L = cv.min - se.min; U = cv.min + se.min
  cv_score_vec = cv_score_mat[,1]
  id.1se = min(which(cv_score_vec <= U & cv_score_vec >= L))
  s_level.1se = s_level_list[id.1se]
  #print(s_level.1se)
  support.1se = sort(grad_sort_idx[1:s_level.1se])
  #print(support.1se)
  TPR.1se = length(intersect(support.1se, supp_true))/length(supp_true)
  FDR.1se = length(setdiff(support.1se, supp_true))/ length(support.1se)
  
  ms_cv_list = list(s_level.min = s_level.min, support.min = support.min, TPR.min = TPR.min, FDR.min = FDR.min,
                    s_level.1se = s_level.1se, support.1se = support.1se, TPR.1se = TPR.1se, FDR.1se = FDR.1se)
  return(ms_cv_list)
}


# ms_select <- function(X, Y, supp_true, scale_X = TRUE){
#     ## variable selection using SIS (marginal screening)
#     ## Input: X, Y: data
#     ##        supp_true: the true support (not used during training)
#     ##        scale_X: whether to scale X before training
#     ## Output: ms_result: list of TPR and FDR 
#   
#     s = length(supp_true)
#     n = dim(X)[1]
#     p = dim(X)[2]
#     if (scale_X) 
#         X = scale(X)
#     Sxx = t(X) %*% X / n
#     Sxy = t(X) %*% Y / n
#     grad = - Sxy
#     grad_sort_idx = order(abs(grad), decreasing = TRUE)
#     TP_ms = rep(0, p)
#     TP_ms[1] = as.numeric(grad_sort_idx[1] %in% supp_true)
#     for (i in c(2:p)){
#         TP_ms[i] = TP_ms[i-1] + as.numeric(grad_sort_idx[i] %in% supp_true)
#     }
#     TPR_ms = TP_ms / s
#     FP_ms = c(1:p) - TP_ms
#     FDR_ms = rep(0, p)
#     for (i in c(1:p)){
#         FDR_ms[i] = FP_ms[i] / i
#     }
#     ms_result = list(TPR=TPR_ms, FDR=FDR_ms)
#     return(ms_result)
# }

scad_cv <- function(X, Y, lambda_list, nfold = 10){
    ## cross validation for SCAD (picasso) using a fixed list of lambdas
    ## Input: X, Y: data
    ##        lambda_list: list of lambda values
    ##        nfold: number of folds for cross validation
    ## Output: list SCADcv, containing: lambda.min, coef_min
    n = dim(X)[1]
    p = dim(X)[2]
  
    cv_score = rep(0, length(lambda_list))
    n0 = floor(n / nfold) * nfold
    ind_sample = matrix(sample(n0), nrow = nfold, ncol = n0 / nfold)
    for (k in 1:length(lambda_list)){
        lambda = lambda_list[k]
        for (j in 1:nfold){
            ind_test = ind_sample[j, ]
            ind_train = setdiff(c(1:n0), ind_test)
            X_test = X[ind_test, ]
            Y_test = Y[ind_test, ]
            X_train = X[ind_train, ]
            Y_train = Y[ind_train, ]
            model_train = picasso(X_train, Y_train, lambda=lambda, 
                                  method="scad", intercept=FALSE)
            beta_train = model_train$beta
            cv_score[k] = cv_score[k] + sum((Y_test - X_test %*% beta_train)^2)/nfold 
        }

        
        
            }
    ind_cv = which(cv_score == min(cv_score))
    lambda.min = lambda_list[ind_cv]
    model_min = picasso(X, Y, lambda=lambda.min, standardize = FALSE,
                        method="scad", intercept=FALSE)
    coef_min = model_min$beta
    SCADcv = list(lambda.min = lambda.min, coef_min = coef_min)
    return(SCADcv)
}

iht_one_cv <- function(X, Y, s_list, eta_list, beta0, nfold = 10, n_cv = 1, maxiter = 1e3, prec = 1e-7){
  ## cross validation for IHT
  ## Input: X, Y: data
  ##        s_list, g_list: list of s and g values
  ##        beta0: initialization
  ##        nfold: number of folds for cross validation
  ##        n_cv: number of cross validations performed
  ## Output: list IHTcv, containing: s.min, g.min, coef_min, cv_score
  n = dim(X)[1]
  p = dim(X)[2]
  cv_score = rep(0, length(s_list))
  n0 = floor(n / nfold) * nfold
  for (ind_cv in 1:n_cv){
    ind_sample = matrix(sample(n0), nrow = nfold, ncol = n0 / nfold)
    for (j in 1:nfold){
      ind_test = ind_sample[j, ]
      ind_train = setdiff(c(1:n0), ind_test)
      X_test = X[ind_test, ]
      Y_test = Y[ind_test, ]
      X_train = X[ind_train, ]
      Y_train = Y[ind_train, ]
      model_iht = iht_one(X_train, Y_train, s_list=s_list, eta_list=eta_list, beta0=beta0, maxiter = maxiter, prec = prec)
      beta_iht = model_iht$beta
      for (k in 1:length(s_list)){
        cv_score[k] = cv_score[k] + sum((Y_test - X_test %*% beta_iht[, k])^2)/nfold/n_cv
      }
    }
  }
  ind_cv = which(cv_score == min(cv_score))
  s_iht.min = s_list[ind_cv]
  #g_iht.min = g_list[ind_cv]
  model_iht = iht_one(X, Y, s_list=c(s_iht.min), eta_list = eta_list, beta0=beta0, maxiter = maxiter, prec = prec)
  coef_min = model_iht$beta
  IHTcv = list(s.min = s_iht.min, eta_list = eta_list, coef_min = coef_min, cv_score = cv_score)
  return(IHTcv)
}




iht_cv <- function(X, Y, s_list, g_list, beta0, nfold = 10, n_cv = 100, maxiter = 1e3, prec = 1e-7){
    ## cross validation for IHT
    ## Input: X, Y: data
    ##        s_list, g_list: list of s and g values
    ##        beta0: initialization
    ##        nfold: number of folds for cross validation
    ##        n_cv: number of cross validations performed
    ## Output: list IHTcv, containing: s.min, g.min, coef_min, cv_score
    n = dim(X)[1]
    p = dim(X)[2]
    cv_score = rep(0, length(s_list))
    n0 = floor(n / nfold) * nfold
    for (ind_cv in 1:n_cv){
    ind_sample = matrix(sample(n0), nrow = nfold, ncol = n0 / nfold)
    for (j in 1:nfold){
        ind_test = ind_sample[j, ]
        ind_train = setdiff(c(1:n0), ind_test)
        X_test = X[ind_test, ]
        Y_test = Y[ind_test, ]
        X_train = X[ind_train, ]
        Y_train = Y[ind_train, ]
        model_iht = iht_cor2(X_train, Y_train, s_list=s_list, g_list=g_list, beta0=beta0, maxiter = maxiter, prec = prec, scale_X = FALSE)
        beta_iht = model_iht$beta
        for (k in 1:length(s_list)){
            cv_score[k] = cv_score[k] + sum((Y_test - X_test %*% beta_iht[, k])^2)/nfold/n_cv
        }
    }
    }
    ind_cv = which(cv_score == min(cv_score))
    s_iht.min = s_list[ind_cv]
    g_iht.min = g_list[ind_cv]
    model_iht = iht_cor2(X, Y, s_list=c(s_iht.min), g_list=c(g_iht.min), beta0=beta0, maxiter = maxiter, prec = prec)
    coef_min = model_iht$beta
    IHTcv = list(s.min = s_iht.min, g.min=g_iht.min, coef_min = coef_min, cv_score = cv_score)
    return(IHTcv)
}

## IHT cv grid

iht_cv_grid <- function(X, Y, s_list, g_list, beta0, nfold = 10, n_cv = 100, maxiter = 1e3, prec = 1e-7){
  ## cross validation for IHT
  ## Input: X, Y: data
  ##        s_list, g_list: list of s and g values
  ##        beta0: initialization
  ##        nfold: number of folds for cross validation
  ##        n_cv: number of cross validations performed
  ## Output: list IHTcv, containing: s.min, g.min, coef_min, cv_score
  #require(doParallel)
  n = dim(X)[1]
  p = dim(X)[2]
  # cv_score = matrix(0, nrow = length(s_list), ncol = length(g_list))
  n0 = floor(n / nfold) * nfold
  
  
  ind_sample = matrix(sample(n0), nrow = nfold, ncol = n0 / nfold)
  
  cv_score = foreach(k1 = 1:length(s_list), .combine = 'rbind')%do%{
  # cat(k1, "th s_hat",'\n')
    val = rep(0, length(g_list))
    for (j in 1:nfold){
      
      ind_test = ind_sample[j, ]
      ind_train = setdiff(c(1:n0), ind_test)
      X_test = X[ind_test, ]
      Y_test = Y[ind_test, ]
      X_train = X[ind_train, ]
      Y_train = Y[ind_train, ]
      model_iht = iht_cor2(X_train, Y_train, s_list=rep(s_list[k1], length(s_list)), g_list=g_list, beta0=beta0, maxiter = maxiter, prec = prec, scale_X = FALSE)
      beta_iht = model_iht$beta
      
      for (k in 1:length(g_list)){
        val[k] = val[k] + sum((Y_test - X_test %*% beta_iht[, k])^2)/nfold/n_cv
      }
      
    }
    val = val
  }
  
  ind_cv = which(cv_score == min(cv_score), arr.ind = TRUE)
  s_iht.min = s_list[ind_cv[1]]
  g_iht.min = g_list[ind_cv[2]]
  model_iht = iht_cor2(X, Y, s_list=c(s_iht.min), g_list=c(g_iht.min), beta0=beta0, maxiter = maxiter, prec = prec)
  coef_min = model_iht$beta
  IHTcv = list(s.min = s_iht.min, g.min=g_iht.min, coef_min = coef_min, cv_score = cv_score)
  return(IHTcv)
}




ms_cv <- function(X, Y, s_list, nfold = 10, n_cv = 100){
    ## cross validation for SIS (marginal screening)
    ## Input: X, Y: data
    ##        s_list: list of s values
    ##        nfold: number of folds for cross validation
    ##        n_cv: number of cross validations performed
    ## Output: list MScv, containing: s.min, coef_min, cv_score
  
    n = dim(X)[1]
    p = dim(X)[2]
    cv_score = rep(0, length(s_list))
    n0 = floor(n / nfold) * nfold
    for (ind_cv in 1:n_cv){
    ind_sample = matrix(sample(n0), nrow = nfold, ncol = n0 / nfold)
    for (j in 1:nfold){
        ind_test = ind_sample[j, ]
        ind_train = setdiff(c(1:n0), ind_test)
        X_test = X[ind_test, ]
        Y_test = Y[ind_test, ]
        X_train = X[ind_train, ]
        Y_train = Y[ind_train, ]
        result_ms = marginal_screening(X_train, Y_train, s_list)
        beta_ms = result_ms$beta
        for (k in 1:length(s_list)){
          cv_score[k] = cv_score[k] + sum((Y_test - X_test %*% beta_ms[, k])^2)/nfold/n_cv
        }
    }
    }
    ind_cv = which(cv_score == min(cv_score))
    s_ms.min = s_list[ind_cv]
    model_ms = marginal_screening(X, Y, s_list = c(s_ms.min))
    coef_min = model_ms$beta
    grad_min = model_ms$grad
    MS_cv = list(s.min = s_ms.min, coef_min = coef_min, grad_min = grad_min, cv_score = cv_score)
    return(MS_cv)
}


