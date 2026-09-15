## Online Supplement to manuscript Counting on count regression: a reexamination of routinely-cited negative binomial specifications

## 07 / 2026 
## 08 / 2026 (add zero truncation)
## 09 / 2026 add countreg; IRRs; and derivatives wrt log alpha or log theta

## -- [on first run] countreg installation (not CRAN)        ####
#install.packages("countreg", repos = "https://zeileis.R-universe.dev")
## Input prep
## Input prep
## -- load functions for auxiliary Poisson regression        ####
source("aux_functions/aux_Pois_fun_ZT.R")   


## -- read data                                              ####
mypath <- "./data"
myData <- read.csv(paste0(mypath, "/TAS_NB_example_ZT.csv", collapse = " "), header = TRUE, row.names=1)



##                                                           ####
## Functions declaration                                     ####
##                                                           ####
## A) NB MLE eqs (allows zero truncation and offset)         ####
## -- own NB PMF                                             ####

## Equations in the paper: 1 (main text) & S1, S18 (Supplementary Material)
## Note:
## -- The pre-built function ?pnbinom uses a different parametrisation

dnegbin_own <-   function(x, p_lambda, p_alpha, zerotrunc = FALSE){
  if(p_alpha == 0){
    theta_nb <- gamma_alpha <- 0
  } else {
    theta_nb <- 1/p_alpha
    gamma_alpha <-gamma(p_alpha)
  }
  sapply(1:length(p_lambda), function(j){
    mu <- p_lambda[j]                                                         # equivalent to p_alpha*a/b
    a <- mu/(mu + p_alpha)
    b <- 1-a                                                                  # equivalent to p_alpha/(p_alpha + mu)
    #bin_coeff <- (gamma(x + p_alpha)/(gamma(x + 1)*gamma_alpha))              
    bin_coeff <- choose((x + p_alpha - 1), x)                                 # should be the same
    pr_x <-  bin_coeff*(a^x)*(b^p_alpha)                                      # should be the same as dnbinom(x, size = p_alpha, p = b)
    if(zerotrunc){
      CDF_0 <- (1 / (1 + theta_nb*mu)^p_alpha)
      pr_x <- pr_x / (1 - CDF_0)                                        
    } 
    return(pr_x)
  })
}


## -- own CDF                                                ####

negbin_CDF <- function(x, p_lambda, p_alpha, zerotrunc = FALSE){
  sum(sapply(0:x, function(k){
    dnegbin_own(k, p_lambda = p_lambda, p_alpha = p_alpha, zerotrunc = zerotrunc)          # pnbinom(x, size = p_alpha, p = b)
  }))
}   



## -- Mean and variance                                      ####

## Equations in the paper: S2 & S19;  S4 & S21 (Supplementary Material)
## Assumes:
## -- p_lambda <- linkingFun(beta_iter = reg_coeff_iter, X = X)  without offset ie., just exp(X %*% beta_iter)

E_NB <- function(p_lambda, p_alpha, rate_param = FALSE, offset_feat = NA, offset_feat_value = NA, zerotrunc = FALSE){
  if(rate_param & (is.na(offset_feat) | length(which(is.na(offset_feat_value)))>0) ){
    stop("Invalid arguments for rate parametrisation")
  }
  if(p_alpha == 0){
    p_theta <- 0
  } else {
    p_theta <- 1/p_alpha
  }
  mu <- p_lambda
  if(rate_param & !is.na(offset_feat)){
    lambda_tilde <- mu * offset_feat_value                                   # mean w/offset =  exp(x'B)*t = exp(x'B + ln t)
    mu <- lambda_tilde 
  }
  if(zerotrunc){
    CDF_0 <- (1 / (1 + p_theta*mu)^p_alpha)
    mu <- mu / (1 - CDF_0)
    mu[which(!is.finite(mu))] <- 0                                           # just in case
  }
  return(mu)
}


Var_NB <- function(p_lambda, p_alpha, rate_param = FALSE, offset_feat = NA, offset_feat_value = NA, zerotrunc = FALSE){
  if(rate_param & (is.na(offset_feat) | length(which(is.na(offset_feat_value)))>0) ){
    stop("Invalid arguments for rate parametrisation")
  }
  if(p_alpha == 0){
    theta_nb <- 0
  } else {
    theta_nb <- 1/p_alpha
  }
  dvar <- sapply(1:length(p_lambda), function(j){
    mu <- E_NB(p_lambda[j], 
               p_alpha, 
               rate_param = rate_param, offset_feat = offset_feat, offset_feat_value = offset_feat_value, zerotrunc = zerotrunc)
    a <- 1 + theta_nb*mu
    c <- mu * a
    if (zerotrunc){
      CDF_0 <- 1/a^p_alpha
      d <- CDF_0 * mu^2
      c - d
    } else {
      c
    }
  })
  return(dvar)
}


## -- Log likelihood                                         ####

## Equations in the paper: 15 (main text) & S23 (Supplementary Material)

loglik_NB <- function(X_data, y_data, p_beta, p_alpha, rate_param = FALSE, offset_feat = NA, offset_feat_value = NA, zerotrunc = FALSE){
  if(p_alpha == 0){
    p_theta <- alpha_lgam <- alpha_log <- 0
  } else {
    p_theta <- 1/p_alpha
    alpha_lgam <- lgamma(p_alpha)
    alpha_log <- log(p_alpha)
  }
  p_lambda <- linkingFun(p_beta, X_data)
  if(rate_param & !is.na(offset_feat)){ 
    p_lambda <- p_lambda*offset_feat_value
  }
  loglambda <- log(p_lambda)
  loglambda[which(!is.finite(loglambda))] <- 0
  a <- lgamma(y_data + p_alpha) - alpha_lgam - lgamma(y_data + 1)
  b <- y_data*loglambda + p_alpha*alpha_log
  c <- (p_alpha + y_data)*log(p_lambda + p_alpha)
  nb_loglik <- sum(a + b - c)
  if(zerotrunc){
    CDF_0 <- (1 / (1 + p_theta*p_lambda)^p_alpha)
    nb_loglik <- nb_loglik - sum(log(1 - CDF_0))  
  }
  return(nb_loglik)
}




## -- Gradient wrt beta                                      ####

## Equations in the paper: S6, S36 & S37 (Supplementary material)
## Assumes:
## -- p_lambda = exp(X %*% beta_iter)

gradient_NB <- function(p_lambda, p_alpha, X, y, rate_param = FALSE, offset_feat = NA, offset_feat_value = NA, zerotrunc = FALSE){
  if(p_alpha == 0){
    p_theta <- 0
  } else {
    p_theta <- 1/p_alpha
  }
  p_mu <- E_NB(p_lambda, 
               p_alpha, 
               rate_param = rate_param, offset_feat = offset_feat, offset_feat_value = offset_feat_value, zerotrunc = zerotrunc) 
  if(rate_param & !is.na(offset_feat)){
    p_lambda <- p_lambda*offset_feat_value
  }
  n_observ <- nrow(X)
  residue_iter <- y - p_lambda      
  stacked_g <- sapply(1:n_observ, function(i){                                       
    #i <- as.numeric(i)
    temp_g_i <- (residue_iter[i] / (1 + p_theta*p_lambda[i]))
    if(zerotrunc){
      zt_adjust_num <-  p_mu[i] 
      zt_adjust_den <- (1 + p_theta*p_lambda[i])^(p_alpha + 1) 
      temp_g_i <- temp_g_i - (zt_adjust_num / zt_adjust_den)
    } 
    temp_g_i*X[i,]
  })
  if(!is.matrix(stacked_g)){ g_gradient <- sum(stacked_g)
  } else { g_gradient <- as.matrix(apply(stacked_g,1,sum)) }
  return(list(gradient = g_gradient, stacked_g = stacked_g,
              pois_residue =  residue_iter)) 
}  





## -- Hessian wrt beta                                       ####

## Equations in the paper: S7; S38 & S39 (Supplementary Material)
## Assumes:
## -- p_lambda = exp(X %*% beta_iter)

hessian_NB <- function(p_lambda, p_alpha, X, y, rate_param = FALSE, offset_feat = NA, offset_feat_value = NA, zerotrunc = FALSE){
  
  ## For debug
  # p_lambda = lambda_iter
  # p_alpha = p_alpha_fix
  #X = X_data
  #y = y_data
  
  
  if(rate_param & !is.na(offset_feat)){
    p_lambda <- p_lambda*offset_feat_value
  }
  n_regressors <- ncol(X)
  n_observ <- nrow(X)
  
  if(p_alpha == 0){
    p_theta <- 0
  } else {
    p_theta <- 1/p_alpha
  }
  Hessian_iter <- matrix(0L,nrow = n_regressors, ncol = n_regressors)
  ## Aternative 1
  phi_sub <- (1 + p_theta*p_lambda)
  H_parent <- ((1 + p_theta*y) / phi_sub^2)*p_lambda
  if(zerotrunc){
    phi_sub_1a <- (1 / phi_sub^(1+p_alpha))
    phi <-  (1/phi_sub^p_alpha) 
    H_zt_A  <- ( p_lambda*phi_sub_1a / (1 - phi) )^2
    H_zt_B1 <- ( p_lambda*phi_sub_1a ) /  (1 - phi)
    #H_zt_B2 <- (p_lambda^2 * (p_theta + 1) * phi_sub^(p_alpha-2)) /  (1 - phi)            # wrong equation in Grogger and Carson!
    H_zt_B2 <- (p_lambda^2 * (p_theta + 1) * (1/phi_sub^(p_alpha+2)) ) /  (1 - phi)
  }
  Hessian_iter <- matrix(0L,nrow = n_regressors, ncol = n_regressors)
  for(i in 1:n_observ){
    x <- X[i,]
    x <- as.numeric(x)
    outprod_X <- (x %o% x) 
    A <- H_parent[i]
    if(zerotrunc){
      B <- H_zt_A[i]
      C <- H_zt_B1[i]
      D <- H_zt_B2[i]
    } else {
      B <- C <- D <- 0
    }
    Hessian_iter <- Hessian_iter - ( A - B + C - D)*outprod_X     
  }
  
  
  ## Alternative 2: based on my paper (uses the mean)
  #p_mu <- E_NB(p_lambda, 
  #             p_alpha, 
  #             rate_param = rate_param, offset_feat = offset_feat, offset_feat_value = offset_feat_value, zerotrunc = zerotrunc)
  # for(i in 1:n_observ){
  #   x <- X[i,]
  #   x <- as.numeric(x)
  #   outprod_X <- (x %o% x) 
  #   negbin_hess_num <- (1 + p_theta*y[i])*p_lambda[i]
  #   negbin_hess_den <- (1 + p_theta*p_lambda[i])^2
  #   H_NB <- (negbin_hess_num / negbin_hess_den) 
  #   if(zerotrunc){
  #     zt_adjust_num_1 <- p_mu[i]
  #     zt_adjust_num_2 <- zt_adjust_num_1 * (p_theta + 1) * p_lambda[i]
  #     zt_adjust_den_1 <- (1 + p_theta*p_lambda[i])^(p_alpha + 1)
  #     zt_adjust_den_2 <- (1 + p_theta*p_lambda[i])^(2 - p_alpha)
  #     zt_adjust  <- ( (zt_adjust_num_1/zt_adjust_den_1) * (1 - zt_adjust_num_1/zt_adjust_den_1) ) - (zt_adjust_num_2/zt_adjust_den_2)
  #     
  #     ## Alternative 3: equivalent, but Grogger's notation
  #     # temp_a <- 1 + p_theta*p_lambda[i]
  #     # temp_b <- 1 - temp_a^(-p_alpha)
  #     # ztrunc_adjust_H_1_num <- p_lambda[i] * temp_a^(-2*(1+p_alpha))
  #     # ztrunc_adjust_H_1_den <- temp_b^2
  #     # ztrunc_adjust_H_2_num <- temp_a^(-(1+p_alpha)) - (p_lambda[i]*(p_theta + 1)*temp_a^(-2+p_alpha))
  #     # ztrunc_adjust_H_2_den <- temp_b
  #     # zt_adjust2 <- ((ztrunc_adjust_H_2_num / ztrunc_adjust_H_2_den) - (ztrunc_adjust_H_1_num / ztrunc_adjust_H_1_den)) * p_lambda[i]
  # 
  #   } else {
  #     zt_adjust  <- 0
  #   }
  #   Hessian_iter <- Hessian_iter - ( (negbin_hess_num / negbin_hess_den) + zt_adjust )*outprod_X     
  # }
  
  
  return(Hessian_iter)
}






## -- Gradient wrt alpha                                     ####

## Equations in the paper: 15, 22 & 23(main text); S24, S25 (Supplementary Material)
## Assumes:
## -- p_lambda = exp(X %*% beta_iter)* offset if rate parametrization


gradient_NB_alpha <- function(y, p_alpha, p_lambda, rate_param = FALSE, offset_feat = NA, offset_feat_value = NA, zerotrunc = FALSE){
  if(!is.matrix(y)) y <- as.matrix(y)
  n_observ <- nrow(y)
  if(p_alpha == 0){
    p_theta <- 0
    alpha_log <-  0
    alpha_digam <- 0
  } else {
    p_theta <- 1/p_alpha
    alpha_log <-  log(p_alpha)
    alpha_digam <-  digamma(p_alpha)
  }
  if(rate_param & !is.na(offset_feat)){
    p_lambda <- p_lambda*offset_feat_value
  }
  ## Option 1
  #a <- digamma(p_alpha + y) - alpha_digam + (alpha_log + 1)
  #b <- log(p_lambda + p_alpha)
  #c <- (p_alpha + y)/(p_lambda + p_alpha)
  ## Option 2: highlights the "residual"
  a <- digamma(p_alpha + y) - alpha_digam
  b <- 1 + p_theta*p_lambda
  c <- (y - p_lambda) / (p_alpha*b)
  if(zerotrunc){
    phi <- 1/b^p_alpha
    dphi <- ( (p_lambda / (p_alpha*b) ) * phi ) - ( log(b) * phi)
    d <- (dphi/(1 - phi))    
  } else {
    d <- 0
  }
  d_dalpha <- sum(a - log(b) - c + d)
  return(d_dalpha)
}




## -- Hessian wrt alpha                                      ####
## Equations in the paper: 16, 25 & 26 (main text); S26, S30 & S31 (Supplementary Material)
## Assumes:
## -- p_lambda = exp(X %*% beta_iter) 

Hessian_NB_alpha <- function(y, p_alpha, p_lambda, rate_param = FALSE, offset_feat = NA, offset_feat_value = NA, zerotrunc = FALSE){
  if(p_alpha == 0){
    p_theta <- 0
    alpha_trigam <- 0
  } else {
    p_theta <- 1/p_alpha
    alpha_trigam <- trigamma(p_alpha)
  }
  if(rate_param & !is.na(offset_feat)){
    p_lambda <- p_lambda*offset_feat_value
  }
  a <- trigamma(y + p_alpha) - alpha_trigam
  b <- (p_lambda^2 + p_alpha*y)
  c <- (p_alpha*(p_alpha + p_lambda)^2)
  if(zerotrunc){ 
    d <- 1 + p_theta*p_lambda
    phi <- 1/d^p_alpha
    dphi <- ( log(d) * phi) - ( (p_lambda / (p_alpha*d) ) * phi )  
    d2phi <-  dphi^2/ phi + ( (p_theta*phi) * (p_lambda / (p_alpha*d))^2 ) 
    e <- (dphi/(1 - phi))^2 +  (d2phi/(1 - phi))
    
    ## equivalent shorter version in eq. S31 Supplementary materials
    f <- p_lambda / (p_alpha + p_lambda)
    g_fun <-  f - log(1 + p_theta*p_lambda) 
    dg_fun <- p_theta * f^2 
    e_short <- (phi / (1 - phi) ) * ( (g_fun^2 / (1 - phi)) + dg_fun)
    
    ## for debug
    #all.equal(e, e_short)
  } else {
    e <- 0
  }
  dL2_dalpha2 <- sum( a + (b/c) + e)
  return( dL2_dalpha2)
}


## -- Gradient wrt theta                                     ####

## Equations in the paper: 6, 8, 24 main text; S32 & S33 Supplementary Materials 
## Assumes:
## -- p_lambda = exp(X %*% beta_iter) * offset

gradient_NB_theta <- function(y, p_alpha, p_lambda, rate_param = FALSE, offset_feat = NA, offset_feat_value = offset_feat_value, zerotrunc = FALSE){
  if(p_alpha == 0){
    p_theta <- 0
    alpha_digam <- 0
  } else {
    p_theta <- 1/p_alpha
    alpha_digam <-  digamma(p_alpha)
  }
  if(rate_param & !is.na(offset_feat)){
    p_lambda_offset <- p_lambda*offset_feat_value                             
    p_lambda <- p_lambda_offset
  }
  a <- (-1/p_theta^2)*(digamma(y + p_alpha) - alpha_digam)
  b <- 1 + p_theta*p_lambda
  c <- (p_alpha^2)*log(b)
  d <- (y - p_lambda) / (p_theta*b)
  if(zerotrunc){
    phi <- 1/b^p_alpha
    dphi <- ( (p_lambda / (p_alpha*b) ) * phi ) - ( log(b) * phi)
    e <- (-1/p_theta^2)*(dphi/(1 - phi))    
  } else {
    e <- 0
  }
  d_dtheta <- sum(a + c + d + e)
  return(d_dtheta)
}

## -- Hessian wrt theta                                      ####

## Equations in the paper: 9, 11, 27 main text; S34 & S35 Supplementary Materials 
## Assumes:
## -- p_lambda <- linkingFun(beta_iter = reg_coeff_iter, X = X)  without offset ie., just exp(X %*% beta_iter)

Hessian_NB_theta <- function(y, p_alpha, p_lambda, rate_param = FALSE, offset_feat = NA, offset_feat_value = NA, zerotrunc = FALSE){
  if(p_alpha == 0){
    p_theta <- 0
    alpha_digam <- 0
  } else {
    p_theta <- 1/p_alpha
    alpha_digam <- digamma(p_alpha)
    alpha_trigam <- trigamma(p_alpha)
  } 
  if(rate_param & !is.na(offset_feat)){
    p_lambda_offset <- p_lambda*offset_feat_value                             
    p_lambda <- p_lambda_offset
  }
  a1 <- (2/p_theta^3)*(digamma(y + p_alpha) - alpha_digam)
  a2 <- (1/p_theta^4)*(trigamma(y + p_alpha ) - alpha_trigam)
  a <-  a1 + a2
  b1 <- 1+p_theta*p_lambda
  b2 <- 2*log(b1)
  b3 <- (p_theta*p_lambda)/b1
  b4 <- p_theta*(1+2*p_theta*p_lambda)/b1^2
  b5 <- y - p_lambda
  b <- 1/p_theta^3*( (b5 * b4) - b3 + b2 )
  if(zerotrunc){
    d <- 1 + p_theta*p_lambda
    phi <- 1/d^p_alpha
    dphi <- ( (p_lambda / (p_alpha*d) ) * phi ) - ( log(d) * phi)
    d2phi <-  dphi^2/ phi + ( (p_theta*phi) * (p_lambda / (p_alpha*d))^2 ) 
    c <- (1/p_theta^4)*(dphi/(1 - phi))^2 +  (1/(1 - phi))*((1/p_theta^4)*d2phi + (2/p_theta^3)*dphi)  
  } else {
    c <- 0
  }
  d2_theta <- sum(a - b + c)
  return(d2_theta)
}


## -- (zero truncated only), x-derivatives                   ####
## Equations in the paper: S43 & S47 (Supplementary Material)
## NOTES: 
## -- the cross-derivatives [alpha, beta] and [beta, alpha] yield the same resul
## -- the equivalent result in theta is -1/theta^2 * the result in alpha
## Assumes:
## -- p_lambda <- linkingFun(beta_iter = reg_coeff_iter, X = X)  without offset ie., just exp(X %*% beta_iter)

xd_NB_alpha <- function(X, y, p_alpha, p_lambda, rate_param = FALSE, offset_feat = NA, offset_feat_value = NA){
  ## Only relevant for the zero truncated case. For the parent specification this is just zero
  
  ## for debug
  #X = X_data 
  #y = y_data
  #p_alpha =  chosen_alpha
  #p_lambda = lambda_iter 
  
  if(p_alpha == 0){
    p_theta <- 0
  } else {
    p_theta <- 1/p_alpha
  }
  if(rate_param & !is.na(offset_feat)){
    p_lambda_offset <- p_lambda*offset_feat_value                             
    p_lambda <- p_lambda_offset
  }
  n_observ<- nrow(X)
  stacked_g<- sapply(1:n_observ, function(i){                                       
    i <- as.numeric(i)
    d <- 1 + p_theta*p_lambda[i]
    phi <- 1/d^p_alpha
    dphi <- ( (p_lambda[i] / (p_alpha*d) ) * phi ) - ( log(d) * phi)
    f <- p_lambda[i] / (p_alpha + p_lambda[i])
    g_fun <-  f - log(1 + p_theta*p_lambda[i]) 
    dg_fun <- p_theta * f^2 
    a <- (p_lambda[i]*y[i])/(p_alpha*(p_alpha + p_lambda[i])^2)
    b <- dg_fun
    c <- (p_lambda[i]/(p_alpha + p_lambda[i])) * (phi/(1 - phi)) * ( (p_lambda[i]/(p_alpha*(p_alpha + p_lambda[i]))) + dphi / (1 - phi) + g_fun )
    temp_g_i  <- a - b - c
    p_alpha*temp_g_i*X[i,]
  })
  if(!is.matrix(stacked_g)){ xd_beta_alpha <- sum(stacked_g)
  } else {xd_beta_alpha <- as.matrix(apply(stacked_g,1,sum)) }
  return( xd_beta_alpha )
}






## -- Hessian wrt log alpha                                  ####

## packages like countreg and VGAM would set z = ln alpha and differentiate the loglikelihood wrt z 
## -- alpha = exp(z)
## -- d alpha d z = exp(z) = alpha 
## -- d^2 alpha d z^2 = exp(z) = alpha 
## -- Scoring: d LL d z = d LL d alpha * exp(z) = d LL d alpha * alpha 
## -- Hessian: d^2LL d z^2 = d^2 LL d alpha^2 * exp(z)^2 + d LL d alpha  * exp(z) = d^2 LL d alpha^2 * alpha^2 + d LL d alpha * alpha 

Hessian_NB_logalpha <- function(y, p_alpha, p_lambda, rate_param = FALSE, offset_feat = NA, offset_feat_value = NA, zerotrunc = FALSE){
  a <- Hessian_NB_alpha(y = y, 
                        p_alpha = p_alpha, 
                        p_lambda = p_lambda, 
                        rate_param = rate_param, offset_feat = offset_feat, offset_feat_value = offset_feat_value, zerotrunc = zerotrunc)
  
  b <- gradient_NB_alpha(y = y, 
                         p_alpha = p_alpha, 
                         p_lambda = p_lambda, 
                         rate_param = rate_param, offset_feat = offset_feat, offset_feat_value = offset_feat_value, zerotrunc = zerotrunc)
  c <- p_alpha^2*a + p_alpha*b
  return(c)
}

## -- Hessian wrt log theta                                  #####

## papers like DOI:10.1080/02664763.2025.2545890 and packages like gamplss would set z = ln theta and differentiate the loglikelihood wrt z 
## -- theta = exp(z)
## -- d theta d z = exp(z) = theta;
## -- d^2 theta d z^2 = exp(z) = theta;
## -- d LL d z = d LL d theta * exp(z) = d LL d theta * theta  
## -- d^2LL d z^2 = d^2 LL d theta^2 * exp(z)^2 + d LL d theta  * exp(z) = d^2 LL d theta^2 * theta^2 + d LL d theta * theta
Hessian_NB_logtheta <- function(y, p_alpha, p_lambda, rate_param = FALSE, offset_feat = NA, offset_feat_value = NA, zerotrunc = FALSE){
  
  ## for debug
  #y=y_data
  #p_alpha = chosen_alpha 
  #p_lambda = lambda_iter
  
  a <- Hessian_NB_theta(y = y, 
                        p_alpha = p_alpha, 
                        p_lambda = p_lambda, 
                        rate_param = rate_param, offset_feat = offset_feat, offset_feat_value = offset_feat_value, zerotrunc = zerotrunc)
  
  b <- gradient_NB_theta(y = y, 
                         p_alpha = p_alpha, 
                         p_lambda = p_lambda, 
                         rate_param = rate_param, offset_feat = offset_feat, offset_feat_value = offset_feat_value, zerotrunc = zerotrunc)
  c <- (1/(p_alpha^2))*a + (1/p_alpha)*b
  return(c)
  
}


## B) Fitting, iterative                                     ####
## -- Half-stepping                                          ####
myHalfStepping_NB <- function(X_data, y_data, beta_iter, Hessian_iter, g_gradient, p_alpha_fix, rate_param = FALSE, offset_feat = NA, offset_feat_value = NA, zerotrunc = FALSE){
  accept_candidate <- FALSE
  converge_flag <- FALSE
  max_iter_halfstep <- 100
  loglik_acceptd <- loglik_NB(X_data, 
                              y_data, 
                              p_beta = beta_iter, 
                              p_alpha = p_alpha_fix,
                              rate_param = rate_param, offset_feat = offset_feat, offset_feat_value = offset_feat_value, zerotrunc = zerotrunc
  )
  iter_halfstep <- 0
  stepsize <- 1
  target_iter <- qr.solve(Hessian_iter, g_gradient, tol = 1e-18)            
  while(!accept_candidate & (iter_halfstep < max_iter_halfstep)){
    temp_delta <- beta_iter - stepsize*target_iter
    loglik_iter <- loglik_NB(X_data, 
                             y_data, 
                             p_beta = temp_delta, 
                             p_alpha = p_alpha_fix,
                             rate_param = rate_param, offset_feat = offset_feat, offset_feat_value = offset_feat_value, zerotrunc = zerotrunc
    )
    criterion_half_stepping <-  loglik_iter > loglik_acceptd
    if(!criterion_half_stepping){ 
      stepsize <- stepsize/2 
      iter_halfstep  <- iter_halfstep + 1
    } else { 
      accept_candidate <- TRUE 
      target_iter <- stepsize*target_iter                                        
      if(iter_halfstep < max_iter_halfstep) {converge_flag <- TRUE}
    }
  } 
  list(target_iter = target_iter, temp_delta = temp_delta, stepsize = stepsize, converge_flag = converge_flag)
}








## -- Newton-Raphson, wrt Beta                               ####
NR_MLE_NB <- function(X_data, y_data, p_alpha_fix, p_beta_init, iter_max = 100, halfstep = TRUE, rate_param = FALSE, offset_feat = NA, offset_feat_value = NA, zerotrunc = FALSE, BHHH = FALSE){
  
  ## FOR DEBUG 
  #iter_max = 100
  #halfstep = FALSE
  #p_alpha_fix = alpha_iter
  #p_beta_init = temp_beta
  
  # if(!zerotrunc){
  #   beta_iter <- p_beta_init
  # } else {
  #   ## From my notes to early version: it seems that initializing w Poisson coeff messes up zero truncated Neg BIn
  #   beta_iter <- matrix(0, nrow = n_regressors, ncol = 1)
  #   ln_target <- log(abs(y_data))
  #   ln_target[which(!is.finite(ln_target))]<-0                                 # in case the target variable takes value 0
  #   #beta_iter[1] <- mean(ln_target)                                           # initialisation from -- thread: https://www2.stat.duke.edu/courses/Spring21/sta440.001/slides/glm-2.html#10
  #   beta_iter[1] <- mean(y_data)                                               # APPARENTLY THIS IS NEEDED FOR NEG BIN ZERO TRUNCATED
  # }
  
  n_regressors <- ncol(X_data)
  beta_iter <- p_beta_init
  iter_count <- 0
  GG_count <- 0
  #tol <- rep(1e-6, n_regressors)
  tol <- 1e-6
  iter_stop <- FALSE
  hessian_is_singular <- FALSE
  lambda_iter <- linkingFun(beta_iter, X_data)
  list_grad <- list()
  while(iter_count <= iter_max & !iter_stop){
    iter_count <- iter_count + 1
    grad <- gradient_NB(p_lambda = lambda_iter,
                        p_alpha = p_alpha_fix, 
                        X = X_data, 
                        y = y_data,
                        rate_param = rate_param, offset_feat = offset_feat, offset_feat_value = offset_feat_value, zerotrunc = zerotrunc
    )           
    g_gradient <- grad$gradient
    
    ## for debug: track gradient - is it going to zero?
    list_grad[[iter_count]] <- g_gradient
    
    Hessian_iter <- hessian_NB(p_lambda = lambda_iter,
                               p_alpha = p_alpha_fix, 
                               X = X_data, 
                               y = y_data,
                               rate_param = rate_param, offset_feat = offset_feat, offset_feat_value = offset_feat_value, zerotrunc = zerotrunc
    )        
    
    ## Warning: Zero truncation tends to produce an Hessian with all positive elements even if the equation is correct
    ## -- see my post https://stats.stackexchange.com/q/676952/513606
    ## As an alternative to the Hessian I consider: BHHH estimator (Outer Product of Gradients)
    ## -- background: see Greene p.561-562 eq.14-17 vs 14-18 and p.1140
    ## -- for a clearer version: Long (1997) page 56
    if(BHHH & zerotrunc & (sum(sign(Hessian_iter)) >= -1*n_regressors^2) ){
      GG_count <- GG_count + 1
      stacked_g <-  grad$stacked_g
      if(length(is.na(nrow(stacked_g)))>0){
        ## business as usual
        G <- t(stacked_g)
      } else {
        G <- as.matrix(stacked_g)                                                 ## we are dealing with odd case in which we are dealing with 1 feature
      }
      GG <- t(G) %*% G
      ## CHECK: equivalent to summing the outer products of gradients
      # outer_grad <- lapply(1:ncol(stacked_g), function(j){
      #   x <- as.vector(stacked_g[,j])
      #   outprod_X  <- outer(x,x)
      # })
      # GG_outgrad <- Reduce("+",outer_grad)                                    # https://stackoverflow.com/questions/11641701/sum-a-list-of-matrices
      # all.equal(unname(GG), GG_outgrad)
      Hessian_iter <- -1*GG                                                     # based on Long p.57
    }
    inv_test <- try(qr.solve(Hessian_iter, g_gradient), silent = TRUE)
    if(inherits(inv_test, "try-error")){
      hessian_is_singular <- TRUE
      iter_stop <- TRUE } 
    if(!hessian_is_singular){
      if(halfstep){
        HalfStep_proc <- myHalfStepping_NB(X_data,
                                           y_data,
                                           beta_iter, 
                                           Hessian_iter, 
                                           g_gradient,
                                           p_alpha_fix,
                                           rate_param = rate_param, offset_feat = offset_feat, offset_feat_value = offset_feat_value, zerotrunc = zerotrunc)
        target_iter <- HalfStep_proc$target_iter
        temp_delta <-  HalfStep_proc$temp_delta
      } else {
        target_iter <- qr.solve(Hessian_iter, g_gradient, tol = 1e-18)
        temp_delta <- beta_iter - target_iter
      }
    } else { warning("Hessian is singular")
      target_iter <- 0 
    }
    if(length(which(sqrt((target_iter)^2) < tol)) < n_regressors){
      beta_iter <- temp_delta 
      lambda_iter <- linkingFun(beta_iter, X_data)
    } else {
      iter_stop <- TRUE
    }
    mean_iter <- E_NB(p_lambda = lambda_iter, 
                      p_alpha = p_alpha_fix, 
                      rate_param = rate_param, 
                      offset_feat = offset_feat,
                      offset_feat_value = offset_feat_value,
                      zerotrunc = zerotrunc)
  }
  return(list(coef_est = beta_iter, lambda = lambda_iter, predict_y = mean_iter, 
              Hessian = Hessian_iter, iter = iter_count, track_gradient = g_gradient, BHHH_instead_times = GG_count  ))
}




## -- Newton-Raphson, wrt alpha                              ####
NR_MLE_NB_alpha <- function(X_data, y_data, p_alpha_init, p_beta_fix, iter_max = 100, rate_param = FALSE, offset_feat = NA, offset_feat_value = NA, zerotrunc = FALSE){
  
  p_lambda_0 <- linkingFun(p_beta_fix, X_data)                                  # lambda = exp(x' beta)
  if(rate_param & !is.na(offset_feat)){
    p_lambda_0 <- p_lambda_0*offset_feat_value
  } 
  p_alpha_0 <- p_alpha_init                                                     # initial guess on alpha 
  iter_count <- 0
  tol <- 1e-6
  iter_stop <- FALSE
  hessian_is_singular <- FALSE
  while(iter_count <= iter_max & !iter_stop){
    iter_count <- iter_count + 1
    if(iter_count == 1){
      p_alpha_iter <-  p_alpha_0
      p_lambda_iter <- p_lambda_0
    } 
    
    g_gradient <-  gradient_NB_alpha(y = y_data,
                                     p_alpha = p_alpha_iter, 
                                     p_lambda = p_lambda_iter,
                                     zerotrunc = zerotrunc
    )           
    Hessian_iter <- Hessian_NB_alpha(y = y_data,
                                     p_alpha = p_alpha_iter, 
                                     p_lambda = p_lambda_iter,
                                     zerotrunc = zerotrunc
    )                 
    inv_test <- !is.finite(g_gradient/Hessian_iter)                              # dealing with scalars now...
    if(inv_test){
      hessian_is_singular <- TRUE
      iter_stop <- TRUE } 
    if(!hessian_is_singular){
      target_iter <- g_gradient/Hessian_iter
    } else { warning("Hessian is singular")
      target_iter <- 0 
    }
    temp_delta <- p_alpha_iter - target_iter
    if(sqrt(target_iter^2) > tol){
      p_alpha_iter <- temp_delta 
      if(p_alpha_iter < 0) p_alpha_iter <- 0
    } else {iter_stop <- TRUE}
  }
  return(list(p_alpha = p_alpha_iter, predict_y = p_lambda_iter, 
              Hessian = Hessian_iter, iter = iter_count  ))
}



## C) Expected (Fisher) Information:                         ####
## -- beta, beta                                             ####

## Equation S8 (non truncated) & S41 (truncated) - Supplementary materials
## Assumes:
## -- p_lambda <- linkingFun(beta_iter = reg_coeff_iter, X = X)  without offset ie., just exp(X %*% beta_iter)

FI_Beta <- function(p_lambda, p_alpha, X, y, rate_param = FALSE, offset_feat = NA, offset_feat_value = NA, zerotrunc = FALSE){
  
  ## DEBUG
  #p_lambda = linkingFun(beta_iter, X_data)   ## For ZT Eq. S41
  #p_alpha = chosen_alpha
  #X = X_data 
  #y = y_data
  
  
  
  if(rate_param & !is.na(offset_feat)){
    p_lambda <- p_lambda*offset_feat_value
  }
  n_regressors <- ncol(X)
  n_observ <- nrow(X)
  
  if(p_alpha == 0){
    p_theta <- 0
  } else {
    p_theta <- 1/p_alpha
  }
  E_Hessian_iter <- matrix(0L,nrow = n_regressors, ncol = n_regressors)
  phi_sub <- (1 + p_theta*p_lambda)
  E_H_parent_noZT <- (p_lambda / phi_sub)                                             # see Eq. S8
  if(zerotrunc){
    phi_sub_1a <- (1 / phi_sub^(1+p_alpha))
    phi <-  (1/phi_sub^p_alpha) 
    H_zt_A  <- ( p_lambda*phi_sub_1a / (1 - phi) )^2
    H_zt_B1 <- ( p_lambda*phi_sub_1a ) /  (1 - phi)
    #H_zt_B2 <- (p_lambda^2 * (p_theta + 1) * phi_sub^(p_alpha-2)) /  (1 - phi)       # wrong equation in Grogger and Carson!
    H_zt_B2 <- (p_lambda^2 * (p_theta + 1) * (1/phi_sub^(p_alpha+2)) ) /  (1 - phi)
  }
  for(i in 1:n_observ){
    x <- X[i,]
    x <- as.numeric(x)
    outprod_X <- (x %o% x) 
    A <- E_H_parent_noZT[i] 
    if(zerotrunc){
      A <- A/(1 - phi[i])                                                               # due to ZT mean
      B <- H_zt_A[i]
      C <- H_zt_B1[i]
      D <- H_zt_B2[i]
    } else {
      B <- C <- D <- 0
    }
    E_Hessian_iter <- E_Hessian_iter - ( A - B + C - D)*outprod_X     
  }
  minus_E_Hessian_iter <- -1*E_Hessian_iter
  return(minus_E_Hessian_iter)
  
  
  
  
  
  ## OLD VERSION GENERATING MIXED SINGS
  
  # if(p_alpha == 0){
  #   p_theta <- 0
  # } else {
  #   p_theta <- 1/p_alpha
  # }
  # if(rate_param & !is.na(offset_feat)){
  #   p_lambda_offset <- p_lambda*offset_feat_value                             
  #   p_lambda <- p_lambda_offset
  # }
  # n_regressors <- ncol(X)
  # n_observ <- nrow(X)
  # FI_Beta_iter <- matrix(0L,nrow = n_regressors, ncol = n_regressors)
  # for(i in 1:n_observ){
  #   x <- X[i,]
  #   x <- as.numeric(x)
  #   outprod_X <- (x %o% x)   
  #   if(zerotrunc){
  #     a <- (1 + p_theta*p_lambda[i])
  #     phi <- 1/a^p_alpha
  #     # b <- (1 + p_theta*E_NB(p_lambda = p_lambda[i], 
  #     #                        p_alpha = p_alpha, 
  #     #                        rate_param = rate_param,
  #     #                        offset_feat = offset_feat,
  #     #                        offset_feat_value = offset_feat_value,
  #     #                        zerotrunc = zerotrunc)) / a^2
  #     
  #     
  #     ## for debug 
  #     b <- (1 + p_theta*y[i])/a^2
  #           
  #     
  #     c <- ( p_lambda[i]*(1/a^(2*(1+p_alpha))) )/ (1-phi)^2
  #     d <- ( 1/(a^(1+p_alpha)) - p_lambda[i]*(1+p_theta)*(1/a^(2+p_alpha)) ) / phi 
  #     
  #     # for debug
  #     FI_Beta_arg <- p_lambda[i] * (-1) * (b - c  + d)
  # 
  #     #FI_Beta_arg <- p_lambda[i] * (b - c  + d)
  #   } else {
  #     ## Equation S8 in supplementary materials
  #     FI_Beta_num <- p_lambda[i]
  #     FI_Beta_den <- 1 + p_theta*p_lambda[i]
  #     FI_Beta_arg <- (FI_Beta_num/FI_Beta_den)
  #   }
  #   FI_Beta_iter <- FI_Beta_iter - (FI_Beta_arg)*outprod_X     
  # }
  # return(FI_Beta_iter)
  
  
}




## -- alpha, alpha                                           ####

## Equation 17, 18 Main text (non truncated); 
## Assumes:
## -- p_lambda <- linkingFun(beta_iter = reg_coeff_iter, X = X)  without offset ie., just exp(X %*% beta_iter)
## -- arbitrary upper bound M
## Notes:
## -- For now zero truncation can be deduced from Eq.S26 and S31 (second derivative) in the supplementary materials.  
## -- Theta must call this function with zerotrunc = FALSE

FI_alpha <- function(p_lambda, p_alpha, X, M =  20, rate_param = FALSE, offset_feat = NA, offset_feat_value = NA, zerotrunc = zerotrunc){
  
  ## Debug
  #p_lambda = lambda_iter                               
  #p_alpha =  chosen_alpha
  #X = X_data 
  
  if(rate_param & !is.na(offset_feat)){
    p_lambda_offset <- p_lambda*offset_feat_value                             
    p_lambda <- p_lambda_offset
  }
  fi_a2 <- sapply(1:length(p_lambda), function(i){
    x_i <- X[i,]
    p_lambda_i <- p_lambda[i]
    fi_a1 <-  sapply(0:M, function(j){
      nb_CDF_j <- negbin_CDF(j, p_lambda = p_lambda_i, p_alpha)
      (1-nb_CDF_j )/(j+p_alpha)^2
    })
    sum(fi_a1) - (p_lambda_i/(p_alpha*(p_alpha + p_lambda_i)))
  })
  if(zerotrunc){
    ## Same "addon" as in second derivative function
    ## Here i use the shorter version in eq. S31 Supplementary materials
    if(p_alpha == 0){
      p_theta = 0
    } else {
      p_theta <- 1/p_alpha
    }
    d <- 1 + p_theta*p_lambda
    phi <- 1/d^p_alpha
    f <- p_lambda / (p_alpha + p_lambda)
    g_fun <-  f - log(1 + p_theta*p_lambda) 
    dg_fun <- p_theta * f^2 
    e_short <- (phi / (1 - phi) ) * ( (g_fun^2 / (1 - phi)) + dg_fun)
    e <- sum(e_short)
  } else {
    e <- 0
  }
  sum(fi_a2) - e
}



## -- theta, theta                                           ####

## Equation 20 and 21 (main text); S45 (Supplementary materials): 
## Assumes:
## -- p_lambda <- linkingFun(beta_iter = reg_coeff_iter, X = X)  without offset ie., just exp(X %*% beta_iter)

FI_theta <- function(p_lambda, p_alpha, X, M =  20, rate_param = FALSE, offset_feat = NA, offset_feat_value = NA, zerotrunc = FALSE){
  
  ## Debug
  #p_lambda = lambda_iter                                
  #p_alpha =  chosen_alpha
  #X = X_data 
  
  
  if(p_alpha == 0){
    p_theta <- 0
  } else {
    p_theta <- 1/p_alpha
  }
  EI_alpha <- FI_alpha(p_lambda = p_lambda,                                
                       p_alpha =  p_alpha,
                       X = X,
                       M = 50,
                       rate_param = rate_param, offset_feat = offset_feat, offset_feat_value = offset_feat_value, 
                       zerotrunc = FALSE                                        ## KEEP IT SET TO FALSE IN THIS CALL
  )
  if(rate_param & !is.na(offset_feat)){
    p_lambda_offset <- p_lambda*offset_feat_value                             
    p_lambda <- p_lambda_offset
  }
  fi_t2 <- sapply(1:length(p_lambda), function(i){
    x_i <- X[i,]
    p_lambda_i <- p_lambda[i]
    fi_t1 <-  sapply(0:M, function(j){
      nb_CDF_j <- negbin_CDF(j, p_lambda = p_lambda_i, p_alpha)
      (1-nb_CDF_j )/(j+p_alpha)
    })
    - sum(fi_t1) + log(1+p_theta*p_lambda_i)
  })
  ##Check main paper: 2/(p_theta^3) * sum(fi_t2)  should be close to 0
  a <- 2/(p_theta^3) * sum(fi_t2)
  b <- 1/(p_theta^4)*EI_alpha                               
  if(zerotrunc){
    d <- 1 + p_theta*p_lambda
    phi <- 1/d^p_alpha
    dphi <- ( (p_lambda / (p_alpha*d) ) * phi ) - ( log(d) * phi)
    d2phi <-  dphi^2/ phi + ( (p_theta*phi) * (p_lambda / (p_alpha*d))^2 ) 
    e <- (dphi/(1 - phi))^2 +  (d2phi/(1 - phi))
    c1 <- 1/(p_theta^4)*(dphi / (1 - phi))^2
    c2 <- 1/(1 - phi)*( 1/(p_theta^4)*d2phi + 2/(p_theta^3)*dphi )
    c <- c1 + c2
    FI_theta_out <- a + b - sum(c)
  } else {
    FI_theta_out <- a + b
    c <- NA
  }
  
  return(list(FI_theta_out = FI_theta_out,
              dLLdalpha_check = a, 
              Lawless = b,
              ZT_term = sum(c)))
}


## -- beta, alpha or theta (= alpha or theta, beta): zero truncated only       ####

## Equation S44 Supplementary Materials
## Assumes:
## -- p_lambda <- linkingFun(beta_iter = reg_coeff_iter, X = X)  without offset ie., just exp(X %*% beta_iter)
## Note:
## -- Only relevant for the zero truncated case. For the parent specification this is just zero
## -- the cross-derivatives [alpha, beta] and [beta, alpha] yield the same result

FI_x_alpha <- function(X, y, p_alpha, p_lambda, rate_param = FALSE, offset_feat = NA, offset_feat_value = NA){
  
  ## Similar to xd_NB_alpha except for the use of the mean function and sign (expectation of minus second derivative)
  
  ## Debug
  #X = X_data
  #y = y_data
  #p_alpha =  chosen_alpha
  #p_lambda = lambda_iter
  
  if(p_alpha == 0){
    p_theta <- 0
  } else {
    p_theta <- 1/p_alpha
  }
  if(rate_param & !is.na(offset_feat)){
    p_lambda_offset <- p_lambda*offset_feat_value                             
    p_lambda <- p_lambda_offset
  }
  n_observ<- nrow(X)
  stacked_g<- sapply(1:n_observ, function(i){                                       
    i <- as.numeric(i)
    d <- 1 + p_theta*p_lambda[i]
    phi <- 1/d^p_alpha
    dphi <- ( (p_lambda[i] / (p_alpha*d) ) * phi ) - ( log(d) * phi)
    f <- p_lambda[i] / (p_alpha + p_lambda[i])
    g_fun <-  f - log(1 + p_theta*p_lambda[i]) 
    dg_fun <- p_theta * f^2 
    y_mean_i <- E_NB(p_lambda = p_lambda[i], 
                     p_alpha = p_alpha, 
                     rate_param = rate_param, 
                     offset_feat = offset_feat,
                     offset_feat_value = offset_feat_value[i],                  
                     zerotrunc = TRUE)
    a <- (p_lambda[i]*y_mean_i)/(p_alpha*(p_alpha + p_lambda[i])^2)
    b <- dg_fun
    c <- (p_lambda[i]/(p_alpha + p_lambda[i])) * (phi/(1 - phi)) * ( (p_lambda[i]/(p_alpha*(p_alpha + p_lambda[i]))) + (dphi / (1 - phi)) + g_fun )
    temp_g_i  <- (a - b - c)
    -1*p_alpha*temp_g_i*X[i,]
  })
  if(!is.matrix(stacked_g)){
    xd_beta_alpha <- sum(stacked_g)
  } else {
    xd_beta_alpha <- as.matrix(apply(stacked_g,1,sum)) 
  }
  return( xd_beta_alpha )
}


FI_x_theta <- function(X, y, p_alpha, p_lambda, rate_param = FALSE, offset_feat = NA, offset_feat_value = NA){
  -1*(p_alpha^2)*FI_x_alpha(X, y, p_alpha, p_lambda, 
                            rate_param = rate_param, offset_feat = offset_feat, offset_feat_value = offset_feat_value)
}

## D) MAIN - Alternating fitting procedure                   ####
NegBReg_altern <- function(myData, target_feat =  target_feat, rate_param = FALSE, offset_feat = NA, zerotrunc = FALSE){
  
  ## 0)Initialise with "one-off" Poisson regression model (self-contained)      ####
  my_Pois_NR <- NR_MLE(myData, target_feat,  
                       #rate_param = rate_param, offset_feat = offset_feat, zerotrunc = zerotrunc, halfstepping = TRUE)
                       rate_param = rate_param, offset_feat = offset_feat, zerotrunc = FALSE, halfstepping = TRUE)
  temp_mean <- my_Pois_NR$predict_y                 
  temp_lambda <- my_Pois_NR$lambda
  temp_beta <- my_Pois_NR$coef_est        
  
  ## Test vs pre-built, rate parametrization case
  # obj_pois <- glm(g_DISRUPTIONS ~ b_IS_BIO + g_CENTRES + g_SUBJECTS,
  #     family = "poisson",
  #     data = myData,
  #     offset = log(d_TIME_offset))
  # obj_pois$coefficients 
  # my_Pois_NR$coef_est
  
  
  ## 1) Alpha, initial                                                          ####
  if(rate_param & !is.na(offset_feat)){
    col_to_remove <- c(which(colnames(myData)==offset_feat),which(colnames(myData)==target_feat))
    regressors_lablels <- colnames(myData)[-col_to_remove]
    offset_feat_value <- as.matrix(myData[, offset_feat])
  } else {
    col_to_remove <-which(colnames(myData)==target_feat)
    regressors_lablels <- colnames(myData)[-col_to_remove]
    offset_feat_value <- NA
  }
  y_data <- as.matrix(myData[,target_feat, drop = FALSE])           
  if(length(regressors_lablels) == 0){                                   ## are we working with the intercept only?
    X_data <-  cbind(1, myData)                                          ## no regressors besides the intercept (to be added)
    X_data <- as.matrix(X_data[,-(1+col_to_remove), drop = FALSE])       ## fixed Jan 2026 - created an issue in stepwise
    colnames(X_data) <- "intercept"
    rownames(X_data) <- rownames(myData)
  } else {
    X_data <- myData[,regressors_lablels]                                ## business as usual
    X_data <- cbind(1, as.matrix(X_data))                                ## remember to always add an intercept
    colnames(X_data) <- c("intercept", regressors_lablels)
  }
  n_regressors <- ncol(X_data) 
  n_observ <- nrow(X_data) 
  degrOfFreed <- n_observ - n_regressors
  
  temp_alpha <- n_observ/sum(((y_data/temp_mean)-1)^2)
  alpha_1 <- NR_MLE_NB_alpha(X_data = X_data, y_data = y_data, p_alpha_init = temp_alpha, p_beta_fix = temp_beta,
                             rate_param = rate_param , offset_feat = offset_feat, offset_feat_value = offset_feat_value, zerotrunc = zerotrunc
  )
  
  ## 2) Alternating procedure                                                   ####
  if(alpha_1$p_alpha < 0){ chosen_alpha <- 0} else { chosen_alpha <- alpha_1$p_alpha }
  alpha_iter_lst  <- list()
  max_iter3 <- big_alpha <- 100 
  tol_gap3 <- 1e-6
  gap3 <- 1                                                               
  hessian_is_singular <- FALSE
  
  ## inherited from MASS::glm.nb
  d1 <- sqrt(2*max(1,degrOfFreed))                                                                                       
  d2 <- 1
  Lm <- loglik_NB(X_data, 
                  y_data, 
                  p_beta = temp_beta, 
                  p_alpha = chosen_alpha,
                  rate_param = rate_param, offset_feat = offset_feat, offset_feat_value = offset_feat_value, zerotrunc = zerotrunc
  )
  Lm0 <- Lm + 2*d1                                                                             
  stopping_criterion <- (abs(Lm0 - Lm)/d1 + abs(gap3)/d2)                                           
  tol_stop <- 1e-4
  iter_count3 <- 1
  alpha_iter_lst[[iter_count3]] <- chosen_alpha 
  
  while ((iter_count3 < max_iter3) && (abs(gap3) > tol_gap3) && (stopping_criterion > tol_stop )) {
    alpha_iter <- as.numeric(alpha_iter_lst[[iter_count3]]) 
    negbin_fit_iter <- NR_MLE_NB(X_data = X_data, y_data = y_data, p_alpha_fix = alpha_iter, p_beta_init = temp_beta,
                                 rate_param = rate_param, offset_feat = offset_feat, offset_feat_value = offset_feat_value, zerotrunc = zerotrunc
    )
    lambda_iter <- negbin_fit_iter$lambda   
    beta_iter <- negbin_fit_iter$coef_est
    temp_beta <- beta_iter 
    mean_iter <- negbin_fit_iter$predict_y
    temp_alpha <- n_observ/sum(((y_data/mean_iter)-1)^2)
    alpha_refined <- NR_MLE_NB_alpha(X_data = X_data, y_data = y_data, p_alpha_init = temp_alpha, p_beta_fix = beta_iter,
                                     rate_param = rate_param , offset_feat = offset_feat, offset_feat_value = offset_feat_value,  zerotrunc = zerotrunc)
    
    if(alpha_refined$p_alpha < 0){ chosen_alpha <- 0} else { chosen_alpha <- alpha_refined$p_alpha }
    gap3 <- chosen_alpha - alpha_iter
    Lm0 <- Lm
    Lm <- loglik_NB(X_data, 
                    y_data, 
                    p_beta = beta_iter, 
                    p_alpha = chosen_alpha,
                    rate_param = rate_param, offset_feat = offset_feat, offset_feat_value = offset_feat_value, zerotrunc = zerotrunc
    )  
    stopping_criterion <- (abs(Lm0 - Lm)/d1 + abs(gap3)/d2) 
    iter_count3 <- iter_count3 + 1
    alpha_iter_lst[[iter_count3]] <- chosen_alpha 
  }
  
  ## 3) Information                                ####
  ## --- Expected Info: beta                       ####
  EIM_beta <- FI_Beta(p_lambda = linkingFun(beta_iter, X_data),   ## For ZT Eq. S41
                      p_alpha = chosen_alpha, 
                      X = X_data, 
                      y = y_data,
                      rate_param = rate_param, offset_feat = offset_feat, offset_feat_value = offset_feat_value, zerotrunc = zerotrunc)
  
  ## --- Expected Info: alpha                      ####
  EIM_alpha <- FI_alpha(p_lambda = lambda_iter,                                 ## For ZT: NO eq at the moment
                        p_alpha =  chosen_alpha, 
                        X_data, 
                        M = 50,
                        rate_param = rate_param, offset_feat = offset_feat, offset_feat_value = offset_feat_value, zerotrunc = zerotrunc)
  
  
  ## --- Expected Info: theta                      ####
  EIM_theta <- FI_theta(p_lambda = lambda_iter,                                 ## For ZT: Eq. S45
                        p_alpha =  chosen_alpha,
                        X = X_data, 
                        M = 50,
                        rate_param = rate_param, offset_feat = offset_feat, offset_feat_value = offset_feat_value, zerotrunc = zerotrunc)
  
  if(zerotrunc){
    ## --- Expected info: x-derivatives                                           ## For ZT: Eq. S
    EIM_offdiag_a <- FI_x_alpha(X = X_data, 
                                y = y_data,
                                p_alpha =  chosen_alpha,
                                p_lambda = lambda_iter, 
                                rate_param = rate_param, offset_feat = offset_feat, offset_feat_value = offset_feat_value)
    
    
    EIM_offdiag <- FI_x_theta(X = X_data, 
                              y = y_data,
                              p_alpha =  chosen_alpha,
                              p_lambda = lambda_iter, 
                              rate_param = rate_param, offset_feat = offset_feat, offset_feat_value = offset_feat_value)
    
    EIM_offdiag_a_log <- chosen_alpha*EIM_offdiag_a                                      # due to chain rule
    EIM_offdiag_log  <- (1/chosen_alpha)*EIM_offdiag                                     # due to chain rule
    
  } else {
    EIM_offdiag_a <- EIM_offdiag <- EIM_offdiag_a_log <- EIM_offdiag_log <- NA
  }
  
  ## --- Expected Info: log theta                  ####
  EIM_logtheta <- 1/chosen_alpha^2 * EIM_theta$FI_theta_out             # The full expression includes the first derivative in theta too, the expectation of which is zero see supplement Proposition S3
  
  
  ## --- Expected Info: log alpha                  ####
  EIM_logalpha <- chosen_alpha^2 * EIM_alpha                            # The full expression includes the first derivative in alpha too, the expectation of which is zero see supplement Proposition S3
  
  
  ## --- Var, and SE (expected info)                            ####
  if(!zerotrunc){
    ## Diagonal Information matrix allows separate inversions
    vcov_NB <- qr.solve(EIM_beta, tol = 1e-18)
    SE_model_NB_expected_beta <- sqrt(diag(vcov_NB))
    SE_model_NB_expected_beta2 <- SE_model_NB_expected_beta3 <- SE_model_NB_expected_beta4 <-  SE_model_NB_expected_beta        # for output reporting
    
    var_alpha_NB <- 1/EIM_alpha
    SE_model_NB_expected_alpha <- sqrt(var_alpha_NB)                               # SE for alpha, expected
    
    var_theta_NB <- 1/EIM_theta$FI_theta_out                                       
    SE_model_NB_expected_theta <- sqrt(var_theta_NB)                               # SE for theta, expected
    
    var_logalpha_NB <- 1/EIM_logalpha
    SE_model_NB_expected_logalpha <- sqrt(var_logalpha_NB)                         # SE for log alpha, expected
    
    var_logtheta_NB <- 1/EIM_logtheta
    SE_model_NB_expected_logtheta <- sqrt(var_logtheta_NB)                         # SE for log theta, expected
    
    
  } else {
    ## Fisher information matrix  no longer diagonal
    EIM <- rbind(cbind(EIM_beta, EIM_offdiag),
                 cbind(t(EIM_offdiag), EIM_theta$FI_theta_out))
    vcov_NB_expect_theta <-  qr.solve(EIM, tol = 1e-18)
    SE_model_NB_expected_all_theta <- sqrt(diag(vcov_NB_expect_theta))
    SE_model_NB_expected_beta <- SE_model_NB_expected_all_theta[1:n_regressors]
    SE_model_NB_expected_theta <- SE_model_NB_expected_all_theta[n_regressors+1]
    
    ## In alpha
    EIM_2 <- rbind(cbind(EIM_beta, EIM_offdiag_a),
                   cbind(t(EIM_offdiag_a), EIM_alpha))
    vcov_NB_expect_alpha <-  qr.solve(EIM_2, tol = 1e-18)
    SE_model_NB_expected_all_alpha <- sqrt(diag(vcov_NB_expect_alpha))
    SE_model_NB_expected_beta2 <- SE_model_NB_expected_all_alpha[1:n_regressors]
    SE_model_NB_expected_alpha <- SE_model_NB_expected_all_alpha[n_regressors+1]
    
    ## In log alpha
    EIM_3 <- rbind(cbind(EIM_beta, EIM_offdiag_a_log),
                   cbind(t(EIM_offdiag_a_log), EIM_logalpha))
    vcov_NB_expect_logalpha <-  qr.solve(EIM_3, tol = 1e-18)
    SE_model_NB_expected_all_logalpha <- sqrt(diag(vcov_NB_expect_logalpha))
    SE_model_NB_expected_beta3 <- SE_model_NB_expected_all_logalpha[1:n_regressors]
    SE_model_NB_expected_logalpha <- SE_model_NB_expected_all_logalpha[n_regressors+1]
    
    ## in log theta
    EIM_4 <- rbind(cbind(EIM_beta, EIM_offdiag_log),
                   cbind(t(EIM_offdiag_log), EIM_logtheta))
    vcov_NB_expect_logtheta <-  qr.solve(EIM_4, tol = 1e-18)
    SE_model_NB_expected_all_logtheta <- sqrt(diag(vcov_NB_expect_logtheta))
    SE_model_NB_expected_beta4 <- SE_model_NB_expected_all_logtheta[1:n_regressors]
    SE_model_NB_expected_logtheta <- SE_model_NB_expected_all_logtheta[n_regressors+1]
    
  } 
  
  
  
  ## Relationship between FI theta and FI_alpha (after https://math.stackexchange.com/a/5145608/1234572)
  FI_theta_Lawless <- chosen_alpha^4*EIM_alpha
  FI_theta_approx <- EIM_theta$FI_theta_out
  FI_theta_approx_Lawless <- EIM_theta$Lawless
  
  
  
  ## --- Observed info                                          ####
  Info_NB_obs_beta <- -1*hessian_NB(p_lambda = lambda_iter, 
                                    p_alpha = chosen_alpha, 
                                    X = X_data, 
                                    y = y_data,
                                    rate_param = rate_param, offset_feat = offset_feat, offset_feat_value = offset_feat_value, zerotrunc = zerotrunc)
  
  
  Info_NB_obs_alpha <- -1*Hessian_NB_alpha(y=y_data, 
                                           p_alpha = chosen_alpha, 
                                           p_lambda = lambda_iter,
                                           rate_param = rate_param, offset_feat = offset_feat, offset_feat_value = offset_feat_value, zerotrunc = zerotrunc)
  
  Info_NB_obs_theta <- -1*Hessian_NB_theta(y=y_data, 
                                           p_alpha = chosen_alpha, 
                                           p_lambda = lambda_iter,
                                           rate_param = rate_param, offset_feat = offset_feat, offset_feat_value = offset_feat_value, zerotrunc = zerotrunc)
  
  
  ## --- Observed info: log alpha; log theta                    ####
  ## Most papers and software would set z = ln theta (or ln alpha) and differentiate the loglikelihood wrt z 
  ## By application of the chain rule
  ## -- d theta d z = exp(z) = theta  (or alpha); and d^2 theta d z^2 = exp(z) too
  ## -- d LL d z = d LL d theta * exp(z) = d LL d theta * theta   [equivalently d LL d alpha * alpha ]
  ## -- d^2LL d z^2 = d^2 LL d theta^2 * exp(z)^2 + d LL d theta  * exp(z) = d^2 LL d theta^2 * theta^2 + d LL d theta * theta  [equivalently d^2 LL d alpha^2 * alpha^2 + d LL d alpha * alpha ] 
  
  Info_NB_obs_logalpha <- -1*Hessian_NB_logalpha(y=y_data, 
                                                 p_alpha = chosen_alpha, 
                                                 p_lambda = lambda_iter,
                                                 rate_param = rate_param, offset_feat = offset_feat, offset_feat_value = offset_feat_value, zerotrunc = zerotrunc)
  
  Info_NB_obs_logtheta <- -1*Hessian_NB_logtheta(y=y_data, 
                                                 p_alpha = chosen_alpha, 
                                                 p_lambda = lambda_iter,
                                                 rate_param = rate_param, offset_feat = offset_feat, offset_feat_value = offset_feat_value, zerotrunc = zerotrunc)
  
  
  ## --- Var, and SE (observed info)                            ####
  if(!zerotrunc){
    ## Diagonal Inormation matrix allow separate inversions
    vcov_NB_obs <-  qr.solve(Info_NB_obs_beta, tol = 1e-18)
    SE_model_NB_obs_beta <- sqrt(diag(vcov_NB_obs))
    SE_model_NB_obs2_beta <- E_model_NB_obs3_beta <- E_model_NB_obs4_beta <- SE_model_NB_obs_beta                               # for output reporting only
    
    var_NB_obs_alpha <- 1/Info_NB_obs_alpha
    SE_model_NB_obs_alpha <- sqrt(var_NB_obs_alpha)
    
    var_NB_obs_theta <- 1/Info_NB_obs_theta
    SE_model_NB_obs_theta <- sqrt(var_NB_obs_theta)
    
    ## SE log alpha can be compared with packages countreg and VGAM:vglm
    var_NB_obs_logalpha <- 1/Info_NB_obs_logalpha
    SE_model_NB_obs_logalpha <- sqrt(var_NB_obs_logalpha)
    
    ## SE log theta can be compared with package gamlss
    var_NB_obs_logtheta <- 1/Info_NB_obs_logtheta 
    SE_model_NB_obs_logtheta <- sqrt(var_NB_obs_logtheta)
    
    
  } else {
    ## form the observed information matrix and invert
    ## G et the off-diagonal elements (x-derivatives) as these are no longer zero 
    OIM_offdiag_alpha <- xd_NB_alpha(X = X_data, 
                                     y = y_data,
                                     p_alpha =  chosen_alpha,
                                     p_lambda = lambda_iter, 
                                     rate_param = rate_param, offset_feat = offset_feat, offset_feat_value = offset_feat_value)
    OIM_offdiag_theta <- -(1/chosen_alpha^2) * OIM_offdiag_alpha
    
    ## log alpha (for comparison with countreg's ztrunc and VGAM posnegbinomial family)
    ## -- z = ln alpha --> alpha = exp(z)
    ## -- d alpha d z = exp(z) = alpha 
    ## -- d2 LL d beta d z = d^2 LL d beta d alpha * exp(z) = d^2 LL d beta d alpha * alpha 
    OIM_offdiag_logalpha <- OIM_offdiag_alpha*chosen_alpha
    
    ## log theta
    ## -- z = ln theta --> theta = exp(z)
    ## -- d theta d z = exp(z) = theta 
    ## -- d2 LL d beta d z = d^2 LL d beta d theta * exp(z) = d^2 LL d beta d theta * theta 
    OIM_offdiag_logtheta <-   OIM_offdiag_theta*(1/chosen_alpha)
    
    
    ## Observed info matrix in theta
    OIM <- rbind(cbind(Info_NB_obs_beta, OIM_offdiag_theta),
                 cbind(t(OIM_offdiag_theta), Info_NB_obs_theta))
    rownames(OIM) <- c(rownames(OIM_offdiag_alpha),"Info_NB_obs_theta") 
    vcov_NB_obs_theta <-  qr.solve(OIM, tol = 1e-18)
    SE_model_NB_obs <- sqrt(diag( vcov_NB_obs_theta))
    SE_model_NB_obs_beta <- SE_model_NB_obs[1:n_regressors]
    SE_model_NB_obs_theta <- SE_model_NB_obs[(n_regressors+1)]
    
    ## Observed info matrix in alpha
    OIM2 <- rbind(cbind(Info_NB_obs_beta, OIM_offdiag_alpha),
                  cbind(t(OIM_offdiag_alpha), Info_NB_obs_alpha))
    rownames(OIM2) <- c(rownames(OIM_offdiag_alpha),"Info_NB_obs_alpha")
    vcov_NB_obs_alpha <-  qr.solve(OIM2, tol = 1e-18)
    SE_model_NB_obs2 <- sqrt(diag( vcov_NB_obs_alpha))
    SE_model_NB_obs2_beta <- SE_model_NB_obs2[1:n_regressors]
    SE_model_NB_obs_alpha <- SE_model_NB_obs2[(n_regressors+1)]
    
    ## Observed info in log alpha
    OIM3 <- rbind(cbind(Info_NB_obs_beta, OIM_offdiag_logalpha),
                  cbind(t(OIM_offdiag_logalpha), Info_NB_obs_logalpha))
    rownames(OIM3) <- c(rownames(OIM_offdiag_logalpha),"Info_NB_obs_logalpha")
    vcov_NB_obs_logalpha <-  qr.solve(OIM3, tol = 1e-18)
    SE_model_NB_obs3 <- sqrt(diag( vcov_NB_obs_logalpha))
    SE_model_NB_obs3_beta <- SE_model_NB_obs3[1:n_regressors]
    SE_model_NB_obs_logalpha <- SE_model_NB_obs3[(n_regressors+1)]
    
    ## Observed info in log theta
    OIM4 <- rbind(cbind(Info_NB_obs_beta, OIM_offdiag_logtheta),
                  cbind(t(OIM_offdiag_logtheta), Info_NB_obs_logtheta))
    rownames(OIM4) <- c(rownames(OIM_offdiag_alpha),"Info_NB_obs_logtheta") 
    vcov_NB_obs_logtheta <-  qr.solve(OIM4, tol = 1e-18)
    SE_model_NB_obs4 <- sqrt(diag( vcov_NB_obs_logtheta))
    SE_model_NB_obs4_beta <- SE_model_NB_obs4[1:n_regressors]
    SE_model_NB_obs_logtheta <- SE_model_NB_obs4[(n_regressors+1)]
    
  }
  
  
  
  
  
  ## 4) out                                                     ####
  return(list(p_alpha_MLE=chosen_alpha, 
              p_lambda_MLE=lambda_iter, 
              p_beta_MLE=beta_iter,
              p_mean = negbin_fit_iter$predict_y,
              log_likelihood=Lm, 
              n_iter=iter_count3,
              negbin_fit_iter$GG_count,
              SE_expected = list(SE_expected_all_theta = list(SE_expected_beta = SE_model_NB_expected_beta,
                                                              SE_expected_theta = SE_model_NB_expected_theta),
                                 SE_expeted_all_alpha = list(SE_expected_beta = SE_model_NB_expected_beta2,      # slightly different for ZT
                                                             SE_expected_alpha = SE_model_NB_expected_alpha),
                                 SE_expected_all_logtheta = list(SE_expected_beta = SE_model_NB_expected_beta4,
                                                                 SE_expected_logtheta = SE_model_NB_expected_logtheta),
                                 SE_expeted_all_logalpha = list(SE_expected_beta = SE_model_NB_expected_beta3,      # slightly different for ZT
                                                                SE_expected_logalpha = SE_model_NB_expected_logalpha)
                                 
              ),
              SE_oberved = list(SE_obs_all_theta = list(SE_obs_beta = SE_model_NB_obs_beta,
                                                        SE_obs_theta = SE_model_NB_obs_theta),
                                SE_obs_all_alpha = list(SE_obs_beta = SE_model_NB_obs2_beta,                          # slightly different
                                                        SE_obs_alpha = SE_model_NB_obs_alpha),
                                SE_obs_all_logtheta = list(SE_obs_beta = SE_model_NB_obs_beta,
                                                           SE_obs_logtheta = SE_model_NB_obs_logtheta),
                                SE_obs_all_logalpha = list(SE_obs_beta = SE_model_NB_obs2_beta,                       # slightly different
                                                           SE_obs_logalpha = SE_model_NB_obs_logalpha)
                                
              ),
              Lawless_FI_theta = list(FI_theta_Lawless = FI_theta_Lawless,
                                      FI_theta_approx = FI_theta_approx,
                                      FI_theta_approx_Lawless = FI_theta_approx_Lawless
              )
  )
  )
}


##                                                           ####
## E) Significance                                           ####
myMLEsignif <- function(my_SE, my_p_est){
  
  modse_zValues <- my_p_est/my_SE
  if(!is.matrix(modse_zValues)){modse_zValues <- as.matrix(modse_zValues)}
  modse_signifc <- 2*pnorm(abs(modse_zValues), lower.tail = FALSE)
  est_summary_table <- cbind.data.frame(round(my_p_est,5),
                                        round(my_SE,5) ,
                                        round(modse_zValues,5), 
                                        round(modse_signifc,5))
  colnames(est_summary_table) <- c("Estimate", "S.E.", "z values", "Pr(>|z|)")
  est_summary_table
}




##                                                           ####
## Run numerical example 1 - baseline                        ####
## a) Set arguments                                          ####
target_feat = "g_DISRUPTIONS"

## no offset, no zt
offset_feat = NA
rate_param = FALSE
zerotrunc = FALSE

## b1) full alternating procedure                            ####
test_NB <- NegBReg_altern(myData = myData, 
                          target_feat =  target_feat,
                          rate_param = rate_param, 
                          offset_feat = offset_feat,
                          zerotrunc = zerotrunc
)

## b2) significance                                          ####
## -- expected, alpha param
test_signif_A11_beta_exp <- myMLEsignif(my_SE = test_NB$SE_expected$SE_expeted_all_alpha$SE_expected_beta, 
                                        my_p_est= test_NB$p_beta_MLE)
test_signif_A21_alpha_exp <- myMLEsignif(my_SE = test_NB$SE_expected$SE_expeted_all_alpha$SE_expected_alpha, 
                                         my_p_est= test_NB$p_alpha_MLE)
tab_signif_A_exp <- rbind(test_signif_A11_beta_exp, test_signif_A21_alpha_exp)
rownames(tab_signif_A_exp)[nrow(tab_signif_A_exp)] <- "alpha"

## -- observed, alpha param
test_signif_A12_beta_obs <- myMLEsignif(my_SE = test_NB$SE_oberved$SE_obs_all_alpha$SE_obs_beta, 
                                        my_p_est= test_NB$p_beta_MLE)
test_signif_A22_alpha_obs <- myMLEsignif(my_SE = test_NB$SE_oberved$SE_obs_all_alpha$SE_obs_alpha, 
                                         my_p_est= test_NB$p_alpha_MLE)
test_signif_A_obs <- rbind(test_signif_A12_beta_obs, test_signif_A22_alpha_obs)
rownames(test_signif_A_obs)[nrow(test_signif_A_obs)] <- "alpha"
colnames(test_signif_A_obs) <- paste0(colnames(test_signif_A_obs),"_obs")

## -- assemble
tab_signif_A <- cbind(tab_signif_A_exp, test_signif_A_obs[,-1])


## -- expected, theta (dispersion param, dp)
test_signif_A11_beta_exp_dp <- myMLEsignif(my_SE = test_NB$SE_expected$SE_expected_all_theta$SE_expected_beta, 
                                           my_p_est= test_NB$p_beta_MLE)
test_signif_A21_theta_exp_dp <- myMLEsignif(my_SE = test_NB$SE_expected$SE_expected_all_theta$SE_expected_theta, 
                                            my_p_est= 1/test_NB$p_alpha_MLE)
tab_signif_A_exp_dp <- rbind(test_signif_A11_beta_exp_dp, test_signif_A21_theta_exp_dp)
rownames(tab_signif_A_exp_dp)[nrow(tab_signif_A_exp_dp)] <- "theta"

## -- observed, theta (dispersion param, dp)
test_signif_A12_beta_obs_dp <- myMLEsignif(my_SE = test_NB$SE_oberved$SE_obs_all_theta$SE_obs_beta, 
                                           my_p_est= test_NB$p_beta_MLE)

test_signif_A22_theta_obs_dp <- myMLEsignif(my_SE = test_NB$SE_oberved$SE_obs_all_theta$SE_obs_theta, 
                                            my_p_est= 1/test_NB$p_alpha_MLE)
tab_signif_A_obs_dp <- rbind(test_signif_A12_beta_obs_dp , test_signif_A22_theta_obs_dp )
rownames(tab_signif_A_obs_dp)[nrow(tab_signif_A_obs_dp)] <- "theta"
colnames(tab_signif_A_obs_dp) <- paste0(colnames(tab_signif_A_obs_dp),"_obs")

## -- assemble
tab_signif_A_dp <- cbind(tab_signif_A_exp_dp, tab_signif_A_obs_dp[,-1])
rownames(tab_signif_A_dp)[nrow(tab_signif_A_dp)] <- "theta"
tab_signif_A_dp_just_theta <- tab_signif_A_dp[nrow(tab_signif_A_dp), ]


## -- expected, log alpha
test_signif_A11_logalpha_exp_beta <- myMLEsignif(my_SE = test_NB$SE_expected$SE_expeted_all_logalpha$SE_expected_beta, 
                                                 my_p_est= test_NB$p_beta_MLE)
test_signif_A21_logalpha_exp <- myMLEsignif(my_SE = test_NB$SE_expected$SE_expeted_all_logalpha$SE_expected_logalpha, 
                                            my_p_est= log(test_NB$p_alpha_MLE))
tab_signif_A_logalpha_exp <- rbind(test_signif_A11_logalpha_exp_beta, test_signif_A21_logalpha_exp)
rownames(tab_signif_A_logalpha_exp)[nrow(tab_signif_A_logalpha_exp)] <- "logalpha"


## -- observed, log alpha
test_signif_A12_logalpha_obs_beta <- myMLEsignif(my_SE = test_NB$SE_oberved$SE_obs_all_logalpha$SE_obs_beta, 
                                                 my_p_est= test_NB$p_beta_MLE)
test_signif_A22_logalpha_obs <- myMLEsignif(my_SE = test_NB$SE_oberved$SE_obs_all_logalpha$SE_obs_logalpha, 
                                            my_p_est= log(test_NB$p_alpha_MLE))

test_signif_A_logalpha_obs <- rbind(test_signif_A12_logalpha_obs_beta, test_signif_A22_logalpha_obs)
rownames(test_signif_A_logalpha_obs)[nrow(test_signif_A_logalpha_obs)] <- "logalpha"
colnames(test_signif_A_logalpha_obs) <- paste0(colnames(test_signif_A_logalpha_obs),"_obs")

## -- assemble
tab_signif_A_logalpha <- cbind(tab_signif_A_logalpha_exp, test_signif_A_logalpha_obs[,-1])
tab_signif_A_logalpha_just_alpha <- tab_signif_A_logalpha[nrow(tab_signif_A_logalpha), ]

## -- Expected, log theta
test_signif_A11_logtheta_exp_beta <- myMLEsignif(my_SE = test_NB$SE_expected$SE_expected_all_logtheta$SE_expected_beta, 
                                                 my_p_est= test_NB$p_beta_MLE)
test_signif_A21_logtheta_exp <- myMLEsignif(my_SE = test_NB$SE_expected$SE_expected_all_logtheta$SE_expected_logtheta, 
                                            my_p_est= log(1/test_NB$p_alpha_MLE))

tab_signif_A_logtheta_exp <- rbind(test_signif_A11_logtheta_exp_beta , test_signif_A21_logtheta_exp)
rownames(tab_signif_A_logtheta_exp )[nrow(tab_signif_A_logtheta_exp)] <- "logtheta"

## -- observed log theta
test_signif_A12_logtheta_obs_beta <- myMLEsignif(my_SE = test_NB$SE_oberved$SE_obs_all_logtheta$SE_obs_beta, 
                                                 my_p_est= test_NB$p_beta_MLE)
test_signif_A22_logtheta_obs <- myMLEsignif(my_SE = test_NB$SE_oberved$SE_obs_all_logtheta$SE_obs_logtheta, 
                                            my_p_est= log(1/test_NB$p_alpha_MLE))

test_signif_A_logtheta_obs <- rbind(test_signif_A12_logtheta_obs_beta, test_signif_A22_logtheta_obs)
rownames(test_signif_A_logtheta_obs)[nrow(test_signif_A_logtheta_obs)] <- "logtheta"
colnames(test_signif_A_logtheta_obs) <- paste0(colnames(test_signif_A_logtheta_obs),"_obs")


## -- assembled tables for each parameterization
tab_signif_A_logtheta <- cbind(tab_signif_A_logtheta_exp, test_signif_A_logtheta_obs[,-1])
tab_signif_A_logtheta_just_theta <- tab_signif_A_logtheta[nrow(tab_signif_A_logtheta), ]




## c) Compare with pre-built                   ####
## --- MASS:glm.nb                             ####
obj_NB <- MASS::glm.nb(g_DISRUPTIONS ~.,
                       data = myData)

obj_NB$coefficients
obj_NB$theta                                                  # in fact this is what we call alpha
1/obj_NB$theta                                                # dispersion parameter, to be compared with GAMLSS                               



## --- countreg:                               ####

obj_NB_cr <- countreg::nbreg(g_DISRUPTIONS ~.,
                             data = myData)

obj_NB_cr$coefficients
obj_NB_cr$coefficients.theta                        # this is ln alpha

summary(obj_NB_cr)[1]
summary(obj_NB_cr)[2]

countreg_A_fulltab_beta <- as.data.frame(summary(obj_NB_cr)[1])
countreg_A_fulltab_logalpha <- as.data.frame(summary(obj_NB_cr)[2])

countreg_A_beta <-  countreg_A_fulltab_beta[,1,drop = FALSE]
countreg_A_beta_SE <-  countreg_A_fulltab_beta[,2,drop = FALSE]
countreg_A_logalpha <- countreg_A_fulltab_logalpha$coefficients.theta.Estimate
countreg_A_logalpha_SE <- countreg_A_fulltab_logalpha$coefficients.theta.Std..Error




## --- GAMLSS                                  ####
obj_NB_gam <- gamlss::gamlss(g_DISRUPTIONS ~.,
                             family = "NBI",
                             data = myData) 

summary(obj_NB_gam)[1:5,]

## extract estimators
gamlss_A_beta <- obj_NB_gam$mu.coefficients
gamlss_A_logtheta <- obj_NB_gam$sigma.coefficients                                 # "as is" output is ln theta
gamlss_A_theta_computed <- exp(obj_NB_gam$sigma.coefficients)
gamlss_A_alpha_computed <- 1/gamlss_A_theta_computed 


## extract standard errors (workaround)
var_gamlss <- vcov(obj_NB_gam)
var_gamlss_coeff <- diag(var_gamlss[1:(nrow(var_gamlss)-1), 1:(ncol(var_gamlss)-1)])
var_gamlss_logtheta <- diag(var_gamlss)[nrow(var_gamlss)]
gamlss_A_beta_SE <- sqrt(var_gamlss_coeff)
gamlss_A_logtheta_SE <- sqrt(var_gamlss_logtheta)



## --- VGAM                                    ####

#library(VGAM)
## Syntax for fitting NB2 as in doi: 10.1111/anzs.12283 section 2
obj_NB_vgam <- VGAM::vglm(g_DISRUPTIONS ~. , 
                          family  = negbinomial,
                          data = myData)

summary(obj_NB_vgam)
coef(obj_NB_vgam, matrix = TRUE)

obj_NB_vgam@coefficients
obj_NB_vgam@predictors

## Where is the dispersion coefficient?
## Based on doi: 10.1111/anzs.12283 equation 1 and 4: they estimate what I call alpha which is what they call k; but it is presented as log k ("eta 2" in eq. 4)
## Additional useful information: lookup "negbinomial" https://cran.r-project.org/web/packages/VGAM/VGAM.pdf
## -- p 606 find argument "nsimEIM": mentions that it is used for computing the diagonal element of the expected information matrix (EIM) corresponding to k based on the simulated Fisher scoring (SFS) algorithm."
## -- P 606 also mentions that the argument "cutoff.prob" i used "to specify how many terms of the infinite series for computing the second diagonal element of the EIM are actually used"
## -- then again p 199 set nsimEIM = NULL to choose the other (EXACT?) algorithm
## -- p 607 mentions the trigamma function and its expectations when computing the EIM with respect to the "size" parameter
## -- p.608 mentions it estimates both parmaters of the neg bin by full maximum likelihood estimation 
## -- p.609 mentions "For negbinomial() the diagonal element of the expected information matrix (EIM) for parameter k involves an infinite series; consequently SFS (see nsimEIM) is used as the backup algorithm
var_vgam <- vcov(obj_NB_vgam)
vgam_A_beta <- obj_NB_vgam@coefficients[c(1,3:length(obj_NB_vgam@coefficients))]
vgam_A_beta_sE <- round(sqrt(diag(var_vgam)[c(1,3:length(obj_NB_vgam@coefficients))]),5)

vgam_A_logalpha <- obj_NB_vgam@coefficients[2]
vgam_A_logalpha_SE <- sqrt(diag(var_vgam)[2])   

vgam_A_alpha_computed <- exp(obj_NB_vgam@coefficients[2])



## --- Table comparative                       ####

tab_fromScratch <- rbind.data.frame(tab_signif_A, 
                                    tab_signif_A_dp_just_theta,
                                    tab_signif_A_logalpha_just_alpha,
                                    tab_signif_A_logtheta_just_theta
)

tab_pre_built_MASS <-  round(summary(obj_NB)$coefficients[,1:2],5)
tab_pre_built_countreg <- round(countreg_A_fulltab_beta[,c(1,2),drop = FALSE], 5)
tab_pre_built_GAMLSS <-  round(cbind.data.frame(gamlss_A_beta, gamlss_A_beta_SE),5)
tab_pre_built_VGAM <- round(cbind.data.frame(vgam_A_beta, vgam_A_beta_sE),5 )

tab_pre_built_alpha_MASS <- round(cbind.data.frame(summary(obj_NB)[17], summary(obj_NB)[18]),5)
tab_pre_built_alpha_countreg <- cbind.data.frame("n.a.", "n.a.")
tab_pre_built_alpha_GAMLSS <- cbind.data.frame("n.a.", "n.a.")
tab_pre_built_alpha_VGAM <- cbind.data.frame("n.a.", "n.a.")

tab_pre_built_theta_MASS <- cbind.data.frame("n.a.", "n.a.")
tab_pre_built_theta_countreg <- cbind.data.frame("n.a.", "n.a.")
tab_pre_built_theta_GAMLSS <- cbind.data.frame("n.a.", "n.a.")
tab_pre_built_theta_VGAM <-  cbind.data.frame("n.a.", "n.a.")

tab_pre_built_logalpha_MASS <- cbind.data.frame("n.a.", "n.a.")
tab_pre_built_logalpha_countreg <- round(countreg_A_fulltab_logalpha[,c(1:2)],5)
tab_pre_built_logalpha_GAMLSS <- cbind.data.frame("n.a.", "n.a.")
tab_pre_built_logalpha_VGAM <-  cbind.data.frame(round(vgam_A_logalpha,5), round(vgam_A_logalpha_SE,5))

tab_pre_built_logtheta_MASS <- cbind.data.frame("n.a.", "n.a.")
tab_pre_built_logtheta_countreg <- cbind.data.frame("n.a.", "n.a.")
tab_pre_built_logtheta_GAMLSS <- cbind.data.frame(round(gamlss_A_logtheta,5), round(gamlss_A_logtheta_SE,5))
tab_pre_built_logtheta_VGAM <-  cbind.data.frame("n.a.", "n.a.")

colnames(tab_pre_built_MASS) <- paste0(c("glmnb_"),colnames(tab_pre_built_MASS))
colnames(tab_pre_built_countreg) <- c("countreg_coeff", "countreg_SE")
colnames(tab_pre_built_GAMLSS) <- c("gamlss_coeff", "gamlss_SE")
colnames(tab_pre_built_VGAM) <- c("vglm_coeff", "vglm_SE")

a <- cbind.data.frame(tab_pre_built_MASS, tab_pre_built_countreg, tab_pre_built_GAMLSS, tab_pre_built_VGAM)
b <- cbind.data.frame(tab_pre_built_alpha_MASS, tab_pre_built_alpha_countreg, tab_pre_built_alpha_GAMLSS, tab_pre_built_alpha_VGAM)
c <- cbind.data.frame(tab_pre_built_theta_MASS, tab_pre_built_theta_countreg, tab_pre_built_theta_GAMLSS, tab_pre_built_theta_VGAM)
d <- cbind.data.frame(tab_pre_built_logalpha_MASS, tab_pre_built_logalpha_countreg, tab_pre_built_logalpha_GAMLSS, tab_pre_built_logalpha_VGAM)
e <- cbind.data.frame(tab_pre_built_logtheta_MASS, tab_pre_built_logtheta_countreg, tab_pre_built_logtheta_GAMLSS, tab_pre_built_logtheta_VGAM)

colnames(e) <- colnames(d) <- colnames(c) <- colnames(b) <- colnames(a)
rownames(b) <- "alpha"
rownames(c) <- "theta"
rownames(d) <- "log alpha"
rownames(e) <- "log theta"

tab_pre_built <- rbind.data.frame(a, b, c, d, e)
myTab4viz <- cbind(tab_fromScratch, tab_pre_built)

## Generate LaTeX export
library(kableExtra)
kbl(myTab4viz, format = "latex")




##                                                           ####
## Run numerical example 2 - offset                          ####
## a) Set arguments                                          ####
target_feat = "g_DISRUPTIONS"

## Offset, no zt
offset_feat = "d_TIME_offset"
rate_param = TRUE
zerotrunc = FALSE


## b1) full alternating procedure                            ####
test_NB_offset <- NegBReg_altern(myData = myData, 
                                 target_feat =  target_feat,
                                 rate_param = rate_param, 
                                 offset_feat = offset_feat,
                                 zerotrunc = zerotrunc
)


test_NB_offset$p_beta_MLE
test_NB_offset$p_alpha_MLE

## b2) Significance                                          ####
## -- expected, alpha param
test_signif_B11_beta_exp <- myMLEsignif(my_SE = test_NB_offset$SE_expected$SE_expeted_all_alpha$SE_expected_beta, 
                                        my_p_est= test_NB_offset$p_beta_MLE)
test_signif_B21_alpha_exp <- myMLEsignif(my_SE = test_NB_offset$SE_expected$SE_expeted_all_alpha$SE_expected_alpha, 
                                         my_p_est= test_NB_offset$p_alpha_MLE)
tab_signif_B_exp <- rbind(test_signif_B11_beta_exp, test_signif_B21_alpha_exp)
rownames(tab_signif_B_exp)[nrow(tab_signif_B_exp)] <- "alpha"

## -- observed, alpha param
test_signif_B12_beta_obs <- myMLEsignif(my_SE = test_NB_offset$SE_oberved$SE_obs_all_alpha$SE_obs_beta, 
                                        my_p_est= test_NB_offset$p_beta_MLE)
test_signif_B22_alpha_obs <- myMLEsignif(my_SE = test_NB_offset$SE_oberved$SE_obs_all_alpha$SE_obs_alpha, 
                                         my_p_est= test_NB_offset$p_alpha_MLE)
test_signif_B_obs <- rbind(test_signif_B12_beta_obs, test_signif_B22_alpha_obs)
rownames(test_signif_B_obs)[nrow(test_signif_B_obs)] <- "alpha"
colnames(test_signif_B_obs) <- paste0(colnames(test_signif_B_obs),"_obs")

## -- assemble
tab_signif_B_alpha <- cbind(tab_signif_B_exp, test_signif_B_obs[,-1])

## -- expected, theta param
test_signif_B11_beta_exp_dp <- myMLEsignif(my_SE = test_NB_offset$SE_expected$SE_expected_all_theta$SE_expected_beta, 
                                           my_p_est= test_NB_offset$p_beta_MLE)
test_signif_B21_theta_exp_dp <- myMLEsignif(my_SE = test_NB_offset$SE_expected$SE_expected_all_theta$SE_expected_theta, 
                                            my_p_est= 1/test_NB_offset$p_alpha_MLE)
tab_signif_B_exp_dp <- rbind(test_signif_B11_beta_exp_dp, test_signif_B21_theta_exp_dp)
rownames(tab_signif_B_exp_dp)[nrow(tab_signif_B_exp_dp)] <- "theta"

## -- observed, theta param
test_signif_B12_beta_obs_dp <- myMLEsignif(my_SE = test_NB_offset$SE_oberved$SE_obs_all_theta$SE_obs_beta, 
                                           my_p_est= test_NB_offset$p_beta_MLE)

test_signif_B22_theta_obs_dp <- myMLEsignif(my_SE = test_NB_offset$SE_oberved$SE_obs_all_theta$SE_obs_theta, 
                                            my_p_est= 1/test_NB_offset$p_alpha_MLE)
tab_signif_B_obs_dp <- rbind(test_signif_B12_beta_obs_dp , test_signif_B22_theta_obs_dp )
rownames(tab_signif_B_obs_dp)[nrow(tab_signif_B_obs_dp)] <- "theta"
colnames(tab_signif_B_obs_dp) <- paste0(colnames(tab_signif_B_obs_dp),"_obs")

## -- assemble
tab_signif_B_dp <- cbind(tab_signif_B_exp_dp, tab_signif_B_obs_dp[,-1])
tab_signif_B_dp_just_theta <- tab_signif_B_dp[nrow(tab_signif_B_dp),]



## -- expected, log alpha
test_signif_B11_logalpha_exp_beta <- myMLEsignif(my_SE = test_NB_offset$SE_expected$SE_expeted_all_logalpha$SE_expected_beta, 
                                                 my_p_est= test_NB_offset$p_beta_MLE)
test_signif_B21_logalpha_exp <- myMLEsignif(my_SE = test_NB_offset$SE_expected$SE_expeted_all_logalpha$SE_expected_logalpha, 
                                            my_p_est= log(test_NB_offset$p_alpha_MLE))
tab_signif_B_logalpha_exp <- rbind(test_signif_B11_logalpha_exp_beta, test_signif_B21_logalpha_exp)
rownames(tab_signif_B_logalpha_exp)[nrow(tab_signif_B_logalpha_exp)] <- "logalpha"


## -- observed, log alpha
test_signif_B12_logalpha_obs_beta <- myMLEsignif(my_SE = test_NB_offset$SE_oberved$SE_obs_all_logalpha$SE_obs_beta, 
                                                 my_p_est= test_NB_offset$p_beta_MLE)
test_signif_B22_logalpha_obs <- myMLEsignif(my_SE = test_NB_offset$SE_oberved$SE_obs_all_logalpha$SE_obs_logalpha, 
                                            my_p_est= log(test_NB_offset$p_alpha_MLE))

test_signif_B_logalpha_obs <- rbind(test_signif_B12_logalpha_obs_beta, test_signif_B22_logalpha_obs)
rownames(test_signif_B_logalpha_obs)[nrow(test_signif_B_logalpha_obs)] <- "logalpha"
colnames(test_signif_B_logalpha_obs) <- paste0(colnames(test_signif_B_logalpha_obs),"_obs")

## -- assemble
tab_signif_B_logalpha <- cbind(tab_signif_B_logalpha_exp, test_signif_B_logalpha_obs[,-1])
tab_signif_B_logalpha_just_alpha <- tab_signif_B_logalpha[nrow(tab_signif_B_logalpha), ]

## -- Expected, log theta
test_signif_B11_logtheta_exp_beta <- myMLEsignif(my_SE = test_NB_offset$SE_expected$SE_expected_all_logtheta$SE_expected_beta, 
                                                 my_p_est= test_NB_offset$p_beta_MLE)
test_signif_B21_logtheta_exp <- myMLEsignif(my_SE = test_NB_offset$SE_expected$SE_expected_all_logtheta$SE_expected_logtheta, 
                                            my_p_est= log(1/test_NB_offset$p_alpha_MLE))

tab_signif_B_logtheta_exp <- rbind(test_signif_B11_logtheta_exp_beta , test_signif_B21_logtheta_exp)
rownames(tab_signif_B_logtheta_exp )[nrow(tab_signif_B_logtheta_exp)] <- "logtheta"

## -- observed log theta
test_signif_B12_logtheta_obs_beta <- myMLEsignif(my_SE = test_NB_offset$SE_oberved$SE_obs_all_logtheta$SE_obs_beta, 
                                                 my_p_est= test_NB_offset$p_beta_MLE)
test_signif_B22_logtheta_obs <- myMLEsignif(my_SE = test_NB_offset$SE_oberved$SE_obs_all_logtheta$SE_obs_logtheta, 
                                            my_p_est= log(1/test_NB_offset$p_alpha_MLE))

test_signif_B_logtheta_obs <- rbind(test_signif_B12_logtheta_obs_beta, test_signif_B22_logtheta_obs)
rownames(test_signif_B_logtheta_obs)[nrow(test_signif_B_logtheta_obs)] <- "logtheta"
colnames(test_signif_B_logtheta_obs) <- paste0(colnames(test_signif_B_logtheta_obs),"_obs")


## -- assemble
tab_signif_B_logtheta <- cbind(tab_signif_B_logtheta_exp, test_signif_B_logtheta_obs[,-1])
tab_signif_B_logtheta_just_theta <- tab_signif_B_logtheta[nrow(tab_signif_B_logtheta), ]





## c) Compare with pre-built                                 ####

## --- MASS:glm.nb                             ####
## Note: offset is specified in a different way compared with glm, via offset()
## -- for an example of how to write an offset in glm.bn: https://stats.stackexchange.com/q/180260/513606
#obj_NB_off <- MASS::glm.nb(g_DISRUPTIONS ~ b_IS_BIO + g_CENTRES + g_SUBJECTS + offset(log(d_TIME_offset)),
obj_NB_off <- MASS::glm.nb(g_DISRUPTIONS ~ g_SUBJECTS + b_IS_BIO+ offset(log(d_TIME_offset)),
                           data = myData)
obj_NB_off$coefficients
obj_NB_off$theta
1/obj_NB_off$theta           # dispersion parameter, to be compared with GAMLSS                               


## --- countreg                                ####

obj_NB_cr_off <- countreg::nbreg(g_DISRUPTIONS ~ g_SUBJECTS + b_IS_BIO,
                                 offset = log(d_TIME_offset),
                                 data = myData)

obj_NB_cr_off$coefficients
obj_NB_cr_off$coefficients.theta                        # this is ln alpha

summary(obj_NB_cr_off)[1]
summary(obj_NB_cr_off)[2]

countreg_B_fulltab_beta <- as.data.frame(summary(obj_NB_cr_off)[1])
countreg_B_fulltab_logalpha <- as.data.frame(summary(obj_NB_cr_off)[2])

countreg_B_beta <-  countreg_B_fulltab_beta[,1,drop = FALSE]
countreg_B_beta_SE <-  countreg_B_fulltab_beta[,2,drop = FALSE]
countreg_B_logalpha <- countreg_B_fulltab_logalpha$coefficients.theta.Estimate
countreg_B_logalpha_SE <- countreg_B_fulltab_logalpha$coefficients.theta.Std..Error


## --- GAMLSS                                  ####
#obj_NB_gam_off <- gamlss::gamlss(g_DISRUPTIONS ~ b_IS_BIO + g_CENTRES + g_SUBJECTS + offset(log(d_TIME_offset)),
obj_NB_gam_off <- gamlss::gamlss(g_DISRUPTIONS ~ g_SUBJECTS + b_IS_BIO + offset(log(d_TIME_offset)),
                                 family = "NBI",
                                 data = myData) 
obj_NB_gam_off$mu.coefficients
exp(obj_NB_gam_off$sigma.coefficients)                                      # dispersion parameter
alpha_gamlss_off <- 1/exp(obj_NB_gam_off$sigma.coefficients)

##
summary(obj_NB_gam_off)

## extract estimators
gamlss_B_beta <- obj_NB_gam_off$mu.coefficients
gamlss_B_logtheta <- obj_NB_gam_off$sigma.coefficients                                 # "as is" output is ln theta
gamlss_B_theta_computed <- exp(obj_NB_gam_off$sigma.coefficients)
gamlss_B_alpha_computed <- 1/gamlss_B_theta_computed 


## extract standard errors (workaround)
var_gamlss_off <- vcov(obj_NB_gam_off)
var_gamlss_coeff_off <- diag(var_gamlss_off[1:(nrow(var_gamlss_off)-1), 1:(ncol(var_gamlss_off)-1)])
var_gamlss_logtheta_off <- diag(var_gamlss_off)[nrow(var_gamlss_off)]
gamlss_B_beta_SE <- sqrt(var_gamlss_coeff_off)
gamlss_B_logtheta_SE <- sqrt(var_gamlss_logtheta_off)









## --- VGAM                                    ####

## Syntax for fitting NB2 as in doi: 10.1111/anzs.12283 section 2
library(VGAM)
# obj_NB_vgam_off <- VGAM::vglm(g_DISRUPTIONS ~ b_IS_BIO + g_CENTRES + g_SUBJECTS + offset(log(d_TIME_offset)), 
#                           family  = negbinomial,
#                           data = myData)

#obj_NB_vgam_off <- VGAM::vglm(g_DISRUPTIONS ~ b_IS_BIO + g_CENTRES + g_SUBJECTS,
obj_NB_vgam_off <- VGAM::vglm(g_DISRUPTIONS ~  g_SUBJECTS + b_IS_BIO,
                              #offset = log(d_TIME_offset),               # this is incorrect according to Thomas
                              offset = cbind(log(d_TIME_offset), 0),      # as per thomas suggestion
                              family  = negbinomial,
                              data = myData)

summary(obj_NB_vgam_off)
coef(obj_NB_vgam_off, matrix = TRUE)

obj_NB_vgam_off@coefficients
obj_NB_vgam_off@predictors

var_vgam_off <- vcov(obj_NB_vgam_off)
vgam_B_beta <- obj_NB_vgam_off@coefficients[c(1,3:length(obj_NB_vgam_off@coefficients))]
vgam_B_beta_sE <- round(sqrt(diag(var_vgam_off)[c(1,3:length(obj_NB_vgam_off@coefficients))]),5)

vgam_B_alpha <- exp(obj_NB_vgam_off@coefficients[2])
vgam_B_logalpha <- (obj_NB_vgam_off@coefficients[2])
vgam_B_logalpha_SE <- sqrt(diag(var_vgam_off)[2])   




## --- Table comparative                       ####
tab_fromScratch_off <- rbind.data.frame(tab_signif_B_alpha,
                                        tab_signif_B_dp_just_theta,
                                        tab_signif_B_logalpha_just_alpha,
                                        tab_signif_B_logtheta_just_theta
)

tab_pre_built_MASS_off  <-  round(summary(obj_NB_off)$coefficients[,1:2],5)
tab_pre_built_countreg_off <- round(countreg_B_fulltab_beta[,c(1,2),drop = FALSE], 5)
tab_pre_built_GAMLSS_off <-  round(cbind.data.frame(gamlss_B_beta, gamlss_B_beta_SE),5)
tab_pre_built_VGAM_off <- round(cbind.data.frame(vgam_B_beta, vgam_B_beta_sE),5 )

tab_pre_built_alpha_MASS_off  <- round(cbind.data.frame(summary(obj_NB_off)[17], summary(obj_NB_off)[18]),5)
tab_pre_built_alpha_countreg_off  <- cbind.data.frame("n.a.", "n.a.")
tab_pre_built_alpha_GAMLSS_off  <- cbind.data.frame("n.a.", "n.a.")
tab_pre_built_alpha_VGAM_off  <- cbind.data.frame("n.a.", "n.a.")

tab_pre_built_theta_MASS_off <- cbind.data.frame("n.a.", "n.a.")
tab_pre_built_theta_countreg_off  <- cbind.data.frame("n.a.", "n.a.")
tab_pre_built_theta_GAMLSS_off <- cbind.data.frame("n.a.", "n.a.")
tab_pre_built_theta_VGAM_off <-  cbind.data.frame("n.a.", "n.a.")

tab_pre_built_logalpha_MASS_off  <- cbind.data.frame("n.a.", "n.a.")
tab_pre_built_logalpha_countreg_off  <- round(countreg_B_fulltab_logalpha[,c(1:2)],5)
tab_pre_built_logalpha_GAMLSS_off  <- cbind.data.frame("n.a.", "n.a.")
tab_pre_built_logalpha_VGAM_off  <-  cbind.data.frame(round(vgam_B_logalpha,5), round(vgam_B_logalpha_SE,5))

tab_pre_built_logtheta_MASS_off  <- cbind.data.frame("n.a.", "n.a.")
tab_pre_built_logtheta_countreg_off  <- cbind.data.frame("n.a.", "n.a.")
tab_pre_built_logtheta_GAMLSS_off  <- cbind.data.frame(round(gamlss_B_logtheta,5), round(gamlss_B_logtheta_SE,5))
tab_pre_built_logtheta_VGAM_off  <-  cbind.data.frame("n.a.", "n.a.")

colnames(tab_pre_built_MASS_off) <- paste0(c("glmnb_"),colnames(tab_pre_built_MASS_off))
colnames(tab_pre_built_countreg_off) <- c("countreg_coeff", "countreg_SE")
colnames(tab_pre_built_GAMLSS_off) <- paste0(c("gamlss_coeff", "gamlss_SE"))
colnames(tab_pre_built_VGAM_off) <- c("vglm_coeff", "vglm_SE")

a <- cbind.data.frame(tab_pre_built_MASS_off, tab_pre_built_countreg_off, tab_pre_built_GAMLSS_off, tab_pre_built_VGAM_off)
b <- cbind.data.frame(tab_pre_built_alpha_MASS_off, tab_pre_built_alpha_countreg_off, tab_pre_built_alpha_GAMLSS_off, tab_pre_built_alpha_VGAM_off)
c <- cbind.data.frame(tab_pre_built_theta_MASS_off, tab_pre_built_theta_countreg_off, tab_pre_built_theta_GAMLSS_off, tab_pre_built_theta_VGAM_off)
d <- cbind.data.frame(tab_pre_built_logalpha_MASS_off, tab_pre_built_logalpha_countreg_off, tab_pre_built_logalpha_GAMLSS_off, tab_pre_built_logalpha_VGAM_off)
e <- cbind.data.frame(tab_pre_built_logtheta_MASS_off, tab_pre_built_logtheta_countreg_off, tab_pre_built_logtheta_GAMLSS_off, tab_pre_built_logtheta_VGAM_off)

colnames(e) <- colnames(d) <- colnames(c) <- colnames(b) <- colnames(a)
rownames(b) <- "alpha"
rownames(c) <- "theta"
rownames(d) <- "log alpha"
rownames(e) <- "log theta"


tab_pre_built_off <- rbind.data.frame(a, b, c, d, e)
myTab4viz_off <- cbind(tab_fromScratch_off , tab_pre_built_off)

## Generate LaTeX export
library(kableExtra)
kbl(myTab4viz_off, format = "latex")


##                                                           ####
## Run numerical example 3 - ZT                              ####
## a) Set arguments                                          ####
target_feat = "g_DISRUPTIONS"

## zero truncation, no offset
offset_feat = NA
rate_param = FALSE
offset_feat_value = NA
zerotrunc = TRUE


## zero truncation and offset
#offset_feat = "d_TIME_offset"
#rate_param = TRUE
#zerotrunc = TRUE

## b1) full alternating procedure                            ####
test_NB_ZT <- NegBReg_altern(myData = myData, 
                             target_feat =  target_feat,
                             rate_param = rate_param, 
                             offset_feat = offset_feat,
                             zerotrunc = zerotrunc
)


test_NB_ZT$p_beta_MLE
test_NB_ZT$p_alpha_MLE


## b2) Significance                                          ####
## -- expected, alpha param
test_signif_C11_beta_exp <- myMLEsignif(my_SE = test_NB_ZT$SE_expected$SE_expeted_all_alpha$SE_expected_beta, 
                                        my_p_est= test_NB_ZT$p_beta_MLE)
test_signif_C21_alpha_exp <- myMLEsignif(my_SE = test_NB_ZT$SE_expected$SE_expeted_all_alpha$SE_expected_alpha, 
                                         my_p_est= test_NB_ZT$p_alpha_MLE)
tab_signif_C_exp <- rbind(test_signif_C11_beta_exp, test_signif_C21_alpha_exp)
rownames(tab_signif_C_exp)[nrow(tab_signif_C_exp)] <- "alpha"

## -- observed, alpha param
test_signif_C12_beta_obs <- myMLEsignif(my_SE = test_NB_ZT$SE_oberved$SE_obs_all_alpha$SE_obs_beta, 
                                        my_p_est= test_NB_ZT$p_beta_MLE)
test_signif_C22_alpha_obs <- myMLEsignif(my_SE = test_NB_ZT$SE_oberved$SE_obs_all_alpha$SE_obs_alpha, 
                                         my_p_est= test_NB_ZT$p_alpha_MLE)
test_signif_C_obs <- rbind(test_signif_C12_beta_obs, test_signif_C22_alpha_obs)
rownames(test_signif_C_obs)[nrow(test_signif_C_obs)] <- "alpha"
colnames(test_signif_C_obs) <- paste0(colnames(test_signif_C_obs),"_obs")

## -- assemble
tab_signif_C_alpha <- cbind(tab_signif_C_exp, test_signif_C_obs[,-1])

## -- expected, theta param
test_signif_C11_beta_exp_dp <- myMLEsignif(my_SE = test_NB_ZT$SE_expected$SE_expected_all_theta$SE_expected_beta, 
                                           my_p_est= test_NB_ZT$p_beta_MLE)
test_signif_C21_theta_exp_dp <- myMLEsignif(my_SE = test_NB_ZT$SE_expected$SE_expected_all_theta$SE_expected_theta, 
                                            my_p_est= 1/test_NB_ZT$p_alpha_MLE)
tab_signif_C_exp_dp <- rbind(test_signif_C11_beta_exp_dp, test_signif_C21_theta_exp_dp)
rownames(tab_signif_C_exp_dp)[nrow(tab_signif_C_exp_dp)] <- "theta"

## -- observed, theta param
test_signif_C12_beta_obs_dp <- myMLEsignif(my_SE = test_NB_ZT$SE_oberved$SE_obs_all_theta$SE_obs_beta, 
                                           my_p_est= test_NB_ZT$p_beta_MLE)

test_signif_C22_theta_obs_dp <- myMLEsignif(my_SE = test_NB_ZT$SE_oberved$SE_obs_all_theta$SE_obs_theta, 
                                            my_p_est= 1/test_NB_ZT$p_alpha_MLE)
tab_signif_C_obs_dp <- rbind(test_signif_C12_beta_obs_dp , test_signif_C22_theta_obs_dp )
rownames(tab_signif_C_obs_dp)[nrow(tab_signif_C_obs_dp)] <- "theta"
colnames(tab_signif_C_obs_dp) <- paste0(colnames(tab_signif_C_obs_dp),"_obs")

## -- assemble
tab_signif_C_dp <- cbind(tab_signif_C_exp_dp, tab_signif_C_obs_dp[,-1])
tab_signif_C_dp_just_theta <- tab_signif_C_dp[nrow(tab_signif_C_dp),]





## -- expected, log alpha
test_signif_C11_logalpha_exp_beta <- myMLEsignif(my_SE = test_NB_ZT$SE_expected$SE_expeted_all_logalpha$SE_expected_beta, 
                                                 my_p_est= test_NB_ZT$p_beta_MLE)
test_signif_C21_logalpha_exp <- myMLEsignif(my_SE = test_NB_ZT$SE_expected$SE_expeted_all_logalpha$SE_expected_logalpha, 
                                            my_p_est= log(test_NB_ZT$p_alpha_MLE))
tab_signif_C_logalpha_exp <- rbind(test_signif_C11_logalpha_exp_beta, test_signif_C21_logalpha_exp)
rownames(tab_signif_C_logalpha_exp)[nrow(tab_signif_C_logalpha_exp)] <- "logalpha"


## -- observed, log alpha
test_signif_C12_logalpha_obs_beta <- myMLEsignif(my_SE = test_NB_ZT$SE_oberved$SE_obs_all_logalpha$SE_obs_beta, 
                                                 my_p_est= test_NB_ZT$p_beta_MLE)
test_signif_C22_logalpha_obs <- myMLEsignif(my_SE = test_NB_ZT$SE_oberved$SE_obs_all_logalpha$SE_obs_logalpha, 
                                            my_p_est= log(test_NB_ZT$p_alpha_MLE))

test_signif_C_logalpha_obs <- rbind(test_signif_C12_logalpha_obs_beta, test_signif_C22_logalpha_obs)
rownames(test_signif_C_logalpha_obs)[nrow(test_signif_C_logalpha_obs)] <- "logalpha"
colnames(test_signif_C_logalpha_obs) <- paste0(colnames(test_signif_C_logalpha_obs),"_obs")

## -- assemble
tab_signif_C_logalpha <- cbind(tab_signif_C_logalpha_exp, test_signif_C_logalpha_obs[,-1])
tab_signif_C_logalpha_just_alpha <- tab_signif_C_logalpha[nrow(tab_signif_C_logalpha), ]

## -- Expected, log theta
test_signif_C11_logtheta_exp_beta <- myMLEsignif(my_SE = test_NB_ZT$SE_expected$SE_expected_all_logtheta$SE_expected_beta, 
                                                 my_p_est= test_NB_ZT$p_beta_MLE)
test_signif_C21_logtheta_exp <- myMLEsignif(my_SE = test_NB_ZT$SE_expected$SE_expected_all_logtheta$SE_expected_logtheta, 
                                            my_p_est= log(1/test_NB_ZT$p_alpha_MLE))

tab_signif_C_logtheta_exp <- rbind(test_signif_C11_logtheta_exp_beta , test_signif_C21_logtheta_exp)
rownames(tab_signif_C_logtheta_exp )[nrow(tab_signif_C_logtheta_exp)] <- "logtheta"

## -- observed log theta
test_signif_C12_logtheta_obs_beta <- myMLEsignif(my_SE = test_NB_ZT$SE_oberved$SE_obs_all_logtheta$SE_obs_beta, 
                                                 my_p_est= test_NB_ZT$p_beta_MLE)
test_signif_C22_logtheta_obs <- myMLEsignif(my_SE = test_NB_ZT$SE_oberved$SE_obs_all_logtheta$SE_obs_logtheta, 
                                            my_p_est= log(1/test_NB_ZT$p_alpha_MLE))

test_signif_C_logtheta_obs <- rbind(test_signif_C12_logtheta_obs_beta, test_signif_C22_logtheta_obs)
rownames(test_signif_C_logtheta_obs)[nrow(test_signif_C_logtheta_obs)] <- "logtheta"
colnames(test_signif_C_logtheta_obs) <- paste0(colnames(test_signif_C_logtheta_obs),"_obs")


## -- assemble
tab_signif_C_logtheta <- cbind(tab_signif_C_logtheta_exp, test_signif_C_logtheta_obs[,-1])
tab_signif_C_logtheta_just_theta <- tab_signif_C_logtheta[nrow(tab_signif_C_logtheta), ]




## c) Compare with pre-built                                 ####
## --- countreg                                ####

obj_NB_cr_ZT <- countreg::zerotrunc(g_DISRUPTIONS ~.,
                                    dist = "negbin",
                                    data= myData)


obj_NB_cr_ZT $coefficients
obj_NB_cr_ZT $coefficients.theta                        # this is ln alpha

summary(obj_NB_cr_ZT )[1]
summary(obj_NB_cr_ZT )[2]

out_countreg_zt <- as.data.frame(summary(obj_NB_cr_ZT )[1])

countreg_C_fulltab_beta <- out_countreg_zt[1:(nrow(out_countreg_zt)-1), 1:2, drop = FALSE]
countreg_C_fulltab_logalpha <- out_countreg_zt[nrow(out_countreg_zt), 1:2, drop = FALSE]

countreg_C_beta <-  countreg_C_fulltab_beta[,1,drop = FALSE]
countreg_C_beta_SE <-  countreg_C_fulltab_beta[,2,drop = FALSE]
countreg_C_logalpha <- countreg_C_fulltab_logalpha$coefficients.theta.Estimate
countreg_C_logalpha_SE <- countreg_C_fulltab_logalpha$coefficients.theta.Std..Error


## --- GAMLSS                                  ####

## As in Hilbe (2011) Ch. 11
library(gamlss.tr)
gen.trun(0, "NBI", type="left", name = "lefttr")
obj_NB_gam_ZT <- gamlss(g_DISRUPTIONS ~.,
                        data = myData,
                        family = "NBIlefttr") 

obj_NB_gam_ZT$mu.coefficients
exp(obj_NB_gam_ZT$sigma.coefficients)                                      # dispersion parameter
alpha_gamlss_ZT <- 1/exp(obj_NB_gam_ZT$sigma.coefficients)

summary(obj_NB_gam_ZT)

## extract estimators
gamlss_C_beta <- obj_NB_gam_ZT$mu.coefficients
gamlss_C_logtheta <- obj_NB_gam_ZT$sigma.coefficients                                 # "as is" output is ln theta
gamlss_C_theta_computed <- exp(obj_NB_gam_ZT$sigma.coefficients)
gamlss_C_alpha_computed <- 1/gamlss_C_theta_computed 


## extract standard errors (workaround)
var_gamlss_ZT <- vcov(obj_NB_gam_ZT)
var_gamlss_coeff_ZT <- diag(var_gamlss_ZT[1:(nrow(var_gamlss_ZT)-1), 1:(ncol(var_gamlss_ZT)-1)])
var_gamlss_logtheta_ZT <- diag(var_gamlss_ZT)[nrow(var_gamlss_ZT)]
gamlss_C_beta_SE <- sqrt(var_gamlss_coeff_ZT)
gamlss_C_logtheta_SE <- sqrt(var_gamlss_logtheta_ZT)








## --- VGAM                                    ####

## see Yee (2020) paper
library(VGAM)
obj_NB_vgam_ZT <- VGAM::vglm(g_DISRUPTIONS ~.,
                             family  = posnegbinomial,
                             data = myData)

vgam_beta_ZT <- obj_NB_vgam_ZT@coefficients[c(1,3:length(obj_NB_vgam_ZT@coefficients))]
vgam_alpha_ZT <- exp(obj_NB_vgam_ZT@coefficients[2])

coef(obj_NB_vgam_ZT, matrix = TRUE)

var_vgam_ZT <- vcov(obj_NB_vgam_ZT)
vgam_C_beta <- obj_NB_vgam_ZT@coefficients[c(1,3:length(obj_NB_vgam_ZT@coefficients))]
vgam_C_beta_sE <- round(sqrt(diag(var_vgam_ZT)[c(1,3:length(obj_NB_vgam_ZT@coefficients))]),5)

vgam_C_alpha <- exp(obj_NB_vgam_ZT@coefficients[2])
vgam_C_logalpha <- (obj_NB_vgam_ZT@coefficients[2])
vgam_C_logalpha_SE <- sqrt(diag(var_vgam_ZT)[2])   

## --- Table comparative                       ####

tab_fromScratch_ZT <- rbind(tab_signif_C_alpha, 
                            tab_signif_C_dp_just_theta,
                            tab_signif_C_logalpha_just_alpha,
                            tab_signif_C_logtheta_just_theta
)

tab_pre_built_MASS_ZT  <-  cbind.data.frame("n.a.", "n.a.")
tab_pre_built_countreg_ZT <- round(countreg_C_fulltab_beta[,c(1,2),drop = FALSE], 5)
tab_pre_built_GAMLSS_ZT <-  round(cbind.data.frame(gamlss_C_beta, gamlss_C_beta_SE),5)
tab_pre_built_VGAM_ZT <- round(cbind.data.frame(vgam_C_beta, vgam_C_beta_sE),5 )

tab_pre_built_alpha_MASS_ZT  <- cbind.data.frame("n.a.", "n.a.")
tab_pre_built_alpha_countreg_ZT <- cbind.data.frame("n.a.", "n.a.")
tab_pre_built_alpha_GAMLSS_ZT  <- cbind.data.frame("n.a.", "n.a.")
tab_pre_built_alpha_VGAM_ZT <- cbind.data.frame("n.a.", "n.a.")

tab_pre_built_theta_MASS_ZT <- cbind.data.frame("n.a.", "n.a.")
tab_pre_built_theta_countreg_ZT <- cbind.data.frame("n.a.", "n.a.")
tab_pre_built_theta_GAMLSS_ZT <- cbind.data.frame("n.a.", "n.a.")
tab_pre_built_theta_VGAM_ZT <-  cbind.data.frame("n.a.", "n.a.")

tab_pre_built_logalpha_MASS_ZT  <- cbind.data.frame("n.a.", "n.a.")
tab_pre_built_logalpha_countreg_ZT  <- round(countreg_C_fulltab_logalpha[,c(1:2)],5)
tab_pre_built_logalpha_GAMLSS_ZT  <- cbind.data.frame("n.a.", "n.a.")
tab_pre_built_logalpha_VGAM_ZT  <-  cbind.data.frame(round(vgam_C_logalpha,5), round(vgam_C_logalpha_SE,5))

tab_pre_built_logtheta_MASS_ZT  <- cbind.data.frame("n.a.", "n.a.")
tab_pre_built_logtheta_countreg_ZT  <- cbind.data.frame("n.a.", "n.a.")
tab_pre_built_logtheta_GAMLSS_ZT  <- cbind.data.frame(round(gamlss_C_logtheta,5), round(gamlss_C_logtheta_SE,5))
tab_pre_built_logtheta_VGAM_ZT  <-  cbind.data.frame("n.a.", "n.a.")

colnames(tab_pre_built_MASS_ZT) <- paste0(c("glmnb_"),colnames(tab_pre_built_MASS_ZT))
colnames(tab_pre_built_countreg_ZT) <- c("countreg_coeff", "countreg_SE")
colnames(tab_pre_built_GAMLSS_ZT) <- paste0(c("gamlss_coeff", "gamlss_SE"))
colnames(tab_pre_built_VGAM_ZT) <- c("vglm_coeff", "vglm_SE")

a <- cbind.data.frame( tab_pre_built_MASS_ZT, tab_pre_built_countreg_ZT, tab_pre_built_GAMLSS_ZT, tab_pre_built_VGAM_ZT)
b <- cbind.data.frame( tab_pre_built_alpha_MASS_ZT, tab_pre_built_alpha_countreg_ZT, tab_pre_built_alpha_GAMLSS_ZT, tab_pre_built_alpha_VGAM_ZT)
c <- cbind.data.frame( tab_pre_built_theta_MASS_ZT, tab_pre_built_theta_countreg_ZT, tab_pre_built_theta_GAMLSS_ZT, tab_pre_built_theta_VGAM_ZT)
d <- cbind.data.frame( tab_pre_built_logalpha_MASS_ZT, tab_pre_built_logalpha_countreg_ZT, tab_pre_built_logalpha_GAMLSS_ZT, tab_pre_built_logalpha_VGAM_ZT )
e <- cbind.data.frame( tab_pre_built_logtheta_MASS_ZT, tab_pre_built_logtheta_countreg_ZT, tab_pre_built_logtheta_GAMLSS_ZT, tab_pre_built_logtheta_VGAM_ZT)

colnames(e) <- colnames(d) <- colnames(c) <- colnames(b) <- colnames(a)
rownames(b) <- "alpha"
rownames(c) <- "theta"
rownames(d) <- "log alpha"
rownames(e) <- "log theta"

tab_pre_built_ZT <- rbind.data.frame(a, b, c, d, e)
myTab4viz_ZT <- cbind(tab_fromScratch_ZT, tab_pre_built_ZT)

## Generate LaTeX export
library(kableExtra)
kbl(myTab4viz_ZT, format = "latex")


##                                                           ####
## Additional metrics (discussion section)                   ####
## --- AIC                                         ####
ll_base <- test_NB$log_likelihood
ll_offs <- test_NB_offset$log_likelihood
ll_ZT <- test_NB_ZT$log_likelihood


AIC_base <-  as.numeric(-2*ll_base+2*(n_regressors))
AIC_offs <- as.numeric(-2*ll_offs+2*(n_regressors - 1))
AIC_ZT <- as.numeric(-2*ll_ZT+2*(n_regressors))

## --- IRR ZT                                      ####

p_lambda <- test_NB_ZT$p_lambda_MLE 
p_alpha <- test_NB_ZT$p_alpha_MLE
p_theta <- 1/p_alpha

if(!zerotrunc){
  something_ZT <- 1
} else {
  ## Grogger and Carson 1991 eq 21 give the derivative that usually underpins IRR for ZTNB but it contains some issues
  ## -- from my calculations (see ZTNB powerpoint deck):  d ln E[Y|Y>0] \ d x = beta * something.
  ## -- something = 1 - E[Y|Y>0]/(1+theta lambda)^(alpha + 1)
  ## -- so: exp(beta) = (d E[Y|Y>0] \ d x)  / something ...is what we call IRR in ZTNB 
  mu <- E_NB(p_lambda, p_alpha, rate_param = rate_param , offset_feat = offset_feat, offset_feat_value = offset_feat_value, zerotrunc = zerotrunc)
  a <- 1 + p_theta*p_lambda
  something_ZT_4_beta <- (1 - mu/a^(p_alpha + 1))                                     # d log EY dx = something_ZT * beta 
  something_ZT_4_expbeta <-  mu/a^(p_alpha + 1)                                         # d log EY dx = exp(beta) * exp (-something_ZT * beta) 
}

cbind.data.frame(myData$g_DISRUPTIONS, test_NB_ZT$p_mean, something_ZT_4_beta, something_ZT_4_expbeta)

IRR_BIO <- exp(tab_signif_C_exp[4,1])
perc_fewer_y <- 1-IRR_BIO
IRR_DUR <- exp(tab_signif_C_exp[2,1])
perc_more_y <-  IRR_DUR - 1
IRR_SUB <- exp(tab_signif_C_exp[3,1])

IRR_BIO_ZT_correction1 <- exp(tab_signif_C_exp[4,1]*something_ZT_4_beta)
IRR_BIO_ZT_correction2 <- IRR_BIO*exp(-something_ZT_4_expbeta*tab_signif_C_exp[4,1])
all.equal(IRR_BIO_ZT_correction1, IRR_BIO_ZT_correction2)

quantile(exp(tab_signif_C_exp[4,1]*something_ZT_4_beta))

IRR_DUR_ZT_correction1 <- exp(tab_signif_C_exp[2,1]*something_ZT_4_beta)
IRR_DUR_ZT_correction2 <- IRR_DUR*exp( (-something_ZT_4_expbeta)*tab_signif_C_exp[2,1])

quantile(exp(tab_signif_C_exp[2,1]*something_ZT_4_beta))

IRR_SUB_ZT_correction1 <- exp(tab_signif_C_exp[3,1]*something_ZT_4_beta)
IRR_SUB_ZT_correction2 <- IRR_SUB*exp( (-something_ZT_4_expbeta)*tab_signif_C_exp[3,1])
quantile(exp( (-something_ZT_4_expbeta)*tab_signif_C_exp[3,1]))
