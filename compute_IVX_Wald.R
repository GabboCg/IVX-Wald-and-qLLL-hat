# Computes the Kostakis, Magdalinos, and Stamatogiannis (2015) IVX-Wald
# statistic for a predictive regression.
# 
# Input
# 
# y    = T-vector of return observations
# X    = T-by-r data matrix of predictor observations
# K    = forecast horizon
# M_n  = bandwith parameter
# beta = scalar in (0,1); specify 'close' to one for best performance
# 
# Output
# 
# A_tilde_IVX_K = r-vector of coefficient estimates
# W_IVX_K       = IVX-Wald statistic
# p_value       = p-value for IVX-Wald statistic
# 
# Reference
# 
# Kostakis, A, T Magdalinos, and MP Stamatogiannis (2015), "Robust
# Econometric Inference for Stock Return Predictability," Review of
# Financial Studies 28, 1506-1553

Compute_IVX_Wald <- function(y, X, K, M_n, beta){

  results_ols_y <- lm(y[2:NROW(y)] ~ X[1:(NROW(X)-1)])
  epsilon_hat <- results_ols_y$resid
  U_hat <- matrix(NaN, nrow = NROW(y)-1, ncol = NCOL(X))
  
  for(j in 1:NCOL(X)){
   
    results_OLS_j <- lm(X[2:NROW(X), j] ~ X[1:(NROW(X)-1), j])
    U_hat[,j] <- results_OLS_j$resid
    
  }
  
  Sigma_hat_ee <- (1/length(epsilon_hat))*(t(epsilon_hat)%*%epsilon_hat)
  Sigma_hat_eu <- (1/length(epsilon_hat))*(t(epsilon_hat)%*%U_hat)
  Sigma_hat_uu <- (1/length(epsilon_hat))*(t(U_hat)%*%U_hat)
  
  Omega_hat_uu <- Sigma_hat_uu
  Omega_hat_eu <- Sigma_hat_eu
  
  if(M_n > 0){
    
    Lambda_hat_uu <- matrix(0, nrow = NCOL(U_hat), ncol = NCOL(U_hat))
    Lambda_hat_ue <- matrix(0, nrow = NCOL(U_hat), ncol = 1)
    
    for(h in 1:M_n){
      
      Lambda_hat_uu <- Lambda_hat_uu + (1-(h/(M_n+1)))*(1/(length(U_hat)-1))*t(U_hat[(h+1):NROW(U_hat)])%*%U_hat[1:(NROW(U_hat)-h)]
      Lambda_hat_ue <- Lambda_hat_ue + (1-(h/(M_n+1)))*(1/(length(U_hat)-1))*t(U_hat[(h+1):NROW(U_hat)])%*%epsilon_hat[1:(NROW(epsilon_hat)-h)]
      
    }
    
    Omega_hat_uu <- Sigma_hat_uu + Lambda_hat_uu + t(Lambda_hat_uu)
    Omega_hat_eu <- Sigma_hat_eu + t(Lambda_hat_ue)
    
  }
  
  R_nz <- 1 - (1/(length(y)-1)^beta)
  d_X <- X[2:NROW(X)] - X[1:(NROW(X)-1)]
  d_X <- rbind(matrix(0, nrow = 1 ,ncol = NCOL(X)), as.matrix(d_X))
  Z_tilde <- matrix(NaN, nrow = NROW(X), ncol = NCOL(X)) # instrument matrix
  Z_tilde[1] <- matrix(0, nrow = 1, ncol = NCOL(X))
  
  for(t in 2:length(X)){
    
    Z_tilde[t] <- R_nz%*%Z_tilde[(t-1)] + d_X[t]
    
  }
  
  y_K <- matrix(NaN, nrow = (NROW(y)-K-1), ncol = 1) 
  X_K <- matrix(NaN, nrow = (NROW(y)-K-1), ncol = NCOL(X)) 
  Z_tilde_K <- matrix(NaN, nrow = (NROW(y)-K-1), ncol = 1) 
  
  for(t in 1:(length(y)-(K-1))){
    
    y_K[t] <- sum(y[t:(t+(K-1))])
    X_K[t] <- sum(X[t:(t+(K-1))])
    Z_tilde_K[t] <- sum(Z_tilde[t:(t+(K-1))])
    
  }
  
  n_K <- length(y_K) - 1
  y_bar_K <- mean(y_K[2:NROW(y_K)])
  x_bar_K <- mean(X_K[1:(NROW(X_K)-1)])
  
  z_tilde_bar_K <- mean(t(Z_tilde_K[1:(NROW(Z_tilde_K)-K)]))
  
  Y_K_under <- y_K[2:NROW(y_K)] - y_bar_K
  X_K_under <- X_K[1:(NROW(X_K)-1)] - kronecker(t(x_bar_K), matrix(1, nrow = n_K, ncol = 1))
  
  Z_tilde_K <- Z_tilde_K[1:(NROW(Z_tilde_K)-1)]
  Z_tilde <- Z_tilde[1:(NROW(Z_tilde)-K)]
  
  A_tilde_IVX_K <- solve(t(X_K_under)%*%Z_tilde)%*%t(Z_tilde)%*%Y_K_under
  
  Omega_hat_FM <- Sigma_hat_ee - Omega_hat_eu*solve(Omega_hat_uu)*t(Omega_hat_eu)
  M_K <- (t(Z_tilde_K)%*%Z_tilde_K)*Sigma_hat_ee - n_K*(z_tilde_bar_K*t(z_tilde_bar_K))*Omega_hat_FM
  Q_H_K <- solve(t(Z_tilde)%*%X_K_under)%*%M_K%*%solve(t(X_K_under)%*%Z_tilde)
  W_IVX_K <- t(A_tilde_IVX_K)%*%solve(Q_H_K)%*%A_tilde_IVX_K
  
  p_value <- 1 - pchisq(W_IVX_K, NCOL(X), lower.tail = TRUE)
  
  results_all <<- list(A_tilde_IVX_K,
                       W_IVX_K,
                       p_value)
  
}