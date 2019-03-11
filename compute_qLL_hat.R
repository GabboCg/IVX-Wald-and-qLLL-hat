# Last modified: 08-18-2015#
#
# Computes the Elliott and Muller (2006) statistic for testing the null
# that beta(t) = beta for all t.
#
# Input
#
# y = T-vector of dependent variable observations
# X = T-by-k data matrix linked with potentially changing coefficients
# Z = T-by-d data matrix linked with constant coefficients (empty matrix
#                                                           %     if all coefficients allowed to change)
# L = lag truncation for Newey-West estimator
#
# Output
#
# qLL_hat = Elliott and Muller (2006) statistic (see Table 1 of paper for
#                                                %           critical values)
#
# Reference
#
# Elliott, G and UK Muller, (2006), "Efficient Tests for General Persistent
# Time Variation in Regression Coefficients," Review of Economic Studies
# 73, 907-940

Compute_qLL_hat <- function(y, X, Z, L){
  
  results <- lm(y ~ X + Z - 1)
  epsilon_hat <- results$resid
  
  X_epsilon_hat <- X*kronecker(matrix(1, nrow = 1, ncol = NCOL(X)), epsilon_hat)
  V_hat_X <- (1/NROW(X))*(t(X_epsilon_hat)%*%X_epsilon_hat)
  
  if(L >= 1){
    for(l in 1:L){
     
      Gamma_hat_l <- (1/NROW(X))*t(X_epsilon_hat[(l+1):(NROW(X_epsilon_hat))])%*%X_epsilon_hat[1:(NROW(X_epsilon_hat)-l)]
      V_hat_X <- V_hat_X + (1-(l/(L+1)))*(Gamma_hat_l + t(Gamma_hat_l))
      
    }
  }
  
  U_hat <- (V_hat_X^(-0.5))%*%t(X_epsilon_hat)
  
  r_bar <- 1 - (10/NROW(X))
  w_hat <- matrix(0, nrow = NCOL(X), ncol = NROW(X))
  
  for(t in 1:NROW(X)){
    if(t==1){
     
      w_hat[,t] <- U_hat[,t]
      
    }else{
      
      w_hat[,t] <- r_bar%*%w_hat[,(t-1)] + (U_hat[,t] - U_hat[,(t-1)])
      
    }
  }
  
  r_bar_trend <- matrix(0, nrow = NROW(y), ncol = 1)
  
  for(ii in 1:NROW(r_bar_trend)){
    
    r_bar_trend[ii] <- r_bar^(ii) 
    
  }
  
  SSR <- matrix(NaN, nrow = NCOL(X), ncol = 1)
  
  for(i in 1:NCOL(X)){
    
    results_i <- lm(t(w_hat)[,i] ~ r_bar_trend - 1)
    SSR[i] <- t(results_i$resid)%*%results_i$resid
    
  }
  
  qLL_hat <- r_bar*sum(SSR) - sum(apply(U_hat^2, 2, sum))
  
  resutl_all <- list(qLL_hat)
  
}