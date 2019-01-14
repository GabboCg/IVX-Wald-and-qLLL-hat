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

IVX_Wald <- function(y, X, K, M_n, beta){
  

  
}