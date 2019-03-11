if(!require("pacman")) install.packages("pacman")
p_load("tidyverse", "xlsx", "sandwich", "lmtest")

library("lpirfs")

GW_dataset <- read.xlsx("data/Returns_short_interest_data.xlsx", 
                        sheetName = "GW variables") 

GW_varset <- as.data.frame(apply(GW_dataset, 2, as.numeric)) %>% 
             filter(yyyymm >= 197212) %>% 
             select(yyyymm, Index, D12, E12, b.m, ntis, tbl, lty, ltr, BAA, 
                    AAA, corpr, infl) %>%
             rename(SP = Index, D12 = D12, E12 = E12, BM = b.m, NTIS = ntis,
                    TBL = tbl, LTY = lty, LTR = ltr, BAA = BAA, AAA = AAA,
                    CORPR = corpr, INFL_lag = infl) %>%  
             mutate(log_DP = log(D12/SP), log_EP = log(E12/SP), 
                    log_DE = log(D12/E12), log_DY = log(D12/lag(SP)),
                    TMS = LTY - TBL, DFY = BAA - AAA, DFR = CORPR - LTR,
                    INFL_lag = lag(INFL_lag)) %>% 
             select(log_DP, log_DY, log_EP, log_DE, BM, NTIS, TBL, LTY, 
                    LTR, TMS, DFY, DFR, INFL_lag) %>% 
             na.omit()

# Stock excess return volatility (annualized) ----------------------------------
RVOL <- as.data.frame(apply(GW_dataset, 2, as.numeric)) %>% 
        filter(yyyymm >= 197201) %>% 
        select(yyyymm, CRSP_SPvw, Rfree) %>% 
        mutate(R_F_lag = lag(Rfree)) %>% 
        rename(SP_R = CRSP_SPvw) %>% 
        select(yyyymm, SP_R, R_F_lag) %>% 
        na.omit() 

RVOL_mat <- matrix(NaN, nrow = nrow(RVOL) - 11, ncol = 1)

for(i in seq_along(RVOL_mat)){
  
  RVOL_mat[i,] <- mean(abs(RVOL$SP_R[i:(i+11)] - RVOL$R_F_lag[i:(i+11)]))
  
}

RVOL_mat <- sqrt(pi/2)*sqrt(12)*RVOL_mat

# merge variables --------------------------------------------------------------
GW_predictor <- cbind(GW_varset, RVOL_mat) %>% 
                rename(RVOL = RVOL_mat) %>% 
                select(log_DP, log_DY, log_EP, log_DE, RVOL, BM, NTIS, TBL, LTY, 
                       LTR, TMS, DFY, DFR, INFL_lag)

# Load equity risk premium data, 1973:01-2014:12 -------------------------------
equity_risk <- data.frame(apply(GW_dataset, 2, as.numeric)) %>% 
               filter(yyyymm >= 197212) %>% 
               select(yyyymm, Rfree, CRSP_SPvw) %>%
               rename(R_SP500 = CRSP_SPvw) %>% 
               mutate(ER = R_SP500 - lag(Rfree), Rfree_lag = lag(Rfree),
                      r = log(1 + R_SP500) - log(1 + Rfree_lag)) %>% 
               na.omit()

# short-interest rate ----------------------------------------------------------
short_interest <- read.xlsx("data/Returns_short_interest_data.xlsx", 
                            sheetName = "Short interest") %>% 
                  select(Date, EWSI) %>% 
                  mutate(log_EWSI = log(EWSI)) %>% 
                  select(Date, log_EWSI)

# Compute log(EWSI) deviation from linear trend --------------------------------
X_linear <- apply(as.matrix(rownames(short_interest)), 2, as.numeric)
result_linear <- lm(short_interest$log_EWSI ~ X_linear)
SII <- scale(result_linear$residual)


# Compute cumulative returns ---------------------------------------------------
h <- c(1, 3, 6, 12)

r_h <- matrix(NaN, nrow = nrow(equity_risk), ncol = NROW(h))

for(i in seq_along(h)){
  for(j in 1:(nrow(equity_risk) - (h[i] - 1))){
  
    r_h[j,i] <- mean(equity_risk$r[j:(j+h[i]-1)])
    
  }
}

# scale predictors -------------------------------------------------------------
GW_predictor_z <- apply(GW_predictor, 2, scale) 

beta_hat <- array(NaN, c(NCOL(GW_predictor), 4, length(h)))

for(i in seq_along(h)){
  for(j in 1:(ncol(GW_predictor))){
   
    if(j <= ncol(GW_predictor)){
      
      X_i_j <- as.matrix(GW_predictor_z[1:(NROW(GW_predictor_z)-h[i]), j])
      Y_i_j <- 100*r_h[2:(NROW(r_h)-(h[i]-1)), i]
      # results_i_j <- newey_west(Y_i_j , X_i_j, h[i])
      
      results_i_j_lm <- lm(100*r_h[2:(NROW(r_h)-(h[i]-1)), i] ~ X_i_j)
      results_i_j_lm_coef <- coeftest(results_i_j_lm, vcov.=NeweyWest(results_i_j_lm, lag=h[i], adjust=FALSE, verbose = TRUE, prewhite = FALSE))
      
      beta_hat[j,,i] <- cbind(results_i_j_lm_coef[2,1],
                              results_i_j_lm_coef[2,3], # coef(summary(results_i_j))[, "t value"][2]
                              NaN, 
                              summary(results_i_j_lm)$r.squared*100)
      
    } 
  }
}

round(beta_hat, 2)

IVX_Wald <- matrix(NaN, nrow = 2, ncol = NROW(h))
beta_hat_SII <- array(NaN, c(1, 4, length(h)))

# IVX Wald and qLL -------------------------------------------------------------

source("Compute_IVX_Wald.R")

for(i in seq_along(h)){
  
  X_i_j <- as.matrix(-SII[1:(NROW(SII)-h[i])])
  Y_i_j <- 100*r_h[2:(NROW(r_h)-(h[i]-1)), i]
  
  results_i_j_lm <- lm(100*r_h[2:(NROW(r_h)-(h[i]-1)), i] ~ X_i_j)
  results_i_j_lm_coef <- coeftest(results_i_j_lm, vcov.=NeweyWest(results_i_j_lm, lag=h[i], adjust=FALSE, verbose = TRUE, prewhite = FALSE))
  
  beta_hat_SII[,,i] <- cbind(results_i_j_lm_coef[2,1],
                          results_i_j_lm_coef[2,3], 
                          NaN, 
                          summary(results_i_j_lm)$r.squared*100)
  
  IVX_results <- Compute_IVX_Wald(equity_risk$r, -SII, h[i], 0, 0.99)
  
  IVX_Wald[1,i] <- IVX_results[[2]]
  IVX_Wald[2,i] <- IVX_results[[3]]
  
}

round(IVX_Wald, 2)

# Compute fixed-regressor wild bootstrap p-values ------------------------------
X_sink <- data.frame(cbind(GW_predictor_z, -SII)) %>% 
          select(-log_DE, -TMS)

results_sink <- lm(r_h[2:NROW(r_h),1] ~ as.matrix(X_sink[1:(NROW(X_sink)-1),]))
epsilon_hat <- as.matrix(results_sink$residuals)
B <- 1000
beta_hat_tstat_star <- array(NaN, c(B, NCOL(GW_predictor) + 1, length(h)))  

saved.seed <- sample(1e9,1)
set.seed(saved.seed)

for(b in 1:B){
  
  print(b)
  
  u_star_b <- as.numeric(rnorm(NROW(r_h)-1, 0,1))
  r_star_b <- rbind(equity_risk$r[1], mean(equity_risk$r) + epsilon_hat*u_star_b)
  r_h_star_b <- matrix(NaN, NROW(equity_risk$r), length(h)) 
  
  for(j in 1:length(h)){
    for(t in 1:(length(equity_risk$r) - (h[j] - 1))){
     
      r_h_star_b[t,j] <- mean(r_star_b[t:(t+(h[j]-1))]) 
      
    }
  }
  
  
  for(j in 1:length(h)){
    for(i in 1:(NCOL(GW_predictor))){
     
      if(i <= NCOL(GW_predictor)){
       
        X_i_j <- as.matrix(GW_predictor_z[1:(NROW(GW_predictor_z)-1-(h[j]-1)), i])
        Y_i_j <- as.vector(100*r_h_star_b[2:(NROW(r_h_star_b)-(h[j]-1)), j])
        
        results_i_j_star_b_lm <- lm(Y_i_j ~ X_i_j)
        results_i_j_star_b <- coeftest(results_i_j_star_b_lm, vcov.=NeweyWest(results_i_j_star_b_lm, lag = h[j], adjust=FALSE, verbose = TRUE, prewhite = FALSE))
    
        beta_hat_tstat_star[b,i,j] <- results_i_j_star_b[2,3] # results_i_j_star_b[[1]][2]/(sqrt(results_i_j_star_b[[2]][2,2])) or coef(summary(results_i_j))[, "t value"][2]
        
      }
    }
  }
}

for(j in 1:length(h)){
  for(i in 1:(NCOL(GW_predictor))){
   
    beta_hat[i, 3, j] <- sum(beta_hat_tstat_star[,i,j] > beta_hat[i, 2, j])/B
     
  }
}

colnames(beta_hat) <- c("Beta", "t_stat", "p_value", "R2")
rownames(beta_hat) <- c("DP","DY","EP","DE","RVOL","BM","NTIS",
                        "TBL","LTY","LTR","TMS","DFY","DFR","INFL")

round(beta_hat, 2)

# Take care of out-of-sample preliminaries -------------------------------------

T <- length(equity_risk$r)
in_sample_end <- 1989
R <- (in_sample_end - 1972)*12 # in-sample period
P <- T - R

FC_PM <- matrix(NaN, nrow = P, ncol = 1)
FC_PR <- array(NaN, c(P, NCOL(GW_predictor), NROW(h)))

# Compute out-of-sample forecasts ----------------------------------------------

for(p in 1:P){
  
  print(p)
  
  FC_PM[p] <- mean(equity_risk$r[1:(R+(p-1))])
  
  for(j in 1:NROW(h)){
    for(i in 1:NCOL(GW_predictor)){
      
      X_i_j_p <- GW_predictor[1:(R+(p-1)-h[j]),i]
      results_i_j_p <- lm(r_h[2:(R+p-h[j]),j] ~ X_i_j_p)
      FC_PR[p,i,j] <- cbind(1, GW_predictor[R+(p-1),i])%*%results_i_j_p$coeff
      
    }
  }
}

# evaluate forecasts -----------------------------------------------------------

R2OS_PR <- array(NaN, c(NCOL(GW_predictor),2,NROW(h)))

for(j in 1:length(h)){
  
  actual_j <- r_h[(R+1):(NROW(r_h)-(h[j]-1)),j]
  u_PM_j <- actual_j - FC_PM[1:(NROW(FC_PM)-(h[j]-1))] 
  u_PR_j <- kronecker(matrix(1, nrow = 1, ncol = NCOL(FC_PR)), actual_j) - FC_PR[1:(NROW(FC_PR)-(h[j]-1)),,j] 
  
  MSFE_PM_j <- mean(u_PM_j^2)
  MSFE_PR_j <- apply(u_PR_j^2, 2, mean)
  
  R2OS_PR_j <- 100*(1-MSFE_PR_j/MSFE_PM_j)
  R2OS_PR[,1,j] = t(R2OS_PR_j)
  
  for(i in 1:NCOL(GW_predictor)){
    
   f_CW_i_j <- u_PM_j^2 - u_PR_j[,i]^2 + (FC_PM[1:(NROW(FC_PM)-(h[j]-1))] - FC_PR[1:(NROW(FC_PR)-(h[j]-1)),i,j])^2
   result_CW_i_j <- lm(f_CW_i_j ~ matrix(1, nrow = NROW(f_CW_i_j)))
   R2OS_PR[i,2,j] <- coeftest(result_CW_i_j, vcov.=NeweyWest(result_CW_i_j, lag=h[j], adjust=FALSE, verbose = TRUE, prewhite = FALSE))[3]
   
  }
}

round(R2OS_PR, 2)
