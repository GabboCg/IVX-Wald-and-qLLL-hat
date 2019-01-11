if(!require("pacman")) install.packages("pacman")
p_load("tidyverse", "xlsx")

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

beta_hat <- array(NaN, c(NCOL(GW_predictor) + 1, 4, length(h)))

for(i in seq_along(h)){
  for(j in 1:(ncol(GW_predictor)+2)){
   
    if(j <= ncol(GW_predictor)){
      
      X_i_j <- GW_predictor_z[1:(NROW(GW_predictor_z)-h[i]), j]
      results_i_j <- lm(100*r_h[2:(NROW(r_h)-(h[i]-1)), i] ~ X_i_j)
      beta_hat[j,,i] <- cbind(results_i_j$coefficients[2],
                              coef(summary(results_i_j))[, "t value"][2],
                              NaN, 
                              summary(results_i_j)$r.squared*100)
      
    } 
    
  }
}




