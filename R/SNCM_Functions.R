#####SNCM Functions#####
library(fdasrvf)
library(caTools)



###Function to produce Diquigiovanni's CP bands###
###Input:
# tr_dat:      Training data as n_tr x nt matrix
# cal_dat:     Calibration data as n_cal x nt matrix
# t_vec:       Time vector
# mean_fn:     string, either "csm" or "Karcher"
# mod_fn:      string, one of "none", "sd", or "optimal"
# alpha:       Significance level
# k_mean:      Karcher mean
####################
###Output:
# Mean function
# Lower bound at each time point
# Upper bound at each time point
# Modulation function
# Calibration scores
sncm_cp<-function(tr_dat, cal_dat, t_vec, mean_fn, mod_fn, alpha, k_mean){
  
  nt<-length(t_vec)
  n_tr<-dim(tr_dat)[1]
  n_cal<-dim(cal_dat)[1]
  k_pct<-ceiling((n_cal+1)*(1-alpha))
  
  ###Compute mean function###
  if(mean_fn=='csm'){
    g_1<-apply(tr_dat,2, mean)
  }
  if(mean_fn=='Karcher'){
    g_1<-k_mean
  }
  
  ###Compute modulation function###
  if(mod_fn=="none"){
    s_1<-rep(1,times = nt)
  }
  if(mod_fn=='sd'){
    s_1<-apply(tr_dat,2,sd)
  }
  if(mod_fn=="optimal"){
    gam_set<-c()
    for(i in 1:n_tr){gam_set[i]<-max(tr_dat[i,]-g_1)}
    gam_p<-ceiling((n_tr+1)*(1-alpha))
    gam<-sort(gam_set)[gam_p]
    H1_idx<-which(gam_set<gam)
    
    if(gam_p>n_tr){
      H_1<-tr_dat
    }else
    {
      H_1<-tr_dat[H1_idx,]  
    }
    s1_num<-c()
    for(i in 1:nt){
      s1_num[i]<-max(abs(H_1[,i]-g_1[i]))
    }
    s1_denom<-trapz(t_vec,s1_num)
    s_1<-s1_num/s1_denom
  }
  
  cal_scores<-rep(NA,times = n_cal)
  for(i in 1:n_cal){
    cal_scores[i]<-max(abs((cal_dat[i,]-g_1)/s_1))
  }
  
  k_s<-sort(cal_scores)[k_pct]
  upr_bnd<-g_1 + k_s*s_1
  lwr_bnd<-g_1 - k_s*s_1
  
  out_list<-list(Mean_fn = g_1, Mod_fn = s_1, Lower_bnd = lwr_bnd, Upper_bnd = upr_bnd, Cal_scores = cal_scores, K_s = k_s, SigLevel = alpha)
  return(out_list)
}

###Function that outputs test score and p-value for a set of test functions using output of sncm_cp###

###Input:
# test_dat:    Test data as n_ts x nt matrix
# cal_out:     Output list from sncm_cp()
# t_vec:       Time vector
####################
###Output:
# vector of p-values
# coverage vector

sncm_cp_pred<-function(ts_dat, cal_out, t_vec){
  
  g_1<-cal_out$Mean_fn
  s_1<-cal_out$Mod_fn
  cal_scores<-cal_out$Cal_scores
  alpha<-cal_out$SigLevel
  
  nt<-length(t_vec)
  n_cal<-length(cal_scores)
  n_ts<-dim(ts_dat)[1]
  
  ts_scores<-rep(NA,n_ts)
  for(i in 1:n_ts){
    ts_scores[i]<-max(abs((ts_dat[i,]-g_1)/s_1))
  }
  
  p_val<-rep(NA, times = n_ts)
  for(i in 1:n_ts){
    count<-(length(which(cal_scores>ts_scores[i])) + runif(1)*(1+length(which(cal_scores==ts_scores[i])))) 
    p_val[i]<-(count)/(n_cal+1)
  }
  covg_vec<-rep(1, times = n_ts)
  covg_vec[p_val<alpha]<-0
  
  ret_list<-list(TestPValue = p_val, TestCoverage = covg_vec)
  return(ret_list)
  
}