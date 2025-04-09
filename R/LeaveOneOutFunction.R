###Leave-one-out implementation for any ICAD method###

##source ICAD functions##
source('R/EFDM_Functions.R')
source('R/SNCM_Functions.R')
source('R/GMD_Functions.R')
source('R/DataGenerationFunctions.R')


###Function for computing calibration scores using EFDM###
###Input:
# dat:           Data set as n x n_t matrix
# t_vec:         Time vector
# method:        ICAD method, one of "EFDM", "SNCM", "GMD", or "GMDM"
# kmean:         If True, compute Karcher mean (use True for EFDM and potentially SNCM)
# efdm_wgt:      Weight for amplitude distance in NCM for efdm_cp
# efdm_iv_dist:  If True, use initial value distance in NCM in efdm_cp
# efdm_iv_wgt:   Weight vector to apply to each component in NCM in the form (amplitude, phase, initial_value) for efdm_cp
# sncm_mean:     Input for mean_fn in sncm_cp
# sncm_mod:      Input for mod_fn in sncm_cp
# sncm_alpha:    Input for alpha in sncm_cp
# gmd_p:         Input for p in functionalCP and functionalCP_m
# gmd_K:         Input for K in functionalCP and functionalCP_m

####################
###Output:
# pvalues for every observation in dat

loo_fn<-function(dat, t_vec, method, kmean = F, efdm_wgt = NULL, efdm_iv_dist = F, efdm_iv_wgt = NULL, sncm_mean = NULL,
                 sncm_mod = NULL, sncm_alpha = NULL, gmd_p = NULL, gmd_K = NULL){
  n<-dim(dat)[1]
  pval_vec<-rep(NA, times = n)
  
  for(i in 1:n){
    dat_use<-dat[-i,]
    ts_obs<-t(dat[i,])
    temp_dat<-tr_cal_split(dat_use)
    tr_dat<-temp_dat$tr_dat
    cal_dat<-temp_dat$cal_dat
    
    if(kmean){
      tw<-time_warping(t(tr_dat),t_vec)
      k_mean<-tw$fmean
    }
    else{
      k_mean<-NULL
    }
    
    if(method=='EFDM'){
      temp_cp<-efdm_cp(tr_dat, cal_dat, t_vec, efdm_wgt, k_mean, efdm_iv_dist, efdm_iv_wgt)
      pval_vec[i]<-efdm_cp_pred(ts_obs, temp_cp, t_vec)
    }
    if(method=='SNCM'){
      temp_cp<-sncm_cp(tr_dat, cal_dat, t_vec, sncm_mean, sncm_mod, sncm_alpha, k_mean)
      pval_vec[i]<-sncm_cp_pred(ts_obs, temp_cp, t_vec)$TestPValue
    }
    if(method=='GMD'){
      temp_cp<-functionalCP(t_vec, t(tr_dat), t(cal_dat), gmd_p, gmd_K)
      pval_vec[i]<-functionalCPpvals(Gmix = temp_cp, test = t(ts_obs))
    }
    if(method=='GMDM'){
      temp_cp<-functionalCP_m(t_vec, t(tr_dat), t(cal_dat), gmd_p, gmd_K)
      pval_vec[i]<-functionalCPpvals_m(Gmix = temp_cp, test = t(ts_obs))
    }
  }
  return(pval_vec)
}






