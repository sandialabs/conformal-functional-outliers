###FDA Outlier Detection Experiment 3###
#This code runs the third simulation experiment#

##source ICAD functions##
source('R/EFDM_Functions.R')
source('R/SNCM_Functions.R')
source('R/GMD_Functions.R')
source('R/DataGenerationFunctions.R')

##define t_vec and alpha##
nt<-250
t_vec<-seq(-8, 8, length.out = nt)
alpha<-0.10

##read in data##
exp4_dat<-readRDS('Data/SimulatedData/exp4_dat.rds')
n<-dim(exp4_dat)[1]

##placeholder matrices for pvalues##
efdm_pval_mat<-matrix(NA, n, dim(exp4_dat)[3])
sncm1_pval_mat<-matrix(NA, n, dim(exp4_dat)[3])
sncm2_pval_mat<-matrix(NA, n, dim(exp4_dat)[3])
gmd1_pval_mat<-matrix(NA, n, dim(exp4_dat)[3])
gmdm1_pval_mat<-matrix(NA, n, dim(exp4_dat)[3])
gmd2_pval_mat<-matrix(NA, n, dim(exp4_dat)[3])
gmdm2_pval_mat<-matrix(NA, n, dim(exp4_dat)[3])
gmd3_pval_mat<-matrix(NA, n, dim(exp4_dat)[3])
gmdm3_pval_mat<-matrix(NA, n, dim(exp4_dat)[3])

for(j in 1:dim(exp4_dat)[3]){
  
  dat<-exp4_dat[,,j]
  
  ##placeholder vectors for pvalues for each data set##
  efdm_pval_vec<-rep(NA, times = n)
  sncm1_pval_vec<-rep(NA, times = n)
  sncm2_pval_vec<-rep(NA, times = n)
  gmd1_pval_vec<-rep(NA, times = n)
  gmdm1_pval_vec<-rep(NA, times = n)
  gmd2_pval_vec<-rep(NA, times = n)
  gmdm2_pval_vec<-rep(NA, times = n)
  gmd3_pval_vec<-rep(NA, times = n)
  gmdm3_pval_vec<-rep(NA, times = n)
  
  for(i in 1:n){
    dat_use<-dat[-i,]
    ts_obs<-t(dat[i,])
    temp_dat<-tr_cal_split(dat_use)
    tr_dat<-temp_dat$tr_dat
    cal_dat<-temp_dat$cal_dat
    
    tw<-time_warping(t(tr_dat),t_vec)
    k_mean<-tw$fmean
    
    #EFDM#
    efdm_tr<-efdm_cp(tr_dat, cal_dat, t_vec, wgt = 0.5, k_mean)
    efdm_pval_vec[i]<-efdm_cp_pred(ts_obs, efdm_tr, t_vec)
    
    #SNCM1#
    sncm1_tr<-sncm_cp(tr_dat, cal_dat, t_vec, mean_fn = 'csm', mod_fn = 'optimal', alpha)
    sncm1_pval_vec[i]<-sncm_cp_pred(ts_obs, sncm1_tr, t_vec)$TestPValue
    
    #SNCM2#
    sncm2_tr<-sncm_cp(tr_dat, cal_dat, t_vec, mean_fn = 'Karcher', mod_fn = 'optimal', alpha, k_mean)
    sncm2_pval_vec[i]<-sncm_cp_pred(ts_obs, sncm1_tr, t_vec)$TestPValue
    
    #GMD1#  
    gmd1_tr<-functionalCP(t_vec, t(tr_dat), t(cal_dat), 10, 1)
    gmd1_pval_vec[i]<-functionalCPpvals(Gmix = gmd1_tr, test = t(ts_obs))
    
    #GMDM1#  
    gmdm1_tr<-functionalCP_m(t_vec, t(tr_dat), t(cal_dat), 10, 1)
    gmdm1_pval_vec[i]<-functionalCPpvals_m(Gmix = gmdm1_tr, test = t(ts_obs))
    
    #GMD2#  
    gmd2_tr<-functionalCP(t_vec, t(tr_dat), t(cal_dat), 5, 2)
    gmd2_pval_vec[i]<-functionalCPpvals(Gmix = gmd2_tr, test = t(ts_obs))
    
    #GMDM2#  
    gmdm2_tr<-functionalCP_m(t_vec, t(tr_dat), t(cal_dat), 5, 2)
    gmdm2_pval_vec[i]<-functionalCPpvals_m(Gmix = gmdm2_tr, test = t(ts_obs))
    
    #GMD3#  
    gmd3_tr<-functionalCP(t_vec, t(tr_dat), t(cal_dat), 4, 3)
    gmd3_pval_vec[i]<-functionalCPpvals(Gmix = gmd3_tr, test = t(ts_obs))
    
    #GMDM2#  
    gmdm3_tr<-functionalCP_m(t_vec, t(tr_dat), t(cal_dat), 4, 3)
    gmdm3_pval_vec[i]<-functionalCPpvals_m(Gmix = gmdm3_tr, test = t(ts_obs))
  }
  
  efdm_pval_mat[,j]<-efdm_pval_vec
  sncm1_pval_mat[,j]<-sncm1_pval_vec
  sncm2_pval_mat[,j]<-sncm2_pval_vec
  gmd1_pval_mat[,j]<-gmd1_pval_vec
  gmdm1_pval_mat[,j]<-gmdm1_pval_vec
  gmd2_pval_mat[,j]<-gmd2_pval_vec
  gmdm2_pval_mat[,j]<-gmdm2_pval_vec
  gmd3_pval_mat[,j]<-gmd3_pval_vec
  gmdm3_pval_mat[,j]<-gmdm3_pval_vec
  
}




