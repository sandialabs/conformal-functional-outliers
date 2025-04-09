###Analysis of Temperature Data###
#This code implements EFDM and leave-one-out EFDM on the temperature data#

##source EFDM and leave-one-out functions##
source('R/EFDM_Functions.R')
source('R/LeaveOneOutFunction.R')
source('R/DataGenerationFunctions.R')

#read in data#
t_vec<-readRDS('Data/TemperatureData/t_vec.rds')        #time vector
t_norm<-seq(0, 1, length.out = length(t_vec))           #transform domain to [0,1]
atl_dat<-readRDS('Data/TemperatureData/atl_dat.rds')    #n x n_t matrix of Atlanta data
dal_dat<-readRDS('Data/TemperatureData/dal_dat.rds')    #n x n_t matrix of Dallas-Love Field data
dfw_dat<-readRDS('Data/TemperatureData/dfw_dat.rds')    #n x n_t matrix of Dallas-Fort Worth data
mia_dat<-readRDS('Data/TemperatureData/mia_dat.rds')    #n x n_t matrix of Miami data
sfo_dat<-readRDS('Data/TemperatureData/sfo_dat.rds')    #n x n_t matrix of San Francisco data

######Standard EFDM analysis for comparing sites#######
#split data sets into training and calibration#
atl_spl<-tr_cal_split(atl_dat)
atl_tr<-atl_spl$tr_dat
atl_cal<-atl_spl$cal_dat

dal_spl<-tr_cal_split(dal_dat)
dal_tr<-dal_spl$tr_dat
dal_cal<-dal_spl$cal_dat

dfw_spl<-tr_cal_split(dfw_dat)
dfw_tr<-dfw_spl$tr_dat
dfw_cal<-dfw_spl$cal_dat

mia_spl<-tr_cal_split(mia_dat)
mia_tr<-mia_spl$tr_dat
mia_cal<-mia_spl$cal_dat

sfo_spl<-tr_cal_split(sfo_dat)
sfo_tr<-sfo_spl$tr_dat
sfo_cal<-sfo_spl$cal_dat

#train on atl, predict on all others#
atl_tw<-time_warping(t(atl_tr), t_norm, parallel = T)
atl_km<-atl_tw$fmean
atl_tr_cp<-efdm_cp(atl_tr, atl_cal, t_norm, k_mean = atl_km, iv_dist = T)
atl_tr_pval<-efdm_cp_pred(rbind(dal_dat, dfw_dat, mia_dat, sfo_dat),atl_tr,cp, t_norm)

#train on dal, predict on all others#
dal_tw<-time_warping(t(dal_tr), t_norm, parallel = T)
dal_km<-dal_tw$fmean
dal_tr_cp<-efdm_cp(dal_tr, dal_cal, t_norm, k_mean = dal_km, iv_dist = T)
dal_tr_pval<-efdm_cp_pred(rbind(atl_dat, dfw_dat, mia_dat, sfo_dat),dal_tr,cp, t_norm)

#train on dfw, predict on all others#
dfw_tw<-time_warping(t(dfw_tr), t_norm, parallel = T)
dfw_km<-dfw_tw$fmean
dfw_tr_cp<-efdm_cp(dfw_tr, dfw_cal, t_norm, k_mean = dfw_km, iv_dist = T)
dfw_tr_pval<-efdm_cp_pred(rbind(atl_dat, dal_dat, mia_dat, sfo_dat),dfw_tr,cp, t_norm)

#train on mia, predict on all others#
mia_tw<-time_warping(t(mia_tr), t_norm, parallel = T)
mia_km<-mia_tw$fmean
mia_tr_cp<-efdm_cp(mia_tr, mia_cal, t_norm, k_mean = mia_km, iv_dist = T)
mia_tr_pval<-efdm_cp_pred(rbind(atl_dat, dal_dat, dfw_dat, sfo_dat),mia_tr,cp, t_norm)

#train on sfo, predict on all others#
sfo_tw<-time_warping(t(sfo_tr), t_norm, parallel = T)
sfo_km<-sfo_tw$fmean
sfo_tr_cp<-efdm_cp(sfo_tr, sfo_cal, t_norm, k_mean = sfo_km, iv_dist = T)
sfo_tr_pval<-efdm_cp_pred(rbind(atl_dat, dal_dat, dfw_dat, mia_dat),sfo_tr,cp, t_norm)

######Leave-one-out analysis for within-site outliers#######
#Within-lot p-values for atl#
atl_loo<-loo_fn(atl_dat, t_norm, method = "EFDM", kmean = T, efdm_wgt = 0.5, efdm_iv_dist = T, efdm_iv_wgt = c(1/3, 1/3, 1/3))

#Within-lot p-values for dal#
dal_loo<-loo_fn(dal_dat, t_norm, method = "EFDM", kmean = T, efdm_wgt = 0.5, efdm_iv_dist = T, efdm_iv_wgt = c(1/3, 1/3, 1/3))

#Within-lot p-values for dal#
dfw_loo<-loo_fn(dfw_dat, t_norm, method = "EFDM", kmean = T, efdm_wgt = 0.5, efdm_iv_dist = T, efdm_iv_wgt = c(1/3, 1/3, 1/3))

#Within-lot p-values for dal#
mia_loo<-loo_fn(mia_dat, t_norm, method = "EFDM", kmean = T, efdm_wgt = 0.5, efdm_iv_dist = T, efdm_iv_wgt = c(1/3, 1/3, 1/3))

#Within-lot p-values for dal#
sfo_loo<-loo_fn(sfo_dat, t_norm, method = "EFDM", kmean = T, efdm_wgt = 0.5, efdm_iv_dist = T, efdm_iv_wgt = c(1/3, 1/3, 1/3))
