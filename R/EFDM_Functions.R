#####EFDM Functions#####
library(fdasrvf)

#Note that all functions herein that take functional data as input assume the shape (number of functions) x (number of time points)

###Function for computing calibration scores using EFDM###
###Input:
# tr_dat:      Training data as n_tr x nt matrix
# cal_dat:     Calibration data as n_cal x nt matrix
# t_vec:       Time vector
# wgt:         Weight for amplitude distance in NCM
# k_mean:      Karcher mean
# iv_dist:     If True, use initial value distance in NCM
# iv_wgt:      Weight vector to apply to each component in NCM in the form (amplitude, phase, initial_value)
####################
###Output:
# Karcher mean
# Calibration scores
efdm_cp<-function(tr_dat, cal_dat, t_vec, wgt = 0.5 , k_mean, iv_dist=F, iv_wgt = c(1/3, 1/3, 1/3)){
  
  nt<-length(t_vec)
  n_tr<-dim(tr_dat)[1]
  n_cal<-dim(cal_dat)[1]
  
  tr_dists<-matrix(NA, n_tr, 2)
  for(i in 1:n_tr){
    temp_tr<-elastic.distance(k_mean, tr_dat[i,], t_vec)
    tr_dists[i,]<-c(temp_tr$Dx, temp_tr$Dy)
  }
  
  if(iv_dist){
    tr_ivd<-abs(tr_dat[,1] - k_mean[1])
    tr_dists<-cbind(tr_dists, tr_ivd)
  }

  cal_dists<-matrix(NA, n_cal, 2)
  for(i in 1:n_cal){
    temp<-elastic.distance(k_mean, cal_dat[i,], t_vec)
    cal_dists[i,]<-c(temp$Dx, temp$Dy)
  }
  if(iv_dist){
    cal_ivd<-abs(cal_dat[,1] - k_mean[1])
    cal_dists<-cbind(cal_dists,cal_ivd)
  }
  
  tr_amp_max<-max(tr_dists[,2])
  tr_ph_max<-max(tr_dists[,1])
  tr_amp_min<-min(tr_dists[,2])
  tr_ph_min<-min(tr_dists[,1])
  if(iv_dist){
    tr_iv_max<-max(tr_dists[,3])
    tr_iv_min<-min(tr_dists[,3])
  }
  
  cal_scd_amp<-(cal_dists[,2] - tr_amp_min)/(tr_amp_max - tr_amp_min)
  cal_scd_ph<-(cal_dists[,1] - tr_ph_min)/(tr_ph_max - tr_ph_min)
  if(iv_dist){
    cal_scd_ivd<-(cal_dists[,3] - tr_iv_min)/(tr_iv_max - tr_iv_min)
  }
  
  cal_scores<-wgt*cal_scd_amp + (1-wgt)*cal_scd_ph
  if(iv_dist){
    cal_scores<-iv_wgt[1]*cal_scd_amp + iv_wgt[2]*cal_scd_ph + iv_wgt[3]*cal_scd_ivd
  }
  
  ret_list<-list(KarcherMean = k_mean, CalibrationScores = cal_scores, TrainDists = tr_dists, AmpPhWeight = wgt, InitialValueDist = iv_dist,
                 InitialValueWeight = iv_wgt)
  return(ret_list)

}

###Function that outputs test score and p-value for a set of test functions using output of efdm_cp###
###Input:
# test_dat:    Test data as n_ts x nt matrix
# cal_out:     Output list from efdm_cp()
# t_vec:       Time vector
####################
###Output:
# vector of p-values
# vector indicating inlier/outlier (1/0)

efdm_cp_pred<-function(ts_dat, cal_out, t_vec){
  
  k_mean<-cal_out$KarcherMean
  cal_scores<-cal_out$CalibrationScores
  tr_dists<-cal_out$TrainDists
  wgt<-cal_out$AmpPhWeight
  iv_dist<-cal_out$InitialValueDist
  iv_wgt<-cal_out$InitialValueWeight
  nt<-length(t_vec)
  n_cal<-length(cal_scores)
  n_ts<-dim(ts_dat)[1]
  
  ts_dists<-matrix(NA, n_ts, 2)
  for(i in 1:n_ts){
    temp_ts<-elastic.distance(k_mean, ts_dat[i,], t_vec)
    ts_dists[i,]<-c(temp_ts$Dx, temp_ts$Dy)
  }
  if(iv_dist){
    ts_ivd<-abs(ts_dat[,1] - k_mean[1])
    ts_dists<-cbind(ts_dists,ts_ivd)
  }
  
  tr_amp_max<-max(tr_dists[,2])
  tr_ph_max<-max(tr_dists[,1])
  tr_amp_min<-min(tr_dists[,2])
  tr_ph_min<-min(tr_dists[,1])
  if(iv_dist){
    tr_ivd_max<-max(tr_dists[,3])
    tr_ivd_min<-min(tr_dists[,3])
  }
  
  ts_scd_amp<-(ts_dists[,2] - tr_amp_min)/(tr_amp_max - tr_amp_min)
  ts_scd_ph<-(ts_dists[,1] - tr_ph_min)/(tr_ph_max - tr_ph_min)
  if(iv_dist){
    ts_scd_ivd<-(ts_dists[,3] - tr_ivd_min)/(tr_ivd_max - tr_ivd_min)
  }
    
  ts_scores<-wgt*ts_scd_amp + (1-wgt)*ts_scd_ph
  if(iv_dist){
    ts_scores<-iv_wgt[1]*ts_scd_amp + iv_wgt[2]*ts_scd_ph + iv_wgt[3]*ts_scd_ivd
  }
  
  p_val<-rep(NA, times = n_ts)
  for(i in 1:n_ts){
    count<-(length(which(cal_scores>ts_scores[i])) + runif(1)*(1+length(which(cal_scores==ts_scores[i])))) 
    p_val[i]<-(count)/(n_cal+1)
  }
  ret<-p_val
  return(ret)
}