###Analysis of Diode Data###
#This code implements EFDM and leave-one-out EFDM on the Zener diode data#

##source EFDM and leave-one-out functions##
source('R/EFDM_Functions.R')
source('R/LeaveOneOutFunction.R')
source('R/DataGenerationFunctions.R')

#read in data#
volt_vec<-readRDS('Data/DiodeData/volt_vec.rds') #vector of voltages over which each device's current is measured
lot1_dat<-readRDS('Data/DiodeData/lot1_dat.rds') #n x n_t matrix of lot 1 data
lot2_dat<-readRDS('Data/DiodeData/lot2_dat.rds') #n x n_t matrix of lot 2 data


#######Standard EFDM Analysis#######
#split lot1 and lot2 data sets into training and calibration sets#
l1_spl<-tr_cal_split(lot1_dat)
l2_spl<-tr_cal_split(lot2_dat)

l1_tr<-l1_spl$tr_dat
l1_cal<-l1_spl$cal_dat

l2_tr<-l2_spl$tr_dat
l2_cal<-l2_spl$cal_dat

#Train on lot 1, predict on lot 2#
l1_tw<-time_warping(t(l1_tr), volt_vec, parallel = T)       #obtain lot 1 Karcher mean
l1_km<-l1_tw$fmean                                        #save lot 1 Karcher mean
l1_tr_cp<-efdm_cp(l1_tr, l1_cal, volt_vec, k_mean = l1_km)  #efdm calibration scores on lot 1
l1tr_l2pval<-efdm_cp_pred(lot2_dat, l1_tr_cp, volt_vec)     #lot 2 pvalues when training on lot 1

#Train on lot 1, predict on lot 2#
l2_tw<-time_warping(t(l2_tr), volt_vec, parallel = T)       #obtain lot 2 Karcher mean
l2_km<-l2_tw$fmean                                        #save lot 2 Karcher mean
l2_tr_cp<-efdm_cp(l2_tr, l2_cal, volt_vec, k_mean = l2_km)  #efdm calibration scores on lot 1
l2tr_l1pval<-efdm_cp_pred(lot1_dat, l2_tr_cp, volt_vec)     #lot 2 pvalues when training on lot 1


#######Leave-one-out EFDM Analysis#######
#Within-lot p-values for Lot 1#
l1_loo<-loo_fn(lot1_dat, volt_vec, method = "EFDM", kmean = T, efdm_wgt = 0.5)

#Within-lot p-values for Lot 2#
l2_loo<-loo_fn(lot2_dat, volt_vec, method = "EFDM", kmean = T, efdm_wgt = 0.5)

