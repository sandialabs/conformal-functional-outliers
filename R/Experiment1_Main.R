###FDA Outlier Detection Experiment 1###
#This code runs the first simulation experiment#

##source ICAD functions#
source('R/EFDM_Functions.R')
source('R/SNCM_Functions.R')
source('R/GMD_Functions.R')


##define t_vec and alpha##
nt<-250
t_vec<-seq(-8, 8, length.out = nt)
alpha<-0.10

##read in data##
tr_list<-readRDS('Data/SimulatedData/tr_list_exp1.rds')
cal_list<-readRDS('Data/SimulatedData/cal_list_exp1.rds')

st_test<-readRDS('Data/SimulatedData/standard_test.rds')
nar_test<-readRDS('Data/SimulatedData/nar_test.rds')
dpk_test<-readRDS('Data/SimulatedData/dpk_test.rds')

full_test<-rbind(st_test, nar_test, dpk_test) ###1:500 - standard; 501:550 - narrow; 551-600 - double peaked###

##Placeholders for results##
efdm_pval_list<-list()
sncm1_pval_list<-list()
sncm2_pval_list<-list()
gmd1_pval_list<-list()
gmd2_pval_list<-list()
gmd3_pval_list<-list()
gmdm_pval_list<-list()


for(j in 1:length(tr_list)){
  
  tr_arr<-tr_list[[j]] 
  cal_arr<-cal_list[[j]] 
  
  ##temporary placeholders##
  efdm_pval<-sncm1_pval<-sncm2_pval<-gmd1_pval<-gmd2_pval<-gmd2_pval<-gmd3_pval<-gmdm_pval<-matrix(NA, dim(full_test)[1], dim(tr_arr)[3])

  #Inner loop:
  for(i in 1:dim(tr_arr)[3]){
    
    #Define tr and calibration data#
    tr_dat<-tr_arr[,,i]
    cal_dat<-cal_arr[,,i]
    
    #compute Karcher mean of training data#
    tw<-time_warping(t(tr_dat),t_vec)
    k_mean<-tw$fmean
    
    #EFDM
    efdm_tr<-efdm_cp(tr_dat, cal_dat, t_vec, wgt = 0.5, k_mean)
    efdm_pred<-efdm_cp_pred(full_test, efdm_tr, t_vec)
    efdm_pval[,i]<-efdm_pred
    
    #SNCM1
    sncm1_tr<-sncm_cp(tr_dat, cal_dat, t_vec, 'csm', 'optimal', alpha, k_mean)
    sncm1_pred<-sncm_cp_pred(full_test, sncm1_tr, t_vec)
    sncm1_pval[,i]<-sncm1_pred$TestPValue
    
    #SNCM2
    sncm2_tr<-sncm_cp(tr_dat, cal_dat, t_vec, 'Karcher', 'optimal', alpha, k_mean)
    sncm2_pred<-sncm_cp_pred(full_test, sncm2_tr, t_vec)
    sncm2_pval[,i]<-sncm2_pred$TestPValue
    
    #GMD1
    gmd1_tr<-functionalCP(t_vec, t(tr_dat), t(cal_dat), 2, 1)
    gmd1_pval[,i]<-functionalCPpvals(Gmix = gmd1_tr, test = t(full_test))
    
    #GMD2
    gmd2_tr<-functionalCP(t_vec, t(tr_dat), t(cal_dat), 3, 2)
    gmd2_pval[,i]<-functionalCPpvals(Gmix = gmd2_tr, test = t(full_test))
    
    #GMD3
    if(dim(tr_dat)[1]<100){next} #skip because not all of them will have GMMs that can be fit
    gmd3_tr<-functionalCP(t_vec, t(tr_dat), t(cal_dat), 5, 2)
    gmd3_pval[,i]<-functionalCPpvals(Gmix = gmd3_tr, test = t(full_test))
    
    #GMDM
    gmdm_tr<-functionalCP_m(input = t_vec, train = t(tr_dat), cal = t(tr_dat), 5, 2)
    gmdm_pval[,i]<-functionalCPpvals_m(Gmix = gmdm_tr, test = t(full_test))
    
    
  }
  
  efdm_pval_list[[j]]<-efdm_pval
  sncm1_pval_list[[j]]<-sncm1_pval
  sncm2_pval_list[[j]]<-sncm2_pval
  gmd1_pval_list[[j]]<-gmd1_pval
  gmd2_pval_list[[j]]<-gmd2_pval
  gmd3_pval_list[[j]]<-gmd3_pval
  gmdm_pval_list[[j]]<-gmdm_pval
  
  
}
