###FDA Outlier Detection Experiment 2###
#This code runs the second simulation experiment#

##Source ICAD functions##
source('R/EFDM_Functions.R')
source('R/SNCM_Functions.R')
source('R/GMD_Functions.R')

##define t_vec and alpha##
nt<-250
t_vec<-seq(-8, 8, length.out = nt)
alpha<-0.10

tr_list<-readRDS('Data/SimulatedData/tr_list_exp1.rds')
cal_list<-readRDS('Data/SimulatedData/cal_list_exp1.rds')
ts_list<-readRDS('Data/SimulatedData/exp2_test.rds')

##define two separate arrays for test data, one for each outlier percentage##
ts_arr1<-ts_list[[1]]
ts_arr2<-ts_list[[2]]

#define placeholders for results#
efdm_pval_list1<-list()
sncm1_pval_list1<-list()
sncm2_pval_list1<-list()
gmd1_pval_list1<-list()
gmd2_pval_list1<-list()
gmd3_pval_list1<-list()
gmdm_pval_list1<-list()

efdm_pval_list2<-list()
sncm1_pval_list2<-list()
sncm2_pval_list2<-list()
gmd1_pval_list2<-list()
gmd2_pval_list2<-list()
gmd3_pval_list2<-list()
gmdm_pval_list2<-list()


#Outer loop:


for(j in 1:length(tr_list)){
  
  tr_arr<-tr_list[[j]] 
  cal_arr<-cal_list[[j]]
  
  efdm_pval1<-matrix(NA, dim(ts_arr1)[1], dim(tr_arr)[3])
  sncm1_pval1<-matrix(NA, dim(ts_arr1)[1], dim(tr_arr)[3])
  sncm2_pval1<-matrix(NA, dim(ts_arr1)[1], dim(tr_arr)[3])
  gmd1_pval1<-matrix(NA, dim(ts_arr1)[1], dim(tr_arr)[3])
  gmd2_pval1<-matrix(NA, dim(ts_arr1)[1], dim(tr_arr)[3])
  gmd3_pval1<-matrix(NA, dim(ts_arr1)[1], dim(tr_arr)[3])
  gmdm_pval1<-matrix(NA, dim(ts_arr1)[1], dim(tr_arr)[3])
  
  efdm_pval2<-matrix(NA, dim(ts_arr2)[1], dim(tr_arr)[3])
  sncm1_pval2<-matrix(NA, dim(ts_arr2)[1], dim(tr_arr)[3])
  sncm2_pval2<-matrix(NA, dim(ts_arr2)[1], dim(tr_arr)[3])
  gmd1_pval2<-matrix(NA, dim(ts_arr2)[1], dim(tr_arr)[3])
  gmd2_pval2<-matrix(NA, dim(ts_arr2)[1], dim(tr_arr)[3])
  gmd3_pval2<-matrix(NA, dim(ts_arr2)[1], dim(tr_arr)[3])
  gmdm_pval2<-matrix(NA, dim(ts_arr2)[1], dim(tr_arr)[3])
  
  #Inner loop:
  for(i in 1:dim(tr_arr)[3]){
    
    #Define tr and calibration data#
    tr_dat<-tr_arr[,,i]
    cal_dat<-cal_arr[,,i]
    ts_dat1<-ts_arr1[,,i]
    ts_dat2<-ts_arr2[,,i]
    
    #compute Karcher mean of training data#
    tw<-time_warping(t(tr_dat),t_vec)
    k_mean<-tw$fmean
    
    #EFDM
    efdm_tr<-efdm_cp(tr_dat, cal_dat, t_vec, wgt = 0.5, k_mean)
    efdm_pred1<-efdm_cp_pred(ts_dat1, efdm_tr, t_vec)
    efdm_pval1[,i]<-efdm_pred1
    efdm_pred2<-efdm_cp_pred(ts_dat2, efdm_tr, t_vec)
    efdm_pval2[,i]<-efdm_pred2
    
    #SNCM1
    sncm1_tr<-sncm_cp(tr_dat, cal_dat, t_vec, 'csm', 'optimal', alpha, k_mean)
    sncm1_pred1<-sncm_cp_pred(ts_dat1, sncm1_tr, t_vec)
    sncm1_pval1[,i]<-sncm1_pred1$TestPValue
    sncm1_pred2<-sncm_cp_pred(ts_dat2, sncm1_tr, t_vec)
    sncm1_pval2[,i]<-sncm1_pred2$TestPValue
    
    #SNCM2
    sncm2_tr<-sncm_cp(tr_dat, cal_dat, t_vec, 'Karcher', 'optimal', alpha, k_mean)
    sncm2_pred1<-sncm_cp_pred(ts_dat1, sncm2_tr, t_vec)
    sncm2_pval1[,i]<-sncm2_pred1$TestPValue
    sncm2_pred2<-sncm_cp_pred(ts_dat2, sncm2_tr, t_vec)
    sncm2_pval2[,i]<-sncm2_pred2$TestPValue

    #GMD1
    gmd1_tr<-functionalCP(t_vec, t(tr_dat), t(cal_dat), 2, 1)
    gmd1_pval1[,i]<-functionalCPpvals(Gmix = gmd1_tr, test = t(ts_dat1))
    gmd1_pval2[,i]<-functionalCPpvals(Gmix = gmd1_tr, test = t(ts_dat2))
    
    #GMD2
    gmd2_tr<-functionalCP(t_vec, t(tr_dat), t(cal_dat), 4, 1)
    gmd2_pval1[,i]<-functionalCPpvals(Gmix = gmd2_tr, test = t(ts_dat1))
    gmd2_pval2[,i]<-functionalCPpvals(Gmix = gmd2_tr, test = t(ts_dat2))
    
    #GMD3
    if(dim(tr_dat)[1]<100){next} #skip because not all of them will have GMMs that can be fit
    gmd3_tr<-functionalCP(t_vec, t(tr_dat), t(cal_dat), 5, 2)
    gmd3_pval1[,i]<-functionalCPpvals(Gmix = gmd3_tr, test = t(ts_dat1))
    gmd3_pval2[,i]<-functionalCPpvals(Gmix = gmd3_tr, test = t(ts_dat2))
    
    #GMDM
    gmdm_tr<-functionalCP_m(input = t_vec, train = t(tr_dat), cal = t(tr_dat), 5, 2)
    gmdm_pval1[,i]<-functionalCPpvals(Gmix = gmdm_tr, test = t(ts_dat1))
    gmdm_pval2[,i]<-functionalCPpvals(Gmix = gmdm_tr, test = t(ts_dat2))
    
    
    
  }
  
  efdm_pval_list1[[j]]<-efdm_pval1
  sncm1_pval_list1[[j]]<-sncm1_pval1
  sncm2_pval_list1[[j]]<-sncm2_pval1
  gmd1_pval_list1[[j]]<-gmd1_pval1
  gmd2_pval_list1[[j]]<-gmd2_pval1
  gmd3_pval_list1[[j]]<-gmd3_pval1
  gmdm_pval_list1[[j]]<-gmdm_pval1
  
  efdm_pval_list2[[j]]<-efdm_pval2
  sncm1_pval_list2[[j]]<-sncm1_pval2
  sncm2_pval_list2[[j]]<-sncm2_pval2
  gmd1_pval_list2[[j]]<-gmd1_pval2
  gmd2_pval_list2[[j]]<-gmd2_pval2
  gmd3_pval_list2[[j]]<-gmd3_pval2
  gmdm_pval_list2[[j]]<-gmdm_pval2
  
}

