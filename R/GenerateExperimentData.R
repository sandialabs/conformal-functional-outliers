#####Generate Data for Simulation Experiments#####
library(abind)
source('R/DataGenerationFunctions.R')

##Use for all simulated data sets##
t_vec<-seq(-8,8, length.out = 250)
set.seed(42)

####Generate simulated data for Experiment 1####
###Test data sets###
##Narrow data##
nar_hph_test<-create_data_set(50, t_vec, ph_sd = 1.75, type = 'nar')

##Double peaked data##
dpk_hph_test<-create_data_set(50, t_vec, ph_sd = 1.75, type = 'dpk')

##Standard test data##
standard_hph_test<-create_data_set(500, t_vec, ph_sd = 1.75, type = 'standard')

###Training/Calibration data sets###
#save as two lists of arrays (one for train and one for calibration) with 3 list elements (one for each value of n) 
#array size: n_tr x nt x 500
#create arrays first by saving matrices into a temporary list, then use abind()#
#n1: n=50, p_sd = 1.75#
tr_list_n1<-list()
cal_list_n1<-list()
for(i in 1:500){
  temp_full_train<-create_data_set(50,t_vec, 1.75, type = 'standard')
  temp_spl<-tr_cal_split(temp_full_train)
  tr_list_n1[[i]]<-temp_spl$tr_dat
  cal_list_n1[[i]]<-temp_spl$cal_dat
}
tr_arr_n1<-abind(tr_list_n1, along = 3)
cal_arr_n1<-abind(cal_list_n1, along = 3)

#n2: n=100, p_sd = 1.75#
tr_list_n2<-list()
cal_list_n2<-list()
for(i in 1:500){
  temp_full_train<-create_data_set(100,t_vec, 1.75, type = 'standard')
  temp_spl<-tr_cal_split(temp_full_train)
  tr_list_n2[[i]]<-temp_spl$tr_dat
  cal_list_n2[[i]]<-temp_spl$cal_dat
}
tr_arr_n2<-abind(tr_list_n2, along = 3)
cal_arr_n2<-abind(cal_list_n2, along = 3)

#n3: n=250, p_sd = 1.75#
tr_list_n3<-list()
cal_list_n3<-list()
for(i in 1:500){
  temp_full_train<-create_data_set(250,t_vec, 1.75, type = 'standard')
  temp_spl<-tr_cal_split(temp_full_train)
  tr_list_n3[[i]]<-temp_spl$tr_dat
  cal_list_n3[[i]]<-temp_spl$cal_dat
}
tr_arr_n3<-abind(tr_list_n3, along = 3)
cal_arr_n3<-abind(cal_list_n3, along = 3)

tr_list_exp1<-list(tr_arr_n1 = tr_arr_n1, tr_arr_n2 = tr_arr_n2, tr_arr_n3 = tr_arr_n3)
cal_list_exp1<-list(cal_arr_n1 = cal_arr_n1, cal_arr_n2 = cal_arr_n2, cal_arr_n3 = cal_arr_n3)


####Generate simulated data for Experiment 2####
###Note that the training/calibration data sets for experiment 2 are the same as in experiment 1###
##Generate test data##
#save as a list of arrays with two elements (one for 5% outliers, one for 10% outliers)
#array size: 100 x nt x 500
#create arrays first by saving matrices into a temporary list, then use abind()#
n_2<-100
p1_in<-n_2*(1-0.05)
p2_in<-n_2*(1-0.10)
p1_out<-n_2-p1_in
p2_out<-n_2-p2_in



#p1: n=100, p1 = 0.05#
test_list_p1<-list()
for(i in 1:500){
  temp_in<-create_data_set(p1_in,t_vec, 1.75, type = 'standard')
  temp_out<-create_data_set(p1_out, t_vec, ph_sd = 1.75, type = 'dpk')
 test_list_p1[[i]]<-rbind(temp_in, temp_out)
}
test_arr_p1<-abind(test_list_p1, along = 3)


#p2: n=100, p2 = 0.10#
test_list_p2<-list()
for(i in 1:500){
  temp_in<-create_data_set(p2_in,t_vec, 1.75, type = 'standard')
  temp_out<-create_data_set(p2_out, t_vec, ph_sd = 1.75, type = 'dpk')
  test_list_p2[[i]]<-rbind(temp_in, temp_out)
}
test_arr_p2<-abind(test_list_p2, along = 3)

##Create final array##
test_list_exp2<-list(test_arr_p1 = test_arr_p1, test_arr_p2 = test_arr_p2)


####Generate simulated data for Experiment 3####
##In this experiment, the test data set is the same as the full train data set##
##The full training sets are randomly selected from the 2 list elements - 5% outliers and 10% outliers, respectively##
p1_rand_idx<-sample(500, 100)
p2_rand_idx<-sample(500, 100)
test_p1_arr<-test_list_exp2[[1]][,,p1_rand_idx]
test_p2_arr<-test_list_exp2[[1]][,,p2_rand_idx]
test_list_exp3<-list(test_p1_arr = test_p1_arr, test_p2_arr = test_p2_arr)
