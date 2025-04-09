#####Data Generation Functions for Simulation Experiments#####

#Standard data generation - for training/calibration data, standard test sets#
standard_data_gen<-function(t_vec, a, p){
  ret<-a*exp(-0.5*(t_vec - p)^2)
  return(ret)
}

#Narrow data generation - for generating narrow test data#
nar_data_gen<-function(t_vec, a, p){
  ret<-a*exp(-2*(t_vec - p)^2)
  return(ret)
}

#double peak data generation - for generating double peak test data#
dpk_data_gen<-function(t_vec, a, p, ph_lev){
  if(ph_lev=='low'){
    ret<-a*exp(-3*(t_vec-0.75 - p)^2/2) + a*exp(-3*(t_vec+0.75 - p)^2/2)
  }
  if(ph_lev=='high'){
    ret<-a*exp(-(t_vec-1.5 - p)^2/2) + a*exp(-(t_vec+1.5 - p)^2/2)
  }
  return(ret)
}

#Function to create a data set#
create_data_set<-function(n, t_vec, ph_sd, type){
  
  nt<-length(t_vec)
  a_i<-rnorm(n, 1, 0.15)
  p_i<-rnorm(n, 0, ph_sd)
  ai_dpk<-rnorm(n,1, 0.10)
  pi_dpk<-rnorm(n,0, 0.15)
  if(ph_sd<1){
    ph_lev<-'low'
  }
  else{
    ph_lev<-'high'
  }
  
  ret_dat<-matrix(NA, n, nt)
  for(i in 1:n){
    
    if(type=='standard'){
      ret_dat[i,]<-standard_data_gen(t_vec, a_i[i], p_i[i])
      next
    }
    if(type=='nar'){
      ret_dat[i,]<-nar_data_gen(t_vec, a_i[i], p_i[i])
      next
    }
    if(type=='dpk'){
      ret_dat[i,]<-dpk_data_gen(t_vec, ai_dpk[i], pi_dpk[i],ph_lev)
      next
    }
  }
  return(ret_dat)
}

#Split full training into train and calibration sets with 2:1 training to calibration ratio#
tr_cal_split<-function(full_tr){
  
  n_full<-dim(full_tr)[1]
  n_tr<-ceiling(2/3*n_full)
  n_cal<-n_full - n_tr
  
  tr_ind<-sample(n_full, n_tr)
  tr_dat<-full_tr[tr_ind,]
  cal_dat<-full_tr[-tr_ind,]
  
  ret_list<-list(tr_dat = tr_dat, cal_dat = cal_dat)
  
}



