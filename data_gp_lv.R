#
# Create synthetic data for GP-LV model
#



# Create synthetic data ---------------------------------------------------

load(file="DATA_RAPID.Rdata",verbose=TRUE)

data_M = DATA[[1]]
names(data_M)
data_M$C[!is.na(data_M$C)] = 0

data = list()
data$t = 0
data$S0 <- 110
data$q = 0.02
data$rf = seq(0.01, 0.03, length.out=8)
data$T = c(0.1, 0.3, 0.5, 0.7, 1, 1.5, 2, 3)
data$K = c(65,75,80,85,90,95,100,105,110,115,120,125,130,135,140,150,155,160,180,190)
data$T_range = range(data$T)
data$K_range = range(data$K)
data$X_grid_unit = cbind_NA(x_to_unit(data$T,data$T_range),x_to_unit(data$K,data$K_range),x_to_unit(data$t,data$T_range))
data$X_unit = as.matrix(expand.grid(x_to_unit(data$T,data$T_range),x_to_unit(data$K,data$K_range),x_to_unit(data$t,data$T_range)))
data$nPoints = length(data$T)*length(data$K)

data$X_grid_unit = data$X_grid_unit[,-3]
data$X_unit = data$X_unit[,-3]

C_map_data = localVolCalls(data$S0,data$rf,data$q,LV_map+10*LV_sd,data$K,data$T,KflatExt=data$S0*Mext) 
IV_map_data = localVolCalls(data$S0,data$rf,data$q,LV_map+10*LV_sd,data$K,data$T,KflatExt=data$S0*Mext,impVol=TRUE) 

data$C = C_map_data
data$C[is.na(data$IV)] = NA
data$IV = IV_map_data 
data$IV[is.na(data$C)] = NA
rownames(data$C) <- data$T*365
colnames(data$C) <- data$K
data$nObs = sum(!is.na(data$C))

#save(data,file="data_gp_lv.Rdata")


DATA <- list()
DATA[[1]] <- data

for(i in 2:10){
  
  data = list()
  data$t = DATA[[i-1]]$t + 1/52
  data$S0 <- DATA[[i-1]]$S0 + 5*rnorm(1)
  print(data$S0)
  data$q = 0.02
  data$rf = seq(0.01, 0.03, length.out=8)
  data$T = DATA[[i-1]]$T - 1/52
  if(i == 4){
    data$T = c(data$T[2:length(data$T)],3)
  }
  print(data$T)
  data$K = c(65,75,80,85,90,95,100,105,110,115,120,125,130,135,140,150,155,160,180,190)
  data$T_range = range(data$T)
  data$K_range = range(data$K)
  data$X_grid_unit = cbind_NA(x_to_unit(data$T,data$T_range),x_to_unit(data$K,data$K_range),x_to_unit(data$t,data$T_range))
  data$X_unit = as.matrix(expand.grid(x_to_unit(data$T,data$T_range),x_to_unit(data$K,data$K_range),x_to_unit(data$t,data$T_range)))
  data$nPoints = length(data$T)*length(data$K)
  LV = f_to_mat(link(f_states[,sample(ncol(f_states),1)]))
  data$C = localVolCalls(data$S0,data$rf,data$q,LV,data$K,data$T,KflatExt=data$S0*Mext) 
  data$IV = localVolCalls(data$S0,data$rf,data$q,LV,data$K,data$T,KflatExt=data$S0*Mext,impVol=TRUE) 
  data$C[is.na(data_M$C)] = NA
  data$IV[is.na(data_M$C)] = NA
  rownames(data$C) <- data$T*365
  colnames(data$C) <- data$K
  data$nObs = sum(!is.na(data$C))
  
  DATA[[i]] <- data
}

#save(DATA,file="data_sequence_gp_lv.Rdata")


