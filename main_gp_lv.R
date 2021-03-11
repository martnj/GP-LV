#
# GP local volatility model
#




# Libraries and functions -------------------------------------------------

require("limSolve")
require("mvtnorm")
source("functions_gp_lv.R")

l.par = list(font.axis=1,font.lab=1,family="serif",cex.lab=1.4,cex.axis=1.4)

plot_LV <- function(f0, f_mu, data){
  
  # Plot a local vol surface
  #
  # GLOBAL: link, l.par
  
  LV = matrix(link(f0 + f_mu), nrow=length(data$T), ncol=length(data$K))
  
  par(l.par)
  zlim = c(0,0.8)
  persp(x=data$K, y=data$T, z=t(LV), theta=35, phi=25, r=sqrt(81),
        col=gray(0.7), border=gray(0.6), axes=T,ticktype="detailed",
        zlim=zlim, zlab='implied volatility', xlab="strike", 
        ylab="maturity")
}


# Set-up model ------------------------------------------------------------

# Define user-specified variables and functions for the model.


# Link-function (positive transform f -> local-vol)
link <- function(f) log(1 + exp(f))
linkInv <- function(sig) log(exp(sig) - 1)

# Log-likelihood function
logL <- function(f0, f_mu, sigma_noise, data){
  
  # Log-likelihood of a local vol surface LV = link(f0 + f_mu) over 'data$T x data$K'.
  # 'f0' is a vector of length = length(data$T)*length(data$K),
  # 'f_mu' and 'sigma_noise' are scalars.
  #
  # If 'data = list(data_1,data_2,...)' => sum individual likelihoods.
  # This requires that 'f0 = c(f0_1,f0_2,...)' of individual f0-vectors.
  #
  # Requires global variables 'Mext' and 'link()'.
  
  # VECTORISE (if 'data' is list of data-lists from several dates)
  if(is.list(data[[1]])){
    res = 0
    idx_stop = 0
    for(i in 1:length(data)){
      idx_start = idx_stop + 1
      idx_stop = idx_start - 1 + data[[i]]$nPoints
      res = res + logL(f0=f0[idx_start:idx_stop],data=data[[i]],f_mu=f_mu,sigma_noise=sigma_noise)
    }
    return(res)
  }
  
  # MAIN FUNCTION
  if(length(f0)!=data$nPoints) stop('length(f0) not matching size TxK of data')
  LV = matrix(link(f0+f_mu),nrow=length(data$T),ncol=length(data$K))
  C_lv = localVolCalls(data$S0,data$rf,data$q,LV,data$K,data$T,KflatExt=data$S0*Mext)
  ll = sum(dnorm(c(data$C),c(C_lv),sigma_noise,log=TRUE),na.rm=TRUE)
  
  return(ll)
}

# Flat extension in moneyness for the FD solver of 'localVolCalls' 
Mext = seq(0.1, 4, by=0.2)

# Hyperpriors: scaled sigmoid Gaussian
ssg_fun <- function(x,y_max,y_min=0) y_min + (y_max-y_min)*1/(1+exp(-x))
ssg_fun_inv <- function(y,y_max,y_min=0) log(y-y_min) - log(y_max-y)
rssg <- function(y_max,y_min=0) ssg_fun(rnorm(length(y_max)),y_max,y_min)
dssg <- function(x,x_max,x_min=0,log=FALSE,z_mu=0,z_sd=1){
  if(!log) return( dnorm(log((x-x_min)/(x_max-x)),mean=z_mu,sd=z_sd)*(x_max-x_min)/((x_max-x)*(x-x_min)) )
  return( log( dnorm(log((x-x_min)/(x_max-x)),mean=z_mu,sd=z_sd)*(x_max-x_min)/((x_max-x)*(x-x_min)) ) )
}

# GP-prior: covariance funcitons
type = 'SE'    # a global variable, one of 'SE', 'Mat12', 'Mat32', 'Mat52'
K_fun <- function(kappa, data, vsdiag=NULL){
  
  # Makes covariance matrix K. Inputs X from 'data$X_grid_unit', 
  # a Cartesian grid of [strikes x maturities] normalised to the unit
  # interval. Optional vector 'vsdiag' scales K vertically. Data
  # can also be list of data-lists, from several dates.
  #
  # Req. global: 'type' in {'SE', 'Mat12', 'Mat32', 'Mat52'}.
  
  if(!is.list(data[[1]])){
    if(is.null(vsdiag)){
      return( covM_SE(X=data$X_grid_unit[,1:(length(kappa)-1)],kappa=kappa,structure='grid',type=type) )
    }else{
      return( diag(vsdiag)%*%covM_SE(X=data$X_grid_unit[,1:(length(kappa)-1)],kappa=kappa,structure='grid',type=type)%*%diag(vsdiag) )
    }
  }else{
    X_unit = c()
    for(i in 1:length(data)){
      X_unit = rbind(X_unit,data[[i]]$X_unit[,1:(length(kappa)-1)])
    }
    if(!is.null(vsdiag)) stop('vertical sclaing not yet implemented if data is list')
    return(covM_SE(X=X_unit,kappa=kappa,type=type))
  }
}
L_fun <- function(kappa, data, vsdiag=NULL){
  
  # Same as 'K_fun' but outputs cholesky decomposition.
  
  if(!is.list(data[[1]])){
    if(is.null(vsdiag)){
      return( covM_SE(X=data$X_grid_unit[,1:(length(kappa)-1)],kappa=kappa,structure='grid',decomp='chol',type=type) )
    }else{
      return( diag(vsdiag)%*%covM_SE(X=data$X_grid_unit[,1:(length(kappa)-1)],kappa=kappa,structure='grid',decomp='chol',type=type) )
    }
  }else{
    X_unit = c()
    for(i in 1:length(data)){
      X_unit = rbind(X_unit,data[[i]]$X_unit[,1:(length(kappa)-1)])
    }
    if(!is.null(vsdiag)) stop('vertical sclaing not yet implemented if data is list')
    return(covM_SE(X=X_unit,kappa=kappa,decomp='chol',type=type))
  }
}
K_mu_cond <- function(kappa, t, f0_prev, chol=FALSE){
  
  # Conditional moments for f0(t)|f0(t-1) = f0_prev with
  # inputs X(t) and X(t-1) from DATA[[t]] and DATA[[t-1]]
  # Given function values: f(t-1) = f0_prev (zero-mean)
  #
  # Req. global: 'DATA', 'type'.
  #
  # OBS: conditional K not yet implemented with vertical scaling.
  
  res = f_cond_mu_K_Xmat(X_star_unit=DATA[[t]]$X_grid_unit,X_obs_unit=DATA[[t-1]]$X_grid_unit,f_obs=f0_prev,kappa=kappa,f_mu=0,type=type)
  if(chol) return( list(K_cond=res$f_star_K,f_mu_cond=res$f_star_mu,L_cond=LfromK(res$f_star_K,method='chol')) )
  else return( list(K_cond=res$f_star_K,f_mu_cond=res$f_star_mu) )
}



# Synthetic data ----------------------------------------------------------

# Load a synthetic data-set with call prices from a single date.


# Load data-set
load(file="data_gp_lv.Rdata", verbose=TRUE)

# Plot implied volatility surface
x = data$K
y = data$T
X = t(data$IV)
zlim = c(0,0.45)
par(l.par)
pmat = persp(x=x,y=y,z=X,theta=35,phi=25,border=NA,col="transparent",ticktype="detailed",axes=T,r=sqrt(81),
             zlab="implied volatility",xlab='strike',ylab="maturity",zlim=zlim)
pointIdx = is.finite(X)
a = rep(x, length(y))[pointIdx]
b = rep(y, each=length(x))[pointIdx]
X3 = matrix(zlim[1], nrow=nrow(X), ncol=ncol(X))
xy3 = trans3d(a, b, X3[pointIdx], pmat)
points(xy3, pch=16, cex=0.5, col=gray(0.5))
xy4 = trans3d(rep(data$S0, length(y)), y, zlim[1], pmat)
for(i in 1:ncol(pointIdx)){
  a = x[pointIdx[,i]]
  b = rep(y[i], length(a))
  xy = trans3d(a, b, X[,i][pointIdx[,i]], pmat)
  points(xy, pch=16, cex=.5, type='o', lwd=1)
}


# MCMC algorithm ----------------------------------------------------------

# Markov chain Monte Carlo sampling of local-vol surface and 
# hyperparameters. 


# Define support of the hyperprior (kappa, f_mu, sigma), with 'kappa = (l_T, l_K, sigma_f)'
kappaMax = c(1,1,1)
kappaMin = c(0.1,0.1,0.1)
lv_mu_max = .5
sigma_noise_max = .75

# Generate initial (kappa,f_mu,sigma)-state:
set.seed(1)
kappa = ssg_fun(rnorm(length(kappaMax)), kappaMax, kappaMin)
xi = rnorm(2) 
f_mu = linkInv( ssg_fun(xi[1], lv_mu_max) )
sigma_noise = ssg_fun(xi[2], sigma_noise_max)

# Generate initial f0-state:
K = K_fun(kappa,data)
L = L_fun(kappa,data)
f0 = L%*%matrix(rnorm(nrow(L)))

# Plot initial LV surface
plot_LV(f0, f_mu, data)

# Compute log-likelihood of current state
llcur = logL(f0, f_mu, sigma_noise, data)

# Number of iterations/states of the Markov chain 
nMH = 5000

# Pre-allocate & monitor:
initial_states = c(f0, kappa, f_mu, sigma_noise, llcur)
states = matrix(nrow=length(initial_states), ncol = nMH)
nIter_f0 <- nIter_kappa <- nIter_likh <- rep(0,nMH)

R <- LR <- NULL
for(i in 1:nMH){
  
  # Ttransition 1: sample f0 given (kappa,f_mu,sigma_noise) with elliptical slice sampling x3.
  L = L_fun(kappa,data)
  ess = ESS_f0(f0,L=L,llcur=llcur,logL=logL,nESS=3,angleRange=pi/2^3,f_mu=f_mu,sigma_noise=sigma_noise,data=data)
  f0 = ess$f0
  llcur = ess$llcur
  
  # Transition 2: sample kappa,f0 given (f_mu,sigma_noise) with surrogate data elliptical slice sampling
  sdess = SDESS(f0=f0,kappa=kappa,kappaMin=kappaMin,kappaMax=kappaMax,logL=logL,llcur=llcur,data=data,K=K,R=R,LR=LR,f_mu=f_mu,sigma_noise=sigma_noise)
  f0 = sdess$f0
  kappa = sdess$kappa
  llcur = sdess$llcur
  K = sdess$K
  R = sdess$R
  LR = sdess$LR
  
  # Transition 3: sample (f_mu,sigma_noise) given (kappa,f0) with elliptical slice sampling x3.
  for(k in 1:3){
    # Set ll-threshold & generate ellipse
    llth <- log(runif(1)) + llcur
    xi_nu <- rnorm(2)
    # Initial proposal & define bracket
    theta <- runif(1,0,2*pi)
    bracket <- c(theta-2*pi,theta)
    nIter = 1
    while(TRUE){
      xi_prop <- xi[1:2]*cos(theta) + xi_nu*sin(theta)
      f_muProp <- linkInv( ssg_fun(xi_prop[1],lv_mu_max) )
      sigma_noiseProp <- ssg_fun(xi_prop[2],sigma_noise_max)
      llprop <- logL(f0,f_muProp,sigma_noiseProp,data)
      # Accept porposal or shrink backet?
      if(llprop > llth){
        xi[1:2] <- xi_prop
        f_mu <- f_muProp
        sigma_noise <- sigma_noiseProp
        llcur <- llprop
        break
      }else{
        nIter = nIter + 1
        bracket[(theta>0)+1] = theta
        theta <- runif(1,bracket[1],bracket[2])
      }
    }
  }
  
  # Store variables
  states[,i] = c(f0,kappa,f_mu,sigma_noise,llcur)
  nIter_f0[i] = ess$nIter
  nIter_kappa[i] = sdess$nIter
  nIter_likh[i] = nIter
  print(i)
}

#save(states,kappaMin,kappaMax,lv_mu_max,sigma_noise_max,type,Mext,link,linkInv,file='mcmc_states_static.Rdata')



# Sequential MCMC algorithm -----------------------------------------------

# Load data (10 dates) and MCMC-results for first date
load(file='mcmc_states_static.Rdata',verbose=TRUE)
load(file="data_sequence_gp_lv.Rdata",verbose=TRUE)

# Extract variables for t=1 from static MCMC sample
idx = seq(1000, ncol(states), length.out=500)
N = data$nPoint
f0_states = states[1:N,idx] 
kappaStates = states[(1:length(kappaMax))+N,idx]
f_muStates = states[length(kappaMax)+1+N,idx]
sigma_noiseStates = states[length(kappaMax)+2+N,idx]
llVec = states[nrow(states),idx]

# List for states and store the t=1 results
nMH = length(idx)
STATES <- list()
STATES[[1]] <- list(f0_states=f0_states,
                    kappaStates=kappaStates,
                    f_muStates=f_muStates,
                    sigma_noiseStates=sigma_noiseStates,
                    llVec=llVec,
                    xi_kappa_mu=rep(0,4),
                    xi_kappa_K=diag(rep(1.5^2,4)),
                    xi_alpha_mu=c(0,0),
                    xi_alpha_K=diag(c(1.5^2,1.5^2)))


# Update support of the hyper-prior; 'kappa = (l_T, l_K, l_t, sigma_f)'
kappaMax[c(1,2,4)] = kappaMax
kappaMax[3] = 1
kappaMin[c(1,2,4)] = kappaMin
kappaMin[3] = 0.1

updateKappa = TRUE
updateAlpha = TRUE
for(t in 2:length(DATA)){
  
  # Pre-allocate for this t
  f0_states = matrix(nrow=DATA[[t]]$nPoints,ncol=nMH)
  kappaStates = matrix(nrow=length(kappaMax),ncol=nMH)
  f_muStates <- sigma_noiseStates <- llVec <- numeric(nMH)
  
  # Update kappa & alpha approximative posterior, from (t-1)-sample 
  if(t==2){
    tmp = ssg_fun_inv(STATES[[t-1]]$kappaStates,kappaMax[-3],kappaMin[-3])
    xi_kappa_mu = numeric(length(kappaMax))
    xi_kappa_mu[-3] = rowMeans(tmp)
    xi_kappa_K = 1.5^2*diag(length(kappaMax))
    xi_kappa_K[-3,-3] = cov(t(tmp))
    xi_kappa_L = t(chol(xi_kappa_K))
    kappa = c( ssg_fun(xi_kappa_L%*%rnorm(length(kappaMax)) + xi_kappa_mu,kappaMax,kappaMin) )
  }else{
    tmp = ssg_fun_inv(STATES[[t-1]]$kappaStates,kappaMax,kappaMin)
    xi_kappa_mu = rowMeans(tmp)
    xi_kappa_K = cov(t(tmp))
    xi_kappa_L = t(chol(xi_kappa_K))
  }
  tmp = ssg_fun_inv(rbind(link(STATES[[t-1]]$f_muStates),STATES[[t-1]]$sigma_noiseStates),c(lv_mu_max,sigma_noise_max))
  xi_alpha_mu = rowMeans(tmp)
  xi_alpha_K = cov(t(tmp))
  xi_alpha_L = t(chol(xi_alpha_K))
  if(t==2){
    alpha = c( ssg_fun(xi_alpha_L%*%rnorm(length(xi_alpha_mu)) + xi_alpha_mu,c(lv_mu_max,sigma_noise_max)) )
    f_mu = linkInv(alpha[1])
    sigma_noise = alpha[2]
  }
  
  # Based on each posterior state f0(t-1)_i, i=1,...,nMH, update {f0(t), kappa, alpha}
  for(i in 1:nMH){
    print(i)
    
    f0_prev = STATES[[t-1]]$f0_states[,i]
    
    #### f0(t) | f0(t-1),kappa,f_mu,sigma_noise ~ ESS ####
    res = K_mu_cond(kappa,t,f0_prev,chol=TRUE)
    K = res$K_cond
    L = res$L_cond   
    f0_mu = res$f_mu_cond 
    
    # Draw initial f0(t)|f0(t-1) and calculate likeihood
    #f0 = f0_prev
    f0 = L%*%matrix(rnorm(length(f0_mu)))
    llcur = logL(f0,f_mu+f0_mu,sigma_noise,data=DATA[[t]])
    
    # Update f0(t)|f0(t-1),L,f0_mu ~ ESS x3
    ess = ESS_f0(f0,L=L,llcur=llcur,logL=logL,nESS=3,data=DATA[[t]],f_mu=f_mu+f0_mu,sigma_noise=sigma_noise)
    f0 = ess$f0  # Don't include predictive mean 'f0_mu' here!
    llcur = ess$llcur
    nEss = ess$nIter
    
    #### kappa,f0(t) | f0(t-1),f_mu,sigma_noise ~ 3 x MH-ESS with fixed 'nu' ####
    for(kk in 1:3){
      nu = forwardsolve(L,f0)
      # ESS for updating kappa: 'xi = ssg_fun_inv(kappa) ~ MVN(xi_kappa_mu,xi_kappa_L)'
      # Current 'xi0' (mean-zero) and ll-threshold:
      xi0 = ssg_fun_inv(kappa,kappaMax,kappaMin) - xi_kappa_mu
      llth = log(runif(1)) + llcur
      # Generate ellipse, initial proposal & define bracket:
      xi0_nu = c( xi_kappa_L%*%rnorm(length(xi0)) )
      theta = runif(1,0,2*pi)
      bracket = c(theta-2*pi,theta)
      nProp_kappa = 1
      while(TRUE){
        xi0_prop <- xi0*cos(theta) + xi0_nu*sin(theta)
        kappaProp <- ssg_fun(xi0_prop + xi_kappa_mu,kappaMax,kappaMin)
        res <- K_mu_cond(kappaProp,t,f0_prev,chol=TRUE)
        L_prop <- res$L_cond
        K <- res$K_cond
        f0_mu_prop = res$f_mu_cond
        f0_prop <- L_prop%*%nu
        llprop <- logL(f0_prop,f_mu+f0_mu_prop,sigma_noise,data=DATA[[t]])
        # Accept porposal or shrink backet?
        if(llprop > llth){
          kappa <- kappaProp
          L <- L_prop
          f0_mu <- f0_mu_prop
          f0 <- f0_prop
          llcur <- llprop
          break
        }else{
          nProp_kappa = nProp_kappa+1
          bracket[(theta>0)+1] = theta
          theta <- runif(1,bracket[1],bracket[2])
        }
      }
    }
    
    #### f0(t) | f0(t-1),kappa,f_mu,sigma_noise ~ ESS ####
    ess = ESS_f0(f0,L=L,llcur=llcur,logL=logL,nESS=3,data=DATA[[t]],f_mu=f_mu+f0_mu,sigma_noise=sigma_noise)
    f0 = ess$f0
    llcur = ess$llcur
    nEss = nEss + ess$nIter
    
    #### f_mu,sigma_noise |kappa,f0(t) ~ ESS ####
    nProp_alpha = 1
    for(kk in 1:3){
      # Current 'xi0', ll-threshold & generate ellipse:
      xi0 = ssg_fun_inv(c(link(f_mu),sigma_noise),y_max=c(lv_mu_max,sigma_noise_max)) - xi_alpha_mu
      llth <- log(runif(1)) + llcur
      xi0_nu <- c( xi_alpha_L%*%rnorm(length(xi0)) )
      # Initial proposal & define bracket:
      theta <- runif(1,0,2*pi)
      bracket <- c(theta-2*pi,theta)
      while(TRUE){
        xi0_prop <- xi0*cos(theta) + xi0_nu*sin(theta)
        f_muProp <- linkInv( ssg_fun(xi0_prop[1] + xi_alpha_mu[1],lv_mu_max) )
        sigma_noiseProp <- ssg_fun(xi0_prop[2] + xi_alpha_mu[2],sigma_noise_max)
        llprop <- logL(f0,f_muProp+f0_mu,sigma_noiseProp,data=DATA[[t]])
        # Accept porposal or shrink backet?
        if(llprop > llth){
          f_mu <- f_muProp
          sigma_noise <- sigma_noiseProp
          llcur <- llprop
          break
        }else{
          nProp_alpha = nProp_alpha+1
          bracket[(theta>0)+1] = theta
          theta <- runif(1,bracket[1],bracket[2])
        }
      }
    }
    
    #### f0(t) | f0(t-1),kappa,f_mu,sigma_noise ~ ESS ####
    ess = ESS_f0(f0,L=L,llcur=llcur,logL=logL,nESS=3,data=DATA[[t]],f_mu=f_mu+f0_mu,sigma_noise=sigma_noise)
    f0 = ess$f0
    llcur = ess$llcur
    nEss = nEss + ess$nIter
    
    # Store
    f0_states[,i] = f0 + f0_mu
    kappaStates[,i] = kappa
    f_muStates[i] = f_mu
    sigma_noiseStates[i] = sigma_noise 
    llVec[i] = llcur
    if(i%%100 == 0) print(paste('Iteration ',i,'(',t,'): #f0-ESS = ',nEss,', #kappa-SDSS = ',nProp_kappa,', #alpha-ESS = ',nProp_alpha,sep=''))
  } 
  
  # Store to STATES
  STATES[[t]] <- list(f0_states=f0_states,
                      kappaStates=kappaStates,
                      f_muStates=f_muStates,
                      sigma_noiseStates=sigma_noiseStates,
                      llVec=llVec,
                      xi_kappa_mu=xi_kappa_mu,
                      xi_kappa_K=xi_kappa_K,
                      xi_alpha_mu=xi_alpha_mu,
                      xi_alpha_K=xi_alpha_K)
}

# Add call prices and implied vols to STATES
for(t in 1:length(STATES)){
  print(t)
  f_states = STATES[[t]]$f0_states + rep(STATES[[t]]$f_muStates,each=nrow(STATES[[t]]$f0_states))
  C_states = matrix(nrow=nrow(f_states),ncol=ncol(f_states))
  IV_states = matrix(nrow=nrow(f_states),ncol=ncol(f_states))
  for(i in 1:ncol(C_states)){
    LV = matrix(link(f_states[,i]),nrow=length(DATA[[t]]$T))
    C_states[,i] = c( localVolCalls(DATA[[t]]$S0,DATA[[t]]$rf,DATA[[t]]$q,LV,DATA[[t]]$K,DATA[[t]]$T,KflatExt=DATA[[t]]$S0*Mext) )
    IV_states[,i] = c( localVolCalls(DATA[[t]]$S0,DATA[[t]]$rf,DATA[[t]]$q,LV,DATA[[t]]$K,DATA[[t]]$T,KflatExt=DATA[[t]]$S0*Mext,impVol=TRUE) )
  }
  
  # Save:
  STATES[[t]]$f_states <- f_states
  STATES[[t]]$C_states <- C_states
  STATES[[t]]$IV_states <- IV_states
}


# save(STATES,kappaMax,kappaMin,lv_mu_max,sigma_noise_max,type,file='mcmc_states_sequential.Rdata')
