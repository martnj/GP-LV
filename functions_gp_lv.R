#
# Functions for GP local volatility model
#



# Generic
rowSD <- function(x) apply(x,1,sd)

plot_IV <- function(data, fname=NA){
  
  # Plot implied volatility surface
  
  x = data$K
  y = data$T
  X = t(data$IV)
  if(!is.na(fname)) png(file=fname)
    zlim = c(0, 0.45)
    par(list(font.axis=1,font.lab=1,family="serif",cex.lab=1.4,cex.axis=1.4))
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
  if(!is.na(fname)) dev.off()
}

# Scaled sigmoid Gaussian: transform, random and density function
ssg_fun <- function(x,y_max,y_min=0) y_min + (y_max-y_min)*1/(1+exp(-x))
ssg_fun_inv <- function(y,y_max,y_min=0) log(y-y_min) - log(y_max-y)
rssg <- function(y_max,y_min=0) ssg_fun(rnorm(length(y_max)),y_max,y_min)
dssg <- function(x,x_max,x_min=0,log=FALSE,z_mu=0,z_sd=1){
  if(!log) return( dnorm(log((x-x_min)/(x_max-x)),mean=z_mu,sd=z_sd)*(x_max-x_min)/((x_max-x)*(x-x_min)) )
  return( log( dnorm(log((x-x_min)/(x_max-x)),mean=z_mu,sd=z_sd)*(x_max-x_min)/((x_max-x)*(x-x_min)) ) )
}


# MCMC samplers
ESS_f0 <- function(f0, L, llcur, logL, nESS=1, angleRange=2*pi,...){
  
  # Eliptical slice sampler for uppdating zero-mean 'f0' 
  # given fixed covariance matrix square-root 'L'.
  #
  # Log-likelihood: 'logL(f0,...)', 'llcur' is likelihood of
  # initial 'f0'.
  #
  # OBS! 'f0', is 'f' with zero mean!
  
  nIter = 1
  for(k in 1:nESS){
    
    # Set ll-threshold & generate ellipse:
    llth <- log(runif(1)) + llcur
    f0nu <- L%*%matrix(rnorm(length(f0)))
    
    # Initial proposal 'theta' and define bracket:
    theta <- runif(1, 0, angleRange)
    bracket <- c(theta - angleRange, theta)
    
    while(TRUE){
      f0prop <- f0*cos(theta) + f0nu*sin(theta)
      llprop <- logL(f0prop,...)
      
      # Accept porposal or shrink backet?
      if(llprop > llth){
        f0 <- f0prop  
        llcur <- llprop
        break
      }else{
        nIter = nIter + 1
        bracket[(theta > 0) + 1] = theta
        theta <- runif(1, bracket[1], bracket[2])
      }
    }
  }
  
  return(list(f0=f0,llcur=llcur,nIter=nIter,L=L))
}

SDESS <- function(f0, kappa, kappaMin, kappaMax, logL, llcur, data, s=0.05,K=NULL,R=NULL,LR=NULL, angleRange=2*pi, xi_mu=0, xi_L=1,...){
  
  # ESS for updating 'kappa' and 'f0' with surrogate data.
  # Log-likelihood: 'logL(f0,data,...)'
  #
  # Req. global:  'K_fun(kappa,data)' 
  #               'ssg_fun(kappa,kappaMax,kappaMin)' 
  #               'ssg_fun_inv(kappa,kappaMax,kappaMin)'
  
  N = length(f0)
  S = s*diag(N)
  
  # 1. Draw surrogate data 'g'
  # 2. Calculate K, R, LR and m for kappa
  # 3. Implied latent variables, 'eta'
  g <- f0 + sqrt(s)*rnorm(N)
  if(is.null(K)) K <- K_fun(kappa,data)
  if(is.null(R)) R <- S - S*solve(S+K,S)
  if(is.null(LR)) LR <- t(chol(R + 1e-5*min(diag(R))*diag(N)))
  m <- 1/s*R%*%g
  eta <- forwardsolve(LR,f0-m)   
  
  # 4. Elliptical slice sampler of 'kappa'
  # Likelihood-threshold & xi-current:
  xi0 <- ssg_fun_inv(kappa,kappaMax,kappaMin) - xi_mu
  th <- log(runif(1)) + llcur + dmvnorm(c(g),rep(0,N),K+S,log=T) 
  
  # Initial xi-proposal & define bracket:
  xi0_nu <- xi_L%*%rnorm(length(kappa))
  theta <- runif(1,0,angleRange)
  bracket <- c(theta-angleRange, theta)
  
  nIter = 1
  while(TRUE){
    # xi-proposal -> kappa-proposal
    xi0_prop <- xi0*cos(theta) + xi0_nu*sin(theta)
    kappaProp <- ssg_fun(xi0_prop+xi_mu,kappaMax,kappaMin)
    
    # K, LR and m for kappa-proposal -> f0-proposal
    K <- K_fun(kappaProp,data)
    R <- S - S*solve(S+K,S)
    m <- 1/s*R%*%g
    eps <- 1e-5*min(diag(R))
    LR <- t(chol(R + eps*diag(N)))
    f0prop <- LR%*%eta + m
    llprop <- logL(f0prop, data=data,...)
    
    # Accept porposal or shrink backet?
    if( (llprop + dmvnorm(c(g),rep(0,N),K+S,log=T)) > th){
      f0 <- f0prop
      kappa <- kappaProp
      llcur <- llprop
      break
    }else{
      nIter = nIter+1
      bracket[(theta>0)+1] = theta
      theta <- runif(1,bracket[1],bracket[2])
    }
  }
  return(list(f0=f0,kappa=kappa,llcur=llcur,nIter=nIter,K=K,R=R,LR=LR))
}


# GP covariance functions
covM_SE <- function(X,Y=NULL,kappa,structure='none',decomp='none',factors=FALSE,type){
  
  # Calculates the covariance matrix 'K = k(X,X)' for
  # squared-exponential kernel with parameters
  # 'kappa = (length_scales,sigma_f)'
  #
  # Allows for X with 'structure' = {'none','duplicated','grid'}
  #
  # Calculate square-root 'K=LL^T' with 'decomp' = {'none','chol','svd','eigen','inv'}
  # where option 'eigen' gives eigenvalue decomposition {Q,d} for 'K=Q diag(d) Q^T
  # and 'inv' gives inverse of 'K'.
  #
  # 'X' is [N x d] matrix with N obs of d-dimensional
  # input, 'x = (x_1,...,x_d)'. Each dimension has 
  # individual scale in 'length_scales'.
  # 
  # For 'X' with structure 'grid', a cartesian product is formed
  # from columns of 'X'. Here 'X[,i]' should contain unique
  # values for dimension 'i', and filled with 'NA' if of 
  # different lengths from other columns.
  # 
  # Optional: covariance matrix 'K = k(X,Y)' when 'Y' is [M x d].
  # Then 'K' will be [N x M] matrix. Decomposition not possible in this case.
  # 
  
  
  
  # Input control:
  if(!is.matrix(X)) X = as.matrix(X)
  if(!is.null(Y) & !is.matrix(Y)) Y = as.matrix(Y)
  kappa = c(kappa)
  
  # Extract 'kappa' -> 'length_scales' and 'sigma_f'
  length_scales = kappa[1:ncol(X)]
  if(length(kappa)>ncol(X)){ 
    sigma_f = kappa[length(kappa)]
  }else{
    sigma_f = 1
  }
  
  if(is.null(Y)){
    # Calculate covariance matrix 'K = k(X,X)'
    if(structure=='none'){
      K = matrix(0,nrow=nrow(X),ncol=nrow(X))
      for(i in 1:ncol(X)){
        if(type=='SE') K = K - 0.5/length_scales[i]^2*as.matrix(dist(X[,i]))^2
        if(type=='Mat12') K = K - 1/length_scales[i]*as.matrix(dist(X[,i]))
        if(type=='Mat32'){
          epnt = sqrt(3)/length_scales[i]*as.matrix(dist(X[,i]))
          K = K - epnt + log(1+epnt)
        }
        if(type=='Mat52'){
          epnt = sqrt(5)/length_scales[i]*as.matrix(dist(X[,i]))
          K = K - epnt + log(1+epnt+epnt^2/3)
        }
      }
      K = sigma_f^2*exp(K)
      if(decomp=='chol'|decomp=='svd'){
        K = LfromK(K,method=decomp)
      }else if(decomp=='eigen'){
        eps = 1e-5*K[1,1]
        res = eigen( K+eps*diag(nrow(K)) )
        K = list(Q=res$vectors,d=res$values)
      }else if(decomp=='inv'){
        eps = 1e-5*K[1,1]
        K = solve( K+eps*diag(nrow(K)) )
      }
    }else if(structure=='duplicated'){
      K = matrix(0,nrow=nrow(X),ncol=nrow(X))
      for(i in 1:ncol(X)){
        if(type=='SE') K = K - 0.5/length_scales[i]^2*dist_duplicated(X[,i])^2
        if(type=='Mat12') K = K - 1/length_scales[i]*dist_duplicated(X[,i])
        if(type=='Mat32'){
          epnt = sqrt(3)/length_scales[i]*dist_duplicated(X[,i])
          K = K - epnt + log(1+epnt)
        }
        if(type=='Mat52'){
          epnt = sqrt(5)/length_scales[i]*dist_duplicated(X[,i])
          K = K - epnt + log(1+epnt+epnt^2/3)
        }
      }
      K = sigma_f^2*exp(K)
      if(decomp=='chol'|decomp=='svd'){
        K = LfromK(K,method=decomp)
      }else if(decomp=='eigen'){
        eps = 1e-5*K[1,1]
        res = eigen( K+eps*diag(nrow(K)) )
        K = list(Q=res$vectors,d=res$values)
      }else if(decomp=='inv'){
        eps = 1e-5*K[1,1]
        K = solve( K+eps*diag(nrow(K)) )
      }
    }else if(structure=='grid'){
      if(decomp=='chol'|decomp=='svd'){
        L = 1
        for(i in 1:ncol(X)){
          L_i = covM_SE(X[!is.na(X[,i]),i],kappa=length_scales[i],decomp=decomp,type=type)
          L = kronecker(L_i , L)
        }
        K = sigma_f*L
      }else if(decomp=='eigen'){
        if(!factors){
          Q = 1
          d = 1
          for(i in 1:ncol(X)){
            res = covM_SE(X[!is.na(X[,i]),i],kappa=length_scales[i],decomp='eigen',type=type)
            Q = kronecker(res$Q , Q)
            d = kronecker_diag(res$d , d)
          }
          K = list(Q=Q,d=sigma_f^2*d)
        }else{
          Q <- list()
          d <- list()
          for(i in 1:ncol(X)){
            res = covM_SE(X[!is.na(X[,i]),i],kappa=length_scales[i],decomp='eigen',type=type)
            Q[[i]] = res$Q
            d[[i]] = sigma_f^(2/ncol(X))*res$d
          }
          K = list(Q_list=Q,d_list=d)
        }
      }else if(decomp=='inv'){
        if(!factors){
          K_inv = 1
          for(i in 1:ncol(X)){
            K_inv = kronecker(covM_SE(X[!is.na(X[,i]),i],kappa=length_scales[i],decomp='inv',type=type) , K_inv)
          }
          K = 1/sigma_f^2*K_inv
        }else{
          K_inv <- list()
          for(i in 1:ncol(X)){
            K_inv[[i]] = 1/sigma_f^(2/ncol(X))*covM_SE(X[!is.na(X[,i]),i],kappa=length_scales[i],decomp='inv',type=type)
          }
          K = list(K_inv_list = K_inv)
        }
      }else{
        K = 1
        for(i in 1:ncol(X)){
          K_i = covM_SE(X[!is.na(X[,i]),i],kappa=length_scales[i],type=type)
          K = kronecker(K_i , K)
        }
        K = sigma_f^2*K
      }
    }else{
      K = NA
    }
  }else{
    
    # Calculate covariance matrix 'K = k(X,Y)'
    # Note: no decomposition possible in this case!
    if(decomp!='none') print('No decomposition possible for K = k(X,Y)!')
    if(ncol(X)!=ncol(Y)) stop(paste('Different dimension: x is',ncol(X),"while y is",ncol(Y),"."))
    if(structure=='none' | structure=='duplicated'){
      K = matrix(0,nrow=nrow(X),ncol=nrow(Y))
      for(i in 1:ncol(X)){
        if(type=='SE') K = K - 0.5/length_scales[i]^2*dist_duplicated(v=X[,i],u=Y[,i])^2
        if(type=='Mat12') K = K - 1/length_scales[i]*dist_duplicated(v=X[,i],u=Y[,i])
        if(type=='Mat32'){
          epnt = sqrt(3)/length_scales[i]*dist_duplicated(v=X[,i],u=Y[,i])
          K = K - epnt + log(1+epnt)
        }
        if(type=='Mat52'){
          epnt = sqrt(5)/length_scales[i]*dist_duplicated(v=X[,i],u=Y[,i])
          K = K - epnt + log(1+epnt+epnt^2/3)
        }
      }
      K = sigma_f^2*exp(K)
    }else if(structure=='grid'){
      K = 1
      for(i in 1:ncol(X)){
        K_i = covM_SE(X=X[!is.na(X[,i]),i],Y=Y[!is.na(Y[,i]),i],kappa=length_scales[i],type=type)
        K = kronecker(K_i , K)
      }
      K = sigma_f^2*K
    }else{
      K = NA
    }
  }
  return(K)
}

dist_duplicated <- function(v,u=NULL,fun=function(x) x){
  #
  # Computes distance matrix 'D = D(v,v)' of vector 'v'
  # and applies element-wise scale/power/function 'fun(D)'. 
  # Efficiency gain when 'v' contains duplicates.
  # Same result as 'fun( as.matrix(dist(v)) )
  #
  # Also for 'cross-distance', 'D = D(v,u)' 
  
  n_v = length(v)
  v_unique = unique(v)
  n_v_unique = length(v_unique)
  if(!is.null(u)){
    n_u = length(u)
    u_unique = unique(u)
    n_u_unique = length(u_unique)
    D = fun( matrix(abs(rep(v_unique,n_u_unique) - rep(u_unique,each=n_v_unique)),nrow=n_v_unique,ncol=n_u_unique) )
    if(n_v_unique==n_v & n_u_unique==n_u) return(D)
    M_v = 1*matrix(rep(v_unique,each=n_v) == rep(v,times=n_v_unique),nrow=n_v)
    M_u = 1*matrix(rep(u_unique,each=n_u) == rep(u,times=n_u_unique),nrow=n_u)
    return(M_v%*%D%*%t(M_u))
  }
  D = fun( as.matrix(dist(v_unique)) )
  if(n_v_unique==n_v) return(D)
  M = 1*matrix(rep(v_unique,each=n_v) == rep(v,times=n_v_unique),nrow=n_v)
  return(M%*%D%*%t(M))
}

LfromK <- function(K,method="svd",eps=NULL){
  if(method=="svd"){
    s = svd(K)
    L = s$u%*%diag(sqrt(s$d))
  }else{
    # For stability of inverse of K
    if(is.null(eps)) eps = 1e-5*K[1,1]
    L = t( chol(K + eps*diag(nrow(K))) )
  }
  return(L)
}

# Prediction equations for GP prior
f_cond_mu_K_Xmat <- function(X_star_unit,X_obs_unit,f_obs,kappa,f_mu=0,L=NULL,v=NULL,vec_output=FALSE,verbose=TRUE,type){
  
  # Predictive equations for f* = f(X*) given obs. f = f(X) i.e.
  # mean m* and cov K* based on kernel of 'type'
  # with 'kappa = (length_scales,sigma_f)'. 
  #
  # 'X_star_unit' is X* scaled to unit hypercube 
  # 'X_obs_unit' is X scaled to unit hypercube
  #
  # 'vec_output=TRUE' outputs c(m*,c(K*),length(m*)) for internal use.
  
  # VECTOR FUNCTION for Gaussian mixture.
  if(is.matrix(kappa)){
    if(ncol(kappa) != ncol(f_obs)) stop('f_obs and kappa are matrices of non-congruent dimensions.')
    #
    # Requires:
    # f_obs = [f_obs_1,f_obs_2,...]  [n X nMCMC] matrix
    # kappa = [kappa_1,kappa_2,...]  [3 X nMCMC] matrix
    # f_mu = [f_mu_1,f_mu_2,...]     nMCMC vector
    
    if(verbose) print('Note: f_star_mu will include mean f_mu. Assumes input f_obs is zero-mean f0.')
    
    myFun <- function(x,n){
      # n = length(f_obs)
      # x = c(f_obs,kappa)
      # y = c(m*,c(K*),length(m*))
      y = f_cond_mu_K_Xmat(X_star_unit=X_star_unit, X_obs_unit=X_obs_unit, f_obs=x[1:n], kappa=x[-(1:n)],f_mu=0,vec_output=TRUE,type=type)
      return(y)
    }
    
    # Calculate moments of Gaussian mixture
    if(length(f_mu)==1) f_mu = rep(f_mu,ncol(kappa))
    X = rbind(f_obs,kappa)
    n = nrow(f_obs)
    f_star_mu = 0
    K_tmp = 0
    for(i in 1:ncol(X)){
      y = myFun(X[,i],n)
      n_star = y[length(y)]
      f_star_mu = f_star_mu + y[1:n_star] + f_mu[i]
      K_tmp = K_tmp + matrix(y[-c(1:n_star,length(y))],nrow=n_star) + matrix(y[1:n_star]+f_mu[i])%*%matrix(y[1:n_star]+f_mu[i],nrow=1)
    }
    f_star_mu = f_star_mu/ncol(X)
    f_star_K = K_tmp/ncol(X) - matrix(f_star_mu)%*%matrix(f_star_mu,nrow=1)
    
    return(list(f_star_mu=f_star_mu,f_star_K=f_star_K,f_mu=mean(f_mu)))
  }
  
  # ORIGINAL FUNCTION
  # Input control:
  kappa = as.numeric(kappa)
  
  # 1. Cholesky of K = K(X,X) + jitter
  if(is.null(L)){
    if(anyNA(X_obs_unit)){ 
      L = covM_SE(X_obs_unit,kappa=kappa,structure='grid',decomp='chol',type=type) 
    }else{ 
      L = covM_SE(X_obs_unit,kappa=kappa,decomp='chol',type=type) 
    }
  }
  
  # 2. Predictive mean and variance
  if(anyNA(X_star_unit)){ 
    K_star = covM_SE(X_star_unit,kappa=kappa,structure='grid',type=type)
  }else{
    K_star = covM_SE(X_star_unit,kappa=kappa,type=type)
  }
  if(anyNA(X_star_unit) & anyNA(X_obs_unit)){ 
    K_star_obs = covM_SE(X_star_unit,X_obs_unit,kappa=kappa,structure='grid',type=type)
  }else{
    if(anyNA(X_obs_unit)){
      if(ncol(X_obs_unit)==2) X_obs_unit = as.matrix(expand.grid(X_obs_unit[!is.na(X_obs_unit[,1]),1],X_obs_unit[!is.na(X_obs_unit[,2]),2]))
      if(ncol(X_obs_unit)==3) X_obs_unit = as.matrix(expand.grid(X_obs_unit[!is.na(X_obs_unit[,1]),1],X_obs_unit[!is.na(X_obs_unit[,2]),2],X_obs_unit[!is.na(X_obs_unit[,3]),3]))
    }
    if(anyNA(X_star_unit)){
      if(ncol(X_star_unit)==2) X_star_unit = as.matrix(expand.grid(X_star_unit[!is.na(X_star_unit[,1]),1],X_star_unit[!is.na(X_star_unit[,2]),2]))
      if(ncol(X_star_unit)==3) X_star_unit = as.matrix(expand.grid(X_star_unit[!is.na(X_star_unit[,1]),1],X_star_unit[!is.na(X_star_unit[,2]),2],X_star_unit[!is.na(X_star_unit[,3]),3]))
    }
    K_star_obs = covM_SE(X_star_unit,X_obs_unit,kappa=kappa,type=type)
  }  
  alpha = forwardsolve(L,f_obs-f_mu)
  alpha = backsolve(t(L),alpha)
  f_star_mu = as.numeric(f_mu + K_star_obs%*%alpha)
  if(is.null(v)) v = forwardsolve(L,t(K_star_obs))
  f_star_K = K_star - t(v)%*%v
  rownames(f_star_K) <- colnames(f_star_K) <- NULL
  
  # Output
  if(vec_output){ return(c(f_star_mu,c(f_star_K),length(f_star_mu))) }
  else{ return(list(f_star_mu=f_star_mu,f_star_K=f_star_K,X_star_unit=X_star_unit,L=L,v=v)) }
}

# Pricing local vol functions

localVolCalls <- function(S0,rf,q,LV,Kgrid,Tgrid,theta=0.5,impVol=FALSE,initialStep="extrapol",Keval=NA,Teval=NA,Cinit=NULL,c0=NULL,KflatExt=NULL,triDiag=TRUE){
  
  #    GRID:
  # Tgrid: 0=Tmin,...TL=Tmax
  # Kgrid: Kmin=K0,...,KI=Kmax
  # LV is (L+1)x(I+1) matrix with local-vol values.
  # theta=1 gives implicit FD, theta=0 gives explicit FD.
  #
  # Cinit: (optional) boundary values, Cinit = [Cleft,Cright]
  # i.e. two stacked (L+1)x1 vectors.
  #
  # Keval/Teval: use to get LV values at a subset of Kgrid x Tgrid.
  # Keval x Teval must be inluded in Kgrid x Tgrid!
  #
  # 'q' and 'rf' can be scalars or vector of same length as 'Tgrid'
  #
  # Update: 'KflatExp' for flat extrapolation in the wings
  
  
  # Extend grid in K-dimension, only if 'c0' and 'Cinit' = NULL
  extendKgrid = FALSE
  extendKgridR = FALSE
  extendKgridFlat = FALSE
  
  if(is.null(Cinit) & is.null(c0)){
    if(!is.null(KflatExt)){
      # Flat extrapolation outside 'Kgrid' from 'KflatExt'
      resExt = extrapolLV(LV=LV,Kin=Kgrid,Kout=KflatExt)
      LV = resExt$LV
      Kgrid = resExt$Kout
      extendKgridFlat = TRUE
    }else{
      # Flat extrapolation 1-step outside 'Kgrid' 
      # Possible since scheme uses LV[,2:(end-1)] 
      if(!any(is.na(LV[,1]))){
        if( (2*Kgrid[1]-Kgrid[2])>=0){
          Kgrid = c(Kgrid[1]-(Kgrid[2]-Kgrid[1]),Kgrid)
          LV = cbind(LV[,1],LV)
          extendKgrid = TRUE
        }
        Kgrid = c(Kgrid, Kgrid[length(Kgrid)]+Kgrid[length(Kgrid)]-Kgrid[length(Kgrid)-1])
        LV = cbind(LV,LV[,ncol(LV)])
        extendKgridR = TRUE
      }
    }
  }
  
  # Define PDE coeffs:
  if(length(q)==1) q = rep(q,length(Tgrid))
  if(length(rf)==1) rf = rep(rf,length(Tgrid))
  alpha <- function(l,i) -q[l+1]
  beta <- function(l,i) -(rf-q)[l+1]*Kgrid[i+1]
  gamma <- function(l,i) 0.5*Kgrid[i+1]^2*LV[l+1,i+1]^2
  
  # Initil value 'c0' at Tmin (only if 'c0 = NULL')
  extendTgrid = FALSE
  if(is.null(c0)){
    c0 = (S0-Kgrid)*(S0>Kgrid)
    
    # IF first maturity is not zero: step to Tmin with implicit 
    # step OR extrapolate LV surface flat:
    if(Tgrid[1]>0){
      if(initialStep == "implicit"){
        TgridTmp = c(0,Tgrid[1])
        Ctmp = matrix(0,nrow=2,ncol=length(Kgrid))
        Ctmp[1,] = c0
        if(is.null(Cinit)){ 
          Ctmp[,1] = S0*exp(-q[1]*TgridTmp) - Kgrid[1]*exp(-rf[1]*TgridTmp)
        }else{ 
          Ctmp[2,1] = Cinit[1,1] 
          Ctmp[2,ncol(Ctmp)] = Cinit[1,2]
        }
        gammaTmp <- function(l,i) 0.5*Kgrid[i+1]^2*LV[l,i+1]^2
        Ctmp = CNnu(Ctmp,Kgrid,TgridTmp,alpha,beta,gammaTmp,theta=1,triDiag=triDiag)
        c0 = Ctmp[2,]
      }else{
        Tgrid = c(0,Tgrid)
        q = c(q[1],q)
        rf = c(rf[1],rf)
        LV = rbind(LV[1,],LV)
        extendTgrid = TRUE
      }
    }
  }
  
  # Set-up initial C matrix
  C = matrix(0,nrow=length(Tgrid),ncol=length(Kgrid))
  C[1,] = c0
  if(is.null(Cinit)){
    C[,1] = S0*exp(-q*Tgrid) - Kgrid[1]*exp(-rf*Tgrid)
  }else{
    # 'Cinit' given as input
    if(extendTgrid){
      C[2:nrow(C),1] = Cinit[,1]
      C[2:nrow(C),ncol(C)] = Cinit[,2]
    }else{
      C[,1] = Cinit[,1]
      C[,ncol(C)] = Cinit[,2]
    }
  } 
  
  # Crank-Nicolson theta-scheme (with non-uniform grid):
  C = CNnu(C = C,Kgrid = Kgrid,Tgrid = Tgrid,alpha = alpha,beta = beta,gamma = gamma,theta = theta,triDiag = triDiag)
  
  # Remove values from extended K-grid and T-grid:
  
  if(extendKgridFlat){
    C = C[,resExt$KinIdx]
    Kgrid = Kgrid[resExt$KinIdx]
  }else{
    if(extendKgrid){
      C = C[,-1]
      Kgrid = Kgrid[-1]
    }
    if(extendKgridR){
      C = C[,-ncol(C)]
      Kgrid = Kgrid[-length(Kgrid)]
    }
  }
  if(extendTgrid){
    Tgrid = Tgrid[-1]
    q = q[-1]
    rf = rf[-1]
    C = C[-1,]
  }
  
  # OUTPUT: values @ 'Teval x Keval'
  if( !is.na(Teval[1]) & !identical(Tgrid,Teval) ){
    idx = findIdx(targetValues=Teval,population=Tgrid)
    C = C[idx,]
    Tgrid = Tgrid[idx]
  }
  if( !is.na(Keval[1]) & !identical(Kgrid,Keval) ){
    idx = findIdx(targetValues=Keval,population=Kgrid)
    if(is.matrix(C)){
      C = C[,idx]
    }else{
      C = C[idx]
    }
    Kgrid = Kgrid[idx]
  }
  if(impVol){
    C = BSimpvol(P=C,S=S0,K=Kgrid,rf=rf,ttm=Tgrid,type="call",q=q)
  }
  return(C)
}

extrapolLV <- function(LV,Kin,Kout,method="flatWings"){
  
  if(method=="flatWings"){
    m = sum(Kout < Kin[1])
    M = sum(Kout > Kin[length(Kin)])
    LV = cbind(matrix(rep(LV[,1],m),nrow=nrow(LV)), LV, matrix(rep(LV[,ncol(LV)],M),nrow=nrow(LV)))
    Kout = c(Kout[Kout < Kin[1]], Kin, Kout[Kout > Kin[length(Kin)]])
    KinIdx = !(Kout < Kin[1] | Kout > Kin[length(Kin)])
    return(list(LV=LV,Kout=Kout,KinIdx=KinIdx))
  }
  return(list(LV=LV,Kout=Kin,KinIdx=!logical(length(Kin))))
}


# Crank-Nicholson solver

CNnu <- function(C,Kgrid,Tgrid,alpha,beta,gamma,theta=0.5,triDiag=TRUE){
  
  # CN scheme that allows for non-uniform Kgrid 
  #
  #    GRID:
  # Tgrid: Tmin=T0,...TL=Tmax
  # Kgrid: Kmin=K0,...,KI=Kmax
  #
  #
  #    INITIAL FUNCTION MATRIX:
  # C is (L+1)x(I+1) matrix with function values.
  # Initial: C[1,] = C(Tmin,K)
  # Boundary: C[,1] = C(T,Kmin) and C[,I+1] = C(T,Kmax)
  # 
  # alpha,beta,gamma are functions of grid-node (l,i)
  #
  
  # Definitions
  L = length(Tgrid)-1
  I = length(Kgrid)-1
  dK <- function(i) Kgrid[i+2] - Kgrid[i+1]
  sigma <- function(i) dK(i-1)/(dK(i)+dK(i-1))
  Sigma <- function(i) 2/(dK(i)+dK(i-1))
  
  if(theta!=0){
    # Matrix ACN:
    dT1 <- function(l) Tgrid[l+1]-Tgrid[l]
    a1 <- function(l,i) beta(l,i)*(1-sigma(i))/dK(i-1) - gamma(l,i)*Sigma(i)/dK(i-1)
    b1 <- function(l,i) 1/(theta*dT1(l))-alpha(l,i)+beta(l,i)*sigma(i)/dK(i)-beta(l,i)*(1-sigma(i))/dK(i-1)+gamma(l,i)*Sigma(i)*(1/dK(i)+1/dK(i-1))
    c1 <- function(l,i) -beta(l,i)*sigma(i)/dK(i) - gamma(l,i)*Sigma(i)/dK(i)
    if(triDiag){
      # Give diagonals of ACN for tri-daig solver (2017-11-09)
      ACNlowDiag <- function(l) theta*dT1(l)*a1(l,1:(I-1))[2:(I-1)]
      ACNtopDiag <- function(l) theta*dT1(l)*c1(l,1:(I-1))[1:(I-2)]
      ACNdiag <- function(l) theta*dT1(l)*b1(l,1:(I-1))
    }else{
      # Original version, gives ACN matrix
      ACN <- function(l){  
        acn = cbind(diag(a1(l,1:(I-1)))[,2:(I-1)],matrix(0,nrow=(I-1),ncol=1)) +
          diag(b1(l,1:(I-1))) + 
          cbind(matrix(0,nrow=(I-1),ncol=1),diag(c1(l,1:(I-1)))[,1:(I-2)])
        return(theta*dT1(l)*acn)
      }
    }
    # Boundary:
    aa <- function(l){
      res = matrix(0,nrow=(I-1),ncol=1)
      res[1,1] = theta*dT1(l)*a1(l,1)*C[l+1,1]
      res[I-1,1] = theta*dT1(l)*c1(l,I-1)*C[l+1,I+1]
      return(res)
    }
  }
  
  if(theta!=1){
    # Matrix BCN:
    dT0 <- function(l) Tgrid[l+2]-Tgrid[l+1]
    a0 <- function(l,i) -beta(l,i)*(1-sigma(i))/dK(i-1) + gamma(l,i)*Sigma(i)/dK(i-1)
    b0 <- function(l,i) 1/((1-theta)*dT0(l))+alpha(l,i)-beta(l,i)*sigma(i)/dK(i)+beta(l,i)*(1-sigma(i))/dK(i-1)-gamma(l,i)*Sigma(i)*(1/dK(i)+1/dK(i-1))
    c0 <- function(l,i) beta(l,i)*sigma(i)/dK(i) + gamma(l,i)*Sigma(i)/dK(i)
    BCN <- function(l){  
      bcn = cbind( diag(a0(l,1:(I-1)))[,2:(I-1)],matrix(0,nrow=(I-1),ncol=1) ) +
        diag( b0(l,1:(I-1)) ) + 
        cbind( matrix(0,nrow=(I-1),ncol=1),diag(c0(l,1:(I-1)))[,1:(I-2)] )
      return((1-theta)*dT0(l)*bcn)
    }
    # Boundary:
    bb <- function(l){
      res = matrix(0,nrow=(I-1),ncol=1)
      res[1,1] = (1-theta)*dT0(l)*a0(l,1)*C[l+1,1]
      res[I-1,1] = (1-theta)*dT0(l)*c0(l,I-1)*C[l+1,I+1]
      return(res)
    }
  }
  
  # Iterate the scheme!
  if(theta==0){
    # Explicit scheme
    for(l in 0:(L-1)){
      C[l+2,2:I] = BCN(l)%*%matrix(C[l+1,2:I],nrow=I-1,ncol=1) + bb(l)
    } 
    return(C)
  }
  if(theta==1){
    # Implicit scheme
    if(triDiag){
      for(l in 0:(L-1)){
        C[l+2,2:I] = Solve.tridiag(ACNlowDiag(l+1),ACNdiag(l+1),ACNtopDiag(l+1),matrix(C[l+1,2:I],nrow=I-1,ncol=1)-aa(l+1))
      } 
    }else{
      for(l in 0:(L-1)){
        C[l+2,2:I] = solve(ACN(l+1),matrix(C[l+1,2:I],nrow=I-1,ncol=1)-aa(l+1))
      } 
    }
    return(C)
  }
  # CN scheme:
  if(triDiag){
    for(l in 0:(L-1)){
      c = matrix(C[l+1,2:I],nrow=I-1,ncol=1)
      C[l+2,2:I] = Solve.tridiag(ACNlowDiag(l+1),ACNdiag(l+1),ACNtopDiag(l+1), BCN(l)%*%c+bb(l)-aa(l+1))
    }
  }else{
    for(l in 0:(L-1)){
      c = matrix(C[l+1,2:I],nrow=I-1,ncol=1)
      C[l+2,2:I] = solve( ACN(l+1),BCN(l)%*%c+bb(l)-aa(l+1) )
    } 
  }
  return(C)
}

# Black-Scholes

BSimpvol <- function(P,S,K,rf,ttm,type="call",q=0,interval=c(-10,10),tol=.Machine$double.eps,maxiter=1000){
  
  # Calculate BS implied volatility with 'uniroot'.
  #
  # Inputs:
  #   'P' scalar / vector / matrix
  #   'S, K, ttm, rf, q' scalars or same as 'P'
  #
  #   If 'P' is matrix, 'K, ttm, rf, q' may
  #   also be vector(s) of same lengt as nrow(P) or ncol(P).
  
  # Vectorised function for non-scalar inputs
  if(length(P)>1){
    if(is.matrix(P)){
      # Case 'P' is matrix: transform to numeric and use 'P' is numeric!
      # First check dimension of K,ttm,rf,q
      if(1<length(K) & length(K)<length(P)){
        if(length(K)==ncol(P)) K = matrix(rep(K,nrow(P)),nrow=nrow(P),byrow=T)
        else K = matrix(rep(K,ncol(P)),nrow=nrow(P))
      }
      if(1<length(ttm) & length(ttm)<length(P)){
        if(length(ttm)==nrow(P)) ttm = matrix(rep(ttm,ncol(P)),nrow=nrow(P))
        else ttm = matrix(rep(ttm,nrow(P)),nrow=nrow(P),byrow=T)
      }
      if(1<length(rf) & length(rf)<length(P)){
        if(length(rf)==nrow(P)) rf = matrix(rep(rf,ncol(P)),nrow=nrow(P))
        else rf = matrix(rep(rf,nrow(P)),nrow=nrow(P),byrow=T)
      }
      if(1<length(q) & length(q)<length(P)){
        if(length(q)==nrow(P)) q = matrix(rep(q,ncol(P)),nrow=nrow(P))
        else q = matrix(rep(q,nrow(P)),nrow=nrow(P),byrow=T)
      }
      impvol = BSimpvol(P=as.numeric(P),S=as.numeric(S),K=as.numeric(K),rf=as.numeric(rf),ttm=as.numeric(ttm),type=type,q=as.numeric(q),interval=interval,tol=tol,maxiter=maxiter)
      impvol = matrix(impvol,nrow=nrow(P))
      rownames(impvol) <- rownames(P)
      colnames(impvol) <- colnames(P)
      return(impvol)
    }else{
      # Case 'P' is numeric
      X = cbind(P,S,K,rf,ttm,q)
      fun <- function(x) BSimpvol(P=x[1],S=x[2],K=x[3],rf=x[4],ttm=x[5],type=type,q=x[6],interval=interval,tol=tol,maxiter=maxiter)
      return( apply(X,1,fun) )
    }
  }
  
  # Original function for scalar inputs:
  if(!is.finite(P)){
    return(NA)
  }
  if(P <= tol){
    return(0)
  }
  objective <- function(sigma) P -  BSformula(S=S,K=K,rf=rf,ttm=ttm,sigma=sigma,type=type,q=q)
  return(uniroot(f=objective,interval=interval,tol=tol,maxiter=maxiter)$root)
}

BSformula <- function(S,K,rf,ttm,sigma,type="call",q=0){
  
  # Black-Scholes formula for call and put option
  #
  # Update: for structure length(ttm) X length(K)
  # i.e. 'sigma' matrix of dimension length(ttm) X length(K) 
  
  if(length(ttm)>1 & length(K)>1) K = matrix(K,nrow=length(ttm),ncol=length(K),byrow=T)
  
  Fwd <- S*exp( (rf-q)*ttm )
  d1 <- ( log(S/K)+(rf-q+0.5*sigma^2)*ttm )/(sigma*sqrt(ttm))
  d2 <- d1 - sigma*sqrt(ttm)
  if(type=="call") return( exp(-rf*ttm)*(Fwd*pnorm(d1)-K*pnorm(d2)) )
  if(type=="call_delta") return( exp(-q*ttm)*pnorm(d1) )
  if(type=="put") return( exp(-rf*ttm)*(K*pnorm(-d2)-Fwd*pnorm(-d1)) )
  if(type=="vega") return( exp(-q*ttm)*S*dnorm(d1)*sqrt(ttm) )
  return(NULL)
}
