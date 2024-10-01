# mvmr_local.R
# Two stage hard thresholding algorithm applied in MVMR individual data when IVs vote for each other
# Select the valid IVs and estimate multiple treatments effects. 
#
# Usage: [VHat,betaHat,ci,betaSD] = mvmr_local(Y,D,Z,X,intercept,alpha,a0,V,S,oracle,V,S)
#
# Y: nx1 outcome vector
# D: nxq treatment matrix
# Z: nxpz candidate IVs
# X: nxpx corvariates
# intercept: whether or not introduce a intercept in linear regression
# alpha: confidence level
# a0: tuning parameter to adjust threshoulding level in first stage; -a0<= beta <= a0
# oracle: whether or not know the the true S and V
# V: true valid IVs

#

# VHat: estimated valid IVs
# betaHat: estimated treatments effects
# ci: alpha-level confidence intervals for treatments effects
# betaSD: standard deviation of estimated treatments effects



mvmr_local <- function(Y,D,Z,X,a0=1,oracle=FALSE,V=NULL,intercept=FALSE,alpha=0.05,vote=FALSE){
  
  #====== Check and Clean Input Type =====================
  # Check Y
  stopifnot(!missing(Y),(is.numeric(Y) || is.logical(Y)),(is.matrix(Y) || is.data.frame(Y)) && ncol(Y) == 1)
  stopifnot(all(!is.na(Y)))
  # Check D
  stopifnot(!missing(D),(is.numeric(D) || is.logical(D)),(is.matrix(D) || is.data.frame(D)) && ncol(D) >= 1)
  stopifnot(all(!is.na(D)))
  # Check Z
  stopifnot(!missing(Z),(is.numeric(Z) || is.logical(Z)),is.matrix(Z))
  stopifnot(all(!is.na(Z)))
  #scale Z
  # Y <- scale(Y,center = TRUE,scale = FALSE)
  # D <- scale(D,center = TRUE,scale = FALSE)
   Z <- scale(Z, center = TRUE, scale = TRUE)
  # Check dimesions
  stopifnot(length(Y) == length(D)/ncol(D), length(Y) == nrow(Z))
  # Check X, if present
  if(!missing(X) && !is.null(X)) {
    stopifnot((is.numeric(X) || is.logical(X)),is.matrix(X) && nrow(X) == nrow(Z))
    stopifnot(all(!is.na(X)))
    W = cbind(Z,X)
  } else {
    W = Z; X = NULL
  }
  if(intercept) {
    W = cbind(W,1)
  }
  # All the other argument
  stopifnot(is.logical(intercept))
  stopifnot(is.numeric(alpha),length(alpha) == 1,alpha <= 1,alpha >= 0)
  stopifnot(is.numeric(a0),length(a0) == 1)
  
  #========= Derive Inputs ==========================
  n=length(Y); pz=ncol(Z); p=ncol(W); q=length(D)/n
  # ITT effects (OLS Estimation)
  gamHat=matrix(nrow =pz, ncol = q)
  EHat=matrix(nrow = n,ncol = q)
  qrW = qr(W)
  GamHat= qr.coef(qrW,Y)[1:pz]
  uHat=qr.resid(qrW,Y)
  for(j in 1:q){
    gamHat[,j]=qr.coef(qrW,D[,j])[1:pz]
    EHat[,j]=qr.resid(qrW,D[,j])
  }
  
  # compute the covariance of W,delta(residual of D regress on Z)
  SigHat=1/n*t(W)%*%W
  ThetaHat=1/(n-p)*t(cbind(uHat,EHat))%*%cbind(uHat,EHat)
  # may change to 1/n
  
  if(oracle==TRUE){
    stopifnot(!missing(S) && !missing(V))
    VHat=V;  BHat=NULL
  }else{
    # Estimate Valid IVs
    if(vote==TRUE){
      SetHats = Vote.VHat(GamHat,gamHat,SigHat,ThetaHat,n,a0)
    }
    else{
      SetHats = Grid.VHat(GamHat,gamHat,SigHat,ThetaHat,n,a0)
    }
    VHat = SetHats$VHat
    BHatIndex = SetHats$BHatIndex
  }
  
  # ====== Obtain divw est, se, and ci==============================
  SigHat.V=solve((solve(SigHat))[VHat,VHat])
  denom=t(gamHat[VHat,])%*%SigHat.V%*%gamHat[VHat,]-length(VHat)*ThetaHat[-1,-1]/n
  betaHat=solve(denom,t(gamHat[VHat,])%*%SigHat.V%*%GamHat[VHat])
  H=as.numeric(sqrt(t(c(1,-betaHat))%*%ThetaHat%*%c(1,-betaHat)))*solve(denom)/n
  sd=sqrt(diag(as.matrix(H)))
  CI=cbind(betaHat-qnorm(1-alpha/2)*sd,betaHat+qnorm(1-alpha/2)*sd) # confidence interval
  
  return(list( VHat=VHat, betaHat=betaHat, CI=CI, betaSD=sd, BHatIndex=BHatIndex))
}



Grid.VHat <- function(GamHat,gamHat,SigHat,ThetaHat,n,a0=1) {
  # Include intercept
  pz=length(GamHat)
  q=length(gamHat)/pz
  p=nrow(SigHat)
  
  
  #=========== estimate V =================================================
  b=seq(-a0,a0,q*a0/sqrt(n))# candidate values
  candidates=rep(list(b), q)
  B=expand.grid(candidates)
  VFlag=matrix(0,nrow=pz,ncol=length(b)^q) # Voting matrix
  SigInv=solve(SigHat)
  for(m in 1:length(b)^q){
    beta_m=as.numeric(B[m,])
    pi_m=GamHat-gamHat%*%beta_m
    sig_m=as.numeric(sqrt(t(c(1,-beta_m))%*%ThetaHat%*%c(1,-beta_m)))
    thresh_m=sig_m*sqrt(a0*log(n)/n)*sqrt(diag(SigInv)[1:pz]) # threshold of pi_bm
    VmHat=which(abs(pi_m)<=thresh_m) # valid set based on bm: V_bm
    VFlag[,m]=1:pz %in% VmHat
  }
  
  
  
  
  
  
  ## voting
  
  maxIV=max(apply(VFlag,2,sum))  #the maximum number of valid IVs
  BHatIndex=which(apply(VFlag,2,sum)==maxIV) #give all couples with the maximum number of valid IVs
  
  
  #select
  VHat=(1:pz)[VFlag[,which.max(apply(VFlag,2,sum))]>0]
  
  
  # Warning check: more than one BHatIndex
  if(length(BHatIndex)>1){
    max_columns <- VFlag[, BHatIndex]
    identical_columns <- apply(max_columns, 2, function(col) all(col == max_columns[, 1]))
    all_identical <- all(identical_columns)
    if(!all_identical){
      warning("VHat Warning: More than one cluster has the largest cluster size. This may be due to identification condition not being met. \n")
      #VHat=(1:pz)[apply(VFlag[,BHatIndex],1,sum)>0]
    }
    }
  
  # Error check: not enough IVs
  if(length(VHat) < q){
    warning("VHat Warning: No enough valid IVs estimated. This may be due to weak IVs or identification condition not being met. Use more robust methods.\n")
    warning("Defaulting to all relevant IVs being valid.\n")
    VHat = 1:pz
  }
  
  return(list(VHat = VHat, BHatIndex=BHatIndex))
}



Vote.VHat <- function(GamHat,gamHat,SigHat,ThetaHat,n,a0=1) {
  # Include intercept
  pz=length(GamHat)
  q=length(gamHat)/pz
  p=nrow(SigHat)
  # 
  # #========= estimate S =================================================
  # gamThresh=sqrt(diag(solve(SigHat))[1:pz] %*% t(diag(ThetaHat[-1,-1])))*sqrt(a0*log(n)/n)
  # #threshold of gam
  # SFlag=(abs(gamHat)>gamThresh)
  # SHat=which(apply(as.matrix(SFlag),1,sum)>0)
  # 
  # # Error check
  # if(length(SHat) < q){
  #   warning("VHat Warning: No enough relevant IVs estimated. This may be due to weak IVs or identification condition not being met. Use more robust methods.\n")
  #   warning("Defaulting to all IVs being relevant.\n")
  #   SHat = 1:pz
  # }
  
  
  #=========== estimate V =================================================
  B=combn(1:pz,q) # B is a q*C(s,q) matrix
  # B=combn(SHat,q)
  
  SigInv=solve(SigHat)
  VFlag=matrix(0,nrow=pz,ncol=ncol(B))
  for(m in 1:ncol(B)){
    bm=B[,m]
    lambda=svd(gamHat[bm,])$d
    SinMIN=min(lambda) # the minmum singular value
    SinMAX=max(lambda) # the maximum singular value
    #SinThresh=(SinMAX+sqrt(a0/2*log(n)/n))*sqrt(a0/2*log(n)/n*log(pz)) # adjust 2
    SinThresh=2*sqrt(a0*log(n)/n*q*log(pz))
    if(SinMIN>SinThresh){
      beta_m=solve(gamHat[bm,],GamHat[bm])
      
      # pi_m=GamHat[SHat]-gamHat[SHat,]%*%beta_m
      # sig_m=as.numeric(sqrt(t(c(1,-beta_m))%*%ThetaHat%*%c(1,-beta_m)))
      # v_m=SigInv[SHat,]-gamHat[SHat,]%*%solve(gamHat[bm,],SigInv[bm,])
      # h=sqrt(diag(v_m%*%SigHat[SHat,SHat]%*%t(v_m)))
      # thresh_m=sig_m*h*sqrt(a0*log(n)/n) ## constant needs a addictive sqrt?
      # VmHat=which(abs(pi_m)<=thresh_m)# valid set based on bm: V_bm
      # VFlag[,m]=SHat %in% VmHat
      
      pi_m=GamHat-gamHat%*%beta_m
      v_m=SigInv[1:pz,]-gamHat%*%solve(gamHat[bm,],SigInv[bm,])
      sig_m=as.numeric(sqrt(t(c(1,-beta_m))%*%ThetaHat%*%c(1,-beta_m)/sum(c(1,-beta_m)^2)))
      h=sqrt(diag(v_m%*%SigHat%*%t(v_m)))
      thresh_m=sig_m*h*sqrt(a0*log(n)/n) ## constant needs a addictive sqrt?
      VmHat=which(abs(pi_m)<=thresh_m)# valid set based on bm: V_bm
      VFlag[,m]=(1:pz) %in% VmHat
    }
  }
  
  
  
  maxIV=max(apply(VFlag,2,sum))  #the maximum number of valid IVs
  BHatIndex=which(apply(VFlag,2,sum)==maxIV) #give all couples with the maximum number of valid IVs
  
  
  ## voting
   VHat=(1:pz)[VFlag[,which.max(apply(VFlag,2,sum))]>0]
   # VHat=SHat[VFlag[,which.max(apply(VFlag,2,sum))]>0]
  
  # Warning check: more than one largest cluster
  if(length(BHatIndex)>1){
    max_columns <- VFlag[, BHatIndex]
    identical_columns <- apply(max_columns, 2, function(col) all(col == max_columns[, 1]))
    all_identical <- all(identical_columns)
    if(!all_identical){
      warning("VHat Warning: More than one cluster has the largest cluster size. This may be due to identification condition not being met. \n")
      #VHat=(1:pz)[apply(VFlag[,BHatIndex],1,sum)>0]
    }
    }
  
  # Error check: not enough IVs
  if(length(VHat) < q){
    warning("VHat Warning: No enough valid IVs estimated. This may be due to weak IVs or identification condition not being met. Use more robust methods.\n")
    warning("Defaulting to all relevant IVs being valid.\n")
    VHat = 1:pz
  }
  
  
  
  return(list(VHat = VHat, BHatIndex=BHatIndex))
}