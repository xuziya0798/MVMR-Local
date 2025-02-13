# mvmr_voting.R
# Select the valid IVs and estimate multiple treatments effects. 
#
# Usage: [VHat,betaHat,ci,betaSD] = mvmr_voting(Y,D,Z,X,intercept,alpha,a0)
#input:
# Y: nx1 outcome vector
# D: nxq treatment matrix
# Z: nxpz candidate IVs
# X: nxpx corvariates
# intercept: whether or not introduce a intercept in linear regression
# alpha: confidence level
# a0: tuning parameter to adjust threshoulding level in first stage; -a0<= beta <= a0
# vote: whether use the voting method

#output:
# VHat: estimated valid IVs
# betaHat: estimated treatments effects
# ci: alpha-level confidence intervals for treatments effects
# betaSD: standard deviation of estimated treatments effects



mvmr_voting <- function(Y,D,Z,X,a0=1,intercept=FALSE,alpha=0.05){
  
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
  
  
  # compute the covariance of W
  SigHat=1/n*t(W)%*%W
  SigInv=solve(SigHat)
  
  #compute Omega: the joint variance of gamhat and Gamhat
  Vi_term <- lapply(1:n, function(i) {
    vi <- c(uHat[i],EHat[i,])   # 提取 v_i 行向量
    vi %*% t(vi)
  })
  
  W_proj <- W %*% SigInv[,1:pz]
  Wi_term <- lapply(1:n, function(i) {
    W_proj[i, ] %*% t(W_proj[i,])
  })
  
  OmegaHat <- (Reduce("+", mapply(function(Vi, Wi) {
    kronecker(Vi, Wi)
  }, Vi_term, Wi_term, SIMPLIFY = FALSE)))/n
  
  
 
  result<-Voting.summary(GamHat,gamHat,SigHat,OmegaHat,n,a0)
  
  
  return(list( VHat=result$VHat, SHat=result$SHat, betaHat=result$betaHat, CI=result$CI, betaSD=result$betaSD))
   

}




Voting.summary <- function(GamHat,gamHat,SigHat,OmegaHat,n,a0=1,alpha=0.05) {
  # Include intercept
  pz=length(GamHat)
  q=length(gamHat)/pz
  p=nrow(SigHat)
  # 
  #========= estimate S =================================================
 
 
  
 
 B=combn(1:pz,q) # B is a q*C(s,q) matrix
 
  
  SigInv=solve(SigHat)
  SFlag=rep(NA,ncol(B))
  
  for(m in 1:ncol(B)){
    bm=B[,m]
    lambda=svd(gamHat[bm,])$d
    SinMIN=min(lambda) # the minmum singular value
    SinMAX=max(lambda) # the maximum singular value

    SinThresh=2*sqrt(a0*log(n)/n*log(pz))
    SFlag[m]=(SinMIN>SinThresh)

  }

  selected_cols <- which(SFlag == 1)
  
  SHat=sort(unique(as.vector(B[, selected_cols])))
  ##Error check
  if(length(SHat) < q){
    warning("VHat Warning: No enough relevant IVs estimated. This may be due to weak IVs or identification condition not being met. Use more robust methods.\n")
    warning("Defaulting to all IVs being relevant.\n")
    SHat = 1:pz
  }

 #=========== estimate V =================================================

  
  VFlag=matrix(0,nrow=length(SHat),ncol=ncol(B))
  for(m in selected_cols){
    bm=B[,m]
    beta_m=solve(gamHat[bm,],GamHat[bm])
    pi_m=GamHat[SHat]-gamHat[SHat,]%*%beta_m
    mat_beta=kronecker(t(c(1,-beta_m)),diag(pz)[SHat,]-gamHat[SHat,]%*%solve(gamHat[bm,],diag(pz)[bm,]))
    thresh_m=sqrt(a0*log(n)/n)*sqrt(diag(mat_beta%*%OmegaHat%*%t(mat_beta))) # threshold of pi_bm
    VmHat=which(abs(pi_m)<=thresh_m)# valid set based on bm: V_bm
    VFlag[,m]=(1:length(SHat)) %in% VmHat
  }

  
  
  maxIV=max(apply(VFlag,2,sum))  #the maximum number of valid IVs
  BHatIndex=which(apply(VFlag,2,sum)==maxIV) #give all couples with the maximum number of valid IVs
  
  
  ## voting
  # VHat=(1:pz)[VFlag[,which.max(apply(VFlag,2,sum))]>0]
  VHat=SHat[VFlag[,which.max(apply(VFlag,2,sum))]>0]
  
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
  
  # ====== Obtain divw est, se, and ci==============================
  
  #2SLS
  SigHat.V=solve((solve(SigHat))[VHat,VHat])
  denom=t(gamHat[VHat,])%*%SigHat.V%*%gamHat[VHat,]
  betaHat=solve(denom,t(gamHat[VHat,])%*%SigHat.V%*%GamHat[VHat])

  #se
  Iv=sort(as.vector(outer(VHat, 0:q, function(j,k) j+k*pz)))
  product=kronecker(t(c(1,-betaHat)),t(gamHat[VHat,])%*%SigHat.V)
  denomInv=solve(denom)
  num=product%*%OmegaHat[Iv,Iv]%*%t(product)
  H=denomInv%*%num%*%denomInv/n
  sd=sqrt(diag(as.matrix(H)))
  CI=cbind(betaHat-qnorm(1-alpha/2)*sd,betaHat+qnorm(1-alpha/2)*sd) # confidence interval
  
  return(list( VHat=VHat, SHat=SHat, betaHat=betaHat, CI=CI, betaSD=sd, BHatIndex=BHatIndex))
  
}
