# This is a simulation for basic settings
# Do the parallel for speeding



library(dplyr)
library(foreach)
library(doParallel)
library(MASS)
library(mvtnorm)
library(ivreg)
source('mvmr_local.R')
source('mvmr_voting.R')


#setwd('Desktop/project/MVMR')
#----------our setting--------------------------
Niter = 1000 # repeate 1000 times
pz=p=10
q=2
n=n.list=4000
c=0.2


RunIter<-function(r,p,q,e0,a0=1,iter){
  library(MASS)
  library(mvtnorm)
  library(ivreg)
  library(lars)
  set.seed(r)
  beta=rep(1,q)
  gam0<-cbind(rnorm(p,1,1),rnorm(p,2,1))
  pi<-rep(1,p)
  if(e0==1){ # majority rule + strong IV
    V<-sample(1:p,floor(0.7*p))
    gam<-gam0
    pi[V]<-0
    }
  else if(e0==2){ # majority rule + weak IV
    V<-sample(1:p,floor(0.7*p))
    gam<-gam0*c
    pi[V]<-0
  }else if(e0==3){ # plurality rule + strong IV
    V<-sample(1:p,floor(0.5*p))
    gam<-gam0; 
    pi[V]<-0
  }else if(e0==4){ # plurality rule + weak IV
    V<-sample(1:p,floor(0.5*p))
    gam<-gam0*c
    pi[V]<-0
  }else if(e0==5){ # validate setting
    pi=0.4*c(rep(1,4),rep(0,6),rep(1,5),rep(0,6))
    V=c(5:10,16:21)
    gam=cbind(c(runif(10,1.5,2.5),rep(0,11)),c(rep(0,10),runif(11,1.5,2.5)))
    beta=c(0.3,0.6)
  }
  re.all<-NULL
  V=sort(V)
  for(n in n.list){
    Gam <- gam* beta + pi
      
      for (it in 1:iter) {
        set.seed(r)
        #generate W,e,E
        if(e0==5){
          Sigma = matrix(nrow =p , ncol = p)
          for(i in 1:p){
            for (j in 1:i) {
              Sigma[i,j]=0.5^(i-j)
              Sigma[j,i]=Sigma[i,j]
            }
          }
          
          Omega0=diag(1+q)*0.5
          Omega0[1,2]=0.25;Omega0[1,3]=0.3
          Omega0=Omega0+t(Omega0)
        }
        else{
           Sigma=diag(p)
           Omega0=matrix(0.4,1+q,1+q)+diag(1+q)*0.6
        }
        mu=rep(0,p);mu_e=rep(0,1+q)
        
        W=mvrnorm(n,mu,Sigma)
        Z=W[,1:p]
        e=mvrnorm(n,mu_e,Omega0); E=e[,1+1:q]
        
        #data generation
        D=Z%*%gam+E
        Y=D%*%beta+Z%*%pi+e[,1]
       
        
        grid.summ<-data.frame(matrix(NA, ncol = 9, nrow = iter))
        vote.summ<-data.frame(matrix(NA, ncol = 9, nrow = iter))
        naive2SLS.summ<-data.frame(matrix(NA, ncol = 9, nrow = iter))
        oracle2SLS.summ<-data.frame(matrix(NA, ncol = 9, nrow = iter))
        naiveTSHT.summ<-data.frame(matrix(NA, ncol = 9, nrow = iter))
        PostAlasso.summ<-data.frame(matrix(NA, ncol = 9, nrow = iter))
        PostAlassoBlock.summ<-data.frame(matrix(NA, ncol = 9, nrow = iter))
   
        grid <- tryCatch(
          {
            res=mvmr_local(Y=Y,D=D,Z=Z,a0=a0)
          },
          error=function(e) {
            return(NULL)# Use default function here
          }
        )
        if(is.null(grid)){
          grid.summ[it,]<- rep(NA,9) 
        }
        else{
          betaHat=res$betaHat; CI=res$CI; SD=res$betaSD; VHat=res$VHat; BHatIndex=res$BHatIndex
          
          AE=abs(betaHat-beta) #abosolute error
          FlagCI=((beta>=CI[,1])+(beta<=CI[,2])==2) # indicator of beta falling in the interval
          grid.summ[it,] <- c(as.numeric(AE[1]),as.numeric(AE[2]),SD[1],SD[2],
                    FlagCI[1],FlagCI[2],all(VHat %in% V),identical(VHat,V),pz-length(VHat))
        }
        vote <- tryCatch(
          {
            res=mvmr_voting(Y=Y,D=D,Z=Z,a0=a0)
          },
          error=function(e) {
            return(NULL)# Use default function here
          }
        )
        if(is.null(vote)){
          vote.summ[it,] <- rep(NA,9)
        }
        else{
          betaHat=res$betaHat; CI=res$CI; SD=res$betaSD; VHat=res$VHat; BHatIndex=res$BHatIndex

          AE=abs(betaHat-beta) #abosolute error
          FlagCI=((beta>=CI[,1])+(beta<=CI[,2])==2) # indicator of beta falling in the interval
          vote.summ[it,]<- c(as.numeric(AE[1]),as.numeric(AE[2]),SD[1],SD[2],
                             FlagCI[1],FlagCI[2],all(VHat %in% V),identical(VHat,V),pz-length(VHat))
        }
        
        #naive 2SLS
        re.naive=ivreg::ivreg(Y ~ D   | Z)
        ci=confint(re.naive)[1+1:q,]
        re.naive=summary(re.naive)
        AE=abs(re.naive$coefficients[1+1:q,1]-beta)
        SD=re.naive$coefficients[1+1:q,2]
        FlagCI=((beta<=(ci)[,2])+(beta>=(ci)[,1])==2)
        naive2SLS.summ[it,]<- c(as.numeric(AE[1]),as.numeric(AE[2]),SD[1],SD[2],
                                FlagCI[1],FlagCI[2],NA,NA,NA)
        #oracle 2SLS
        re.oracle=ivreg::ivreg(Y ~ D+Z[,setdiff(1:p,V)]   | Z)
        ci=confint(re.oracle)[1+1:q,]
        re.oracle=summary(re.oracle)
        AE=abs(re.oracle$coefficients[1+1:q,1]-beta)
        SD=re.oracle$coefficients[1+1:q,2]
        FlagCI=((beta<=(ci)[,2])+(beta>=(ci)[,1])==2)
        oracle2SLS.summ[it,]<- c(as.numeric(AE[1]),as.numeric(AE[2]),SD[1],SD[2],
                           FlagCI[1],FlagCI[2],1,1,NA)
        #TSHT
        for (j in 1:q){
          res=TSHT(Y=Y,D=as.matrix(D[,j]),Z=Z)
          if(length(res$betaHat)>1){
            AE[j]=abs(res$betaHat$MaxClique1-beta[j])
            ci[j,]=res$ci$MaxClique1
            SD[j]=res$beta.sdHat$MaxClique1
          }else{
            AE[j]=abs(res$betaHat-beta[j])
            ci[j,]=res$ci
            SD[j]=res$beta.sdHat
          }
        }
        FlagCI=((beta<=(ci)[,2])+(beta>=(ci)[,1])==2)
        naiveTSHT.summ[it,]<- c(as.numeric(AE[1]),as.numeric(AE[2]),SD[1],SD[2],
                                 FlagCI[1],FlagCI[2],NA,NA,NA)

        #Post-Alasso
        res.a <- tryCatch(
          {
            res=MVadap.dt(Y=Y,D=D,Z=Z)
          },
          error=function(e) {
            return(NULL)# Use default function here
          }
        )
        if(is.null(res.a)){
          PostAlasso.summ[it,] <- rep(NA,9)
        }
        else{
          AE=abs(res.a$betaa.dt[1:q]-beta)
          ci=res.a$ci_dt
          SD=res.a$sd_dt
          FlagCI=((beta<=(ci)[,2])+(beta>=(ci)[,1])==2)
          inV=as.numeric(gsub("Z", "", res.a$invalid.dt))
          VHat=setdiff(1:pz,inV)
          PostAlasso.summ[it,] <- c(as.numeric(AE[1]),as.numeric(AE[2]),SD[1],SD[2],
                                    FlagCI[1],FlagCI[2],all(VHat %in% V),identical(VHat,V),length(inV))
        }



        #Post-Alasso-Block
        res.b <- tryCatch(
          {
            res=MVadap.dtblock(Y=Y,D=D,Z=Z)
          },
          error=function(e) {
            return(NULL)# Use default function here
          }
        )
        if(is.null(res.b)){
          PostAlassoBlock.summ[it,] <- rep(NA,9)
        }
        else{
          betaHat=res.b$betaHat; CI=res.b$CI; SD=res.b$betaSD; VHat=res.b$VHat; BHatIndex=res.b$BHatIndex

          AE=abs(res.b$betaa.dt[1:q]-beta)
          ci=res.b$ci_dt
          SD=res.b$sd_dt
          FlagCI=((beta<=(ci)[,2])+(beta>=(ci)[,1])==2)
          inV=as.numeric(gsub("Z", "", res.b$invalid.dt))
          VHat=setdiff(1:pz,inV)
          PostAlassoBlock.summ[it,] <- c(as.numeric(AE[1]),as.numeric(AE[2]),SD[1],SD[2],
                                         FlagCI[1],FlagCI[2],all(VHat %in% V),identical(VHat,V),length(inV))
        }


        }
      
     
  
      result.mat <- rbind(
        colMeans(grid.summ),
        colMeans(vote.summ),
        colMeans(naive2SLS.summ),
        colMeans(oracle2SLS.summ),
        colMeans(naiveTSHT.summ),
        colMeans(PostAlasso.summ),
        colMeans(PostAlassoBlock.summ)
        )
      colnames(result.mat)<-c('AE1','AE2','SD1','SD2','Flag1CI','Flag2CI','AllValid','Identical','ninV')
      result.df <- data.frame(e=rep(e0,7), beta = rep(beta[1],7), n=rep(n,7), 
                              method = c('grid','vote','naive 2SLS','oracle 2SLS','naive TSHT','PostAlasso','PostAlassoBlock'))
      result.df <- cbind(result.df, result.mat)
      re.all<-rbind(re.all, result.df)
  }
  return(re.all)
}



# for scenario 1,2,3,4

clnum<-detectCores() 
cl <- makeCluster(getOption("cl.cores", clnum))
registerDoParallel(cl)

res1 <- foreach(r = 1:Niter, .combine = 'rbind') %dopar% RunIter(r=r,p=p,q=q,e0=1,iter=1)
res2 <- foreach(r = 1:Niter, .combine = 'rbind') %dopar% RunIter(r=r,p=p,q=q,e0=2,iter=1)
res3 <- foreach(r = 1:Niter, .combine = 'rbind') %dopar% RunIter(r=r,p=p,q=q,e0=3,iter=1)
res4 <- foreach(r = 1:Niter, .combine = 'rbind') %dopar% RunIter(r=r,p=p,q=q,e0=4,iter=1)
#-------validation setting: the majority rule---------------
pz=p=21
n.list=c(1000,4000)
res5 <- foreach(r = 1:Niter, .combine = 'rbind') %dopar% RunIter(r=r,p=p,q=q,e0=5,iter=1)

stopCluster(cl)

options(digits=3)

 result<- rbind(res1,res2,res3,res4,res5)

 

result <- result %>%
  group_by(e, method, beta[1], n) %>%
  summarize(
    AE1 = mean(na.omit(AE1)),
    AE2 = mean(na.omit(AE2)),
    SD1=mean(na.omit(SD1)),
    SD2=mean(na.omit(SD2)),
    Flag1CI = mean(na.omit(Flag1CI)),
    Flag2CI = mean(na.omit(Flag2CI)),
    AllValid=mean(na.omit(AllValid)),
    Identical=mean(na.omit(Identical)),
    ninV=mean(na.omit(ninV))
  )
write.csv(result,"mvmr-summ.csv",row.names = FALSE)



