cormatrix2=function(x,nv){
  CC.int=matrix(x,nv,nv,byrow=TRUE)
  CC.int=CC.int+t(CC.int)
  diag(CC.int)<-1
  return(CC.int)
}
#five different options
##1. The three generating functions (genfun) options##
#1A. Time invariant
invariant<-function(nt,MaxAbsValue)  #nt is the sample size, and MaxAbsValue is the maximum absolute value of the function.
{genfun=rep(NA,(nt)) #Creating the function first with NAs (missing values).
genfun=rep(MaxAbsValue,(nt)) #Here the actual invariant function is created.
return(genfun)
}

#1B. Linear
linear<-function(nt,MaxAbsValue)  #nt is the sample size, and MaxAbsValue is the maximum absolute value of the function.
{genfun=rep(NA,(nt)) #Creating the function first with NAs (missing values).
genfun=seq(0,MaxAbsValue,length.out=(nt))  #Here the actual linear function is created.

return(genfun)
}

#1C. Sine
sine<-function(nt,MaxAbsValue) #nt is the sample size, and MaxAbsValue is the maximum absolute value of the function.
{genfun=rep(NA,(nt)) #Creating the function first with NAs (missing values).
tt=1:(nt) #Defining a time parameter in order to create the sine function
genfun=MaxAbsValue*sin(2*pi*tt/(nt)) #Here the actual sine function is created.
return(genfun)
}


choose.coefaint<-function(FUN,nt,MaxAbsValue,nv){
  FUNchoose=c(invariant,linear,sine)
  aint=matrix(NA,(nt),nv)
  for(i in 1:nv){
    aint[,i]=FUNchoose[[FUN[i]]](nt,MaxAbsValue[i])}
  return(aint)
}

#in order to be able to start stationairity check
stat.check=function(nt,rho){
  WAAR=c()
  duizend=rep(NA,nt)
  for (t in 1:(nt)){
    duizend[t]=ifelse (max(abs(eigen(matrix(rho[t,],nv,nv))$values))<1,1,0)
  }
  WAAR=all(duizend==1)
  return(WAAR)
}


choose.coefaint<-function(FUN,nt,MaxAbsValue,nv){
  FUNchoose=c(invariant,linear,sine)
  aint=matrix(NA,nt,nv)
  for(i in 1:nv){
    aint[,i]=FUNchoose[[FUN[i]]](nt,MaxAbsValue[i])}
  return(aint)
}




choose.coefrho<-function(FUNr,nt,MaxAbsValue,nv){
  FUNchoose=c(invariant,linear,sine)
  rho=matrix(NA,nt,nv*nv)
  for(i in 1:(nv*nv)){
    rho[,i]=FUNchoose[[FUNr[i]]](nt,MaxAbsValue[i])
  }

  s.check=stat.check(nt=nt,rho=rho)
  conv=FALSE
  while (!conv){
    if (s.check==TRUE) {conv=TRUE}
    else{
      for(i in 1:(nv*nv)){
        rho[,i]=FUNchoose[[FUNr[i]]](nt,MaxAbsValue[i])
      }
      s.check=stat.check(nt=nt,rho)
    }}
  return(rho)
}



tvvarSIM<-function(aint,rho,nt,nv){
  aint=aint
  sigma2 <- cov.mat<-cormatrix2(0.1,nv)
  y=matrix(0,nt,nv)
  vecsigma2=matrix(sigma2,nv*nv,1)
  pphi=kronecker(matrix(rho[1,],nv,nv,byrow=T),matrix(rho[1,],nv,nv,byrow=T))
  matrix(solve(diag(nv*nv)-pphi)%*%vecsigma2,nv,nv)

  y[1,]<-rmvnorm(1,mean=solve(diag(nv)-matrix(rho[1,],nv,nv,byrow=T))%*%matrix(aint[1,],nv,1),sigma=matrix(solve(diag(nv*nv)-pphi)%*%vecsigma2,nv,nv))

  for (t in 2:nt){
    y[t,]=matrix(aint[t,],nv,1)+matrix(rho[t-1,],nv,nv,byrow=T)%*%y[t-1,]+matrix(rmvnorm(1,sigma=sigma2),nv,1) #t=>1  # u should be a martingal difference noise

  }

  colnames(y)=paste("y",1:nv,sep="")
  return(list(y=y,aint=aint,rho=rho))
}
