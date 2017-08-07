


tvvarGAM <- function(data,
                     nv,
                     nt,
                     nb=10,
                     SIMdata,
                     simulated=TRUE,
                     plot=TRUE,
                     estimates=FALSE){ #data is the data,nt is the number of time points, nv=number of variables and k is the number of knots

  if (estimates==TRUE)
  {
    Results_GAM<-array(NA,c(c(nv+1,nv),nt,3))


    #nt=100 # 100 500
    tt=1:nt

    #%##########################################%###
    ####  Part 1: creating the data #############
    #%##########################################%##
    #The functions in this part are created in the file (source code for simulation article 3.R)



    y<-data


    #%##########################################%###
    ####  Part 2: ESTIMATING GAM##### #############
    #%##########################################%##

    for (ii in 1:nv){

      mod<-tvvarGAM(y,nv,nt,nb,SIMdata = FALSE,plot=FALSE)[[ii]]

      mat_dat<-matrix(c(tt,rep(rep(1,nt),nv)),length(tt),nv+1)
      coln_Data<-paste("y",1:nv,"L",sep="")
      coln_Data_full<-c("tt",coln_Data)
      colnames(mat_dat)<-coln_Data_full
      newd<-as.data.frame(mat_dat)
      Xp=predict(mod,newd,type="lpmatrix",seWithMean = TRUE)
      kdim=dim(Xp)[2]/c(nv+1)
      newpre=predict(mod,new.data=newd,type="terms",se=TRUE)
      Results_GAM[1,ii,1:nt,2]<-Xp[,1:kdim]%*% coef(mod)[1:kdim]#basis functions intercept!

      Numbrep=5000
      modr<-mvrnorm(Numbrep,coef(mod),mod$Vp+diag((nv+1)*kdim)*10^(-30))

      #The confidence intervals
      int.ci<-matrix(NA,nt,Numbrep)
      for (m in 1:Numbrep){
        int.ci[,m]<-  Xp[,1:kdim]%*%modr[m,1:kdim]
      }

      Results_GAM[1,ii,1:nt,1]<-apply(int.ci,1,quantile,c(.975))
      Results_GAM[1,ii,1:nt,3]<-apply(int.ci,1,quantile,c(.025))


      for (j in 1:nv){
        Results_GAM[j+1,ii,1:nt,2]<-Xp[,(j*kdim+1):((j+1)*kdim)]%*% coef(mod)[(j*kdim+1):((j+1)*kdim)]

        #The confidence intervals
        phi.ci<-matrix(NA,nt,Numbrep)
        for (m in 1:Numbrep){
          phi.ci[,m]<-Xp[,(j*kdim+1):((j+1)*kdim)]%*%modr[m,(j*kdim+1):((j+1)*kdim)]
        }
        Results_GAM[j+1,ii,1:nt,1]<-apply(phi.ci,1,quantile,c(.975))
        Results_GAM[j+1,ii,1:nt,3]<-apply(phi.ci,1,quantile,c(.025))
      }


      outlist<-list('Estimate'=Results_GAM[, , , 2],'CI_low'=Results_GAM[, , , 3],'CI_high'=Results_GAM[, , , 1])

    }
    return(Results_GAM=outlist)
  }

  else {




  tt=1:nt

  Data1=matrix(0,nt,(nv*2),byrow=T)
  for (h in (nv+1):(nv*2)){
    Data1[,(h-nv)]=data[,(h-nv)]# data wordt hier gelagged

    Data1[,h]=c(NA,data[1:(nt-1),(h-nv)])# data wordt hier gelagged
  }

  Data1=as.data.frame(Data1)
  colnames(Data1)=c(paste("y",1:nv,sep=""),paste("y",1:nv,"L",sep=""))

  coln=colnames(Data1)[1:nv]
  colnL=colnames(Data1)[(nv+1):(nv*2)]
  allcol2=c()
  for(i in 1:nv){allcol2[i]=paste("s(tt,by=",colnL[i],",k=nb",")",sep="")}
  allcol3=paste(allcol2,collapse="+")


  model=list()
  for (j in 1:nv){
    ff <- as.formula(paste(coln[j]," ~ ","s(tt,k=nb)","+",allcol3))
    model[[j]]<-gam(ff,data=Data1,seWithMean=TRUE)

  }




  if(plot==TRUE & simulated==TRUE){
  par(mfrow=c(nv,(nv+1)))
  tt=1:nt
  for (i in 1:nv){
    mod<-model[[i]]
    k=0

    for ( j in 1:(nv+1))
      if(j==1)
      {plot.gam(mod,seWithMean = TRUE,select=1,rug=F,shift=coef(mod)[1],ylab="intercept")
        lines(tt,SIMdata$aint[,i],col="red")
      }
    else {plot.gam(mod,seWithMean = TRUE,select=j,rug=F,ylim=c(-1,1))
      k=1+k

      lines(tt,SIMdata$rho[,k],col="red")

    }
  }
  }

  else if (plot==TRUE) {  par(mfrow=c(nv,(nv+1)))
    for (i in 1:nv){
      mod<-model[[i]]
      k=0

      for ( j in 1:(nv+1))
        if(j==1)
        {plot.gam(mod,seWithMean = TRUE,select=1,rug=F,shift=coef(mod)[1],ylab="intercept")
               }
      else {plot.gam(mod,seWithMean = TRUE,select=j,rug=F,ylim=c(-1,1))
        k=1+k


      }
    }
  }

  else {
    return(model)

  }

}
}


