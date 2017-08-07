



tvvarDATA <- function(data,
                    nv,
                    nt,
                    nb
                    )



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

