


tvvarGAM <- function(data, # the n x p data matrix
                     nb = 10,
                     consec,
                     SIMdata,
                     simulated = FALSE,
                     plot = FALSE,
                     estimates = TRUE,
                     tvvar=TRUE,
                     tresholding=TRUE,
                     pbar){ #data is the data,nt is the number of time points, nv=number of variables and k is the number of knots



  # --------- Fill in defaults ---------

  if(missing(pbar)) pbar <- TRUE
  if(missing(consec)) consec <- 1:nrow(data) # if not provided, assume all measurements are subsequent

  # --------- Compute Aux Variables ---------

  nt <- nrow(data)
  nv <- ncol(data)
  tt=1:nt


  # --------- Case I: estimates = TRUE ---------

  if(estimates==TRUE) {


    Results_GAM<-array(NA,c(c(nv+1,nv),nt,3))


    #%##########################################%###
    ####  Part 1: creating the data #############
    #%##########################################%##
    #The functions in this part are created in the file (source code for simulation article 3.R)

    y <- data

    #%##########################################%###
    ####  Part 2: ESTIMATING GAM##### #############
    #%##########################################%##

    mod_all <- tvvarGAM(data = y,
                        nb = nb,
                        SIMdata = FALSE,
                        simulated = FALSE,
                        plot = FALSE,
                        estimates = FALSE,
                        pbar = TRUE)

    if(tresholding==TRUE){
    for (ii in 1:nv){

      mod <- mod_all[[ii]]

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

      #-----------------Tresholding---------------



      int.ci<-matrix(NA,nt,Numbrep)
      for (m in 1:Numbrep){
        int.ci[,m]<-  Xp[,1:kdim]%*%modr[m,1:kdim]
      }




      Results_GAM[1,ii,1:nt,1]<-apply(int.ci,1,quantile,c(.975))
      Results_GAM[1,ii,1:nt,3]<-apply(int.ci,1,quantile,c(.025))

      for(x in 1:nt){
        ifelse((Results_GAM[1,ii,x,1]< 0 &&  Results_GAM[1,ii,x,3] > 0) || ( Results_GAM[1,ii,x,1] > 0 &&  Results_GAM[1,ii,x,3] < 0)==TRUE,Results_GAM[1,ii,x,2]<-0,Results_GAM[1,ii,x,2])

      }


      for (j in 1:nv){
        Results_GAM[j+1,ii,1:nt,2]<-Xp[,(j*kdim+1):((j+1)*kdim)]%*% coef(mod)[(j*kdim+1):((j+1)*kdim)]

        #The confidence intervals
        phi.ci<-matrix(NA,nt,Numbrep)
        for (m in 1:Numbrep){
          phi.ci[,m]<-Xp[,(j*kdim+1):((j+1)*kdim)]%*%modr[m,(j*kdim+1):((j+1)*kdim)]
        }
        Results_GAM[j+1,ii,1:nt,1]<-apply(phi.ci,1,quantile,c(.975))
        Results_GAM[j+1,ii,1:nt,3]<-apply(phi.ci,1,quantile,c(.025))

        for(x in 1:nt){
          ifelse((Results_GAM[j+1,ii,x,1]< 0 &&  Results_GAM[j+1,ii,x,3] > 0) || ( Results_GAM[j+1,ii,x,1] > 0 &&  Results_GAM[j+1,ii,x,3] < 0)==TRUE,Results_GAM[j+1,ii,x,2]<-0,Results_GAM[j+1,ii,x,2])

        }


        }



      outlist<-list('Estimate'=Results_GAM[, , , 2],'CI_low'=Results_GAM[, , , 3],'CI_high'=Results_GAM[, , , 1])

    }


    return(Results_GAM=outlist)
  }





  #---------non tresholding ------------------------
    else{
      for (ii in 1:nv){

      mod <- mod_all[[ii]]

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




    # --------- Case I: estimates = FALSE ---------

   } else {




    # Use lagData() from the mgm package

    lagD_obj <- mgm:::lagData(data,
                              lags = 1,
                              consec = consec)


    # Back to Laura's variable names:
    Data2 <- cbind(lagD_obj$data_response, lagD_obj$l_data_lags[[1]])
    Data2 <- rbind(rep(NA, ncol(Data2)), Data2) # to make Laura's code below work
    Data2 <- as.data.frame(Data2)
    coln<-colnames(data)
    colnL=paste(coln,"L",sep="")
    colnames(Data2)=c(coln,colnL)
    Data1 <- Data2




    # --------- Case: tvvar = TRUE ---------
    if(tvvar==TRUE) {
    allcol2=c()
    for(i in 1:nv){allcol2[i]=paste("s(tt,by=",colnL[i],",k=nb",")",sep="")}
    allcol3=paste(allcol2,collapse="+")


    # Progress bar
    if(pbar==TRUE) pb <- txtProgressBar(min = 0, max=nv, initial=0, char="-", style = 3)

    model=list()
    for (j in 1:nv){

      ff <- as.formula(paste(coln[j]," ~ ","s(tt,k=nb)","+",allcol3))
      model[[j]]<-gam(ff,data=Data1,seWithMean=TRUE)

      # Update Progress Bar
      if(pbar==TRUE) setTxtProgressBar(pb, j)


    }
    }

    # --------- Case: tvvar = FALSE and thus a standard VAR is estimated ---------

    else{
      allcol2=c()
      for(i in 1:nv){allcol2[i]=paste(colnL[i],sep="")}
      allcol3=paste(allcol2,collapse="+")


      # Progress bar
      if(pbar==TRUE) pb <- txtProgressBar(min = 0, max=nv, initial=0, char="-", style = 3)

      model=list()
      for (j in 1:nv){

        ff <- as.formula(paste(coln[j]," ~ ","tt","+",allcol3))
        model[[j]]<-gam(ff,data=Data1,seWithMean=TRUE)

        # Update Progress Bar
        if(pbar==TRUE) setTxtProgressBar(pb, j)


      }




    }









    # --------- Plotting Stuff ---------


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
    }  else if (plot==TRUE) {  par(mfrow=c(nv,(nv+1)))
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
    } else {


      # Return Estimates

      return(model)

    }

  }
}


