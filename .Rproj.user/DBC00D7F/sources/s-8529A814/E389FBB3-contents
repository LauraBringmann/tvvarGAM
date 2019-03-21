



tvvarGAM <- function(data, # the n x p data matrix
                     nb = 10,
                     consec,
                     SIMdata=NULL,
                     simulated = FALSE,
                     plot = FALSE,
                     scale = FALSE,
                     beepvar,
                     dayvar,
                     estimates = FALSE,
                     tvvarOpt = "TVVAR",
                     thresholding = FALSE,
                     pbar){


  #---------- Input check ---------

  if(!(is.numeric(data))) stop('Object data has to be numeric.')
  ifelse(simulated==TRUE,ifelse(is.null(SIMdata)==FALSE,ifelse(plot==TRUE,TRUE,TRUE),stop("SIMdata had to be defined")),TRUE)
  if(tvvarOpt=="VAR" & estimates==TRUE & plot==FALSE) stop("You cannot get the estimates of a VAR model with this argument, please put estimates=FALSE")
  if(tvvarOpt=="VAR" & plot==TRUE) stop("With this function you cannot get a plot of the parameters of a VAR model, please put plot=FALSE")

  # --------- Fill in defaults ---------

  if(missing(pbar)) pbar <- TRUE
  if(missing(consec)) consec <- NULL
  if(missing(beepvar)) beepvar <- NULL
  if(missing(dayvar)) dayvar <- NULL

  # --------- Compute Aux Variables ---------

  # ----- Compute consec argument -----


  # Input checks (can only specify consec OR beepvar and dayvar)

  if(!is.null(consec) & !is.null(beepvar)) stop("Please specify the consecutiveness of measurements either via consec, or via dayvar and beepvar")
  if(!is.null(consec) & !is.null(dayvar)) stop("Please specify the consecutiveness of measurements either via consec, or via dayvar and beepvar")

  if(!is.null(dayvar)) if(is.null(beepvar)) stop("Argument beepvar not specified.")
  if(!is.null(beepvar)) if(is.null(dayvar)) stop("Argument dayvar not specified.")

  if(!is.null(beepvar) & !is.null(dayvar)) {

    consec <- mgm:::beepday2consec(beepvar = beepvar,
                                   dayvar = dayvar)

  }
  else{consec <- 1:nrow(data)}
  # --------- Compute Aux Variables ---------

  nt <- nrow(data)
  nv <- ncol(data)
  if(nv==1){if(is.numeric(data)==FALSE) stop('Object data has to be numeric')}
  else{if(all(apply(data[,],2,is.numeric))==FALSE) stop('Object data has to be numeric')}


  tt=1:nt

  # Define colnames, if not provided with data
  if(is.null(colnames(data))) colnames(data) <- paste0("X", 1:nv)

  coln<-colnames(data)#Defining the colnames
  if(is.null(coln)) stop("Colnames need to be defined")
  colnL=paste(coln,"L",sep="") #And the lagged colnames

  call <- list(data= data,nb=nb,consec=consec,
               SIMdata=SIMdata,
               simulated=simulated,
               plot=plot,
               scale=scale,
               estimates=estimates,
               tvvarOpt=tvvarOpt,
               thresholding=thresholding)

  #%##########################################%###
  ####  Part 1: creating the data #############
  #%##########################################%##
  #The functions in this part are created in the file (source code for simulation article 3.R)

  y <- data

  #%##########################################%###
  ####  Part 2: ESTIMATING GAM##### #############
  #%##########################################%##

  mod_all <- tvvarDATA(data = y,
                       nb = nb,
                       pbar = pbar,scale=scale)$model

  # --------- Case: estimates = TRUE ---------

  if(estimates==TRUE | plot==TRUE) {


    if(plot==TRUE & simulated==TRUE){
      par(mfrow=c(nv,(nv+1)))
      tt=1:nt
      k=0
      for (i in 1:nv){
        mod<-mod_all[[i]]


        for ( j in 1:(nv+1))

          if(j==1)
          {plot.gam(mod,seWithMean = TRUE,select=1,rug=F,shift=coef(mod)[1],xlab="Time",ylab=paste("Intercept of variable",coln[i],sep=""))
            lines(tt,SIMdata$aint[,i],col="red")
          }
        else {plot.gam(mod,seWithMean = TRUE,select=j,rug=F,ylim=c(-1,1),xlab="Time",ylab=paste(coln[i]," is regressed on ",colnL[j-1],sep=""))
          k=1+k

          lines(tt,SIMdata$rho[,k],col="red")

        }
      }

    }  else if (plot==TRUE) {  par(mfrow=c(nv,(nv+1)))
      for (i in 1:nv){
        mod<-mod_all[[i]]


        for ( j in 1:(nv+1))
          if(j==1)
          {plot.gam(mod,seWithMean = TRUE,select=1,rug=F,shift=coef(mod)[1],ylab=paste("Intercept of variable",coln[i],sep=""))
          }
        else {plot.gam(mod,seWithMean = TRUE,select=j,rug=F,ylim=c(-1,1),xlab="Time",ylab=paste(coln[i]," is regressed on ",colnL[j-1],sep=""))


        }
      }

    }



    if(estimates==TRUE){

      Results_GAM<-array(NA,c(c(nv+1,nv),nt,3))




      #-----------------thresholding---------------
      if(thresholding==TRUE){
        for (ii in 1:nv){

          mod <- mod_all[[ii]]

          mat_dat<-matrix(c(tt,rep(rep(1,nt),nv)),length(tt),nv+1)
          coln_Data<-colnL
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



          Results<-list('Estimate'=Results_GAM[, , , 2],'CI_low'=Results_GAM[, , , 3],'CI_high'=Results_GAM[, , , 1])


          outlist<-list(call=call,Results_GAM=Results,model=mod_all)

        }


        return(outlist)
      }





      #---------non thresholding ------------------------
      else{
        for (ii in 1:nv){

          mod <- mod_all[[ii]]

          mat_dat<-matrix(c(tt,rep(rep(1,nt),nv)),length(tt),nv+1)
          coln_Data<-colnL
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


          Results<-list('Estimate'=Results_GAM[, , , 2],'CI_low'=Results_GAM[, , , 3],'CI_high'=Results_GAM[, , , 1])


          outlist<-list(call=call,Results_GAM=Results,model=mod_all)


        }
        return(outlist)


      }


    }




    # Return Estimates





  }
  else{
    outlist<-list(call=call,model=mod_all)
    return(outlist)
  }

}

