tvvarDATA <- function(data, # the n x p data matrix
                     nb = 10,
                     consec,
                     beepvar,
                     dayvar,
                     scale=FALSE,
                     tvvarOpt="TVVAR",
                     pbar){ #data is the data,nt is the number of time points, nv=number of variables and k is the number of knots

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

  } # if: specification of consecutiveness via beepvar and dayvar



  nt <- nrow(data)
  nv <- ncol(data)
  tt=1:nt

  # Standardize or scale your data
  if(scale==TRUE){data<-scale(data) }
  # Use lagData() from the mgm package

lagD_obj <- mgm:::lagData(data,
                          lags = 1,
                          consec = consec)


# Back to the variable names:
Data2 <- cbind(lagD_obj$data_response, lagD_obj$l_data_lags[[1]])
Data2 <- rbind(rep(NA, ncol(Data2)), Data2) # to make the code below work
Data2 <- as.data.frame(Data2)
coln<-colnames(data)
colnL=paste(coln,"L",sep="")
colnames(Data2)=c(coln,colnL)
Data1 <- Data2




# --------- Case: tvvarOpt = TVVAR ---------
if(tvvarOpt=="TVVAR") {
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
  outlist<-list(model=model)




  return(outlist)

}

# --------- Case: tvvarOpt = "VAR" and thus a standard VAR is estimated ---------

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



  outlist<-list(model=model)




  return(outlist)

}
}




