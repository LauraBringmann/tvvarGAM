\name{tvvarGAM}
\alias{tvvarGAM}



\title{
Time-varying vector autoregressive models using GAM.
}
\description{
This is the main function of the package tvvarGAM. With this tvvarGAM one can estimate a smooth function for all variables simultaneously, including the intercepts, using generalized additive models (GAM). The output is an object of the mgcv package (the package used to estimate GAM's). In this way, one can easily estimate and plot time-varying smooth functions for all coefficients in a TV-VAR model.
}


\usage{
tvvarGAM(data, nb = 10, consec, SIMdata = NULL,
         simulated = FALSE, plot = FALSE, estimates = FALSE,
         tvvarOpt="TVVAR", tresholding = FALSE, pbar = TRUE,
         ...)
}

\arguments{
  \item{data}{
     n x p data matrix,
}
  \item{nb}{
     A scalar indicating the number of basis functions used for each smooth function. The default setting is 10, just as in the mgcv package.
}

\item{plot}{
     Plots the smooth functions of your model. The default is FALSE.
}

\item{simulated}{
    Indicates if your data was simulated. The default is FALSE.
}

\item{SIMdata}{
    Optional argument. In case your data is simulated you have to feed the simulated data using the tvvarSIM function.
}

\item{estimates}{
    Optional argument. In case you would like the exact parameter values of the smooth functions you can put Estimates=TRUE. This will take much longer than just the standard GAM function. The default is set to FALSE.
}

\item{tvvarOpt}{If tvvarOpt='TVVAR' a time-varying vector autoregressive model is estimated. Other option is 'VAR', which estimates the standard vector autoregressive model using ordinary least squares. Default is 'TVVAR'.}

\item{thresholding}{If the 95\% credible interval contains a zero the edge-weight is put to zero. Default is FALSE.}

\item{pbar}{If pbar = TRUE a progress bar is shown. Defaults to pbar = TRUE.}


}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The value that are returned is a list as desrcribed in \code{\link{gamObject}}.

In case you want the parameters of the smooth functions directly, use Estimates=TRUE.
This returns a  LIST
 \item{Estimate }{This gives the estimates of all smooth functions. The first argument is the coefficient number of the independent variable starting with the intercept (e.g., interecept, y1L, y2L). The second argument indicates which dependent variable it is (e.g., y1, y2). The third argument is indicative of the number of time points.}
\item{CI_low}{This gives the credible interval (lower limit) of all smooth functions. The first argument is the coefficient number of the independent variable starting with the intercept (e.g., interecept, y1L, y2L). The second argument indicates which dependent variable it is (e.g., y1, y2). The third argument is indicative of the number of time points.}
\item{CI_high}{This gives the credible interval (upper limit) of all smooth functions. The first argument is the coefficient number of the independent variable starting with the intercept (e.g., interecept, y1L, y2L). The second argument indicates which dependent variable it is (e.g., y1, y2). The third argument is indicative of the number of time points.}





}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Laura Bringmann <l.f.bringmann@rug.nl>
Jonas Haslbeck <jonashaslbeck@gmail.com>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also \code{\link{mgcv}}
}
\examples{
#Simulate data
#Example: 1 variable
set.seed(12345)
nt<-100
nv<-1
rho<-choose.coefrho(FUN=c(2),nt=nt,nv=nv,MaxAbsValue = c(0.3)) #linear function with 1 parameter (autoregressiveeffect)
aint<-choose.coefaint(FUN=c(3),nt=nt,nv=nv,MaxAbsValue = c(1)) #cosine function

SIMdata<-tvvarSIM(aint,rho,nt=nt,nv=nv)
DATA<-tvvarSIM(aint,rho,nt=nt,nv=nv)$y

GAM<-tvvarGAM(data=DATA,plot=TRUE)


#Example: 2 variables
set.seed(12345)
nt<-1000
nv<-2
rho<-choose.coefrho(FUN=c(2,2,2,2),nt=nt,nv=nv,MaxAbsValue = c(0.3,0.4,0.1,0.5)) #linear function with 4 parameters (autoregressive and crosslagged effects)
aint<-choose.coefaint(FUN=c(3,3),nt=nt,nv=nv,MaxAbsValue = c(1,1)) #cosine function

SIMdata<-tvvarSIM(aint,rho,nt=nt,nv=nv)
DATA<-tvvarSIM(aint,rho,nt=nt,nv=nv)$y

#Estimate and plot

GAM<-tvvarGAM(data=DATA,plot=FALSE,estimates=FALSE)

#Looking at the object shows you the formula that has been used
GAM$model[[1]] # this shows you the input for the first variable
summary(GAM$model[[1]]) # This is the summary statistics for the first variable


#If all coefficients in the model are estimated as smooth variables you can use the following plotting function

tvvarGAM(data=DATA,plot=TRUE,simulated=TRUE,SIMdata=SIMdata) # in case it is simulated data
tvvarGAM(data=DATA,plot=TRUE,simulated=FALSE)# In case it is not simulated data


#If you want directly the parameters of the estimated coefficients you can do the following:
GAM<-tvvarGAM(data=DATA,estimates=TRUE,plot=FALSE)


#In this case, you can just plot the smooth functions like this and adjust the plot as you like:
#To plot the results of y2 ~ y1L (y2 regressed on y1L)
       plot.ts(GAM$Results_GAM$Estimate[2,2,],ylab="y2 regressed on y1L",ylim=c(-1,1))
       lines(GAM$Results_GAM$CI_low[2,2,],lty=2)
       lines(GAM$Results_GAM$CI_high[2,2,],lty=2)
#To plot the results of y2 ~ intercept (intercept of y2 )
       plot.ts(GAM$Results_GAM$Estimate[1,2,],ylab="intercept y2",ylim=c(-1,1))
       lines(GAM$Results_GAM$CI_low[1,2,],lty=2)
       lines(GAM$Results_GAM$CI_high[1,2,],lty=2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
