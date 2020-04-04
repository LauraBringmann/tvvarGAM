#' @import stats graphics utils mgcv mvtnorm


#' @title Fit the xxx model (gam)
#'
#' @description \code{tvvarGAM} xxx.
#'
#' @param data An \eqn{(nt\times nv)}{(nt x nv)} data matrix, or an object of
#'   class 'tvvarSIM'.
#' @param nb xxx (default = 10).
#' @param consec xxx (default = \code{NULL}).
#' @param scale xxx (default = \code{FALSE}).
#' @param beepvar xxx (default = \code{NULL}).
#' @param dayvar xxx (default = \code{NULL}).
#' @param tvvarOpt xxx (default = \code{"TVVAR"}).
#' @param thresholding xxx (default = \code{FALSE}).
#' @param pbar xxx (default = \code{TRUE}).
#'
#' @return The function returns a list (an object of class \code{tvvarGAM}) with 3
#'   elements:
#'   \item{call}{xxx.}
#'   \item{Results_GAM}{xxx.}
#'   \item{model}{xxx.}
#'
#' @section Details: xxx.
#'
#' @references
#' \insertRef{RobertsLaughlin1996}{GGUM}
#'
#' \insertRef{Robertsetal2000}{GGUM}
#'
#' @author Laura Bringmann, \email{l.f.bringmann@rug.nl}
#'
#' @examples
#' # Example 1 - xxx
#'
#' # Example 2 - xxx
#' @export
tvvarGAM <- function(data         = NULL,     # An (nt x nv) data matrix *or* an object of class 'tvvarSIM'
                     nb           = 10,
                     consec       = NULL,
                     scale        = FALSE,
                     beepvar      = NULL,
                     dayvar       = NULL,
                     tvvarOpt     = "TVVAR",
                     thresholding = FALSE,
                     pbar         = TRUE)
{
  #---------- Input check ---------
  ifelse (is.null(data),
          stop("Parameter 'data' is empty! \n  Either supply a data matrix or a simulated data object of class 'tvvarSIM'."),
          ifelse(class(data) == "tvvarSIM",
                 {
                   SIMdata   <- data
                   data      <- data$y
                   simulated <- TRUE
                 },
                 ifelse(is.numeric(data),
                        simulated <- FALSE,
                        stop("Parameter 'data' is empty! \n  Either supply a data matrix or a simulated data object of class 'tvvarSIM'.")
                 )
          )
  )

  # ----- Compute consec argument -----
  ifelse (is.null(consec),
          ifelse (is.null(beepvar) || is.null(dayvar),
                  ifelse (is.null(beepvar) && is.null(dayvar),
                          consec <- 1:nrow(data),
                          stop("Parameter 'consec' was not provided; only 'dayvar' or 'beepvar' was provided.\n  In such cases, provide BOTH 'dayvar' and 'beepvar'.")),
                  consec <- beepday2consec(beepvar = beepvar, dayvar  = dayvar)),
          if (!is.null(beepvar) || !is.null(dayvar))
            stop("Please specify the consecutiveness of measurements either via consec, OR via dayvar and beepvar.")
  )

  # --------- Compute Aux Variables ---------
  nt <- nrow(data)
  nv <- ncol(data)
  tt <- 1:nt

  # Define colnames, if not provided with data:
  if (is.null(colnames(data))) colnames(data) <- paste0("X", 1:nv)
  coln  <- colnames(data)
  # The lagged colnames:
  colnL <- paste0(coln, "L")

  call <- list(data         = if (simulated) SIMdata else data,
               nb           = nb,
               consec       = consec,
               simulated    = simulated,
               beepvar      = beepvar,
               dayvar       = dayvar,
               scale        = scale,
               tvvarOpt     = tvvarOpt,
               thresholding = thresholding)

  # --------- Estimating GAM ---------
  mod_all <- tvvarDATA(data     = data,
                       tvvarOpt = tvvarOpt,
                       nb       = nb,
                       pbar     = pbar,
                       scale    = scale,
                       consec   = consec)$model

  # --------- Retrieving results ---------
  Results_GAM <- array(NA, c(nv+1, nv, nt, 3))

  estimates     <- lapply(1:nv,      function(x) plot(mod_all[[x]], select = "None", n = nt))
  estimates.fit <- lapply(estimates, function(x) sapply(1:(nv+1), function(y) x[[y]]$fit))
  estimates.se  <- lapply(estimates, function(x) sapply(1:(nv+1), function(y) x[[y]]$se))
  estimates.int <- lapply(1:nv,      function(x) cbind(rep(coef(mod_all[[x]])[1], nt), matrix(0, nt, nv)))

  for (ii in 1:nv)
  {
    Results_GAM[, ii, , 1] <- t(estimates.int[[ii]] + estimates.fit[[ii]] + estimates.se[[ii]])
    Results_GAM[, ii, , 2] <- t(estimates.int[[ii]] + estimates.fit[[ii]])
    Results_GAM[, ii, , 3] <- t(estimates.int[[ii]] + estimates.fit[[ii]] - estimates.se[[ii]])
  }

  if (thresholding)
  {
    tmp.sgn                <- sign(Results_GAM[, ii, , 1] * Results_GAM[, ii, , 3]) > 0
    Results_GAM[, ii, , 2] <- Results_GAM[, ii, , 2] * (tmp.sgn * 1)
  }

  Results <- list('Estimate' = Results_GAM[, , , 2],
                  'CI_low'   = Results_GAM[, , , 3],
                  'CI_high'  = Results_GAM[, , , 1])

  outlist <- list(call        = call,
                  Results_GAM = Results,
                  model       = mod_all)
  class(outlist) <- "tvvarGAM"
  return(outlist)
}

# Define plot() method for class "tvvarGAM":
#' @export
plot.tvvarGAM <- function(x,   # object of class 'tvvarGAM'
                          ...)
{
  ifelse(class(x$call$data) == "tvvarSIM",
         data <- x$call$data$y,
         data <- x$call$data)

  coln  <- colnames(data)
  colnL <- paste0(coln, "L") # the lagged colnames
  nv    <- ncol(data)
  nt    <- nrow(data)
  tt    <- 1:nt

  par(mfrow = c(nv, nv+1),
      oma = c(2, 2, .25, .25),
      mar = c(2, 2, 1, 1),
      mgp = c(2, 1, 0),
      xpd = NA)

  for (i in 1:nv) {
    mod <- x$model[[i]]

    for (j in 1:(nv+1)) {
      plot.gam(mod,
               seWithMean = TRUE,
               select     = j,
               rug        = FALSE,
               ylim       = if (j == 1) NULL else c(-1, 1),
               shift      = if (j == 1) coef(mod)[1] else 0,
               xlab       = "Time",
               ylab       = if (j == 1) paste0("Intercept of ", coln[i]) else paste0(coln[i], " ~ ", colnL[j-1]),
               bty        = "n"
      )

      if (x$call$simulated)
      {
        if (j == 1) lines(tt, x$call$data$aint[, i], col = "red") else lines(tt, x$call$data$rho[, (i-1)*nv + (j-1)], col = "red")
      }
    }
  }
}

# Define summary() method for class "tvvarGAM":
#' @export
summary.tvvarGAM <- function(object,
                             ...)
{
  object[[3]]
}
