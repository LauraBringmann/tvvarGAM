# Accessory functions for tvvarGAM() and tvvarDATA().

# Generate an easy-going correlation matrix:
cormatrix <- function(
  nv,   # number of variables
  min,  # minimal cor value
  max   # maximum cor value
)
{
  tmp                 <- matrix(runif(nv * nv, min, max), nv, nv)
  tmp[lower.tri(tmp)] <- t(tmp)[lower.tri(tmp)]
  diag(tmp)           <- 1
  return(tmp)
}

# Three different options to generate data:
# A. Time invariant:
invariant <- function(
  nt,          # number of time points
  MaxAbsValue  # maximum absolute value of the function
)
{
  rep(MaxAbsValue, nt)
}
# B. Linear:
linear <- function(
  nt,          # number of time points
  MaxAbsValue  # maximum absolute value of the function
)
{
  seq(0, MaxAbsValue, length.out = nt)
}
# C. Sine:
sine <- function(
  nt,          # number of time points
  MaxAbsValue # maximum absolute value of the function
)
{
  tt <- seq(0, nt, length.out = nt)
  MaxAbsValue * sin(2 * pi * tt / nt)
}

# Here you can choose the function to simulate a VAR or TV-VAR process:
choose.coefaint <- function(
  nt,         # number of time points
  nv,         # number of variables
  FUN,        # vector of length nv with 1s=invariant, 2s=linear, 3s=sine (one per variable)
  MaxAbsValue # vector of length nv with maximum absolute value of the function (one per variable)
)
{
  FUNchoose <- c(invariant, linear, sine)
  aint      <- matrix(NA, nt, nv)
  for (i in 1:nv)
  {
    aint[, i] <- FUNchoose[[FUN[i]]](nt, MaxAbsValue[i])
  }
  return(aint)
}

# In order to be able to start stationarity check:
stat.check <- function(
  nt, # number of time points
  nv, # number of variables
  rho # ??
)
{
  WAAR <- c()
  tmp  <- rep(NA, nt)
  for (t in 1:nt)
  {
    tmp[t] <- ifelse (max(abs(eigen(matrix(rho[t, ], nv, nv))$values)) < 1, 1, 0)
  }
  WAAR <- all(tmp == 1)
  return(WAAR)
}

# Here you can choose a function for the auto and cross effects:
choose.coefrho <- function(
  nt,         # number of time points
  nv,         # number of variables
  FUN,        # vector of length (nv*nv) with 1s=invariant, 2s=linear, 3s=sine (one per variable)
  MaxAbsValue # vector of length (nv*nv) with maximum absolute value of the function (one per variable)
)
{
  FUNchoose <- c(invariant, linear, sine)
  rho       <- matrix(NA, nt, nv * nv)

  # Proceed while stationarity is not met:
  conv    <- FALSE
  while (!conv)
  {
    for (i in 1:(nv*nv))
    {
      rho[, i] <- FUNchoose[[FUN[i]]](nt, MaxAbsValue[i])
    }
    s.check <- stat.check(nt, nv, rho)
    if (s.check) conv <- TRUE
  }

  return(rho)
}

# Here is the actual simulation function:
#' @title xxx
#'
#' @description \code{tvvarSIM} xxx.
#'
#' @param nt Number of time points.
#' @param nv Number of variables.
#' @param FUN.aint Vector of length \eqn{nv}{nv} with 1s=invariant, 2s=linear,
#'   3s=sine. Default is a vector of 1s.
#' @param max.aint Vector of length \eqn{nv}{nv} with maximum value for the
#'   intercept. Default is a vector of 1s.
#' @param FUN.rho Vector of length \eqn{nv^2}{nv^2} with 1s=invariant,
#'   2s=linear, 3s=sine. Default is a vector of 1s.
#' @param max.rho Vector of length \eqn{nv^2}{nv^2} with maximum value for the
#'   off-diagonal correlations. Default is a vector of .4s.
#' @param min.cor Minimum correlation value. Default is .1.
#' @param max.cor Maximum correlation value. Default is .5.
#'
#' @return The function returns a list (an object of class \code{tvvarSIM})
#'   with 4 elements:
#'   \item{y}{xxx.}
#'   \item{aint}{xxx.}
#'   \item{rho}{xxx.}
#'   \item{sigma}{xxx.}
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
tvvarSIM <- function(
  nt,                        # number of time points
  nv,                        # number of variables
  FUN.aint = rep(1, nv),     # vector of length nv with 1s=invariant, 2s=linear, 3s=sine
  max.aint = rep(1, nv),     # vector of length nv with maximum value for the intercept
  FUN.rho  = rep(1 , nv*nv), # vector of length (nv*nv) with 1s=invariant, 2s=linear, 3s=sine
  max.rho  = rep(.4, nv*nv), # vector of length (nv*nv) with maximum for the off-diagonal correlations
  min.cor  = .1,             # minimum correlation value
  max.cor  = .5              # maximum correlation value
)
{
  # FUN.aint <- if (is.null(FUN.aint)) rep(1, nv)

  aint     <- choose.coefaint(nt, nv, FUN.aint, max.aint)
  rho      <- choose.coefrho (nt, nv, FUN.rho,  max.rho )
  sigma    <- cormatrix(nv, min.cor, max.cor)
  vecsigma <- matrix(sigma, ncol = 1)
  pphi     <- kronecker(matrix(rho[1, ], nv, nv, byrow = T),
                        matrix(rho[1, ], nv, nv, byrow = T))

  y      <- matrix(NA, nt, nv)
  y[1, ] <- rmvnorm(1,
                    mean  = solve(diag(nv) - matrix(rho[1, ], nv, nv, byrow = TRUE)) %*% matrix(aint[1, ], nv, 1),
                    sigma = matrix(solve(diag(nv * nv) - pphi) %*% vecsigma, nv, nv))
  for (t in 2:nt)
  {
    y[t, ] <- matrix(aint[t, ], nv, 1) + matrix(rho[t-1, ], nv, nv, byrow = T) %*% y[t-1, ] +
      matrix(rmvnorm(1, sigma = sigma), nv, 1) # t=>1  # u should be a martingale difference noise
  }

  colnames(y) <- paste("y", 1:nv, sep = "")

  tmp        <- list(y     = y,
                     aint  = aint,
                     rho   = rho,
                     sigma = sigma #JT# yes or no ?
  )
  class(tmp) <- "tvvarSIM"
  return(tmp)
}
