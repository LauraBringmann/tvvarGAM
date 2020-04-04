#' @title xxx
#'
#' @description \code{tvvarDATA} xxx.
#'
#' @param data An \eqn{(nt\times nv)}{(nt x nv)} data matrix.
#' @param tvvarOpt xxx (default = \code{"TVVAR"}).
#' @param nb xxx (default = 10).
#' @param pbar xxx (default = \code{TRUE}).
#' @param scale xxx (default = \code{FALSE}).
#' @param consec xxx (default = \code{NULL}).
#'
#' @return The function returns a list with 1 element:
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
tvvarDATA <- function(data,               # An (nt x nv) data matrix
                      tvvarOpt = "TVVAR",
                      nb       = 10,
                      pbar     = TRUE,
                      scale    = FALSE,
                      consec   = NULL
)
{
  # --------- Compute Aux Variables ---------
  nt <- nrow(data)
  nv <- ncol(data)

  # Standardize or scale your data:
  if (scale == TRUE) data <- scale(data)

  # Use lagData() from the mgm package:
  lagD_obj <- lagData(data,
                      lags   = 1,
                      consec = consec)

  # Preprocess data for gam():
  data.gam <- data.frame(lagD_obj$data_response, lagD_obj$l_data_lags[[1]])
  data.gam <- rbind(NA, data.gam)
  coln     <- colnames(data)
  colnL    <- paste0(coln, "L")
  colnames(data.gam) <- c(coln, colnL)
  tt       <- 1:nrow(data.gam)

  # --------- Fit GAM ---------
  switch(tvvarOpt,
         "TVVAR" =
           {
             formula.right <- paste0("s(tt,by=", colnL, ",k=nb",")", collapse = " + ")
             formula.right <- paste0("s(tt,k=nb) + ", formula.right)
           },

         "VAR"   =
           {
             formula.right <- paste0(colnL, collapse = " + ")
             formula.right <- paste0("tt + ", formula.right)
           },

         "VAR.mixed" =
           {
             # placeholder for future extension
           })

  # Initialize progress bar:
  if (pbar == TRUE) pb <- txtProgressBar(min = 0, max = nv, initial = 0, char = "-", style = 3)

  model <- list()
  for (j in 1:nv) {
    formula.full  <- as.formula(paste(coln[j], " ~ ", formula.right))
    model[[j]]    <- gam(formula.full, data = data.gam, seWithMean = TRUE)

    # Update progress bar:
    if(pbar == TRUE) setTxtProgressBar(pb, j)
  }

  outlist <- list(model = model)
  return(outlist)
}
