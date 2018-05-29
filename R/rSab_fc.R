#' Gibbs update for additive effects covariance
#' 
#' @usage rSab_fc(a,b,Sab0=NULL,eta0=NULL)
#' @param a vector of row random effects
#' @param b vector of row random effects
#' @param Sab0 prior (inverse) scale matrix for the prior distribution
#' @param eta0 prior degrees of freedom for the prior distribution
#' @author Peter Hoff
#' @export rSab_fc
#'
rSab_fc<-function(a,b,Sab0=NULL,eta0=NULL) 
{
  if(is.null(Sab0)){ Sab0<-diag(2) } 
  if(is.null(eta0)){ eta0<-4 }

  solve(rwish(solve(eta0*Sab0+crossprod(cbind(a, b))), eta0+length(a)))
}



