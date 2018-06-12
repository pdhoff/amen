#' Gibbs update for multiplicative effects covariance
#' 
#' @usage rSuv_fc(U,V, Suv0=NULL,kappa0=NULL)
#' @param U matrix of row random effects
#' @param V matrix of row random effects
#' @param Suv0 prior (inverse) scale matrix for the prior distribution
#' @param kappa0 prior degrees of freedom for the prior distribution
#' @author Peter Hoff
#' @export rSuv_fc
#'
rSuv_fc<-function(U,V,Suv0=NULL,kappa0=NULL) 
{
  if(is.null(Suv0)){ Suv0<-diag(2*ncol(U)) }
  if(is.null(kappa0)){ kappa0<-2+nrow(Suv0) }

  UV<-cbind(U,V)
  solve(rwish(solve(kappa0*Suv0+t(UV)%*%UV),nrow(U)+kappa0))
}



