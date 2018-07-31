#' log density for GBME models
#'
#' Calculation of the log conditional density of the 
#' latent AMEN matrix \code{Z} given observed data \code{Y}. 
#' 
#' @param Z n X n latent relational matrix following an AMEN model
#' @param Y n X n observed relational matrix 
#' @param llYZ a vectorizable function taking two arguments, y and z. See details below. 
#' @param EZ n X n mean matrix for \code{Z} based on AMEN model (including additive effects) 
#' @param rho dyadic correlation in AMEN model for \code{Z} 
#' @param s2 residual variance in AMEN model for \code{Z}
#' 
#' @return a symmetric matrix where entry i,j is proportional 
#' to the log conditional bivariate density of \code{z[i,j],z[j,i]}. 
#' 
#' @details This function is used for updating dyadic pairs of 
#' the latent variable matrix \code{Z} based on \code{Y} and 
#' an AMEN model for \code{Z}. The function \code{llYZ} specifies 
#' the log likelihood for each single \code{z[i,j]} based on 
#' \code{y[i,j]}, that is, \code{llYZ} gives the log probability 
#' density (or mass function) of \code{y[i,j]} given \code{z[i,j]}. 
#' 
#' @examples 
#' ## For (overdispersed) Poisson regression, use
#' llYZ<-function(y,z){ dpois(y,z,log=TRUE) } 
#' 
#' @author Peter Hoff
#' @export ldZgbme
ldZgbme<-function(Z,Y,llYZ,EZ,rho,s2=1)
{
  ## from log p(z|ez,rho,s2)
  c<-.5*(1/sqrt(1+rho) + 1/sqrt(1-rho))/sqrt(s2)
  d<-.5*(1/sqrt(1+rho) - 1/sqrt(1-rho))/sqrt(s2)
  E<-Z-EZ
  lpZ<- -.5*(c*E + d*t(E))^2 ; lpZ<-lpZ+t(lpZ) ; diag(lpZ)<-diag(lpZ)/2

  ## from log p(y|z)
  llZ<-llYZ(Y,Z)
  llZ<-llZ+t(llZ) ; llZ[is.na(Y)]<-0  

  lpZ+llZ
}

