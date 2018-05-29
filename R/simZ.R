#' Simulate Z given its expectation and covariance
#' 
#' Simulate Z given its expectation and covariance
#' 
#' 
#' @usage simZ(EZ, rho, s2 = 1)
#' @param EZ expected value of Z
#' @param rho dyadic correlation
#' @param s2 dyadic variance
#' @return a simulated value of Z
#' @author Peter Hoff
#' @export simZ
simZ <-
function(EZ,rho,s2=1)
{ 
  c<-(sqrt(1+rho) + sqrt(1-rho))/2
  d<-(sqrt(1+rho) - sqrt(1-rho))/2
  EC<-matrix(rnorm(length(EZ)),nrow(EZ),nrow(EZ))
  EC<- sqrt(s2)*( c*EC + d*t(EC) )
  EZ+EC    
} 

