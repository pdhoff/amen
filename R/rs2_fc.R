#' Gibbs update for dyadic variance
#' 
#' Gibbs update for dyadic variance
#' 
#' 
#' @usage rs2_fc(E, rho,nu0=NULL,s20=NULL)
#' @param E square residual relational matrix
#' @param rho current value of rho 
#' @param nu0 prior degrees of freedom 
#' @param s20 prior estimate of s2
#' @return a new value of s2
#' @author Peter Hoff
#' @export rs2_fc
rs2_fc <-
function(E,rho,nu0=NULL,s20=NULL)
{ 
    if(is.null(nu0)){ nu0<-1 } 
    if(is.null(s20)){ s20<-1 } 

    H<-mhalf( solve(matrix(c(1,rho,rho,1),2,2)) )
    EM<-cbind(E[upper.tri(E)],t(E)[upper.tri(E)] ) %*%H 
    ED<-diag(E)/sqrt(1+rho)
    1/rgamma(1, (length(EM)+length(ED)+nu0)/2 , (sum(EM^2)+sum(ED^2)+nu0*s20)/2 )
}
