#' Simulate a tobit relational matrix
#' 
#' Simulates a tobit relational matrix
#' 
#' 
#' @usage simY_tob(EY, rho, s2)
#' @param EY square matrix giving the expected value of the relational matrix
#' @param rho dyadic correlation
#' @param s2 dyadic variance
#' @return a square matrix
#' @author Peter Hoff
#' @export simY_tob
simY_tob <- function(EY,rho,s2)
{
  YS<-simZ(EY,rho,s2) 
  diag(YS)<-NA  
  YS[YS<=0]<-0 
  YS
}
