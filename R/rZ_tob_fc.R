#' Simulate Z based on a tobit model
#' 
#' Simulates a random latent matrix Z given its expectation, dyadic correlation
#' and a nonnegative relational matrix Y
#' 
#' 
#' @usage rZ_tob_fc(Z, EZ,rho,s2,Y)
#' @param Z a square matrix, the current value of Z
#' @param EZ expected value of Z
#' @param rho dyadic correlation
#' @param s2 dyadic variance
#' @param Y square relational matrix with nonnegative entries
#' @return a square matrix, the new value of Z 
#' @author Peter Hoff
#' @export rZ_tob_fc
rZ_tob_fc<-function(Z,EZ,rho,s2,Y)
{ 
  # simulates Z under the contraints
  # (1)  Y[i,j]>0 => Z[i,j]=Y[i,j]
  # (2)  Y[i,j]=0 => Z[i,j]<0
  
  sz<-sqrt( s2*(1-rho^2) )
  ut<-upper.tri(EZ)
  lt<-lower.tri(EZ)
  
  for(tri in 1:2)
  { 
    if(tri==1){ up<-ut & Y<=0 }
    if(tri==2){ up<-lt & Y<=0 }
    
    ez<- EZ[up] + rho*( t(Z)[up]  - t(EZ)[up] )
    zup<-ez+sz*qnorm(runif(sum(up),0,pnorm(-ez/sz)))
    zerr<-which(abs(zup)==Inf)
    if(length(zerr)>0){ zup[zerr]<-(Z[up])[zerr] }
    Z[up]<-zup
  }
  Z 
}

