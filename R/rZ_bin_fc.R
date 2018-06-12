#' Simulate Z based on a probit model
#' 
#' Simulates a random latent matrix Z given its expectation, dyadic correlation
#' and a binary relational matrix Y
#' 
#' 
#' @usage rZ_bin_fc(Z, EZ, rho, Y)
#' @param Z a square matrix, the current value of Z
#' @param EZ expected value of Z
#' @param rho dyadic correlation
#' @param Y square binary relational matrix
#' @return a square matrix , the new value of Z
#' @author Peter Hoff
#' @export rZ_bin_fc
rZ_bin_fc <-
function(Z,EZ,rho,Y)
{ 
  # simulates Z under the contraints
  # (1)  Y[i,j]=1   => Z[i,j]>0
  # (2)  Y[i,j]=0   => Z[i,j]<0
  
  sz<-sqrt(1-rho^2)
  ut<-upper.tri(EZ)
  lt<-lower.tri(EZ)
  
  Y[is.na(Y)]<- -1
  for(y in c((-1):1))
  { 
    lb<-c(-Inf,-Inf,0)[y+2] ; ub<-c(Inf,0,Inf)[y+2]
    
    for(tri in 1:2)
    { 
      if(tri==1){ up<-ut & Y==y }
      if(tri==2){ up<-lt & Y==y }
      
      ez<- EZ[up] + rho*( t(Z)[up]  - t(EZ)[up] )
      zup<-ez+sz*qnorm(runif(sum(up),pnorm((lb-ez)/sz),pnorm((ub-ez)/sz)))
      zerr<-which(abs(zup)==Inf)
      if(length(zerr)>0){ zup[zerr]<-(Z[up])[zerr] }
      Z[up]<-zup
    }
  }
  
  ##
  c<-(sqrt(1+rho) + sqrt(1-rho))/2
  d<-(sqrt(1+rho) - sqrt(1-rho))/2
  E<-matrix(rnorm(nrow(Y)^2),nrow(Y),nrow(Y))
  ZP<-EZ + c*E + d*t(E) 
  A<-( (Y== -1) | ( sign(ZP) == sign(Y-.5)) ) ; diag(A)<-TRUE
  A<-A & t(A)
  Z[A]<-ZP[A]
  ##

  ## this line now redundant because of previous chunk
  diag(Z)<-rnorm(nrow(Z),diag(EZ),sqrt(1+rho))
  ##

  Z
}

