#' Plot results of an AME object
#' 
#' A set of plots summarizing the MCMC routine for an AME fit, as well as some
#' posterior predictive checks.
#' 
#' 
#' @param x the result of fitting an AME model
#' @param ... additional parameters (not used)
#' @return a series of plots
#' @author Peter Hoff
#' @S3method plot ame
plot.ame <-
function(x, ...)
{  
  fit<-x 
  require(amen) 

  gof<-nrow(fit$GOF)>1
  p<-ncol(fit$BETA)

  if(!gof | p==0 )
  {
    par(mfrow=c(1+2*gof,2),mar=c(3,3,1,1),mgp=c(1.75,0.75,0))
  }
  if(gof & p>0 )
  {
    par(mar=c(3,3,1,1),mgp=c(1.75,0.75,0))
    layout(matrix(c(1,3,5,1,3,5,2,4,6,2,4,7),3,4)  )
  }



  mVC<-apply(fit$VC,2,median)
  matplot(fit$VC,type="l",lty=1,ylab="VC")
  abline(h=mVC,col=1:length(mVC) )

  if(p>0) 
  { 
    mBETA<-apply(fit$BETA,2,median)
    matplot(fit$BETA,type="l",lty=1,col=1:p,ylab="BETA")
    abline(h=mBETA,col=1:p )
    abline(h=0,col="gray")
  }

  if(gof)
  { 
    ## ame 
    if(length(dim(fit$GOF))==2)
    {
      for(k in 1:5)
      {
        hist(fit$GOF[-1,k],xlim=range(fit$GOF[,k]),main="",prob=TRUE,
             xlab=colnames(fit$GOF)[k],col="lightblue",ylab="",yaxt="n")
             abline(v=fit$GOF[1,k],col="red")
      }
    } 
    ## ame_rep
    if(length(dim(fit$GOF))==3)  
    { 
      DG<-sweep( fit$GOF[,,-1,drop=FALSE],c(1,2),fit$GOF[,,1],"-")
      DG<-zapsmall(DG)
      for(k in 1:5)
      {
        boxplot(t(DG[k,,]),col="lightblue",ylab=dimnames(fit$GOF)[[1]][k])
        abline(h=0,col="gray")
      }
    } 
  } 

}

