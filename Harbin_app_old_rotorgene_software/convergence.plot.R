# Function to check convergence of boundary values

convergence.plot<-function(datavec, plotmain){

  # datavec: vector of normalised GOI values
  nGOIvals<-length(datavec)
  nbounds<-7
  boundvals<-quantile(datavec,probs=c(0.2,0.4,0.6,0.8),na.rm=TRUE)
  
  d<-density(datavec,na.rm=TRUE)
  
  histplot_f <- function(){
    plot(d,main=plotmain,xlab="Normalised GOI value",xlim=c(0,max(d$x)))
    lines(d,lwd=2)
    polygon(d,col="darkgrey")
    for(i in 1:nbounds){
      abline(v=boundvals[i],lty=3)
    }}
    
    histplot <- histplot_f()
    return(histplot)
    
  }
