# Function to check convergence of boundary values

convergence.check<-function(datavec){

  # datavec: vector of normalised GOI values
  nGOIvals<-length(datavec)
  nbounds<-7
  boundvals<-quantile(datavec,probs=c(0.2,0.4,0.6,0.8),na.rm=TRUE)

  return(boundvals)
}
