#'@export
extract.pvalues<-function(x)
  {

  temp<-sapply(x,function(x) x[,"p-value"])
  return(t(temp))
  }
