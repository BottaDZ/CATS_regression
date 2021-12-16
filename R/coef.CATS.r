#'@export
coef.CATS<-function(m)
{

  temp<-sapply(m$model,coef)
  return(t(temp))
}

