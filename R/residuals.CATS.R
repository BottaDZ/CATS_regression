#' Dunn-Smith residuals of CATS model
#'
#' Calculates Dunn-Smith residuals of a fitted CATS model.
#'
#' Advantage of these residulas is that if model was specified correctly, they
#' follow standard normal distribution.
#'
#'@param m a fitted CATS model
#'@return matrix of residulas
#'@export
residuals.CATS<-function(m)
{
  nb.resid.CATS<-function(x)
  {
    theta<-x$family$getTheta(TRUE)
    pred<-predict(x,type="response")
    out<-pnbinom(x$model$Y-1,mu=pred,size=theta)+runif(length(x$model$Y))*dnbinom(x$model$Y,mu=pred,size=theta)
    out<-qnorm(pmax(pmin(out,1-1E-15),1E-15))
    return(out)
  }

  binom.resid.CATS<-function(x,n)
  {
    Y<-x$model[[1]]
    if (n>1) Y<-Y[,1]
    pred<-predict(x,type="response")
    out<-pbinom(Y-1,size=n,prob=pred)+runif(length(Y))*dnbinom(Y,size=n,prob=pred)
    out<-qnorm(pmax(pmin(out,1-1E-15),1E-15))
    return(out)
  }

  gaussian.resid.CATS<-function(x)
  {
    out<-residuals(x,type="response")
    out<-out/sd(out)
    return(out)
  }

  tweedie.resid.CATS<-function(x)
  {
    phi.hat <- x$deviance/sum(x$prior.weights)
    mu.hat <- fitted(x)
    power.hat <- x$family$getTheta(TRUE)
    out<-qnorm(ptweedie(x$model$Y, mu=mu.hat, phi=phi.hat, power=power.hat))
  }


  if (m$params$family=="poisson")
  {
    pred<-predict.CATS(m)
    out<-ppois(m$abund-1,pred)+runif(length(m$abund))*dpois(m$abund,pred)
    out<-qnorm(pmax(pmin(out,1-1E-15),1E-15))
  }

  if (m$params$family=="binomial") out<-t(sapply(m$model,nb.resid.CATS,n=m$n))

  if (m$params$family=="tw") out<-t(sapply(m$model,tweedie.resid.CATS))

  if (m$params$family=="nb") out<-t(sapply(m$model,nb.resid.CATS))

  if (m$params$family=="gaussian") out<-t(sapply(m$model,gaussian.resid.CATS))

  return(as.vector(out))
}
