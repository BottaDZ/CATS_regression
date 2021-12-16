#'Diagnostic plots
#'
#'It creates two diagnotic plots: (1) QQ-plot of residuals can be used for checking
#'normality of residuals. (2) Residuals vs predicted abundances allows checking
#'constant dispersion (homoscedasticity) of residuals.
#'
#'@param m a fitted CATS model
#'@param ... additional \code{\link[base]{graphical parameters}} passed to plot function
#'@export
plot.CATS<-function(m,...)
  {
  res<-as.vector(residuals(m))
  pred<-as.vector(predict(m))
  par(mfrow=c(1,2))
  qqnorm(res)
  qqline(res)
  plot(res~pred,xlab="predicted abundaces",ylab="Dunn-Smith residuals",...)
  }
