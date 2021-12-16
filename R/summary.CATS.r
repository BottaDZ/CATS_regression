#' Predicted abundances
#'
#' It gets the abundances predicted by the fitted CATS model
#'

#'@param m a model fitted by CATSregression function
#'@return predicted values in a numeric vector

summary.CATS<-function(x)
{
  lapply(x$model,summary)
}

