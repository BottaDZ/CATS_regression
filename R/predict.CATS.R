#' Predicted abundances
#'
#' It gets predicted abundances in a matrix form
#'
#'@param m a fitted CATS model
#'@return matrix of predicted abundances (i.e. their expected value)
#'@export

predict.CATS<-function(m)
{
  out<-t(sapply(m$model,predict,type="response"))
  return(out)
}
