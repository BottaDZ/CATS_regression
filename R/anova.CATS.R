#'Randomization-based Analysis of Deviance for CATS model
#'
#'Calculates separate analysis of deviance tables for each plot (sample). P-values
#'are calculated by randomazation (reshuffling) of trait values. Correlations
#'among traits are conserved during randomization
#'
#'@param m a fitted CATS model
#'@param n.random number of randomizations
#'@return list of deviance tables
#'@export
anova.CATS<-function(m,n.random=49)
  {
  out<-vector(mode="list",length=length(m$model))
  Chi.sq<-vector(mode="list",length=length(m$model))
  n.sp<-nrow(m$params$trait)
  n.site<-length(m$model)
  if (inherits(m$model[[1]],"gam"))
    {

    for (j in 1:n.site)
      {
      temp<-anova(m$model[[j]])
      Chi.sq[[j]]<-c(temp$pTerms.table[,"Chi.sq"], temp$s.table[,"Chi.sq"])
      }

    for (i in 1:n.random)
      {
      random.trait<-m$params$trait[sample(1:n.sp),]
      m.random<-CATSregression(m$params$abund,random.trait,family=m$params$family,
                               prior=m$params$prior,formul=m$params$formula.short,
                               est.prior = m$params$est.prior)
      for (j in 1:n.site)
        {
        temp<-anova(m.random$model[[j]])
        Chi.sq[[j]]<-rbind(Chi.sq[[j]],temp$pTerms.table[,"Chi.sq"], temp$s.table[,"Chi.sq"])
        }
      }

    }
   else
    {
    }
  for (j in 1:n.site)
    {
    temp<-Chi.sq[[j]]
    out[[j]]<-cbind(temp[1,],colMeans(temp[-1,]),(rowSums(t(temp[-1,])>=temp[1,])+1)/(n.random+1))
    colnames(out[[j]])<-c("observed LR","random LR","p-value")
    }
  names(out)<-names(m$model)
  return(out)
  }
