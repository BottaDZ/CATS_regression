#' Cheking CATS model by parametric bootstrap
#'
#'This function checks if the fitted model can describe the real assembly process.
#'
#'The parametric bootstrap approach is applied: random samples are  generated
#'assuming that the fitted model describes accurately the real data. Then fit of
#'the original model to these parametric bootstrap data is measured by Kullback-Leibler
#'R-squared. If the model is accurate, the R-squared of real data will not differ
#'significantly from distribution R-squared of bootstrap samples.
#'
#'@param m a fitted CATS model
#'@param nrand number of bootstrap samples
#'@param draw.plot Is a diagnostic plot drawn?
#'@param save.R2 Should matrix of R-squared values be get as output?
#'
#'@return if \code{save.R2==T} matrix of R-squared values in bootsrap samples
#'(plots in colunms, bootstrap samples in rows); otherwise \code{NULL}
#'@export
checkCATS<-function(m,nrand=100,draw.plot=TRUE,save.R2=FALSE)
  {
  distrib<-m$params$family
  cat
  n.site<-nrow(m$params$abund) # number of plots
  n.sp<-ncol(m$params$abund) # number of species
  R.sq<-matrix(NA,nrow=nrand,ncol=n.site)
  obs<-matrix(NA,nrow=n.site,ncol=n.sp)
  loglikelihood<-matrix(NA,nrow=n.site,ncol=n.sp)
  loglikelihood0<-matrix(NA,nrow=n.site,ncol=n.sp)

  for (i in 1:nrand)
    {
    for (j in 1:n.site)
      {
      if (distrib=="poisson")
        {
        obs[j,]<-rpois(n=n.sp,lambda=m$fitted[j,])
        loglikelihood[j,]<-dpois(x=obs[j,],lambda=m$fitted[j,],log=T)
        loglikelihood0[j,]<-dpois(x=obs[j,],lambda=mean(m$fitted[j,]),log=T)
        }
      if (distrib=="nb")
        {
        theta<-m$model[[j]]$family$getTheta(TRUE)
        obs[j,]<-rnbinom(n=n.sp,mu=m$fitted[j,],size=theta)
        loglikelihood[j,]<-dnbinom(x=obs[j,],mu=m$fitted[j,],size=theta,log=T)
        loglikelihood0[j,]<-dnbinom(x=obs[j,],mu=mean(m$fitted[j,]),size=theta,log=T)
        }
      if (distrib=="binomial")
        {
        obs[j,]<-rbinom(n=n.sp,prob=m$fitted,size=m$n)
        loglikelihood[j,]<-dbinom(x=obs[j,],prob=m$fitted[j,],size=m$params$n,log=T)
        loglikelihood0[j,]<-dbinom(x=obs[j,],prob=mean(m$fitted[j,]),size=m$params$n,log=T)
        }
      if (distrib=="gaussian")
        {
        sigm<-sd(residuals(m$model[[j]]))
        obs[j,]<-rnorm(n=n.sp,mean=m$fitted,sd=sigm)
        loglikelihood[j,]<-dnorm(x=obs[j,],mean=m$fitted[j,],sd=sigm,log=T)
        loglikelihood0[j,]<-dnorm(x=obs[j,],mean=mean(m$fitted[j,]),sd=sigm,log=T)
        }
      if (distrib=="tw")
        {
        phi.hat <- x$deviance/sum(x$prior.weights)
        power.hat <- x$family$getTheta(TRUE)
        obs[j,]<-rtweedie(n=n.sp,mu=m$fitted[j,],phi=phi.hat,power=power.hat)
        loglikelihood[j,]<-log(dtweedie(y=obs[j,],mu=m$fitted[j,],phi=phi.hat,power=power.hat))
        loglikelihood0[j,]<-log(dtweedie(y=obs[j,],mu=mean(m$fitted[j,]),phi=phi.hat,power=power.hat))
        }


      }
    R.sq[i,]<-1-(rowSums(loglikelihood)/rowSums(loglikelihood0))
    }


  #cat(length(rep(1:length(m$R.sq),each=nrand)))
  if (draw.plot)
     {
     w<-rbind(R.sq,m$R.sq)
     colnames(w)<-seq(1,ncol(R.sq))
     boxplot(w,xlab="plot ID",ylab=expression(paste("R"^"2")), ylim=c(0,1))
     points(1:ncol(R.sq),m$R.sq,pch=19,col="red",cex=1.5)
     # abline(h=0,col="red")
     }
  if (save.R2) return(R.sq)
  }

