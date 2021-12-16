#' Fitting Community Assembly via Trait Selection (CATS) model
#'
#' This function fits CATS model by generalized version of the procedure described by Warton et al. (2015)
#'
#' This is the main function of the package. It fits the CATS model. As Warton et al.
#' (2015) has shown fitting GLM with Poisson distribution results in the same parameter
#' estimates as fitting model by the original maximum entropy approach. Advantage
#' of this approach is its higher flexibility: it allows changing distribution
#' (e.g. using negative binomial instead of Poisson for over-dispersed counts) and
#' modelling non-linear relationships by GAM instead of GLM.
#'
#' In the recent version six distributions are supported:
#' \itemize{
#'    \item Poisson (family="poisson")
#'    \item binomial (family="binomial")
#'    \item negative binomial (family="nb")
#'    \item Tweedie with 1<power<2 (family="tw")
#'    \item normal (family="gaussian")
#'         }
#'
#'\bold{Specifying model formula}
#' The fitted model should be specified either in parameter \code{formul} or in
#' \code{complexity}. The later allow three options:
#' \itemize{
#'    \item linear using all traits as predictor; using Poisson distribution it equivalent to
#'    fitting the original CATS model
#'    \item quadratic it includes all pairwise interactions (including quadratic terms)
#'    into the model
#'    \item gam it applied smoothing for each trait separatelyi
#'         }
#'
#'
#'
#'@param abund matrix or data frame of abundances (species should be in columns)
#'@param trait vector of single trait or data frame of several traits
#'@param trait.name Optional name of trait, if \code{trait} is a vector. If missing, the name will be "trait1"
#'@param prior vector of meta-community level relative abundances
#'@param est.prior Should meta-community level relative abundances be estimated from abundances
#'@param complexity Predefined model formulas (see Details)
#'@param formul Right-side of the formula to be fitted
#'@param family Distribution of abundances
#'@param log.link Should log-link be used for binomial data instead of the cannonical logit?
#'@param n number of trials in the case of binomial data
#'@return \code{CATSregression} returns an object of class CATS, which is a list containing the following components
#' \itemize{
#'         \item \code{params} list of input parameters
#'         \item \code{abund} observed abundances
#'         \item \code{model} list of objects of class gam or logbin
#'         \item \code{fitted} fitted abundances
#'         \item \code{R.sq} vector of plotwise goodness-of-fits measured by Kullback-Leibler R-squared
#'         \item \code{R.sq.adj} vector of plotwise goodness-of-fits measured by adjusted Kullback-Leibler R-squared
#'         \item \code{overall.R.sq} goodness-of-fit for the whole dataset measured by Kullback-Leibler R-squared
#'         \item \code{R.sq.adj} goodness-of-fit for the whole dataset measured by of adjusted Kullback-Leibler R-squared
#'         }
#'@references	Warton, D. I.; Shipley, B. & Hastie, T. (2015) CATS regression –
#'a model‐based approach to studying trait‐based community assembly
#'\emph{Methods in Ecology and Evolution} \bold{6}(4): 389-398.
#'\url{http://dx.doi.org/10.1111/2041-210X.12280}
#'@export




CATSregression<-function(abund, trait,prior=NULL, est.prior=FALSE,
                         complexity=c("linear","quadratic","gam"),
                         family=c("poisson","binomial","tw","nb",
                                  "gaussian"),
                         log.link=TRUE,
                         trait.name=NULL,formul=NULL,n=NULL,...)
{
  family<-match.arg(family)
  complexity=match.arg(complexity)
  if (family=="binomial" & log.link==T) require(logbin) else require(mgcv)


  if (!(is.data.frame(abund)|is.matrix(abund))) stop("Abundance should be matrix or data frame!")
  n.site = nrow(abund)  # number of sites
  n.sp = ncol(abund)  # number of species
  if (is.null(rownames(abund)))
    rownames(abund)<-paste("site",1:n.site,sep="")
  if (is.null(colnames(abund)))
    colnames(abund)<-paste("sp",1:n.sp,sep="")
  Y = as.vector(abund)

  if (is.vector(trait))
  {
    trait<-data.frame(trait1=trait)
    if (!is.null(trait.name))
    {
      colnames(trait)<-trait.name[1]
    }
  }
  if (!is.data.frame(trait)) stop("Trait should be vector or data frame!")
  if (!n.sp==nrow(trait)) stop("Abundance and trait data differ in number of species")

  if (log.link==F)
    if ((!is.null(prior))|est.prior)
      if (family=="binomial")
        warning("Effect of meta-community level processes may estimated incorrectly by this settings")


  if (!is.null(n))
    if (!(family=="binomial"))
      n<-NULL

  if (family=="binomial")
    if (is.null(n)) n<-1

  if (!is.null(n))
  {
    if (any(n-Y<0)) stop("Number of trials (n) should not be lower then number of successes!")
    if (!is.vector(n) | length(n)>1) stop("Parameter n should be a single integer!")
    if (n>1) Y<-cbind(Y,n-Y)
  }
  W<-Y
  for (i in 1:ncol(trait)) W<-cbind(W,trait[rep(1:n.sp,each=n.site),i])

  if (est.prior)
    prior<-colSums(abund)/sum(abund)

  if (!is.null(prior))
  {
    if (!is.vector(prior)) stop("Prior should be a vector")
    if (!n.sp==length(prior)) stop("Length of prior differs from number of species")
    prior[prior==0]<-1E-10
    W<- cbind(W,rep(prior,each=n.site))
  }

  W<-as.data.frame(W)

  W$plotID<-as.factor(rep(rownames(abund),n.sp))



  if ( ncol(as.matrix(Y))==2 ) {
    if (is.null(prior)) colnames(W)<-c("Y","n-Y",colnames(trait),"plotID") else colnames(W)<-c("Y","n-Y",colnames(trait),"prior","plotID")
  } else {
    if (is.null(prior)) colnames(W)<-c("Y",colnames(trait),"plotID") else colnames(W)<-c("Y",colnames(trait),"prior","plotID")
  }

  res<-list()
  res$params<-list()
  res$params$abund<-abund
  res$params$trait<-trait
  res$params$prior<-prior
  res$params$est.prior<-est.prior
  res$params$trait.name<-trait.name
  res$params$family=family
  res$params$n<-n


  if (is.null(formul))
    if (log.link & family=="binomial")
      {
      if (complexity=="linear") formul<-paste(colnames(trait),collapse="+")
      if (complexity=="quadratic")
        {
        poly.trait<-as.data.frame(poly(as.matrix(trait),degree=2),)
        ColumnNames<-colnames(poly.trait)
        TraitNames<-colnames(trait)
        for (i in 1:ncol(poly.trait))
          {
          temp<-as.numeric(unlist(strsplit(colnames(poly.trait)[i],split="[.]")))
          temp2<-""
          for (j in 1:length(temp))
            if (temp[j]>0)
              {
              if (temp2=="") temp2<-paste(temp2,TraitNames[j],temp[j],sep="")
              else temp2<-paste(temp2,"_",TraitNames[j],temp[j],sep="")
              }
          ColumnNames[i]<-temp2
          }
        colnames(poly.trait)<-ColumnNames
        W<-cbind(W,poly.trait)
        formul<-paste(colnames(poly.trait),collapse="+")
        }
      if (complexity=="gam") formul<-paste("B(",colnames(trait),")",collapse="+",sep="")
      }
    else
      {
      if (complexity=="linear") formul<-paste(colnames(trait),collapse="+")
      if (complexity=="quadratic")
        formul<-paste("poly(",paste(colnames(trait),collapse=","),",degree=2)",sep="")
      if (complexity=="gam")
        formul<-paste("s(",colnames(trait),")", collapse="+",sep="")
      }

  res$params$formula.short<-formul
  if (!is.null(prior)) formul<-paste(formul,"+offset(log(prior))")
  if ( ncol(as.matrix(Y))==2 ) formul<-paste("cbind(Y,n-Y) ~ ",formul)
  else formul<-paste("Y ~ ",formul)
  res$params$formula.full<-formul


  mod<-vector(mode="list",length=n.site)
  names(mod)<-rownames(abund)
  res$fitted<-vector()
  res$R.sq<-rep(NA,n.site)
  res$R.sq.adj<-rep(NA,n.site)
  overall.dev0<-rep(NA,n.site)

  for (i in 1:n.site)
    {
    if (log.link & family=="binomial")
      {
      if (complexity=="gam") m<-do.call(logbin.smooth,c(list(formula=as.formula(formul),data=W,
                                                               subset=(W$plotID==rownames(abund)[i])),
                                                          list(...)))
      else m<-do.call(logbin,c(list(formula=as.formula(formul),data=W,
                                      subset=(W$plotID==rownames(abund)[i])),
                                 list(...)))
      }
    else
      m<-do.call(gam,c(list(as.formula(formul),data=W,
                           family=family,method="REML",
                           subset=(W$plotID==rownames(abund)[i])),
                           list(...)))
    mod[[i]]<-m
    res$fitted<-rbind(res$fitted,predict(m,type="response"))

    if ( ncol(as.matrix(Y))==2 ) formul0<-"cbind(Y,n-Y) ~ 1"
    else formul0<-"Y ~ 1"

    if (log.link & family=="binomial")
      m0<-do.call(logbin,c(list(as.formula(formul0),data=W,
                                  subset=(W$plotID==rownames(abund)[i])),
                                  list(...)))
    else
      {if (family=="nb")
        m0<-do.call(gam,c(list(as.formula(formul0),data=W,
                               family=negbin(theta=m$family$getTheta(TRUE)),method="REML",
                               subset=(W$plotID==rownames(abund)[i])),
                               list(...)))
      else
        m0<-do.call(gam,c(list(as.formula(formul0),data=W,
                               family=family,method="REML",
                               subset=(W$plotID==rownames(abund)[i])),
                          list(...)))
      }
    dev<-deviance(m)
    dev0<-deviance(m0)
    overall.dev0[i]<-dev0
    k<-extractAIC(m)[1]
    res$R.sq[i]<-(dev0-dev)/dev0
    res$R.sq.adj[i]<-(dev0-dev-k)/dev0
    }

  names(res$R.sq)<-row.names(abund)
  names(res$R.sq.adj)<-row.names(abund)
  names(mod)<-row.names(abund)
  res$abund<-abund
  res$model<-mod
  res$overall.R.sq<-sum(res$R.sq*overall.dev0)/sum(overall.dev0)
  res$overall.R.sq.adj<-sum(res$R.sq.adj*overall.dev0)/sum(overall.dev0)

  class(res)<-'CATS'
  return(res)
}
