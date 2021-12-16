\name{CATS}
\alias{CATS}

\title{Description of \code{CATS} class}

\description{
CATS object is the output of CATSregression function and imput for other functions in the package
}
\value{
\item{params}{list of input parameters}
\item{W}{abundances, traits, plot ID and offset in vectorized form}
\item{model}{Single object or list of objects (if method=="multiple.models") of class \code{gam} or \code{logbin}}
\item{fitted}{fitted abundances}
\item{R.sq}{Kullback-Leibler R-squared}
\item{R.sq.adj}{adjusted Kullback-Leibler R-squared}
}

