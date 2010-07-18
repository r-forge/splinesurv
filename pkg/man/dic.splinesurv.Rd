\name{dic.splinesurv}
\alias{dic.splinesurv}
\title{Deviance Information Criterion (DIC) for splinesurv fits}
\description{ Computes a Monte Carlo estimate of DIC for a fitted \code{splinesurv} object.}
\usage{
dic.splinesurv(object)
}
\arguments{
    \item{object}{An object of class \code{splinesurv}.}
}
\details{
    DIC is computed by estimating the mean deviance of the parameters in the chain, as well as the deviance at the posterior mean.
}
\value{
    A list with component \code{dic} and several others that may be used for diagnostics.
}
\seealso{
\code{\link{splinesurv}}
}
\keyword{survival}
