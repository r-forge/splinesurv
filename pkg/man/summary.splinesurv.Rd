\name{summary.splinesurv}
\alias{summary.splinesurv}
\alias{print.summary.splinesurv}
\title{Summary method for splinesurv objects}
\description{Prints a summary of a fit from \code{\link{splinesurv}}, including
    posterior means and quantiles of the posterior distribution for the regression
    parameters.}
\usage{
\method{summary}{splinesurv}(object,quantiles=c(.025,.975),...)
}
\arguments{
    \item{object}{an object of type \code{\link{splinesurv}}.}
    \item{quantiles}{a vector of quantiles at which the posterior distribution
        should be printed. By default prints the boundary points of a two-sided
        95\% credible interval.}
    \item{...}{additional parameters for \code{\link{print}}.}
}
\value{
    An item of type \code{summary.splinesurv} with components \code{call, coef,
    iter, burnin, nknots.l, nknots.u, ord.l, ord.u, knotspacing.l, knotspacing.u}
    as documented in \code{\link{splinesurv}}, as well as:
    \item{knotrange.l, knotrange.u}{vectors giving the boundary points of
        the knots for the baseline hazard and frailty density splines respectively.}
    \item{quantiles}{matrix with the posterior mean and desired quantiles of
        the regression parameters.}
}
\seealso{
\code{\link{splinesurv}}
}
\keyword{survival}