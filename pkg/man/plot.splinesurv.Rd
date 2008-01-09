\name{plot.splinesurv}
\alias{plot.splinesurv}
\title{Plots for splinesurv objects}
\description{Function to plot the the baseline hazard, frailty density, and
    posterior density of regression parameters for objects of type \code{\link{splinesurv}}.
} 
\usage{
\method{plot}{splinesurv}(x,which=c("base","frail","coef","all"),...)
}
\arguments{
    \item{x}{an object of type \code{\link{splinesurv}}.}
    \item{which}{a string specifying which plot to produce. Setting \code{"base"}
        plots the baseline hazard (default), \code{"frail"} plots the frailty density,
        \code{"coef"} plots the posterior density of the coefficients. Setting
        \code{"all"} plots all plots at once.}
    \item{...}{additional parameters passed on to \code{\link{plot}}}
}
\keyword{survival}