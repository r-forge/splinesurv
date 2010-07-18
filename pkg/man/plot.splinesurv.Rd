\name{plot.splinesurv}
\alias{plot.splinesurv}
\title{Plots for splinesurv objects}
\description{Function to plot the the baseline hazard, frailty density, and
    posterior density of regression parameters for objects of type
\code{\link{splinesurv}}, with optional pointwise credible bands.
} 
\usage{
\method{plot}{splinesurv}(x, which = c("hazard", "survival", "frailty", "coef", "all"), 
    newdata = NULL, iter=NULL, fn = mean, marginal = c("none", "mc", "numerical"), 
    plotknots = TRUE, npoints = 100, npost = 100, alpha = .05, legend = NULL,
    lty = 1, col = 2, lwd = 2, lty.knots = 1, col.knots = 8, lwd.knots = 1, xlab = NULL, 
    ylab = NULL, main = NULL, xlim = NULL, ylim = NULL, tk = FALSE, ...)
}
\arguments{
    \item{x}{an object of type \code{\link{splinesurv}}.}
    \item{which}{a string specifying which plot to produce. Setting \code{"hazard"}
        plots the baseline hazard (default), \code{"survival"} plots the survivor function, \code{"frailty"} plots the frailty density,
        \code{"coef"} plots the posterior density of the coefficients. Setting
        \code{"all"} plots all plots at once.}
    \item{newdata}{ a data frame in the same format as the original data containing data for prediction of the hazard or survival function. If it has more than one row, multiple functions will be plotted, but credible bands will not.}
    \item{iter}{integer, indicating which iteration of the MCMC chain to plot. If \code{NULL} (default), the posterior means are used}
    \item{plotknots}{boolean indicating whether frailty spline knots should be plotted as vertical lines}
    \item{fn}{function to be applied to the posterior samples. Defaults to 
        \code{mean}, but can be set to \code{median} or any other summary function.}
    \item{marginal}{ string specifying whether to plot marginal survivor functions.
        If \code{"none"}, conditional survivor functions are plotted. If \code{"mc"},
        marginal survivor functions are computed by Monte Carlo integration, and
        if \code{"numerical"}, by numerical integration. Ignored if 
        \code{which} is not \code{"survival"}} 
    \item{npoints}{number of points to be used for drawing}
    \item{npost}{number of samples from the posterior to be used for estimating
        the curve}
    \item{alpha}{level for plotting pointwise confidence bands. Setting
        \code{alpha=NULL} omits bands.}
    \item{legend}{legend text, if desired}
    \item{lty}{line type for primary curve}
    \item{col}{color for primary curve}
    \item{lwd}{line width for primary curve}
    \item{lty.knots}{line type for knots}
    \item{col.knots}{color of the knot lines}
    \item{lwd.knots}{line width of the knot lines}
    \item{xlab}{x-axis label}
    \item{ylab}{y-axis label}
    \item{main}{main plot title}
    \item{xlim}{x-axis limits}
    \item{ylim}{y-axis limits}
    \item{tk}{boolean indicator whether to use the tkrplot GUI}
    \item{...}{additional parameters passed on to \code{\link{plot}}}
}
\keyword{survival}
