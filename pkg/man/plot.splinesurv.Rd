\name{plot.splinesurv}
\alias{plot.splinesurv}
\title{Plots for splinesurv objects}
\description{Function to plot the the baseline hazard, frailty density, and
    posterior density of regression parameters for objects of type \code{\link{splinesurv}}.
} 
\usage{
\method{plot}{splinesurv}(x, which = c("hazard", "survival", "frailty", "coef", "all"), 
newdata = NULL, iter=NULL, plotknots = TRUE, npoints = 100, legend = NULL, lty = 1, 
col = 2, lwd = 2, lty.knots = 1, col.knots = 8, lwd.knots = 1, xlab = NULL, 
ylab = NULL, main = NULL, xlim = NULL, ylim = NULL, ...)
}
\arguments{
    \item{x}{an object of type \code{\link{splinesurv}}.}
    \item{which}{a string specifying which plot to produce. Setting \code{"hazard"}
        plots the baseline hazard (default), \code{"survival"} plots the survivor function, \code{"frailty"} plots the frailty density,
        \code{"coef"} plots the posterior density of the coefficients. Setting
        \code{"all"} plots all plots at once.}
    \item{newdata}{ a data frame in the same format as the original data containing data for prediction.}
    \item{iter}{integer, indicating which iteration of the MCMC chain to plot. If \code{NULL} (default), the posterior means are used}
    \item{plotknots}{boolean indicating whether frailty spline knots should be plotted as vertical lines}
    \item{npoints}{number of points to be used for drawing}
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
    \item{...}{additional parameters passed on to \code{\link{plot}}}
}
\keyword{survival}
