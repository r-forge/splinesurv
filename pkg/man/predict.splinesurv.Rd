\name{predict.splinesurv}
\alias{predict.splinesurv}
\title{Prediction for splinesurv objects}
\description{Function to predict the the baseline hazard, frailty density, 
    linear predictor, or risk for objects of type \code{\link{splinesurv}}.
} 
\usage{
\method{predict}{splinesurv}(object, type = c("hazard","survival","lp","risk","frailty"),
        x = NULL, marginal = c("none", "mc", "numerical"), 
        newdata = NULL, iter = NULL, fn = mean,
        alpha = NULL, npost = 100, ...)
}
\arguments{
    \item{object}{an object of type \code{\link{splinesurv}}.}
    \item{type}{a string specifying the type of prediction. Setting
    \code{"hazard"} predicts the baseline hazard (default), \code{"survival"} predicts the survivor function, \code{"lp"} predicts
    the linear predictor, \code{"risk"} produces risk estimates, and
    \code{"frailty"} predicts the frailty density.}
    \item{x}{ a vector containing the values at which prediction is desired. If
    \code{NULL}, a sequence of length 100 over the data range is used.}
    \item{marginal}{ string specifying whether to plot marginal survivor functions.
        If \code{"none"}, conditional survivor functions are plotted. If \code{"mc"},
        marginal survivor functions are computed by Monte Carlo integration, and
        if \code{"numerical"}, by numerical integration. Ignored if 
        \code{type} is not \code{"survival"}} 
    \item{newdata}{ a data frame in the same format as the original data
    containing data for prediction. If cluster membership is not set to
    \code{NA}, the frailty is included in predictions of the hazard or 
    survival.}
    \item{iter}{integer, indicating which iteration of the MCMC chain to use. If
    \code{NULL} (default), the posterior means are used}
    \item{npost}{number of samples from the posterior to be used for estimating
        the curve}
    \item{fn}{function to be applied to the posterior samples. Defaults to 
    \code{mean}, but can be set to \code{median} or any other summary function.}
    \item{alpha}{level of pointwise confidence bands. Setting \code{alpha=NULL}
    omits bands.}
    \item{...}{additional parameters currently ignored}
}
\keyword{survival}
