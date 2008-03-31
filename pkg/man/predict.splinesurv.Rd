\name{predict.splinesurv}
\alias{predict.splinesurv}
\title{Prediction for splinesurv objects}
\description{Function to predict the the baseline hazard, frailty density, 
    linear predictor, or risk for objects of type \code{\link{splinesurv}}.
} 
\usage{
\method{predict}{splinesurv}(object,type=c("hazard","lp","risk","frailty"),     x=NULL,newdata=NULL,iter=NULL,...)
}
\arguments{
    \item{object}{an object of type \code{\link{splinesurv}}.}
    \item{type}{a string specifying the type of prediction. Setting \code{"hazard"}
        predicts the baseline hazard (default), \code{"lp"} predicts the linear predictor, \code{"risk"} produces risk estimates, and \code{"frailty"} predicts the frailty density.}
    \item{x}{ a vector containing the values at which prediction is desired. If \code{NULL}, a sequence of length 100 over the data range is used.}
    \item{newdata}{ a data frame in the same format as the original data containing data for prediction.}
    \item{iter}{integer, indicating which iteration of the MCMC chain to use. If \code{NULL} (default), the posterior means are used}
    \item{...}{additional parameters currently ignored}
}
\keyword{survival}
