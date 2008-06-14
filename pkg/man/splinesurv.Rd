\name{splinesurv}
\alias{splinesurv}
\alias{splinesurv.formula}
\alias{splinesurv.agdata}
\alias{splinesurv.data.frame}
\alias{print.splinesurv}
\alias{splinesurv.object}
\title{Bayesian survival analysis with B-spline baseline hazard and random
    effects density}
\description{Allows the fitting of proportional hazards survival models to
    possibly clustered data using Bayesian methods. The function follows a MCMC
    method to sample from the posterior distribution of the regression
    parameters, frailties, and parameters specifying the hazard and frailty
    distribution.
    
    The baseline hazard and random effects density are modeled as convex
    combinations of a parametric component (for example, a gamma frailty, or a
    weibull baseline hazard) and a nonparametric component modeled as a
    penalized B-splines, with the penalty depending on either the integrated
    squared second derivative of the curve, or the sum of squared second
    differences in the spline parameters. 
}
\synopsis{\method{splinesurv}{formula}(formula, data = parent.frame(), ...)}
\usage{
\method{splinesurv}{formula}(formula, data = parent.frame(), verbose = 3,
           hazard = NULL, frailty = NULL, control = NULL, coda = FALSE, 
           initial = NULL, usec = TRUE, ...)
}
\arguments{
    \item{formula}{a formula object, similar to \code{\link{coxph}}. The
    response to the left of the \code{~} should be a survival object generated
    by \code{\link{Surv}}. The right side may contain a \code{cluster(x)} term
    if \code{x} is the variable that indicates cluster membership.
    }
    \item{data}{a \code{data.frame} with columns corresponding to the terms in
    the formula.}
    \item{verbose}{an integer from 0 to 5 that determines the quantity of output
    printed to the screen. Setting \code{verbose=0} is completely silent.}
    \item{hazard}{a list containing parameters defining the baseline hazard,
        with the following optional components. For any component that is not set,
        the default is used.
        \describe{
        \item{\code{type}}{a string that determines whether the hazard is
            parametric, nonparametric, or has both components. Takes on values
            \code{"spline"}, \code{"parametric"} or \code{"both"}, with
            \code{"spline"} as default.}
        \item{\code{spline.adaptive}}{logical, determines whether adaptive knot
            selection should be used. Defaults to \code{TRUE}.} 
        \item{\code{spline.ord}}{an integer greater than 1, giving the order
            of the spline defining the baseline hazard. Defaults to
            \code{spline.ord=4}, corresponding to cubic B-splines.}
        \item{\code{spline.nknots}}{an integer giving the number of interior
            knots used in the baseline hazard spline. If \code{NULL}, the number
            of knots is chosen automatically, which is the default. If
            \code{spline.adaptive=TRUE}, this only serves as the initial number
            of knots.}
        \item{\code{spline.nknots.prior}}{a string specifying the type of prior
            on the number of knots to use. Options are \code{"poisson"},
            \code{"geometric"}, \code{"poissonmix"}, \code{"negbin"} and 
            \code{"power"}. }  
        \item{\code{spline.nknots.hyper}}{the parameter(s) for the prior on
            the number of knots. For the Poisson prior, the mean, for the
            Geometric, the probability, for the Poisson-mixture, two means,
            and for the Power prior, the exponent.}
        \item{\code{spline.knotspacing}}{string that determines the way that
            knots are automatically chosen, if applicable. Possible values are
            \code{"quantile"} to place knots at quantiles of observed event
            times, or \code{"equal"} to equally space knots over the range of
            observed times. Defaults to \code{"equal"}. If
            \code{spline.adaptive=TRUE}, this affects selection of candidate
            knot positions.}
        \item{\code{spline.knots}}{a vector of length
            \code{spline.nknots+2*spline.ord-2} giving the positions of all the
            knots used in the baseline hazard spline (including boundary knots).
            If \code{NULL}, knots are chosen automatically.  Defaults to
            \code{NULL}.}
        \item{\code{spline.penalty}}{a string giving the type of penalty to
            be used for the baseline hazard splines. Possible values are
            \code{"none"} for a standard Gaussian prior (default),
            \code{"2diff"}, for a penalty on the second differences, or
            \code{"2deriv"} for a penalty on the integrated squared second
            derivative. }
        \item{\code{spline.hyper}}{a vector of length 2, containing
            hyperparameters for the spline weights. These control the smoothness
            of the spline when \code{spline.adaptive=FALSE}. Defaults to
            \code{c(0.01,0.01)}.} 
        \item{\code{param.dist}}{a string giving the desired parametric form of
            the hazard. Currently supported are \code{"exponential"} and
            \code{"weibull"}.}
        \item{\code{...}}{other parameters can be set as well, including
            tuning parameters, prior variances, hyperparameters, initial values
            for weights, etc, but this should not be attempted by someone
            unfamiliar with the source code.}
        }
    }
    \item{frailty}{a list containing parameters defining the frailty density,
        analogous to \code{hazard}, except for the following 
        differences:
        \describe{
        \item{\code{param.dist}}{a string giving the desired parametric form
            of the frailty. Currently supported are \code{"gamma"} and
            \code{"lognormal"}.}
        }
    }
    \item{control}{a list containing control parameters for the MCMC and
    optimization, with the following optional components. For any component that
    is not set, the default is used.
        \describe{
        \item{\code{burnin}}{an integer giving the number of iterations
            discarded as burn-in. Default is \code{burnin=500}.}
        \item{\code{maxiter}}{an integer giving the total number of MCMC
            iterations, must be greater than \code{burnin}. Default is
            \code{maxiter=1000}.}
        \item{\code{thin}}{an integer giving the number of iterations to be
            discarded between every recorded iteration. That is, if \code{iter=10},
            a total of \code{maxiter*thin} iterations will be run, but only
            \code{maxiter} will be recorded.}
        \item{\code{tun.auto}}{logical, determines whether to adaptively
            calibrate the tuning parameters during the burn-in period to achieve
            an acceptance rate of 0.25.  Defaults to \code{TRUE}.}
        \item{\code{tun.int}}{an integer giving the interval in iterations
            between calibration of the tuning parameters, if
            \code{tun.auto=TRUE}.  Default is \code{100}.}
        }
    }
    \item{initial}{a list containing initial values for the chain. Not
        implemented and currently ignored.}
    \item{coda}{a logical variable indicating whether the \code{coda} package
        should be used to return the parameter history as \code{mcmc} objects.
        Defaults to \code{FALSE}.}
    \item{usec}{a logical variable, determines whether fast C code should be
        used. Defaults to \code{TRUE}, and should only be set to \code{FALSE} 
        for debugging.}
    \item{...}{additional parameters (currently ignored).}
}
\value{
    An object of class \code{splinesurv}, with the following components
    \item{call}{the original call to the model-fitting function}
    \item{posterior.mean}{a list containing the posterior means of all
        parameters, with the following components 
    \describe{
	    \item{\code{coefficients}}{regression coefficients}
	    \item{\code{frailty}}{frailty estimates}
        \item{\code{hazard.spline.par}}{parameters defining the spline component
        (if applicable)}
        \item{\code{hazard.param.par}}{parameters defining the parametric
        component (if applicable)}
        \item{\code{hazard.weight}}{weight of the spline component (if
        applicable)}
        \item{\code{frailty.spline.par}}{parameters defining the spline
        component (if applicable)}
        \item{\code{frailty.param.par}}{parameters defining the parametric
        component (if applicable)}
        \item{\code{frailty.weight}}{weight of the spline component (if
        applicable)}
        }
    }
    \item{history}{a list containing the parameter history of the MCMC
    iterations, with the following components, either as matrices or vectors of
    length \code{maxiter}:
        \describe{
        \item{\code{frailty}}{matrix of frailty estimates at each iteration.}
        \item{\code{coefficients}}{matrix of regression parameter estimates
            at each iteration.}
        \item{\code{hazard.spline.par}}{matrix of parameters defining the
            baseline hazard spline.}
        \item{\code{hazard.spline.knots}}{matrix of knots defining the
            baseline hazard spline.}
        \item{\code{frailty.spline.par}}{matrix of parameters defining the
            frailty density spline.}
        \item{\code{frailty.spline.knots}}{matrix of knots defining the
            frailty density spline.}
        \item{\code{frailty.spline.fvar}}{vector of estimated variances of the
            frailty density spline.}
        \item{\code{priorvar}}{matrix of prior variance parameters for the
            regression parameters, hazard spline, frailty density spline, hazard
            parameters, frailty parameters, hazard weight, frailty weight, in
            that order.}
        \item{\code{accept}}{binary vectors of Metropolis-Hastings
            acceptance indicators for the regression parameters, hazard spline,
            frailty density spline, hazard parameters, frailty parameters,
            hazard weight, frailty weight, and frailties respectively.}
        }
        If in the input, \code{coda=TRUE}, then the returned objects are of
        class \code{mcmc} instead of \code{matrix}.
    }
    \item{frailty}{a list analogous to \code{frailty} in the input, additionally
    containing the number and placement of knots, if not specified in the input}
    \item{hazard}{a list analogous to \code{hazard} in the input, additionally
    containing the number and placement of knots, if not specified in the input}
    \item{control}{a list analogous to \code{control} in the input, containing
    the control parameters used in the procedure, even those not set explicitly
    in the call.}
}
\references{E. Sharef, D. Ruppert and R. Strawderman. \dQuote{A Bayesian
approach to clustered survival analysis with nonparametric baseline hazard and
frailty distributions}, \emph{in preparation}.
}
\seealso{
\code{\link{summary.splinesurv}}, \code{\link{plot.splinesurv}}, \code{\link{coxph}}
}
\examples{
    ## Generate a small survival data set:
    s <- sim.sample(m = 10, Ji = rep(10, 10))
    agdata <- s$agdata
    
    ## Run a (very) short MCMC chain
    fit <- splinesurv(Surv(time, delta) ~ Z + cluster(i), data = agdata, control = list(maxiter = 50, burnin = 25))
    
    ## Run another chain,  with a Weibull hazard and linear B-spline frailties
    # with fixed knots.
    fit2 <- splinesurv(Surv(time, delta) ~ Z + cluster(i), data = agdata, control = list(maxiter = 50, burnin = 25), hazard = list(type = "parametric", param.dist = "weibull"), frailty = list(type = "spline", spline.ord = 2, spline.adaptive = FALSE))
    ## View summaries and plots of the fits
    summary(fit)
    plot(fit, "all")
    summary(fit2)
    plot(fit2, "all")
}
\keyword{survival}
