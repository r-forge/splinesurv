\name{splinesurv-package}
\alias{splinesurv-package}
\docType{package}
\title{
Nonparametric bayesian survival analysis
}
\description{
Utilities for nonparametric Bayesian analysis of clustered
survival data. The baseline hazard function and frailty density
are modeled using penalized B-splines, and may optionally have
a parametric component.
}
\details{
    The core function is \code{\link{splinesurv}}, and methods for printing and 
    plotting are provided as well. See the \code{\link{splinesurv}} documentation
    for details on model specification and fitting procedure.

    \tabular{ll}{
    Package: \tab splinesurv\cr
    Version: \tab 0.3-4\cr
    Depends: \tab R (>= 2.5.0), MASS, survival\cr
    Suggests: \tab coda\cr
    License: \tab GPL (>= 2)\cr
    Packaged: \tab Wed Dec 26 18:47:45 2007; esharef\cr
    Built: \tab R 2.6.0; ; 2007-12-26 18:47:51; unix\cr
}

Index:
\preformatted{
\link{plot.splinesurv}         Plots for splinesurv objects
\link{predict.splinesurv}         Prediction for splinesurv objects
\link{splinesurv}              Bayesian survival analysis with B-spline
                        baseline hazard and random effects density
\link{summary.splinesurv}      Summary method for splinesurv objects
}
}
\author{
Emmanuel Sharef <ess28@cornell.edu>
Robert L. Strawderman <rls54@cornell.edu>
David Ruppert <dr24@cornell.edu>
}

\keyword{ survival }

