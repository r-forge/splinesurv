#{{{ #Header
# Todo: (low) rug plot of frailties and event times
# Todo: (low) precompute posterior curves
# Todo: (low) allow exporting and importing initial values
# Todo: (low) allow plotting spline and parametric component separately
# Todo: (low) better comments (using Doxygen maybe)
# Todo: (low) Needs input checking

#****h* /00main
#  NAME
#    00main --- main fitting routines
#  FUNCTION
#   The splinesurv package contains 
#   utilities for nonparametric Bayesian analysis of clustered 
#   survival data. The baseline hazard function and frailty density are modeled
#   using penalized B-splines. Options include adaptive knot selection and the 
#   inclusion of a parametric component.
#
#   The most important function is splinesurv.agdata, which does most
#   of the model fitting. After initializing, the method either continues in
#   R (if the option usec=TRUE is set), callling the various mh.* procedures,
#   or calls SplineSurvMainLoop in the C compiled code, which does the same
#   thing, only faster. The R implemented routines are thus primarily for
#   debugging.
#   
#   Source code modules are organized as follows: initRoutine contains functions
#   used to initialize the algorithm. After initialization, RFitting contains
#   the routines for running estimation entirely within R, and CFitting contains
#   analogous functions written in C. These two Fitting modules work independently
#   of one another, and give the same result, although the CFitting routines of
#   course run much more quickly. The RFitting routines thus serve primarily as
#   a sanity check, since they are easier to understand. Each of the fitting modules
#   has a number of submodules, see their descriptions for detail. The routines use data
#   structures described in 01structures.
#   
#   The S3Methods module contains user-callable functions, many of which are
#   visible when the package is installed. The simSurvival module contains tools
#   to generate simulated survival data, used in conducting simulation studies.
#
#   The call graph below is rather small and unhelpful. See a larger but 
#   equally unhelpful pdf version here:
#   href:callgraph.pdf
#
#   |exec dot -Tpdf -o callgraph.pdf ../../man/calls2.dot
#   |dotfile ../../man/calls2.dot
# USAGE
#   For usage instructions, see the R package documentation
# CONTENTS
#    splinesurv.agdata --- main estimation function
#    SplineSurvMainLoop --- main loop in C
# AUTHOR
#  Emmanuel Sharef
#*******

#****h* /01structures
# NAME
#    01structures --- data structures used
# FUNCTION
#    The R and C implementations of this routine use particular data structures
#    to contain curves, regression information and estimation history, which are
#    documented here
# CONTENTS
#    CCurve --- structure to store a curve in C
#    CHistory --- structure to store MCMC history in C
#    CRegression --- structure to store regression information in C
#    RCurve --- structure to store a curve in R
#    RHistory --- structure to store MCMC history in R
#    RRegression --- structure to store regression information in R
#******

#****h* /RFitting
# NAME
#   RFitting --- model fitting routines in R
# FUNCTION
#   Functions used to sample the MCMC chain from within R, without calling the C
#   main loop. Some of these may make calls to C code occasionally, but the model-
#   fitting on the whole is done in R. This gives the same results as C, but it's
#   easier to debug and understand.
# 
#   This module is organized into several submodules. MetropolisHastings contains
#   the functions that consitute the main MCMC loop. These rely on likelihood computations
#   in the makeLikelihood module. Updating and evaluating B-splines and curves is
#   handled in the curveUpdate module, which relies on C functions accessible via the
#   CWrappers module. Other utlities related to B-splines are in splineUtils, and 
#   miscellaneous utilities are in miscUtils.
# CONTENTS
#   curveUpdate --- Update B-spline curves
#   CWrappers --- Wrappers for functions written in C
#   makeLikelihood --- Compute likelihoods of parameters
#   MetropolisHastings --- MH and RJMCMC steps
#   miscUtils --- miscellaneous
#   splineUtils --- Utilities for evaluating splines and related integrals
#   ZZdebug --- debugging functions
#******


#****h* RFitting/splineUtils
# NAME
#     splineUtils --- Utilities for evaluating splines and related integrals
# FUNCTION
#     Evaluate B-splines and integrals over B-splines, either in R or fast
#     C code.
# CONTENTS
#    evalBinte --- compute the integrals of each spline basis function
#    evalCinte --- compute partial integrals over the spline basis
#    evalEinte --- compute the expected distance from 1 of a B-spline
#    frailtysplinefvar --- compute the spline component frailty variance
#    ki --- addressing of spline indices
#    makesplinebasis --- construct B-spline basis functions
#    mysplineDesign --- works like splineDesign()
#    nBsmom --- compute the N-th moment of a B-spline basis function
#    nknotsPrior --- evaluate the prior on the number of knots
#    plotspline --- plot a spline given a set of knots and parameters
#    splineconv --- compute the convolution of two splines
#    splinederivint --- compute the convolution of the derivatives of two spline bases
#******

#****h* RFitting/miscUtils
# NAME
#   miscUtils --- miscellaneous
# FUNCTION
#   Miscellanous useful utilities
# CONTENTS
#    accrate.predict.lm --- predicted acceptance rate distance from 25%
#    haspar --- check if a curve has a parametric component
#    hasspline --- check if a curve has a spline component
#    makeoutputcurve --- construct the curve to be returned in the output
#    mdiag  --- replacement for diag()
#    repairfrailtypar --- fix frailty parameters
#    submean --- compute the mean of a subset of a vector
#*******

#****h* /simSurvival
# NAME
#     simSurvival --- Simulate survival data
# FUNCTION
#     Generate simulated clustered survival data with arbitrary baseline hazard
#     and frailty distributions. The main user-visible function here is sim.sample,
#     which calls other simulation routines.
# CONTENTS
#    bs.survfn --- compute survivor function for B-spline hazard
#    dnegbin --- negative binomial distribution
#    generateevents --- Generate random event times
#    generaterandom --- Generate random numbers
#    makeagdata --- convert simulated data into Anderson-Gill format
#    MYmvrnorm  --- multivariate normal random variates
#    rinvgamma --- generate inverse-gamma random variates
#    sim.sample --- main simulation function
#******

#****h* /initRoutine
# NAME  
#   initRoutine --- Initialize the splinesurv routine
# FUNCTION
#   Compute initial values for all parameters of interest and allocate the
#   necessary storage. The module CinitRoutine contains functions implemented
#   in C that are used during the initialization only.
# CONTENTS
#    CinitRoutine --- additional C functions used for initialization only
#    fitparametric --- fit a parametric component to a curve
#    inithistory --- initialize the history structure
#    inverthessian --- invert a Hessian matrix
#    makePenalty.2deriv --- compute a penalty matrix on second derivatives of B-splines
#    makePenalty.2diff --- compute a penalty matrix on second differences
#    makeknots --- make knots for a curve with a spline component
#    makepenalty --- construct a penalty matrix
#    nknotsPriorMean --- compute the prior mean of the number of knots of a curve
#    numHess.par  --- compute a numerical Hessian for the parametric component
#******

#****h* RFitting/curveUpdate
# NAME
#   curveUpdate --- Update B-spline curves
# FUNCTION
#   Re-evaluate B-spline curves, usually with new parameters, different
#   knots or component weights. The routines in R are relatively slow and
#   easier to understand, the ones in C are designed to be faster.
# CONTENTS
#    evalparametric --- evaluate the parametric component of a curve
#    evalspline --- evaluate the spline component of a curve
#    updatecurvex --- change observations of a curve
#    updatehistory --- update RHistory structure
#    updateparametric --- update parametric component parameters
#    updateregression --- update regression coefficients
#    updatespline --- update the curve's spline parameters
#    weightcurve --- re-weight the spline and parametric components of a curve
#******

#****h* RFitting/makeLikelihood
# NAME
#   makeLikelihood --- Compute likelihoods of parameters
# FUNCTION
#   Compute parameter likelihoods for use in Metropolis-Hastings steps
# CONTENTS
#    mkgr.spline.haz --- gradient of hazard spline parameters
#    mkhess.coef --- hessian of regression coefficients
#    mklik.coef --- likelihood of regression coefficients
#    mklik.frail --- likelihood of frailty for cluster i
#    mklik.param.frail --- likelihood of parametric component for frailty
#    mklik.param.haz --- likelihood of parametric component parameters for hazard
#    mklik.spline.frail --- likelihood of frailty spline parameters
#    mklik.spline.frail.init --- likelihood of frailty spline parameters (initialization)
#    mklik.spline.haz --- likelihood of hazard spline parameters
#    mklik.weight.frail --- likelihood of weight of spline component for frailty
#    mklik.weight.haz --- likelihood of weight of spline component for hazard
#    smoothpen --- compute the smoothness penalty for a curve
#******

#****h* RFitting/MetropolisHastings
# NAME
#   MetropolisHastings --- MH and RJMCMC steps
# FUNCTION
#   Execute each type of MH and Reversible-Jump step needed for the chain
# CONTENTS
#    acceptreject --- accept of reject a M-H step
#    mh --- prototype Metropolis-Hastings
#    mh.bdm --- RJMCMC for Birth-death-move steps
#    mh.coef --- MH for regression coefficients
#    mh.frail --- MH for frailties
#    mh.frailty.param --- MH for frailty parametric component parameters
#    mh.frailty.spline --- MH for frailty spline parameters
#    mh.hazard.param --- MH for hazard parametric component parameters
#    mh.hazard.spline --- MH for hazard spline parameters
#    mh.weight --- MH for spline component weight for either hazard or frailty
#    updatepostvar.curve --- update the prior variance for a curve
#    updatepostvar.coef --- update prior variance for regression coefficients
#******

#****h* RFitting/CWrappers
# NAME 
#   CWrappers --- Wrappers for functions written in C
# FUNCTION
#   Make it easy to call functions written in C from inside R
# CONTENTS
#    cevalBinte --- wrapper for the C implementation of evalBinte
#    cevalCinte --- wrapper for the C implementation of cevalCinte
#    cevalEinte --- wrapper for the C implementation of cevalEinte
#    cmakePenalty.2deriv --- wrapper for cMakePenalty2diff
#    cmkgr.spline.haz --- wrapper for cInitGrHazSpline
#    cmklik.spline.frail --- wrapper for cInitLikFrailSpline
#    cmklik.spline.haz --- spline hazard likelihood in C wrapper
#    csplinedesign --- wrapper for csplinedesign
#******

#****h* /S3Methods
# NAME
#   S3Methods --- Methods for S3 classes
# FUNCTION
#   Define the handling for the splinesurv and splinesurv.summary classes.
#
#   The main user-facing function is splinesurv, which controls method-dispatch
#   and generally dispatches to splinesurv.formula, if called with a formula
#   argument. The latter is responsible for creating an input data frame in the
#   format required by splinesurv.agdata, which does the fitting.
#
#   Other routines are available to summarize and plot splinesurv fitted objects.
# CONTENTS
#    plot.splinesurv --- plot method for splinesurv objects
#    post.fvar --- compute posterior frailty variance
#    predict.splinesurv --- prediction method for splinesurv objects
#    print.splinesurv --- print a splinesurv object
#    print.summary.splinesurv --- print summary.splinesurv object
#    printcurvesummary --- summarize a curve for summary.splinesurv
#    splinesurv --- method dispatch for splinesurv
#    splinesurv.data.frame --- splinesurv method for data frames
#    splinesurv.formula --- formula interface for splinesurv
#    splinesurvtkplot --- plot the curve using tcltk
#    summary.splinesurv --- creates an object of class summary.splinesurv
#******
    
#****h* RFitting/ZZdebug
# NAME
#   ZZdebug --- debugging functions
# FUNCTIONS
#   Functions whose sole purpose is debugging, that have not been removed from the codebase
# CONTENTS
#    rmkgr.spline.haz --- R reimplementation of cInitGrHazSpline
#    rmklik.spline.haz --- R re-implementation of the likelihood function in C
#******

#****d* 01structures/RCurve
# NAME
#    RCurve --- structure to store a curve in R
#  FUNCTION
#    This type of structure contains the information about a curve, either the hazard or frailty.
#    This includes all parameters for spline and parametric components, weights, knot positions,
#    etc, basis functions, tuning parameters, etc. It is the R analogue of the CCurve structure.
#
#    In R, this is implemented as a simple list. Most of the components are defined at
#    initialization, others are added during the procedure.
#
#    The components are:
#
#        type                    "spline", "parametric" or "both"
#        hasspline               boolean, whether a splien component is included
#        haspar                  boolean, whether a parametric component is included
#        name                    "hazard" or "frailty"
#        spline.ord              order of the spline
#        spline.adaptive         boolean, whether to use adaptive knot selection
#        spline.knotspacing      "equal", "quantile", type of knot distribution 
#        spline.nknots           number of knots
#        spline.nknots.prior     prior on the number of knots
#        spline.nknots.hyper     hyperparameters for the number of knots
#        spline.ncandknots       number of candidate knots
#        spline.candknots        candidate knots
#        spline.maxoccknots      maximum number of occupied candidate knots
#        spline.bdmconst         constant for birth-death-move step
#        spline.knots            vector of knots
#        spline.par              spline parameters
#        spline.min              minimum value of a spline parameter
#        spline.penalty          type of smoothness penalty to use ("none","2diff","2deriv")
#        spline.penaltyfactor    scaling of the penalty
#        spline.priorvar         prior variance for the spline penalty
#        spline.hyper            hyperparameters for the prior variance
#        spline.basis            B-spline basis
#        spline.basisint         integrals of each B-spline basis function
#        spline.basiscum         cumulative integrals of each basis function (hazard only)
#        spline.basisexp         expectation of each spline basis function (frailty only)
#        spline.candcov          covariance matrix used to generate candidate parameter sets
#        spline.candcholcov      cholesky factorization of spline.candcov
#        spline.fvar             frailty variance (frailty only)
#        spline.fixedind         index of the spline component that remains fixed (frailty only)
#        spline.tun              tuning parameter for the spline parameters
#        spline.accept           acceptance rate of spline parameters
#        param.dist              parametric distribution type
#        param.par               parameters for the distribution
#        param.priorvar          prior variance for the parametric parameters
#        param.hyper             hyperparameters for the variance
#        param.tun               tuning parameter
#        param.accept            acceptance rate of parametric parameters
#        weight                  weight of spline component
#        weight.priorvar         variance of the weight
#        weight.hyper            hyperparameters of the weight
#        weight.tun              tuning parameter for the weight
#        weight.accept           acceptance rate of the weight
#        x                       observations (times for hazard, frailties for frailty curve)
#        spline.y                spline component evaluated at x
#        param.y                 parametric component evaluated at x
#        y                       both components weighted
#        spline.ycum             cumulative hazard evaluated at x using spline component
#        param.ycum              cumulative hazard evaluated at y using parametric component
#        ycum                    cumulative hazard
#******

#****d* 01structures/RRegression
# NAME
#    RRegression --- structure to store regression information in R
#  FUNCTION
#    This type of structure contains the information about the regression components.
#    It is the R analogue of the CRegression structure.
#
#    In R, this is implemented as a simple list. Most of the components are defined at
#    initialization, others are added during the procedure.
#
#    The components are:
#        m               number of clusters
#        Ji              vector of cluster sizes
#        cluster         cluster index for each subject
#        status          vector of status indicators
#        covariates      matrix of covariates for each subject
#        coefficients    regression coefficient estimates
#        candcov         covariance matrix for candidate generation
#        lp              vector of linear predictors
#        priorvar        prior variance of regression coefficients
#        hyper           hyperparameters for prior variance
#        tun             tuning parameters
#        accept          most recent acceptance indicator
#******

#****d* 01structures/RHistory
# NAME
#    RHistory --- structure to store MCMC history in R
#  FUNCTION
#    This type of structure contains the state of the MCMC procedure at each iteration.
#    It is the R analogue of the CHistory structure.
#
#    In R, this is implemented as a simple list. Most of the components are defined at
#    initialization, others are added during the procedure. Each of the component has
#    maxiter rows, and as many columns as needed to store the vector at each iteration.
#    For spline parameters and knots, enough storage is allocated to store the maximum
#    number of knots permitted by the adaptive selection procedure (if used). Unused storage
#    is filled with -Inf.
#
#    See inithistory for the allocation procedure.
#
#       frailty                 frailty estimates for each cluster
#       coefficients            regression coefficient estimates
#       loglik                  full log-likelihood
#       hazard.spline.par       hazard spline parameters
#       hazard.spline.knots     hazard spline knots
#       frailty.spline.par      frailty spline parameters
#       frailty.spline.knots    frailty spline knots
#       frailty.spline.fvar     frailty spline variance
#       hazard.param.par        hazard parametric component parameters
#       frailty.param.par       frailty parametric component parameters
#       hazard.weight           hazard spline component weight
#       frailty.weight          frailty spline component weight
#       priorvar                matrix of prior variances 
#       accept                  matrix of acceptance rates
#******
#}}}

{{{ #Simulation
##############################################################
# \section{Simulation} Generate simulated data
##############################################################

# Todo: diagnostic only, remove eventually
#****f* splineUtils/plotspline
#  NAME
#    plotspline --- plot a spline given a set of knots and parameters
#  FUNCTION
#    Create a plot of a spline curve, optionally plotting knots and component
#    splines
#  INPUTS
#    knots          vector of length K+Q containing a set of knot positions, with
#                   the Q border knots repeated at each end.
#    theta          vector of length K with spline parameters
#    ord            spline order Q
#    npoints        number of points at which to evaluate the spline
#    plotknots      whether to plot the knots as vertical lines
#    plotmean       whether to plot the spline mean as a vertical line
#    plotsplines    whether to plot the component spline basis functions
#    norm           whether to normalize the spline (as for the frailty)
#    ...            other parameters for plot
#  OUTPUTS
#    none
#  SYNOPSIS
plotspline <- function(knots, theta, ord, npoints = 1000, plotknots = T, plotmean = F,
        plotsplines = F, norm = F, xlim = NULL, ylim = NULL, col = "red", ...)
#  SOURCE
#
{
    knots <- knots[knots>-Inf]
    theta <- theta[theta>-Inf]
    if(is.null(xlim)) xlim = range(knots)
    # Generate point set for evaluation
    x = seq(from = xlim[1], to = xlim[2], length = npoints)
    dx = diff(xlim) / npoints
    #   Compute spline basis
    spl <- csplinedesign(knots, x = x, ord = ord)
    if(norm){
        Bint <- cevalBinte(knots, ord)
        for(i in 1:dim(spl)[1]) spl[i, ] <- spl[i, ] / Bint
    }
    #   Evaluate the spline at the set of points x
    splsum <- spl%*%exp(theta) 
    if(norm) splsum <- splsum / sum(exp(theta))
    ymax <- max(splsum)
    if(is.null(ylim)) ylim <- c(0, ymax)
    #   Plot basis functions
    if(plotsplines){
        matplot(x, spl / max(spl) * ylim[2] / 2, type = "l", lty = 2, xlim = xlim,
                ylim = ylim, ...)
        lines(x, splsum, col = col, lwd = 2, lty = 1)
    }else{
        plot(x, splsum, col = col, type = "l", lwd = 2, lty = 1, xlim = xlim,
                ylim = ylim, ...)
    }
    if(plotknots) abline(v = knots, col = "grey")
    # Compute and plot the mean
    if(plotmean){
        Bint <- colSums(spl) / dx
        spl.norm <- spl%*%mdiag(1 / Bint)
        Eint <- rep(0, dim(spl)[2])
        for(i in 1:dim(spl)[2]) Eint[i] <- sum(spl.norm[, i] * x) / dx
        E <- Eint%*%exp(theta) / sum(exp(theta))
        abline(v = E, col = "grey", lwd = 2)
    }
}
#************ plotspline 

# my MYmvrnorm for testing (to make sure normal variates in C code are the same)
# Todo: Remove this function eventually
#****f* simSurvival/MYmvrnorm
#  NAME
#    MYmvrnorm  --- multivariate normal random variates
#  FUNCTION
#    Generate multivariate normal variables. This is a plug-in replacement
#    for mvrnorm from MASS, but it uses the Cholesky factorization method,
#    so that the numbers from the equivalent C routine are identical.
#  INPUTS
#    n      number of variables to generate
#    mu     mean vector
#    Signma covariance matrix
#  OUTPUTS
#    A matrix (or vector) of multivariate normal numbers.
#  SYNOPSIS
MYmvrnorm <- function (n = 1, mu, Sigma) 
#  SOURCE
#
{
    l <- length(mu)
    candnorm <- matrix(rnorm(n * l), nrow = n)
    out <- candnorm%*%chol(Sigma, pivot = TRUE)
    out <- out + rep(mu, rep(n, l)) 
    if(n == 1) out <- as.numeric(out)
    out
}
#************ MYmvrnorm 

#****f* simSurvival/generaterandom
#  NAME
#    generaterandom --- Generate random numbers
#  FUNCTION
#    Generate random numbers from any of the following distributions:
#     fixed         not random, fixed at a certain value
#     weibull       Weibull, parametrized by rate and shape
#     gamma         Gamma, parametrized by mean and variance
#     normal        Normal, parametrized by mean and variance
#     lognormal     Lognormal, parametrized by mean and variance
#     normmix       Mixture of normals, parametrized by means, variances and weights
#     lognormmix    Mixture of lognormals, parametrized by means, variances and weights
#     unifmix       Mixture of uniforms, parametrized by weights and bounds
#  INPUTS
#     n         number of variates to generate 
#     type      string, one of the types above
#     params    a list with parameters, which are different for each type. See the code
#               for each type
#  OUTPUTS
#     out       n random numbers with the specified distribution
#  SYNOPSIS
generaterandom <- function(n, type, params)
#  SOURCE
#
{
    if(!(type%in%c("fixed", "weibull", "gamma", "normal", "lognormal", "normmix",
            "lognormmix", "unifmix"))) stop("Invalid distribution type")
    # fixed at params$value
    if(type == "fixed"){
        if(!("value"%in%names(params))) stop("Parameter value not specified for type fixed")
        return(rep(params$value, n))
    }
    # weibull with rate params$lambda0 and shape params$gamweib
    if(type == "weibull"){
        if(!all(c("lambda0", "gamweib")%in%names(params))) 
                stop("Parameters lambda0, gamweib not specified for type weibull")
        lambda0 <- params$lambda0
        gamweib <- params$gamweib
        return(rweibull(n, shape = gamweib, scale = lambda0^(-1 / gamweib)))
    }
    # gamma with mean params$mu and variance params$sigma2
    if(type == "gamma"){
        if(!all(c("mu", "sigma2")%in%names(params)))
                stop("Parameters mu, sigma2 not specified for type gamma")
        mu <- params$mu
        sigma2 <- params$sigma2
        return(rgamma(n, shape = mu^2 / sigma2, scale = sigma2 / mu))
    }
    # normal with mean params$mu and variance params$sigma2
    if(type == "normal"){
        if(!all(c("mu", "sigma2")%in%names(params))) 
                stop("Parameters mu, sigma2 not specified for type normal")
        mu <- params$mu
        sigma2 <- params$sigma2
        return(rnorm(n, mean = mu, sd = sqrt(sigma2)))
    }
    # lognormal with mean params$mu and variance params$sigma2
    if(type == "lognormal"){
        if(!all(c("mu", "sigma2")%in%names(params)))
                stop("Parameters mu, sigma2 not specified for type lognormal")
        mu <- params$mu
        sigma2 <- params$sigma2
        sigma2prime <- log(1 + sigma2 / mu^2)        
        muprime <- log(mu) - 1 / 2 * sigma2prime
        return(rlnorm(n, meanlog = muprime, sdlog = sqrt(sigma2prime)))
    }
    # normal mixture with weights params$w, means params$mu and variances params$sigma2
    if(type == "normmix"){
        if(!all(c("mu", "sigma2", "w")%in%names(params)))
                stop("Parameters mu, sigma2, w not specified for type normmix")
        mu <- params$mu
        sigma2 <- params$sigma2
        w <- params$w / sum(params$w)
        if(length(w) == 1) w = rep(w, length(mu))
        if(length(mu) != length(sigma2) | length(mu) != length(w) |
                length(mu) < 2) stop("Bad parameter lengths for type normmix")
        out <- MYmvrnorm(n, mu, mdiag(sigma2)) 
        return(t(out)[findInterval(runif(n), cumsum(w)) + 1 + 0:(n - 1) * length(w)])
    }
    # lognormal mixture with weights params$w, means params$mu and variances params$sigma2
    if(type == "lognormmix"){
        if(!all(c("mu", "sigma2", "w")%in%names(params)))
                stop("Parameters mu, sigma2, w not specified for type lognormmix")
        mu <- params$mu
        sigma2 <- params$sigma2
        w <- params$w / sum(params$w)
        if(length(w) == 1) w = rep(w, length(mu))
        if(length(mu) != length(sigma2) | length(mu) != length(w) |
                length(mu) < 2) stop("Bad parameter lengths for type lognormmix")
        sigma2prime <- log(1 + sigma2 / mu^2)        
        muprime <- log(mu) - 1 / 2 * sigma2prime
        out <- MYmvrnorm(n, muprime, mdiag(sigma2prime)) 
        return(exp(t(out)[findInterval(runif(n), cumsum(w)) + 1 + 0:(n - 1) * length(w)]))   
    }
    # uniform mixture with weights params$w, bounds params$bounds
    if(type == "unifmix"){
        if(!all(c("w", "bounds")%in%names(params)))
                stop("Parameters w, bounds not specified for type unifmix")
        w <- params$w / sum(params$w)
        bounds <- matrix(params$bounds, ncol = 2, byrow = TRUE)
        out <- rep(0, n)
        for(i in 1:n){
            which <- sum(runif(1) > cumsum(w)) + 1
            out[i] <- runif(1, bounds[which, 1], bounds[which, 2])
        }
        return(out)
    }
}
#************ generaterandom 

#****f* simSurvival/generateevents
#  NAME
#    generateevents --- Generate random event times
#  FUNCTION
#    Generate a set of random event times for a sample, given its size,
#    a single covariate, frailties and regression coefficients.
#    Many baseline hazard specifications are supported, see code comments.
#  INPUTS
#    m      number of clusters
#    Ji     vector of cluster sizes
#    beta   a single true regression coefficient
#    Ui     length m vector of "true" frailties
#    Zij    length sum(Ji) vector of covariates
#    type   type of baseline hazard specification
#    params parameters for the baseline hazard
#  OUTPUTS
#    Tij    length sum(Ji) vector of event times
#  SYNOPSIS
generateevents <- function(m, Ji, beta, Ui, Zij, type, params)
#  SOURCE
#
{
    Tij <- rep(0, sum(Ji))
    Uij <- rep(Ui, Ji)
    #  Weibull baseline hazard, baseline params$lambda0 and shape params$gamweib
    if(type == "weibull"){
        if(!all(c("lambda0", "gamweib")%in%names(params))) stop("Parameters lambda0, gamweib, w not specified for type weibull")
        lambda0 <- params$lambda0
        gamweib <- params$gamweib
        for(ind in 1:sum(Ji)){
            Tij[ind] <- ((-exp(-beta * Zij[ind]) * log(1 - runif(1)) / 
                (lambda0 * Uij[ind]))^(1 / gamweib))
        }
    }
    #  stepfunction baseline hazard, with breakpoints at params$breaks and hazards
    #  at params$haz
    if(type == "stepfunction"){
        breaks <- params$breaks
        haz <- params$haz
        if(length(haz) != length(breaks) + 1) stop("Step function params: haz should be one longer than breaks")
        for(ind in 1:sum(Ji)){
            accept <- FALSE
            Tijprop <- 0
            maxhaz <- max(haz) * Uij[ind] * exp(beta * Zij[ind])
            #  Use accept-reject sampling
            while(!accept){
                Tijprop <- Tijprop -1 / maxhaz * log(runif(1))
                thishaz <- haz[findInterval(Tijprop, breaks) + 1] * Uij[ind] * 
                    exp(beta * Zij[ind])
                if(maxhaz * runif(1) < thishaz){
                    Tij[ind] <- Tijprop
                    accept <- TRUE
                }
            }
        }
    }
    #  B-spline baseline hazard. params$b is an object returned by bs(),
    #  and params$w is a set of weights
    if(type == "bspline"){
        b <- params$b
        w <- params$w
        rbound <- attr(b, "Boundary.knots")[2]
        survfn <- bs.survfn(b, w, 1000, rbound)
        if(survfn>.25) warning(paste("Baseline survival over B - spline support is high:", 
            format(survfn, digits = 3)))
        for(ind in 1:sum(Ji)){
            accept <- FALSE
            Tijprop <- 0
            # Accept-reject sampling
            while(!accept){
                maxhaz <- max(w) * Uij[ind] * exp(beta * Zij[ind])
                Tijprop <- Tijprop -1 / maxhaz * log(runif(1))
                thishaz <- predict(b, min(Tijprop, rbound - 1e-5))%*%w * Uij[ind] * 
                    exp(beta * Zij[ind])
                if(maxhaz * runif(1) < thishaz){
                    Tij[ind] <- Tijprop
                    accept <- TRUE
                }
            }
        }
    }
    return(Tij)
}
#************ generateevents 

#****f* simSurvival/bs.survfn
#  NAME
#    bs.survfn --- compute survivor function for B-spline hazard
#  FUNCTION
#   Given a B-spline baseline hazard, either compute the survivor
#   function at a set of times, or at a single time.
#  INPUTS
#   b       B-spline basis in the form of bs() output
#   w       weights of each basis function
#   n       number of intervals to use
#   t       time at which to compute survival (whole function if NULL)
#  OUTPUTS
#   If t = NULL, the survivor function at n times, else at time t.
#  SYNOPSIS
bs.survfn <- function(b, w, n, t = NULL)
#  SOURCE
#
{
    rbound <- attr(b, "Boundary.knots")[2]
    x <- seq(from = 0, to = rbound, length = n)
    h <- (predict(b, x)%*%w)
    if(is.null(t)){
        H <- cumsum(h * rbound / n)
        S <- exp(-H)
        return(S)
    }else{
        Ht <- sum(h[x < t] * rbound / n)
        St <- exp(-Ht)
        return(St)
    }
}
#************ bs.survfn 

#****f* simSurvival/makeagdata
#  NAME
#    makeagdata --- convert simulated data into Anderson-Gill format
#  FUNCTION
#    Create an Anderson-Gill formatted data frame given a set of
#    event times and covariates.
#  INPUTS
#    m      number of clusters
#    Ji     cluster size
#    Tij    event times simulated by generateevents
#    deltaij censoring indicators
#    Zij    covariates
#  OUTPUTS
#    agdata Data frame with columns
#       i       cluster indicator
#       j       subject ID
#       time    event time
#       delta   event indicator
#       ...     covariate columns
#  SYNOPSIS
makeagdata <- function(m, Ji, Tij, deltaij, Zij)
#  SOURCE
#
{
    agdata <- data.frame(
        i = rep(1:m, Ji),
        j = c(sapply(Ji, function(x) return(1:x)), recursive = TRUE),
        time = Tij,
        delta = deltaij
    )
    agdata <- cbind(agdata, Zij)
    return(agdata)  
}
#************ makeagdata 

#****f* simSurvival/sim.sample
#  NAME
#    sim.sample --- main simulation function
#  FUNCTION
#    Generate a simulated sample of survival data
#  INPUTS
#    m      number of clusters
#    Ji     cluster size
#    params a list with multiple optional components. If any
#           of these is not given, defaults are used.
#       beta        "true" regression coefficient
#       haz.type    type of baseline hazard (see generateevents)
#       haz.params  parameters for baseline hazard
#       frail.type  type of frailty distribution (see generaterandom)
#       frail.params    parameters for frailty
#       Z.type      type of covariate
#       Z.params    params for covariate
#       C.type      type of censoring process
#       C.params    params for censoring
#  OUTPUTS
#    agdata     An A-G data frame containing a simulated sample
#    Ui         "true" frailties for the sample
#    params     parameters used to generate the sample
#  SYNOPSIS
sim.sample <- function(m = 10, Ji = rep(5, 10), params = NULL)
#  SOURCE
#
{
    if(length(Ji) == 1) Ji <- rep(Ji, m)
    if(length(Ji) != m) stop("Length of Ji does not match m")
    params.in <- params
    #   Default parameters
    params.default <- list(
        beta = 1,
        haz.type = "weibull",
        haz.params = list(lambda0 = 1, gamweib = 1.8),
        frail.type = "lognormal",
        frail.params = list(mu = 1, sigma2=.25),
        Z.type = "normal",
        Z.params = list(mu = 0, sigma2 = 1),
        C.type = "weibull",
        C.params = list(lambda0 = 1, gamweib = 1.8)
    )
    params <- params.default
    if(!is.null(params.in)){
        for(n in names(params.in)) eval(parse(text = paste("params$", n, " <- params.in$", 
            n, sep = "")))
        if(params.in$haz.type == "bspline" & is.null(params.in$haz.params)){
            tmax = 3;
            N = 4
            b <- bs(0, knots = seq(from = 0, to = tmax, length = N + 1)[ - (N + 1)], 
                Boundary.knots = c(0, 2), degree = 3)
            w <- c(.3, .2, .4, .2, .6, .8, 1) * 4
            params$haz.params <- list(b = b, w = w)
        }
    }
    # Make covariates
    Zij <- generaterandom(sum(Ji), params$Z.type, params$Z.params)
    # Make censoring
    Cij <- generaterandom(sum(Ji), params$C.type, params$C.params)
    if(!is.null(params$C.max)) Cij[Cij > params$C.max] <- params$C.max
    # Make frailties
    Ui <- generaterandom(m, params$frail.type, params$frail.params)
    # Make event times
    Tij <- generateevents(m, Ji, params$beta, Ui, Zij, params$haz.type, params$haz.params)
    # Make event indicators
    deltaij <- as.numeric(Tij < Cij)
    # apply censoring
    Tij[deltaij == 0] <- Cij[deltaij == 0]
    
    # make output data frame
    agdata <- makeagdata(m, Ji, Tij, deltaij, data.frame(Z = Zij))
    return(list(agdata = agdata, Ui = Ui, params = params))
}
#************ sim.sample 

}}}

{{{ #Utility
##############################################################
# \section{Utility} Utility and Extraction functions
##############################################################

#****f* miscUtils/hasspline
#  NAME
#    hasspline --- check if a curve has a spline component
#  INPUTS
#    curve      an RCurve object
#  OUTPUTS
#    boolean indicating whether the curve has a spline component
#  SYNOPSIS
hasspline <- function(curve) return(curve$type == "spline" | curve$type == "both")
#************ hasspline 

#****f* miscUtils/haspar
#  NAME
#    haspar --- check if a curve has a parametric component
#  INPUTS
#    curve      an RCurve object
#  OUTPUTS
#    boolean indicating whether the curve has a parametric component
#  SYNOPSIS
haspar <- function(curve) return(curve$type == "parametric" | curve$type == "both")
#************ haspar 

#****f* miscUtils/mdiag
#  NAME
#    mdiag  --- replacement for diag()
#  FUNCTION
#    works like diag(), except for length 1 vector inputs, in which case it
#    returns a 1x1 matrix
#  INPUTS
#    x      a vector or matrix
#  OUTPUTS
#    if x is a vector, a matrix with x as its diagonal
#    if x is a matrix, a vector containing its diagonal
#  SYNOPSIS
mdiag <- function(x) 
#  SOURCE
#
if(is.vector(x) && length(x) == 1 && x < 1) return(matrix(x, 1, 1)) else return(diag(x))
#************ mdiag 

#****f* miscUtils/repairfrailtypar
#  NAME
#    repairfrailtypar --- fix frailty parameters
#  FUNCTION
#   For identifiability, it is necessary to fix one of the frailty spline
#   parameters at 0, and not include it in optimization and estimation
#   routines. This function adds a 0 to a vector at a specified index.
#  INPUTS
#   par     a vector
#   ind     integer index
#  OUTPUTS
#   The input vector par, with a 0 inserted at its ind entry
#  SYNOPSIS
repairfrailtypar <- function(par, ind)
#  SOURCE
#
{
    if(ind == 1) return(c(0, par))
    if(ind == length(par)) return(c(par, 0))
    return(c(par[1:(ind - 1)], 0, par[ind:length(par)]))
}
#************ repairfrailtypar 

#****f* initRoutine/inverthessian
#  NAME
#    inverthessian --- invert a Hessian matrix
#  FUNCTION
#   The Hessian matrices returned by optim() are not always positive-definite
#   and invertible. This function adjusts the diagonal of an input matrix until
#   it is invertible. This gives an acceptable covariance matrix for subsequent use.
#  INPUTS
#    hess   a matrix
#  OUTPUTS
#    Sigma  the inverse of hess, or of something very close to it
#  SYNOPSIS
inverthessian <- function(hess){
#  SOURCE
#
    K <- dim(hess)[1]
    # Try to invert the Hessian
    Sigma <- try(solve(-hess), silent = TRUE)
    d <- 10
    # modify the diagonal until it is invertible
    while(inherits(Sigma, "try-error")){ 
        Sigma <- try(solve(-(hess - 10^(-d) * mdiag(K))), silent = TRUE);
        d <- d - 1
    }
    # modify the diagonal until it is positive definite
    while(!all(eigen(Sigma)$val > 0)){ 
        Sigma <- solve(-(hess - 10^(-d) * mdiag(K))) ;
        d <- d - 1  
    }
    return(Sigma)
}
#************ inverthessian 

#****f* initRoutine/numHess.par
#  NAME
#    numHess.par  --- compute a numerical Hessian for the parametric component
#  FUNCTION
#    Since there are so many parametric component specifications and parametrizations,
#    it is easiest to compute the Hessian of the parameters numerically. Given a set
#    of parameters and a likelihood function, this routine computes the Hessian of
#    the function at the given parameters.
#
#    The method used is basic finite differences
#  INPUTS
#    param.par  parameters of the parametric component
#    fun        likelihood function for the parametric component
#    eps        precision to use for numerical Hessian
#  OUTPUTS
#  SYNOPSIS
numHess.par <- function(param.par, fun, eps = 1e-5, ...)
#  SOURCE
#
{
    # Utility function to compute numerical derivatives
    numDer.par <- function(param.par, fun, eps = 1e-5, ...)
    {        
        lik1 <- fun(param.par, ...)
        nd <- rep(0, length(param.par))
        #   take finite differences along the set of parameters
        for(i in 1:length(nd))
        {
            param.par2 <- param.par
            param.par2[i] <- param.par2[i] + eps
            lik2 <- fun(param.par2, ...)
            nd[i] <- (lik2 - lik1) / eps
        }
        return(nd)
    }
    # allocate storage for numerical Hessian
    nh <- matrix(0, length(param.par), length(param.par))
    #   base gradient
    gr1 <- numDer.par(param.par, fun, eps, ...)
    #   compute gradients by finite differences
    for(i in 1:length(param.par))
    {
        param.par2 <- param.par
        param.par2[i] <- param.par2[i] + eps
        gr2 <- numDer.par(param.par2, fun, eps, ...)
        nh[i, ] <- (gr2 - gr1) / eps
    }
    return(nh)
}
#************ numHess.par 

#****f* initRoutine/inithistory
#  NAME
#    inithistory --- initialize the history structure
#  FUNCTION
#    The history structure keeps track of the posterior samples from the
#    chain. It contains components for all the parameters that change
#    in the course of the chain. The function updatehistory is called
#    at the end of each iteration to update this structure.
#  INPUTS
#    hazard     a hazard RCurve
#    frailty    a frailty RCurve
#    regression a RRegression structure
#    control    an RControl structure
#  OUTPUTS
#    an RHistory structure
#  SYNOPSIS
inithistory <- function(hazard, frailty, regression, control)
#  SOURCE
#
{
    history <- NULL
    maxiter <- control$maxiter
    history$frailty <- matrix(0, maxiter, length(frailty$x))
    history$coefficients <- matrix(0, maxiter, length(regression$coefficients))
    history$loglik <- matrix(0, maxiter, 1)
    # ssorate for spline knots and parameters
    if(hazard$hasspline) {
        history$hazard.spline.par <- matrix(-Inf, maxiter, 
            hazard$spline.maxoccknots + hazard$spline.ord)
        history$hazard.spline.knots <- matrix(-Inf, maxiter, 
            hazard$spline.maxoccknots + 2 * hazard$spline.ord)
    }
    if(frailty$hasspline) {
        history$frailty.spline.par <- matrix(-Inf, maxiter, 
            frailty$spline.maxoccknots + frailty$spline.ord)
        history$frailty.spline.knots <- matrix(-Inf, maxiter, 
            frailty$spline.maxoccknots + 2 * frailty$spline.ord)
        history$frailty.spline.fvar <- matrix(0, maxiter, 1)
    }
    # storage for parametric parameters
    if(hazard$haspar) history$hazard.param.par <- matrix(0, 
        maxiter, length(hazard$param.par))
    if(frailty$haspar) history$frailty.param.par <- matrix(0, 
        maxiter, length(frailty$param.par))
    # storage for weights
    if(hazard$hasspline & hazard$haspar) history$hazard.weight <- matrix(0, maxiter, 1)
    if(frailty$hasspline & frailty$haspar) history$frailty.weight <- matrix(0, maxiter, 1)
    # storage for prior variances
    history$priorvar <- matrix(0, maxiter, 7); 
    colnames(history$priorvar) <- c("coefficients", "hazard.spline", "frailty.spline",
            "hazard.param", "frailty.param", "hazard.weight", "frailty.weight")
    # storage for acceptance history
    history$accept <- matrix(0, maxiter, 8)
    colnames(history$accept) <- c(colnames(history$priorvar), "frailty")
    # update with first iteration
    history <- updatehistory(history, 1, hazard, frailty, regression)
    return(history)
}
#************ inithistory 

#****f* curveUpdate/updatehistory
#  NAME
#    updatehistory --- update RHistory structure
#  FUNCTION
#    Store the state of the chain at the current iteration in the RHistory
#    structure.
#  INPUTS
#    history        an RHistory structure
#    i              integer, current iteration
#    hazard         a hazard RCurve
#    frailty        a frailty RCurve
#    regression     a RRegression structure
#  OUTPUTS
#    history        updated Rhistory
#  SYNOPSIS
updatehistory <- function(history, i, hazard, frailty, regression)
#  SOURCE
#
{
    # store frailties and coefficients
    history$frailty[i, ] <- frailty$x
    history$coefficients[i, ] <- regression$coefficients
    history$loglik[i, ] <- regression$loglik
    # store spline knots and parameters
    if(hazard$hasspline) {
        history$hazard.spline.par[i, 1:length(hazard$spline.par)] <- hazard$spline.par
        history$hazard.spline.knots[i, 1:length(hazard$spline.knots)] <- hazard$spline.knots
    }
    if(frailty$hasspline) {
        history$frailty.spline.par[i, 1:length(frailty$spline.par)] <- frailty$spline.par
        history$frailty.spline.knots[i, 1:length(frailty$spline.knots)] <- frailty$spline.knots
        history$frailty.spline.fvar[i] <- frailty$spline.fvar
    }
    # store parametric components and weights
    if(hazard$haspar) history$hazard.param.par[i, ] <- hazard$param.par
    if(frailty$haspar) history$frailty.param.par[i, ] <- frailty$param.par
    if(hazard$hasspline & hazard$haspar) history$hazard.weight[i] <- hazard$weight
    if(frailty$hasspline & frailty$haspar) history$frailty.weight[i] <- frailty$weight
    # store prior variance terms
    history$priorvar[i, ] <- c(regression$priorvar, hazard$spline.priorvar, frailty$spline.priorvar,
        hazard$param.priorvar, frailty$param.priorvar, hazard$weight.priorvar, frailty$weight.priorvar)
    # store acceptance history
    history$accept[i, ] <- c(regression$accept, hazard$spline.accept, frailty$spline.accept,
        hazard$param.accept, frailty$param.accept, hazard$weight.accept, frailty$weight.accept, frailty$accept)

    return(history)
}
#************ updatehistory 

#****f* simSurvival/rinvgamma
#  NAME
#    rinvgamma --- generate inverse-gamma random variates
#  FUNCTION
#    Generate inverse-gamma random numbers
#  INPUTS
#    n      number of variates to generate
#    shape  shape parameter
#    scale  scale parameter
#  OUTPUTS
#    n inverse-gamma random numbers
#  SYNOPSIS
rinvgamma <- function(n, shape, scale = 1)
#  SOURCE
#
{

    return(1 / rgamma(n, shape, rate = scale))
} 
#************ rinvgamma 

#****f* simSurvival/dnegbin
#  NAME
#    dnegbin --- negative binomial distribution
#  FUNCTION
#    Compute the negative binomial distribution function. This is
#    for use as a prior on the number of knots for adaptive selection.
#  INPUTS
#    x  the values at which the distribution should be computed
#    r  number of trials
#    p  probability of success
#  OUTPUTS
#    The NB(r,p) distribution evaluated at x
#  SYNOPSIS
dnegbin <- function(x, r, p) 
#  SOURCE
#
{
    out <- gamma(x + r) / factorial(x) / gamma(r) * p^r * (1 - p)^x
    out
}

#************ dnegbin 

#****f* miscUtils/accrate.predict.lm
#  NAME
#    accrate.predict.lm --- predicted acceptance rate distance from 25%
#  FUNCTION
#   This function is used to automatically select the tuning parameters.
#  INPUTS
#    m      a model returned by lm()
#    x      a tuning parameter
#  OUTPUTS
#    the predicted squared difference between the model m evaluated at x
#    and 25%
#  SYNOPSIS
accrate.predict.lm <- function(x, m) (predict(m, data.frame(x = x))-.25)^2
#************ accrate.predict.lm 

#****f* miscUtils/submean
#  NAME
#    submean --- compute the mean of a subset of a vector
#  FUNCTION
#    For a vector or matrix x and a subset, compute the mean of the subset of values
#    of x.
#  INPUTS
#    x          a vector or matrix
#    subset     a subset, either as a range or as logical
#    f          a function to be applied (mean by default)
#  OUTPUTS
#    The function f applied to a subset of x
#  SYNOPSIS
submean <- function(x, subset, f = mean) 
#  SOURCE
#
{
    if(is.null(x)) return(NULL)
    if(!is.null(dim(x))) return(apply(x[subset, ,drop = F], 2, f))
    else return(f(x[subset]))
}
#************ submean 

#****f* miscUtils/makeoutputcurve
#  NAME
#    makeoutputcurve --- construct the curve to be returned in the output
#  FUNCTION
#    The RCurve construct contains many values that change with each iteration.
#    This function cleans the Rcurve and retains only a subset for output.
#  INPUTS
#    curve      an RCurve, either hazard or frailty
#  OUTPUTS
#    outcurve   a structure containing only a subset of the components of the input
#  SYNOPSIS
makeoutputcurve <- function(curve)
#  SOURCE
#
{
    outcurve <- list(name = curve$name,
                    type = curve$type,
                    spline.adaptive = curve$spline.adaptive,
                    spline.nknots = curve$spline.nknots,
                    spline.nknots.prior = curve$spline.nknots.prior,
                    spline.maxoccknots = curve$spline.maxoccknots,
                    spline.nknots.hyper = curve$spline.nknots.hyper,
                    spline.knotspacing = curve$spline.knotspacing,
                    spline.ord = curve$spline.ord,
                    spline.norm = curve$spline.norm,
                    spline.penalty = curve$spline.penalty,
                    spline.penaltyfactor = curve$spline.penaltyfactor,
                    spline.hyper = curve$spline.hyper,
                    param.dist = curve$param.dist,
                    param.hyper = curve$param.hyper,
                    weight.hyper = curve$weight.hyper
                )
    return(outcurve)   
}
#************ makeoutputcurve 

#****f* initRoutine/nknotsPriorMean
#  NAME
#    nknotsPriorMean --- compute the prior mean of the number of knots of a curve
#  FUNCTION
#    If adaptive selection is used, the algorithm initializes the number of spline
#    knots at the prior mean. This function computes the prior mean number of
#    spline knots for different priors
#  INPUTS
#    curve  an RCurve structure
#  OUTPUTS
#    the prior mean number of knots, for type spline.nknots.prior and parameters
#    spline.nknots.hyper
#  SYNOPSIS
nknotsPriorMean <- function(curve)
#  SOURCE
#
{
    if(curve$spline.nknots.prior == "poisson")
        return(curve$spline.nknots.hyper)
    if(curve$spline.nknots.prior == "geometric")
        return(round(1 / curve$spline.nknots.hyper))
    if(curve$spline.nknots.prior == "poissonmix")
        return(round(mean(curve$spline.nknots.hyper)))
    if(curve$spline.nknots.prior == "negbin")
        return(round(weighted.mean(1:curve$spline.maxoccknots, 
            dnegbin(1:curve$spline.maxoccknots, curve$spline.nknots.hyper[1],
            curve$spline.nknots.hyper[2]))))
    if(curve$spline.nknots.prior == "power")
        return(round(weighted.mean(1:curve$spline.maxoccknots,
            (1:curve$spline.maxoccknots)^curve$spline.nknots.hyper)))
}
#************ nknotsPriorMean 

#****f* splineUtils/nknotsPrior
#  NAME
#    nknotsPrior --- evaluate the prior on the number of knots
#  FUNCTION
#    In the course of RJMCMC for selecting the number of knots, the prior probability
#    of a given number of knots must be evaluated.
#  INPUTS
#    x      number of knots at which to evaluate the prior
#    curve  an RCurve structure, with spline.nknots.prior and spline.nknots.hyper components
#  OUTPUTS
#  SYNOPSIS
nknotsPrior <- function(x, curve)
#  SOURCE
#
{
    if(curve$spline.nknots.prior == "poisson")
        return(dpois(x, curve$spline.nknots.hyper))
    if(curve$spline.nknots.prior == "geometric")
        return(dgeom(x, curve$spline.nknots.hyper))
    if(curve$spline.nknots.prior == "poissonmix")
        return(mean(dpois(x, curve$spline.nknots.hyper)))
    if(curve$spline.nknots.prior == "negbin")
        return(dnegbin(x, curve$spline.nknots.hyper[1], curve$spline.nknots.hyper[2]))
    if(curve$spline.nknots.prior == "power")
        return(x^curve$spline.nknots.hyper)
}
#************ nknotsPrior 

}}}

{{{ #Initialize
##############################################################
# \section{Initialize} Initialize knots, penalties, etc for spline components
##############################################################

#****f* initRoutine/makeknots
#  NAME
#    makeknots --- make knots for a curve with a spline component
#  FUNCTION
#    Automatically initialize the set of spline knots if they are not given
#    in the input. Also initializes candidate knots for adaptive knot selection.
#  INPUTS
#    curve      an RCurve structure
#    x          a set of data points to be used for constructing knots
#    bounds     optional boundary knots (length 2 vector)
#  OUTPUTS
#    the input RCurve, with additional spline.knots and spline.candknots components.
#  SYNOPSIS
makeknots <- function(curve, x, bounds = NULL)
#  SOURCE
#
{
    if(!curve$hasspline) return(curve)
    #   extract needed curve components
    BUF <- 0.01     # boundary buffer 
    knots <- curve$spline.knots;
    nknots <- curve$spline.nknots;
    ncandknots <- curve$spline.ncandknots;
    knotspacing <- curve$spline.knotspacing;
    adaptive <- curve$spline.adaptive
    ord <- curve$spline.ord
    candknots <- NULL
    K <- ord + nknots
    if(is.null(bounds)) {   # this version requires boundary knots to be given
        browser()
    }
    if(adaptive) nintknots <- ncandknots else nintknots <- nknots
    if(is.null(knots)){
        # distribute knots and candidate knots as quantiles of the data
        if(knotspacing == "quantile"){
            ibounds <- c(min(x), max(x))
            lrep <- ord; rrep <- ord
            if(ibounds[1] == bounds[1]) {nintknots <- nintknots + 1; lrep <- ord - 1}
            if(ibounds[2] == bounds[2]) {nintknots <- nintknots + 1; rrep <- ord - 1}
            candknots <- quantile(unique(x), seq(from = BUF, to = 1 - BUF, length = nintknots))
            # select the occupied knots as a random subset of the candidate knots
            occknots <- sort(sample(1:nintknots, nknots))
            knots <- candknots[occknots]
            candknots <- c(rep(bounds[1], lrep), candknots, rep(bounds[2], rrep))
            knots <- c(rep(bounds[1], lrep), knots, rep(bounds[2], rrep))
            attr(candknots, "occupied") <- c(rep(2, lrep), (1:nintknots)%in%occknots, 
                rep(2, rrep))
        }
        # distribute knots and candidate knots equally over the data range
        if(knotspacing == "equal"){
            dbounds <- diff(bounds)
            # distribute candidate knots equally
            candknots <- seq(from = bounds[1] + BUF * dbounds, to = bounds[2] - BUF * dbounds,
                length = nintknots + 2)
            candknots <- candknots[ - 1];candknots <- candknots[ - length(candknots)]
            occknots <- sort(sample(1:nintknots, nknots))
            # select the occupied knots as a random subset of the candidate knots
            knots <- candknots[occknots]
            knots <- c(rep(bounds[1], ord), knots, rep(bounds[2], ord))
            candknots <- c(rep(bounds[1], ord), candknots, rep(bounds[2], ord))
            attr(candknots, "occupied") <- c(rep(2, ord), (1:nintknots)%in%occknots, rep(2, ord))
        }
        # half of the knots are equally distributed, the other half are quantiles
        if(knotspacing == "mixed"){
            dbounds <- diff(bounds)
            # quantile candknots
            candknots1 <- quantile(unique(x), seq(from = BUF, to = 1 - BUF,
                length = floor(nintknots / 2)))
            candknots2 <- seq(from = bounds[1] + BUF * dbounds, to = bounds[2] - BUF * dbounds,
                length = ceiling(nintknots / 2))
            candknots <- sort(sample(unique(c(candknots1, candknots2)), nintknots))
            occknots <- sort(sample(1:nintknots, nknots))
            knots <- candknots[occknots]
            knots <- c(rep(bounds[1], ord), knots, rep(bounds[2], ord))
            candknots <- c(rep(bounds[1], ord), candknots, rep(bounds[2], ord))
            attr(candknots, "occupied") <- c(rep(2, ord), (1:nintknots)%in%occknots,
                rep(2, ord))
        }
    }
    # attributes for the knots object
    attr(knots, "boundary") <- bounds
    attr(knots, "index") <- seq(from=-(ord - 1), length = length(knots), by = 1)
    attr(knots, "order") <- ord
    curve$spline.knots <- knots
    curve$spline.candknots <- candknots
    return(curve)
}
#************ makeknots 

#****f* splineUtils/makesplinebasis
#  NAME
#    makesplinebasis --- construct B-spline basis functions
#  FUNCTION
#    Compute the spline.basis, spline.basiscum and spline.basisexp components of an
#    RCurve with a spline component.
#  INPUTS
#    curve      an RCurve structure
#    quick      if TRUE, only the basis itself is computed, not the other two components
#    usec       boolean, whether to use C code for fast computation
#  OUTPUTS
#    curve      an RCurve with updated spline.basis, spline.basisint,
#               spline.basiscum and spline.basisexp components
#  SYNOPSIS
makesplinebasis <- function(curve, quick = FALSE, usec = TRUE)
#  SOURCE
#
{
    if(!curve$hasspline) return(curve)
    knots <- curve$spline.knots; ord <- curve$spline.ord; x <- curve$x
    # Evaluate the spline basis at the set of x values of the curve
    if(usec) B <- csplinedesign(knots, x = x, ord = ord) else 
        B <- splineDesign(knots, x = x, ord = ord)
    # Compute the integral of each of the basis functions
    if(usec) Bint <- cevalBinte(knots, ord) else Bint <- evalBinte(knots, ord)
    if(curve$spline.norm) for(i in 1:dim(B)[1]) B[i, ] <- B[i, ] / Bint
    if(!quick) { 
        if(curve$name == "hazard")  {
            #  compute cumulative integrals of each basis function to the x values
            if (usec) C <- cevalCinte(knots, ord, x, Bint)
            else C <- evalCinte(knots, ord, x, Bint)
        }
        else C <- NULL
        if(curve$name == "frailty") {
            # compute one minus the expectation over each of the basis functions
            if(usec) E <- cevalEinte(knots, ord)
            else E <- evalEinte(knots, ord)
        }
        else E <- NULL
        curve$spline.basiscum <- C
        curve$spline.basisexp <- E
    }
    curve$spline.basisint <- Bint
    curve$spline.basis <- B
    return(curve)
}
#************ makesplinebasis 

#****f* splineUtils/ki
#  NAME
#    ki --- addressing of spline indices
#  FUNCTION
#    allows addressing spline parameters by their "true" indices (which can be negative)
#    rather than by their vector indices
#  SYNOPSIS
ki <- function(knots, j) return(j + attr(knots, "ord"))
#************ ki 

#****f* splineUtils/evalBinte
#  NAME
#    evalBinte --- compute the integrals of each spline basis function
#  FUNCTION
#    The set of spline basis functions are defined entirely by the set of knots
#    and the order of the spline. This function computes the integral of each of 
#    those basis functions, which is used for normalizing the frailty basis splines
#    and for computing cumulative integrals for the hazard.
#  INPUTS
#    knots  a set of spline knots as produced by makeknots
#    ord    integer order of the spline
#  OUTPUTS
#    a vector containing the integral over each spline basis function
#  SYNOPSIS
evalBinte <- function(knots, ord)
#  SOURCE
#
{
    K <- sum(knots > attr(knots, "b")[1] & knots < attr(knots, "b")[2])
    Binte <- rep(0, K + ord)
    for(j in 1:(K + ord)){
        Binte[j] <- (knots[ki(knots, j)] - knots[ki(knots, j - ord)]) / ord
    }
    return(Binte)
}
#************ evalBinte 

#****f* splineUtils/evalCinte
#  NAME
#    evalCinte --- compute partial integrals over the spline basis
#  FUNCTION
#    In order to compute the cumulative baseline hazard, the integral of each spline
#    basis function from 0 to each observation must be computed.
#  INPUTS
#    knots      a set of spline knots as output by makeknots
#    ord        integer spline order
#    obs        vector of observations at which the integrals should be evaluated
#    Binte      basis function integrals produced by evalBinte
#  OUTPUTS
#    a matrix of size length(obs) x N (where N is the number of basis functions)
#    containing in entry (i,j) the integral of basis function j from 0 to obs[i]
#  SYNOPSIS
evalCinte <- function(knots, ord, obs, Binte)
#  SOURCE
#
{
    K <- sum(knots > attr(knots, "b")[1] & knots < attr(knots, "b")[2])
    Cinte <- matrix(0, length(obs), K + ord)
    # Compute a spline basis of order ord+1
    knots2 <- c(attr(knots, "b")[1], knots, attr(knots, "b")[2])
    attr(knots2, "i") <- c(min(attr(knots, "i")) - 1, attr(knots, "i"), 
        max(attr(knots, "i")) + 1)
    attr(knots2, "b") <- attr(knots, "b")
    Bordp1 <- splineDesign(knots2, x = obs, outer.ok = TRUE, ord = ord + 1)
    # compute the integrals
    for(i in 1:length(obs)){
        for(j in 1:(K + ord)){
            # If obs is greater than the rightmost support point, return the full integral
            if(obs[i] >= knots[ki(knots, j)]) Cinte[i, j] <- Binte[j]
            # otherwise use the formula for the partial integral
            if(obs[i] < knots[ki(knots, j)] & obs[i] >= knots[ki(knots, j - ord)]) 
                Cinte[i, j] <- Binte[j] * sum(Bordp1[i, (j + 1):(K + ord + 1)])
        }
    }
    return(Cinte)
}
#************ evalCinte 

#****f* splineUtils/evalEinte
#  NAME
#    evalEinte --- compute the expected distance from 1 of a B-spline
#  FUNCTION
#    Since the frailty density must have mean 1, it is important to be able to 
#    compute the difference between 1 and the frailty density mean. This function
#    calculates the expected values of a set of normalized B-spline basis functions
#    and subtracts them from 1. The difference between 1 and the frailty mean is then
#       curve$spline.basisexp%*%exp(curve$spline.par)
#    where the former is the output of this function and the latter is the set of
#    spline weights.
#  INPUTS
#    knots      a set of spline knots as output by makeknots
#    ord        integer spline order
#  OUTPUTS
#    a vector containing one minus the expectation of each of the normalized basis functions
#  SYNOPSIS
evalEinte <- function(knots, ord)
#  SOURCE
#
{
    K <- sum(knots > attr(knots, "b")[1] & knots < attr(knots, "b")[2])
    Einte <- rep(0, K + ord)
    for(j in 1:(K + ord)){
        # Einte[j] contains the 1st moment of the j-th spline of order ord defined
        # on a given set of knots
        Einte[j] <- nBsmom(1, ord, j, knots)
    }
    return(1 - Einte)
}
#************ evalEinte 

#****f* splineUtils/nBsmom
#  NAME
#    nBsmom --- compute the N-th moment of a B-spline basis function
#  FUNCTION
#    Computes the N-th moment of the j-th B-spline basis function of order ord defined
#    on the given set of knots.
#  INPUTS
#    N       moment to compute
#    ord     order of the spline
#    j       which basis function to compute the moment for
#    knots   set of knots to use
#  OUTPUTS
#    double, containing the N-th B-spline moment
#  SYNOPSIS
nBsmom <- function(N, ord, j, knots)
#  SOURCE
#
{
    lknot <- knots[ki(knots, j - ord)]
    rknot <- knots[ki(knots, j)]
    if(lknot == rknot) return(rknot^N)
    # This is a magic recursive formula that can be computed easily using integration
    # by parts
    return(ord / (rknot - lknot) / (N + 1) * (-nBsmom(N + 1, ord - 1, j - 1, knots) + 
        nBsmom(N + 1, ord - 1, j, knots)))
}
#************ nBsmom 

#****f* splineUtils/splineconv
#  NAME
#    splineconv --- compute the convolution of two splines
#  FUNCTION
#    This is needed to construct a penalty on the integrated squared second derivative.
#    This function computes the convolution of two spline basis functions, that is,
#       int_0^infty ( x^k B_1(x) * B_2(x) ) dx
#    where the splines may be of different orders, but are defined on the same set of knots.
#  INPUTS
#    k      the power of x used in the convolution
#    n1     order of the first spline
#    j1     largest knot of the first spline
#    n2     order of the second spline
#    j2     largest knot of the second spline
#    knots  set of knots on which the splines are defined
#  OUTPUTS
#    The k-th order convolution of the splines defined by the input parameters
#  SYNOPSIS
splineconv <- function(k, n1, j1, n2, j2, knots)
#  SOURCE
#
{
    # if the splines don't overlap, the convolution is 0
    if(j1 - n1 >= j2 | j2 - n2 >= j1) return(0)
    # if both splines are first-order, the convolution is trivial
    if(n1 == 1 & n2 == 1){
        out <- 1 / (k + 1) * (knots[ki(knots, j1)]^(k + 1) - knots[ki(knots, j1 - 1)]^(k + 1))
        return(out)
    }
    # By symmetry, we can assume that n1>n2 wlog. If this is not the case, switch them around
    if(n2 > n1){
        n3 <- n1; n1 <- n2; n2 <- n3
        j3 <- j1; j1 <- j2; j2 <- j3
    }
    # use a magic formula that can be derived by integration by parts
    out <- 0
    denom1 <- knots[ki(knots, j1 - 1)] - knots[ki(knots, j1 - n1)]
    denom2 <- knots[ki(knots, j1)] - knots[ki(knots, j1 - n1 + 1)]
    if(denom1 > 0){
        out <- out + 1 / denom1 * splineconv(k + 1, n1 - 1, j1 - 1, n2, j2, knots)
        out <- out - knots[ki(knots, j1 - n1)] / denom1 * 
            splineconv(k, n1 - 1, j1 - 1, n2, j2, knots)
    }
    if(denom2 > 0){
        out <- out + knots[ki(knots, j1)] / denom2 * 
            splineconv(k, n1 - 1, j1, n2, j2, knots)
        out <- out - 1 / denom2 * splineconv(k + 1, n1 - 1, j1, n2, j2, knots)
    }
    return(out)
}
#************ splineconv 

#****f* splineUtils/splinederivint
#  NAME
#    splinederivint --- compute the convolution of the derivatives of two spline bases
#  FUNCTION
#    Used to compute the penalty matrix for the integrated squared second derivative. This
#    routine computes the integral from 0 to infinity of the l1 derivative and the l2
#    derivative of the j1 and j2-th splines of order n1 and n2 defined on a set of knots.
#    For instance, in order to compute the [j1,j2] entry of the penalty matrix,
#    the function makePenalty.2deriv calls
#       out[j1, j2] <- splinederivint(2, ord, j1, 2, ord, j2, knots)
#  INPUTS
#    l1     derivative of the first spline
#    n1     order of the first spline
#    j1     index of the first spline
#    l2     derivative of the second spline
#    n2     order of the second spline
#    j2     index of the second spline
#    knots  set of knots on which the splines are defined
#  OUTPUTS
#    a double containing the convolution of the derivatives of two basis functions
#  SYNOPSIS
splinederivint <- function(l1, n1, j1, l2, n2, j2, knots)
#  SOURCE
#
{
    # if they don't overlap, the integral is 0
    if(j1 - n1 >= j2 | j2 - n2 >= j1) return(0)
    # For the 0th derivatives, use the regular convolution method
    if(l1 == 0 & l2 == 0) return(splineconv(0, n2, j1, n2, j2, knots))
    # symmetry
    if(l2 > l1){
        l3 <- l1; l1 <- l2; l2 <- l3
        n3 <- n1; n1 <- n2; n2 <- n3
        j3 <- j1; j1 <- j2; j2 <- j3
    }
    out <- 0
    # recursive method to step-by-step reduce this problem to a regular convolution problem
    denom1 <- knots[ki(knots, j1 - 1)] - knots[ki(knots, j1 - n1)]
    denom2 <- knots[ki(knots, j1)] - knots[ki(knots, j1 - n1 + 1)]
    if(denom1 > 0) out <- out + (n1 - 1) / denom1 *
        splinederivint(l1 - 1, n1 - 1, j1 - 1, l2, n2, j2, knots)
    if(denom2 > 0) out <- out - (n1 - 1) / denom2 *
        splinederivint(l1 - 1, n1 - 1, j1, l2, n2, j2, knots)
    return(out)
}
#************ splinederivint 

#****f* splineUtils/mysplineDesign
#  NAME
#    mysplineDesign --- works like splineDesign()
#  FUNCTION
#    a plug-in replacement for splineDesign() in the splines package. It calls the
#    same internal code, but skips the input error checking and reformatting. This should
#    only be called with known-good input.
#  INPUTS
#    knots  knots in the format required by splineDesign
#    x      set of observations at which the basis should be evaluated
#    ord    spline order
#  OUTPUTS
#    a spline basis matrix in which entry (i,j) contains the j-th basis function
#    evaluated at x[i]
#  SYNOPSIS
mysplineDesign <- function (knots, x, ord = 4) 
#  SOURCE
#
{
    nk <- length(knots)
    x <- as.numeric(x)
    derivs <- integer(1)
    # call the internal function used by splineDesign()
    temp <- .Call("spline_basis", knots, ord, x, derivs, PACKAGE = "splines")
    ncoef <- nk - ord
    design <- rep(0, ncoef)
    jj <- 1:ord + attr(temp, "Offsets")
    design[jj] <- temp
    dim(design) <- c(1, ncoef)
    design
}
#************ mysplineDesign 

#****f* initRoutine/makepenalty
#  NAME
#    makepenalty --- construct a penalty matrix
#  FUNCTION
#    Construct a penalty matrix for use with penalized spline fitting. Options are
#    a penalty on the squared second differences of the spline parameters, or a penalty
#    on the integrated squared second derivative.
#  INPUTS
#    curve      an RCurve structure
#    usec       boolean, whether to use fast C code
#  OUTPUTS
#    curve      the input curve, with spline.penaltymatrix component updated
#  SYNOPSIS
makepenalty <- function(curve, usec = TRUE)
#  SOURCE
#
{
    if(!curve$hasspline) return(curve)
    penalty <- curve$spline.penalty
    ord <- curve$spline.ord; nknots <- curve$spline.nknots; knots <- curve$spline.knots
    # second difference penalty
    if(penalty == "2diff") P <- makePenalty.2diff(ord + nknots)
    # second derivative penalty
    if(penalty == "2deriv" | penalty == "log2deriv") {
        if(usec) P <- cmakePenalty.2deriv(ord, knots)
        else  P <- makePenalty.2deriv(ord, knots)
        # adjust for normalized B-splines for frailties
        if(curve$spline.norm){
            Bint <- curve$spline.basisint
            P <- P / (Bint%*%t(Bint))
        }
    }
    if(penalty == "none") P <- 0
    curve$spline.penaltymatrix <- P
    return(curve)
}
#************ makepenalty 

#****f* initRoutine/makePenalty.2diff
#  NAME
#    makePenalty.2diff --- compute a penalty matrix on second differences
#  FUNCTION
#    This computes a matrix P such that 
#       x %*% P %*% x
#    is the sum of squared second differences of x.
#  INPUTS
#    K      the desired size of the matrix
#  OUTPUTS
#    P      a K x K matrix that generates squared second differences.
#  SYNOPSIS
makePenalty.2diff <- function(K){
#  SOURCE
#
    D <- matrix(0, K - 2, K)
    for(i in 1:(K - 2)){
        D[i, i] <- 1
        D[i, i + 1] <- -2
        D[i, i + 2] <- 1
    }
    P <- t(D)%*%D
    return(P)
}
#************ makePenalty.2diff 

#****f* initRoutine/makePenalty.2deriv
#  NAME
#    makePenalty.2deriv --- compute a penalty matrix on second derivatives of B-splines
#  FUNCTION
#    This function computes a matrix P such that
#       exp(x) %*% P %*% exp(x)
#    is the integral of the squared second derivative of a B-spline with a given set of
#    knots and component weights exp(x)
#  INPUTS
#    ord    order of the B-spline
#    knots  set of basis knots
#  OUTPUTS
#    P      a K x K matrix that penalizes the integrated squared second derivative
#  SYNOPSIS
makePenalty.2deriv <- function(ord, knots)
#  SOURCE
#
{
    #   compute the number of spline components K
    nspline <- sum(knots > attr(knots, "b")[1])
    out <- matrix(0, nspline, nspline)
    for(j1 in 1:nspline){
        for(j2 in j1:nspline){
            # compute convolutions of second derivatives for each pair of basis functions
            out[j1, j2] <- splinederivint(2, ord, j1, 2, ord, j2, knots)
        }
    }
    # the matrix is symmetric
    for(j1 in 2:nspline){
        for(j2 in 1:(j1 - 1)) out[j1, j2] <- out[j2, j1]
    }
    return(out)
}
#************ makePenalty.2deriv 

#****f* makeLikelihood/smoothpen
#  NAME
#    smoothpen --- compute the smoothness penalty for a curve
#  FUNCTION
#    Penalize lack of smoothness of a B-spline curve as part of a likelihood computation
#  INPUTS
#    curve      an RCurve structure
#    der        the derivative of the smoothness penalty to compute
#  OUTPUTS
#    value of the smoothness penalty
#  SYNOPSIS
smoothpen <- function(curve, der = 0)
#  SOURCE
#
{
    type <- curve$spline.penalty
    name <- curve$name
    theta <- curve$spline.par
    # extract the penalty matrix
    P <- curve$spline.penaltymatrix
    sigma2 <- curve$spline.priorvar
    if(der >= 2) stop("second derivative not implemented")
    # second difference penalty
    if(type == "2diff"){
        if(der == 0) return(max( t(theta)%*%P%*%theta / (2 * sigma2), 0))
        if(der == 1) return( P%*%theta /sigma2)
    }
    # second derivative penalty
    if(type == "2deriv"){
        et <- exp(theta)
        if(der == 0) return(max( t(et)%*%P%*%et / (2 * sigma2), 0))
        if(der == 1) return( mdiag(et)%*%P%*%et / sigma2 )
    }
    # penalty on the log 2nd derivative
    if(type == "log2deriv"){
        et <- exp(theta)
        ePe <- as.numeric(t(et)%*%P%*%et) 
        if(der == 0) return(max(log(ePe + 1)/ (2 * sigma2), 0))
        if(der == 1) return( mdiag(et)%*%P%*%et / sigma2 /(ePe + 1))
    }
    # gaussian "penalty", which isn't really a smoothness penalty at all.
    if(type == "none"){
        if(der == 0) return(theta%*%theta / (2 * sigma2))
        if(der == 1) return( theta / sigma2)
    }
}
#************ smoothpen 


}}}

{{{ #CurveUpdate
##############################################################
# \section{CurveUpdate} Curve updating routines
##############################################################

#****f* initRoutine/fitparametric
#  NAME
#    fitparametric --- fit a parametric component to a curve
#  FUNCTION
#    Given a curve with a parametric component, compute a good set of initial
#    values for the parametric component parameters. The parametric component
#    distributions are parametrized in a way that allows Gaussian priors on the
#    parameters, and this function incorporates that.
#
#    See also evalparametric for the parametrization used.
#  INPUTS
#    curve      an RCurve structure with a parametric component
#    x          a set of data points used for estimation
#  OUTPUTS
#    curve      updated curve with curve$param.par holding initial values
#  SYNOPSIS
fitparametric <- function(curve, x)
#  SOURCE
#
{
    name <- curve$name
    dist <- curve$param.dist
    if(dist == "none") return(curve)
    # frailty curve
    if(name == "frailty")
    {
        # compute the variance of the frailties and set the parameter
        # according to the parametrization selected.
        Ui <- x
        if(dist == "gamma")
            par <- log(var(Ui))
        if(dist == "lognormal")
        {
            varu <- var(Ui)
            par <- log(log(varu + 1))
        }
        curve$param.par <- par
        curve$x <- Ui
    }
    # hazard curve
    if(name == "hazard")
    {
        agdata <- x
        varnames <- colnames(agdata)[ - (1:4)]
        qvarnames <- paste("`", varnames, "`", sep = "")
        # use survreg to fit a parametric component to the hazard
        # and transform the estimated parameters to the parametrization
        fit <- survreg(as.formula(paste("Surv(time, delta)~", 
            paste(qvarnames, collapse = " + "))), data = agdata, dist = dist)
        if(dist == "exponential"){
            par <- log(fit$icoef[1])
        }
        if(dist == "weibull"){
            lambda <- exp(-fit$icoef[1])
            gamma <- 1 / exp(fit$icoef[2])
            par <- c(log(lambda), log(gamma))
        }
        if(dist == "lognormal"){
            par <- c(fit$icoef[1], fit$icoef[2])
        }
        names(par) <- NULL
        curve$param.par <- par
        curve$x <- agdata$time
    }
    # Evaluate the curve at the parameter values chosen.
    curve <- evalparametric(curve)
    return(curve)
}
#************ fitparametric 

#****f* curveUpdate/evalparametric
#  NAME
#    evalparametric --- evaluate the parametric component of a curve
#  FUNCTION
#    Evaluate the parametric component of a curve, either at all observations
#    or at a single observation.
#  INPUTS
#    curve      an RCurve structure
#    i          the index of the observation that should be evaluated (0=all)
#  OUTPUTS
#    curve      the input curve, with x[i] reevaluated at curve$param.par
#  SYNOPSIS
evalparametric <- function(curve, i = 0)
#  SOURCE
#
{
    if(!curve$haspar) return(curve)
    if(i == 0) ind <- 1:length(curve$x) else ind <- i
    name <- curve$name
    dist <- curve$param.dist
    if(dist == "none") return(curve)
    # extract parameters and values at which to evaluate
    par <- curve$param.par
    x <- curve$x[ind]
    if(name == "hazard"){
        if(dist == "exponential"){
            # exponential components are parametrized by their log-baseline
            lambda <- exp(par)
            y <- rep(lambda, length(x))
            ycum <- x * lambda
        }
        if(dist == "weibull"){
            # weibull components are parametrized by their log-baseline and log-scale
            lambda <- exp(par[1])
            alpha <- exp(par[2])
             y <- alpha * lambda * x^(alpha - 1)
             ycum <- lambda * x^alpha
        }
        if(dist == "lognormal")
            stop("lognormal distribution currently not fully supported")
    }
    if(name == "frailty"){
        ycum <- NULL
        if(dist == "gamma"){
            # gamma components are parametrized by minus their log-shape
            alpha <- exp(-par)
            y <- dgamma(x, shape = alpha, rate = alpha)
        }
        if(dist == "lognormal"){
            # lognormal components are parametrized by their log-variance
            alpha <- exp(par)
            y <- exp(-(log(x) + alpha / 2)^2 / (2 * alpha)) / (x * sqrt(2 * pi * alpha))
        }
    }
    curve$param.y[ind] <- y 
    curve$param.ycum[ind] <- ycum
    if(curve$hasspline) {
         # reweight the curve if it has a spline component
         curve <- weightcurve(curve, i)
    }else{
        curve$y[ind] <- curve$param.y[ind]
        curve$ycum[ind] <- curve$param.ycum[ind]
    }  
    return(curve)
}
#************ evalparametric 

#****f* curveUpdate/evalspline
#  NAME
#    evalspline --- evaluate the spline component of a curve
#  FUNCTION
#    Evaluate the spline component of a curve with a new set of component weights, either
#    at all observations, or at a single index.
#  INPUTS
#    curve      an RCurve structure
#    i          index of the observation that should be evaluated (0=all) 
#    quick      for hazard curve, whether the cumulative basis integrals should be computed
#  OUTPUTS
#    curve      the curve, with x[i] evaluated at the new curve$spline.par parameters
#  SYNOPSIS
evalspline <- function(curve, i = 0, quick = FALSE)
#  SOURCE
#
{
    if(!curve$hasspline) return(curve)
    # extract observations and parameters
    if(i == 0) {
		ind <- 1:length(curve$x)
		curve$spline.y <- NULL
	}else{ind <- i}
    spline.par <- curve$spline.par
    #   normalize the frailty spline parameters
    if(curve$name == "frailty") spline.par <- spline.par - log(sum(exp(spline.par)))
    # Evaluate the spline component
    curve$spline.y[ind] <- drop(curve$spline.basis[ind, ,drop = FALSE]%*%exp(spline.par))
    if(curve$name == "hazard" & !quick)
        # for the hazard, evaluate the spline integrals
        curve$spline.ycum[ind] <- drop(curve$spline.basiscum[ind, ,drop = FALSE]%*%exp(spline.par))
    else
        curve$spline.ycum <- NULL
    if(curve$haspar) {
        # Reweight the curve
        curve <- weightcurve(curve, i)
    }else{
        curve$y[ind] <- curve$spline.y[ind]
        curve$ycum[ind] <- curve$spline.ycum[ind]
    }
    return(curve)
}
#************ evalspline 

#****f* splineUtils/frailtysplinefvar
#  NAME
#    frailtysplinefvar --- compute the spline component frailty variance
#  FUNCTION
#    Compute the variance of the frailty density specified by a curve's spline component.
#  INPUTS
#    frailty    an RCurve structure
#  OUTPUTS
#    fvar       variance of the distribution specified by the frailty curve    
#  SYNOPSIS
frailtysplinefvar <- function(frailty)
#  SOURCE
#
{
    # compute the second moment of the frailty
    moment2 = 1 - cevalEinte(frailty$spline.knots, frailty$spline.ord, 2)
    moment2 <- drop(moment2)%*%drop(exp(frailty$spline.par))
    # normalize
    moment2 <- moment2 / sum(exp(frailty$spline.par))
    # subtract the mean
    fvar <- moment2 - 1
    return(fvar)
}
#************ frailtysplinefvar 

#****f* curveUpdate/updatespline
#  NAME
#    updatespline --- update the curve's spline parameters
#  FUNCTION
#    Update a curve with new spline parameters. This simply consists of copying the parameters
#    into the curve structure and re-evaluating it using evalspline.
#  INPUTS
#    curve      an RCurve structure
#    spline.par a new set of parameters for the spline component
#  OUTPUTS
#    curve      the input curve with the spline component updated to the new parameters.
#  SYNOPSIS
updatespline <- function(curve, spline.par)
#  SOURCE
#
{
    if(!curve$hasspline) return(curve)
    # For the frailty, make sure that the parameter at the fixed index is 0
    if(curve$name == "frailty") spline.par <- spline.par - spline.par[curve$spline.fixedind]
    # copy new parameters into the curve
    curve$spline.par <- spline.par
    # evaluate the curve
    curve <- evalspline(curve)
    return(curve)
}
#************ updatespline 

#****f* curveUpdate/updatecurvex
#  NAME
#    updatecurvex --- change observations of a curve
#  FUNCTION
#    Change the set of "x" values of a curve, that is, the set of points at which the
#    curve is evaluated. This is called by mh.frail, when the frailties are updated.
#    This function should be called if curve$x[i] was changed.
#  INPUTS
#    curve      an RCurve structure
#    i          the index that was updated
#  OUTPUTS
#    curve      the input curve, with the basis evaluated at x[i]
#  SYNOPSIS
updatecurvex <- function(curve, i)
#  SOURCE
#
{
    if(curve$name != "frailty") stop("Only frailty bases can be updated")
    if(curve$hasspline){
        knots <- curve$spline.knots; ord <- curve$spline.ord; x <- curve$x[i]
        # recompute the basis at x[i]
        curve$spline.basis[i, ] <- mysplineDesign(knots, x, ord) / curve$spline.basisint 
        # evaluate the spline component
        curve <- evalspline(curve, i)
    }
    # evaluate the parametric component
    curve <- evalparametric(curve, i)
    return(curve)   
}
#************ updatecurvex 

#****f* curveUpdate/updateparametric
#  NAME
#    updateparametric --- update parametric component parameters
#  FUNCTION
#    Update a curve to use a new set of parameters for the parametric component.
#  INPUTS
#    curve      an RCurve structure
#    param.par  parameters for the parametric component
#  OUTPUTS
#    curve      the input curve, re-evaluated at the new parameters
#  SYNOPSIS
updateparametric <- function(curve, param.par)
#  SOURCE
#
{
    if(!curve$haspar) return(curve)
    # copy the new parameters into the curve and re-evaluate
    curve$param.par <- param.par
    curve <- evalparametric(curve)
    return(curve)
}
#************ updateparametric 

#****f* curveUpdate/weightcurve
#  NAME
#    weightcurve --- re-weight the spline and parametric components of a curve
#  FUNCTION
#    Re-evaluate a curve if the relative weight of the spline and parametric component
#    has changed, or if either of the component values has changed. Can evaluate either
#    the entire curve or only a single index.
#  INPUTS
#    curve      an RCurve structure
#    i          the index at which to evaluate (0=all)
#  OUTPUTS
#    curve      the input curve, with curve$y and curve$ycum updated
#  SYNOPSIS
weightcurve <- function(curve, i = 0)
#  SOURCE
#
{
    if(i == 0) ind <- 1:length(curve$x) else ind <- i
    # reweigh curve$y
    curve$y[ind] <- curve$weight * curve$spline.y[ind] + (1 - curve$weight) *
        curve$param.y[ind]
    # for the hazard, reweigh curve$ycum
    if(curve$name == "hazard")
        curve$ycum[ind] <- curve$weight * curve$spline.ycum[ind] + (1 - curve$weight) *
            curve$param.ycum[ind]
    return(curve)
}
#************ weightcurve 

#****f* curveUpdate/updateregression
#  NAME
#    updateregression --- update regression coefficients
#  FUNCTION
#     Recompute the RRegression structure for a new set of coefficients
#  INPUTS
#     regression        an RRegression structure
#     coef              new set of regression coefficients
#  OUTPUTS
#     regression        RRegression structure with the lp component updated
#  SYNOPSIS
updateregression <- function(regression, coef)
#  SOURCE
#
{
    regression$coefficients <- coef
    regression$lp <- drop(regression$covariates%*%regression$coefficients)
    return(regression)
}
#************ updateregression 

}}}

{{{ #Likelihood
##############################################################
# \section{Likelihood} Likelihoods, gradients and hessians
##############################################################

#****f* makeLikelihood/mklik.coef
#  NAME
#    mklik.coef --- likelihood of regression coefficients
#  FUNCTION
#    Compute the likelihood of regression coefficients, given everything else.
#  INPUTS
#    coef       regression coefficients to be evaluated
#    hazard     RCurve for hazard
#    frailty    RCurve for frailty
#    regression RRegression structure
#  OUTPUTS
#    lik    loglikelihood of coef
#  SYNOPSIS
mklik.coef <- function(coef, hazard, frailty, regression)
#  SOURCE
#
{
    status <- regression$status
    regression <- updateregression(regression, coef)
    lp <- regression$lp
    frailrep <- rep(frailty$x, regression$Ji)
    # point process likelihood
    lik <- status%*%lp
    lik <- lik - sum(frailrep * hazard$ycum * exp(lp))
    # prior penalty
    lik <- lik - sum(coef^2) / (2 * regression$priorvar)
    return(as.numeric(lik))
}
#************ mklik.coef 

#****f* makeLikelihood/mkhess.coef
#  NAME
#    mkhess.coef --- hessian of regression coefficients
#  FUNCTION
#    Compute the hessian of regression coefficients, given everything else.
#  INPUTS
#    coef       regression coefficients to be evaluated
#    hazard     RCurve for hazard
#    frailty    RCurve for frailty
#    regression RRegression structure
#  OUTPUTS
#    hess    hessian of coef
#  SYNOPSIS
mkhess.coef <- function(coef, hazard, frailty, regression)
#  SOURCE
#
{
    status <- regression$status
    regression <- updateregression(regression, coef)
    lp <- regression$lp
    frailrep <- rep(frailty$x, regression$Ji)
    Z <- regression$covariates
    hess <- -t(Z)%*%(rep(frailrep * exp(lp) * hazard$ycum, dim(Z)[2]) * Z)
    hess <- hess - mdiag(rep(1, length(coef))) / regression$priorvar
    return(hess)
}
#************ mkhess.coef 

#****f* makeLikelihood/mklik.frail
#  NAME
#    mklik.frail --- likelihood of frailty for cluster i
#  FUNCTION
#    Compute the loglikeilhood of the frailty of cluster i
#  INPUTS
#    i          index of the frailty. The frailty value is stored in frailty$x[i]
#    hazard     RCurve for hazard
#    frailty    RCurve for frailty
#    regression RRegression structure
#  OUTPUTS
#    lik    loglikelihood of frailty$x[i]
#  SYNOPSIS
mklik.frail <- function(i, hazard, frailty, regression)
#  SOURCE
#
{
    ind <- which(regression$cluster == i)
    Ui <- frailty$x[i]
    status <- regression$status[ind]
    lp <- regression$lp[ind]
    cumhaz <- hazard$ycum[ind]
    lik <- log(frailty$y[i])
    lik <- lik + sum(status * log(Ui))
    lik <- lik - sum(Ui * cumhaz * exp(lp))
    return(lik)   
}
#************ mklik.frail 

#****f* makeLikelihood/mklik.spline.haz
#  NAME
#    mklik.spline.haz --- likelihood of hazard spline parameters
#  FUNCTION
#    Compute the loglikelihood of parameters spline.par for the spline component of
#    the hazard curve.
#  INPUTS
#    spline.par a vector of parameters for each of the spline basis functions
#    hazard     RCurve for hazard
#    frailty    RCurve for frailty
#    regression RRegression structure
#  OUTPUTS
#    lik    loglikelihood of spline.par
#  SYNOPSIS
mklik.spline.haz <- function(spline.par, hazard, frailty, regression)
#  SOURCE
#
{
    if(!hazard$hasspline) return(0)
    if(any(is.na(spline.par))) return(-Inf)
    status <- regression$status
    lp <- regression$lp
    hazard <- updatespline(hazard, spline.par)
    frailrep <- rep(frailty$x, regression$Ji)
    # point process likelihood
    lik <- sum(status * log(hazard$y)) 
    lik <- lik - sum(frailrep * hazard$ycum * exp(lp))
    # smoothness penalty
    lik <- lik - hazard$spline.penaltyfactor * smoothpen(hazard, 0)
    # penalize parameters that are too small
    lik <- lik - sum(ifelse(spline.par< hazard$spline.min, 
        (spline.par - hazard$spline.min)^2, 0))    
    lik <- as.numeric(lik)
    return(lik)
}
#************ mklik.spline.haz 

#****f* makeLikelihood/mkgr.spline.haz
#  NAME
#    mkgr.spline.haz --- gradient of hazard spline parameters
#  FUNCTION
#    Compute the gradient of parameters spline.par for the spline components of a
#    hazard curve. This is only called during initialization, to give gradients for
#    use by optim().
#  INPUTS
#    spline.par a vector of parameters for each of the spline basis functions
#    hazard     RCurve for hazard
#    frailty    RCurve for frailty
#    regression RRegression structure
#  OUTPUTS
#    gr    gradient of the loglikelihood of spline.par
#  SYNOPSIS
mkgr.spline.haz <- function(spline.par, hazard, frailty, regression)
#  SOURCE
#
{
    if(!hazard$hasspline) return(rep(0, length(spline.par)))
    status <- regression$status
    lp <- regression$lp
    hazard <- updatespline(hazard, spline.par)
    frailrep <- rep(frailty$x, regression$Ji)
    Det <- diag(exp(hazard$spline.par))
    gr <- Det%*%t(hazard$spline.basis)%*%(status / hazard$y)
    gr <- gr - Det%*%t(hazard$spline.basiscum)%*%(frailrep * exp(lp))
    gr <- hazard$weight * gr
    gr <- gr - hazard$spline.penaltyfactor * smoothpen(hazard, 1)
    gr <- gr - ifelse(spline.par< hazard$spline.min, 2 * (spline.par - hazard$spline.min), 0)  
    gr <- as.numeric(gr)
    return(gr)
}
#************ mkgr.spline.haz 

#****f* makeLikelihood/mklik.spline.frail.init
#  NAME
#    mklik.spline.frail.init --- likelihood of frailty spline parameters (initialization)
#  FUNCTION
#    This is a variant of mklik.spline.frail for use during initialization. It differs in
#    that spline.par does not contain the fixed index, and hence repairfrailtypar must be
#    called to add it back in. Moreover, it produces a frailty mean close to 1 by
#    penalizing the distance.
#  INPUTS
#    spline.par a vector of parameters for each of the spline basis functions
#    hazard     RCurve for hazard
#    frailty    RCurve for frailty
#    regression RRegression structure
#  OUTPUTS
#    lik    loglikelihood of spline.par
#  SYNOPSIS
mklik.spline.frail.init <- function(spline.par, hazard, frailty, regression)
#  SOURCE
#
{
    if(any(is.na(spline.par))) return(-Inf)
    spline.par <- repairfrailtypar(spline.par, frailty$spline.fixedind)
    frailty <- updatespline(frailty, spline.par)
    M <- frailty$spline.meanpenalty
    # base likeihood and smoothness penalty
    lik <- sum(log(frailty$y)) - frailty$spline.penaltyfactor * smoothpen(frailty, 0)
    # penalty for being far from mean 1
    lik <- lik - M * (frailty$spline.basisexp%*%exp(frailty$spline.par))^2
    # penalty for having too small parameters
    lik <- lik - sum(ifelse(spline.par< frailty$spline.min, 
        (spline.par - frailty$spline.min)^2, 0))    
    # do not allow too large parameters
    if(any(spline.par > 20)) lik<- -Inf #needed for numerics
    lik <- as.numeric(lik)
    return(lik)
}
#************ mklik.spline.frail.init 

#****f* makeLikelihood/mklik.spline.frail
#  NAME
#    mklik.spline.frail --- likelihood of frailty spline parameters
#  FUNCTION
#    Compute loglikelihood of spline.par for a frailty spline curve.
#  INPUTS
#    spline.par a vector of parameters for each of the spline basis functions
#    hazard     RCurve for hazard
#    frailty    RCurve for frailty
#    regression RRegression structure
#  OUTPUTS
#    lik        loglikelihood of spline.par
#  SYNOPSIS
mklik.spline.frail <- function(spline.par, hazard, frailty, regression)
#  SOURCE
#
{
    if(!frailty$hasspline) return(0)
    if(any(is.na(spline.par))) return(-Inf)
    frailty <- updatespline(frailty, spline.par)
    M <- frailty$spline.meanpenalty
    lik <- sum(log(frailty$y))
    lik <- lik - frailty$spline.penaltyfactor * smoothpen(frailty, 0)
    lik <- lik - sum(ifelse(spline.par< frailty$spline.min,
        (spline.par - frailty$spline.min)^2, 0))    
    if(any(spline.par > 20)) lik<- -Inf #needed for numerics
    lik <- as.numeric(lik)
    return(lik)
}
#************ mklik.spline.frail 

#****f* makeLikelihood/mklik.param.haz
#  NAME
#    mklik.param.haz --- likelihood of parametric component parameters for hazard
#  FUNCTION
#    Compute loglikelihood of parametric component parameters for the hazard curve
#  INPUTS
#    param.par  a vector of parameters for each of the parametric component
#    hazard     RCurve for hazard
#    frailty    RCurve for frailty
#    regression RRegression structure
#  OUTPUTS
#    lik        loglikelihood of param.par
#  SYNOPSIS
mklik.param.haz <- function(param.par, hazard, frailty, regression)
#  SOURCE
#
{
    if(!hazard$haspar) return(0)
    # update parametric component
    hazard <- updateparametric(hazard, param.par)
    status <- regression$status
    lp <- regression$lp
    frailrep <- rep(frailty$x, regression$Ji)
    # likelihood computation
    lik <- sum(status * log(hazard$y) - frailrep * hazard$ycum * exp(lp))
    lik <- lik - sum(param.par^2) / (2 * hazard$param.priorvar)
    return(lik)
}
#************ mklik.param.haz 

#****f* makeLikelihood/mklik.param.frail
#  NAME
#    mklik.param.frail --- likelihood of parametric component for frailty
#  FUNCTION
#   Compute loglikelihood of parametric component parameters for the frailty curve
#  INPUTS
#    param.par  a vector of parameters for each of the parametric component
#    hazard     RCurve for hazard
#    frailty    RCurve for frailty
#    regression RRegression structure
#  OUTPUTS
#    lik        loglikelihood of param.par
#  SYNOPSIS
mklik.param.frail <- function(param.par, hazard, frailty, regression)
#  SOURCE
#
{
    if(!frailty$haspar) return(0)
    frailty <- updateparametric(frailty, param.par)
    lik <- sum(log(frailty$y)) - sum(param.par^2) / (2 * frailty$param.priorvar)
    return(lik)
}
#************ mklik.param.frail 

#****f* makeLikelihood/mklik.weight.haz
#  NAME
#    mklik.weight.haz --- likelihood of weight of spline component for hazard
#  FUNCTION
#   Compute loglikelihood of relative weight of spline and parametric components
#   if the hazard curve has both.
#  INPUTS
#    weight     weight of the spline component (0<=weight<=1)
#    hazard     RCurve for hazard
#    frailty    RCurve for frailty
#    regression RRegression structure
#  OUTPUTS
#    lik        loglikelihood of weight
#  SYNOPSIS
mklik.weight.haz <- function(weight, hazard, frailty, regression)
#  SOURCE
#
{
    hazard$weight <- weight
    hazard <- weightcurve(hazard)
    status <- regression$status
    lp <- regression$lp
    frailrep <- rep(frailty$x, regression$Ji)
    # point process likelihood
    lik <- sum(status * log(hazard$y) - frailrep * hazard$ycum * exp(lp))
    # prior on the weight
    lik <- lik + (hazard$weight.hyper[1] - 1) * log(hazard$weight) +
        (hazard$weight.hyper[2] - 1) * log(1 - hazard$weight)
    return(lik) 
}
#************ mklik.weight.haz 

#****f* makeLikelihood/mklik.weight.frail
#  NAME
#    mklik.weight.frail --- likelihood of weight of spline component for frailty
#  FUNCTION
#   Compute loglikelihood of relative weight of spline and parametric components
#   if the frailty curve has both.
#  INPUTS
#    weight     weight of the spline component (0<=weight<=1)
#    hazard     RCurve for hazard
#    frailty    RCurve for frailty
#    regression RRegression structure
#  OUTPUTS
#    lik        loglikelihood of weight
#  SYNOPSIS
mklik.weight.frail <- function(weight, hazard, frailty, regression)
#  SOURCE
#
{
    frailty$weight <- weight
    frailty <- weightcurve(frailty)
    # base likelihood
    lik <- sum(log(frailty$y))
    # prior on the weight
    lik <- lik + (frailty$weight.hyper[1] - 1) * log(frailty$weight) +
        (frailty$weight.hyper[2] - 1) * log(1 - frailty$weight)
    return(lik) 
}
#************ mklik.weight.frail 

}}}

{{{ #Metropolis
##############################################################
# \section{Metropolis} Metropolis - Hastings
##############################################################

#****f* MetropolisHastings/acceptreject
#  NAME
#    acceptreject --- accept of reject a M-H step
#  FUNCTION
#    Handles accept-reject portion of every Metropolis-Hastings step. Given a
#    base likelihood and a candidate loglikelihood, generates a random number and
#    decides whether to accept or reject the step.
#  INPUTS
#    baselik    loglikelihood of the base model
#    candlik    loglikelihood of the candidate model
#    ratio      additional multiplier (e.g. prior ratio)
#  OUTPUTS
#    acc        boolean, whether to accept or reject the jump
#  SYNOPSIS
acceptreject <- function(baselik, candlik, ratio = 1)
#  SOURCE
#
{
    if(is.nan(candlik)) candlik <- -Inf
    # compute acceptance probability
    r <- exp(candlik - baselik) * ratio
    p.accept <- min(1, r)
    if(is.nan(p.accept)) p.accept <- 0
    # generate uniform and accept or reject
    if(runif(1) < p.accept){
        return(TRUE)
    }else{
        return(FALSE)
    }
}
#************ acceptreject 

#****f* MetropolisHastings/mh
#  NAME
#    mh --- prototype Metropolis-Hastings
#  FUNCTION
#    Metropolis-Hastings step for general parameters and likelihood functions,
#    for which the candidates are generated as multivariate normal.
#  INPUTS
#    par        base model parameters
#    fun        likelihood function to use
#    candcov    covariance matrix for candidate generation
#    tun        tuning parameter for candidate generation
#  OUTPUTS
#    par        new parameters after MH step
#    acc        boolean, whether the step was accepted
#  SYNOPSIS
mh <- function(par, fun, candcov, tun, ...)
#  SOURCE
#
{
    # base likelihood
    baselik <- fun(par, ...)
    # generate candidate
    cand <- MYmvrnorm(1, par, candcov * tun)
    # candidate likelihood
    candlik <- fun(cand, ...)
    # accept-reject and update parameters
    acc <- acceptreject(baselik, candlik)
    if(acc) out <- cand else out <- par
    return(list(par = out, acc = acc))    
}
#************ mh 

#****f* MetropolisHastings/mh.frail
#  NAME
#    mh.frail --- MH for frailties
#  FUNCTION
#    Update the frailty estimates by Metropolis-Hastings
#  INPUTS
#    hazard     RCurve for hazard
#    frailty    RCurve for frailty
#    regression RRegression structure
#  OUTPUTS
#    frailty    Rcurve with updated frailty values
#  SYNOPSIS
mh.frail <- function(hazard, frailty, regression)
#  SOURCE
#
{
    acc <- rep(0, regression$m)
    # update each of the frailties separately
    for(i in 1:regression$m){
        u <- frailty$x[i]
        v <- frailty$tun
        # generate candidates from gamma distribution
        cand <- rgamma(1, shape = u^2 / v, scale = v / u)
        # check if the candidate is out of bounds or NaN
        if(is.nan(cand) || (frailty$hasspline && 
          (cand > attr(frailty$spline.knots, "b")[2] | 
            cand < attr(frailty$spline.knots, "b")[1]))) next;
        # choose another fraity to compensate for the change in the mean
        j <- i
        while(j == i) j <- floor(runif(1, 1, length(frailty$x) + 1))
        # adjust the second candidate to make sure the mean remains 1
        candj <- u + frailty$x[j] - cand
        # check that candj is in bounds as well.
        if(is.nan(candj) || (frailty$hasspline && 
          (candj > attr(frailty$spline.knots, "b")[2] | 
            candj < attr(frailty$spline.knots, "b")[1]))) next;

        # base likelihood
        baselik <- mklik.frail(i, hazard, frailty, regression) +
            mklik.frail(j, hazard, frailty, regression)
        temp <- frailty
        temp$x[i] <- cand
        temp$x[j] <- candj
        temp <- updatecurvex(temp, i)
        temp <- updatecurvex(temp, j)
        # candidate likelihoood
        candlik <- mklik.frail(i, hazard, temp, regression) +
            mklik.frail(j, hazard, temp, regression)

        # transition u->cand
        puc <- suppressWarnings(dgamma(cand, shape = u^2 / v, scale = v / u)) 
        # transition cand->u
        pcu <- suppressWarnings(dgamma(u, shape = cand^2 / v, scale = v / cand)) 

        acc[i] <- acceptreject(baselik, candlik, pcu / puc)
        if(acc[i]) frailty <- temp
    }
    frailty$accept <- mean(acc)
    return(frailty)
}
#************ mh.frail 

#****f* MetropolisHastings/mh.frailty.spline
#  NAME
#    mh.frailty.spline --- MH for frailty spline parameters
#  FUNCTION
#    Metropolis-Hastings steps for spline parameters for frailty
#  INPUTS
#    hazard     RCurve for hazard
#    frailty    RCurve for frailty
#    regression RRegression structure
#  OUTPUTS
#    frailty    Rcurve with updated frailty parameters
#  SYNOPSIS
mh.frailty.spline <- function(hazard, frailty, regression)
#  SOURCE
#
{
    if(!frailty$hasspline) return(frailty)
    sumacc <- 0
    nj <- length(frailty$spline.par)
    ord2 <- round(frailty$spline.ord) / 2
    # base likelihood
    baselik <- mklik.spline.frail(frailty$spline.par, hazard, frailty, regression)
    # update each parameter separately
    for(j in 1:nj){
        # get position j and position k which will be used to keep the mean at 1
        cand <- frailty$spline.par
        k <- j
        Ej <- frailty$spline.basisexp[j]
        while(j == k | k < ord2 | k > nj - ord2) k <- floor(runif(1, 1, nj + 1))
        Ek <- frailty$spline.basisexp[k]
        cand[j] <- cand[j] + frailty$spline.tun * rnorm(1, 0, frailty$spline.candsd[j])
        newmean <- Ej * (exp(cand[j]) - exp(frailty$spline.par[j]))
        newinner <- exp(frailty$spline.par[k]) - newmean / Ek
        if(newinner <= 0) next
        candk = log(newinner)
        cand[k] = candk
        # candidate likelihood
        candlik <- mklik.spline.frail(cand, hazard, frailty, regression)
        thisacc <- acceptreject(baselik, candlik, 1)
        if(thisacc){
            baselik <- candlik
            frailty <- updatespline(frailty, cand)
            frailty$spline.fvar <- frailtysplinefvar(frailty)
        }
        sumacc <- sumacc + thisacc
    }
    frailty$spline.accept <- sumacc / nj
    return(frailty)
}
#************ mh.frailty.spline 

#****f* MetropolisHastings/mh.hazard.spline
#  NAME
#    mh.hazard.spline --- MH for hazard spline parameters
#  FUNCTION
#    Metropolis-Hastings steps for spline parameters for hazard
#  INPUTS
#    hazard     RCurve for hazard
#    frailty    RCurve for frailty
#    regression RRegression structure
#  OUTPUTS
#    hazard     Rcurve with updated hazard parameters
#  SYNOPSIS
mh.hazard.spline <- function(hazard, frailty, regression)
#  SOURCE
#
{
    if(!hazard$hasspline) return(hazard)
    sumacc <- 0
    nj <- length(hazard$spline.par)
    cand <- rep(0, nj)
    for(j in 1:nj) cand[j] <- hazard$spline.par[j] + hazard$spline.tun *
        rnorm(1, 0, hazard$spline.candsd[j])
    baselik <- mklik.spline.haz(hazard$spline.par, hazard, frailty, regression)
    for(j in 1:nj){
        thiscand <- hazard$spline.par
        thiscand[j] <- cand[j]
        candlik <- mklik.spline.haz(thiscand, hazard, frailty, regression)
        thisacc <- acceptreject(baselik, candlik, 1)
        if(thisacc){
            baselik <- candlik
            hazard <- updatespline(hazard, thiscand)
        }
        sumacc <- sumacc + thisacc
    }
    hazard$spline.accept <- sumacc / nj
    return(hazard)
}
#************ mh.hazard.spline 

#****f* MetropolisHastings/mh.frailty.param
#  NAME
#    mh.frailty.param --- MH for frailty parametric component parameters
#  FUNCTION
#    Metropolis-Hastings steps for parametric component parameters for frailty
#  INPUTS
#    hazard     RCurve for hazard
#    frailty    RCurve for frailty
#    regression RRegression structure
#  OUTPUTS
#    frailty     Rcurve with updated frailty parameters
#  SYNOPSIS
mh.frailty.param <- function(hazard, frailty, regression)
#  SOURCE
#
{
    if(!frailty$haspar) return(frailty)
    mhout <- mh(frailty$param.par, mklik.param.frail, frailty$param.candcov,
        frailty$param.tun, hazard = hazard, frailty = frailty, regression = regression)
    if(mhout$acc) frailty <- updateparametric(frailty, mhout$par)
    frailty$param.accept <- mhout$acc
    return(frailty)
}
#************ mh.frailty.param 

#****f* MetropolisHastings/mh.hazard.param
#  NAME
#    mh.hazard.param --- MH for hazard parametric component parameters
#  FUNCTION
#    Metropolis-Hastings steps for parametric component parameters for hazard
#  INPUTS
#    hazard     RCurve for hazard
#    frailty    RCurve for frailty
#    regression RRegression structure
#  OUTPUTS
#    hazard     Rcurve with updated hazard parameters
#  SYNOPSIS
mh.hazard.param <- function(hazard, frailty, regression)
#  SOURCE
#
{
    if(!hazard$haspar) return(hazard)
    mhout <- mh(hazard$param.par, mklik.param.haz, hazard$param.candcov, hazard$param.tun,
        hazard = hazard, frailty = frailty, regression = regression)
    if(mhout$acc) hazard <- updateparametric(hazard, mhout$par)
    hazard$param.accept <- mhout$acc
    return(hazard)
}
#************ mh.hazard.param 

#****f* MetropolisHastings/mh.coef
#  NAME
#    mh.coef --- MH for regression coefficients
#  FUNCTION
#    Metropolis-Hastings step for regression coefficients
#  INPUTS
#    hazard     RCurve for hazard
#    frailty    RCurve for frailty
#    regression RRegression structure
#  OUTPUTS
#    regression   RRegression with updated coefficients
#  SYNOPSIS
mh.coef <- function(hazard, frailty, regression)
#  SOURCE
#
{
    mhout <- mh(regression$coefficients, mklik.coef, regression$candcov, regression$tun,
        hazard = hazard, frailty = frailty, regression = regression)
    if(mhout$acc) regression <- updateregression(regression, mhout$par)
    regression$accept <- mhout$acc
    return(regression)
}
#************ mh.coef 

#****f* MetropolisHastings/mh.weight
#  NAME
#    mh.weight --- MH for spline component weight for either hazard or frailty
#  FUNCTION
#    Metropolis-Hastings step for relative weight of spline and parametric components.
#    Works for either hazard or frailty curve, depending on the setting of "which".
#  INPUTS
#    which      string, can be either "hazard" or "frailty"
#    hazard     RCurve for hazard
#    frailty    RCurve for frailty
#    regression RRegression structure
#  OUTPUTS
#    curve      updated RCurve for either hazard or frailty
#  SYNOPSIS
mh.weight <- function(which, hazard, frailty, regression)
#  SOURCE
#
{
    # get the curve and likelihood function for the value of "which"
    which <- match.arg(which, c("hazard", "frailty"))
    if(which == "frailty"){
    	curve <- frailty
        fun <- mklik.weight.frail
    }
    if(which == "hazard"){
    	curve <- hazard
        fun <- mklik.weight.haz
    }
    # generate candidate weight as beta
    if(!curve$haspar | !curve$hasspline) return(curve)
    w <- min(max(curve$weight, .01), .99)
    v <- curve$weight.tun
    alpha <- w * (w * (1 - w) / v - 1)
    beta <- (1 - w) / w * alpha
    cand <- rbeta(1, alpha, beta)
    if(is.nan(cand)){
    	curve$weight.accept <- FALSE;
        return(curve)
    }
    # compute transition ratio
    alphac <- cand * (cand * (1 - cand) / v - 1)
    betac <- (1 - cand) / cand * alphac
    baselik <- fun(w, hazard, frailty, regression)
    candlik <- fun(cand, hazard, frailty, regression)
    puc <- suppressWarnings(dbeta(cand, alpha, beta))
    pcu <- suppressWarnings(dbeta(w, alphac, betac))
    acc <- acceptreject(baselik, candlik, pcu / puc)
    if(acc){
        curve$weight <- cand
        curve <- weightcurve(curve)
    }
    curve$weight.accept <- acc
    return(curve)
}
#************ mh.weight 

#****f* MetropolisHastings/mh.bdm
#  NAME
#    mh.bdm --- RJMCMC for Birth-death-move steps
#  FUNCTION
#    This function handles the reversible-jump MCMC steps for adding, removing, or moving
#    knots in the adaptive knot selection procedure. It can work with either the hazard
#    or frailty curve.
#  INPUTS
#    which      string, can be either "hazard" or "frailty"
#    hazard     RCurve for hazard
#    frailty    RCurve for frailty
#    regression RRegression structure
#  OUTPUTS
#    curve      updated RCurve for either hazard or frailty
#  SYNOPSIS
mh.bdm <- function(which, hazard, frailty, regression)
#  SOURCE
#
{
    # get all the needed components
    if(which == "hazard") curve <- hazard
    if(which == "frailty") curve <- frailty
    if(!curve$hasspline) return(curve);
    knots <- curve$spline.knots
    ord <- curve$spline.ord
    params <- curve$spline.par
    nknots <- curve$spline.nknots
    candknots <- curve$spline.candknots
    ncandknots <- curve$spline.ncandknots
    occ <- attr(candknots, "occupied")
    occind <- which(occ == 1)
    # evaluate the prior probability of k knots, k+1 knots, and k-1 knots
    pk <- nknotsPrior(nknots, curve)
    pkp1 <- if(nknots < curve$spline.maxoccknots) nknotsPrior(nknots + 1, curve) else 0
    pkm1 <- if(nknots > 1) nknotsPrior(nknots - 1, curve) else 0
    # compute probability of birth, death and move steps
    pb <- curve$spline.bdmconst * min(1, pkp1 / pk) # P(birth)
    pd <- curve$spline.bdmconst * min(1, pkm1 / pk) # P(death)
    pm <- max(0, 1 - pb - pd) # P(move)
    u <- runif(1)
    if(curve$name == "hazard") 
        baselik <- mklik.spline.haz(curve$spline.par, hazard, frailty, regression)
    if(curve$name == "frailty") 
        baselik <- mklik.spline.frail(curve$spline.par, hazard, frailty, regression)
    if(u < pd){
        # Death step
        j <- sample(1:nknots, 1) + ord # index of dying knot
        x <- knots[j]
            #cat("Remove knot ", j, " at ", knots[j], "(occ", occ[occind[j - ord]], ")\n")
        knots2 <- knots
        params2 <- params
        knots2 <- knots2[ - j] # remove the knot
        # update parameters symetrically to birth step
        params2 <- params2[ - (j - 1)]
        if(ord > 2) for(j2 in (j - ord + 1):(j - 2)){
             r2 <- (x - knots2[j2]) / (knots2[j2 + ord - 1] - knots2[j2])
             inner <- 1 / r2 *exp(params2[j2]) - (1 - r2) / r2 * exp(params2[j2 - 1])
             if(inner > 0) params2[j2] <- log(inner) else params2[j2] <- curve$spline.min
        }
        attr(knots2, "boundary") <- attr(knots, "boundary")
        attr(knots2, "order") <- attr(knots, "order")
        attr(knots2, "index") <- 1:(nknots - 1) - ord
        # create a temporary curve that is identical to the original, except for
        # the spline.knots and spline.params components
        temp <- curve
        temp$spline.knots <- knots2
        temp$spline.nknots <- nknots - 1
        occ2 <- occ; occ2[occind[j - ord]] <- 0
        attr(temp$spline.candknots, "occupied") <- occ2
        # Compute the spline basis of the temp curve
        temp <- makesplinebasis(temp)
        # Jacobian of the transformation
        J <- exp(params2[j - 1]) / (exp(params[j - 1]) - exp(params[j - 2]))
        if(ord > 2) for(j2 in (j - ord + 2):(j - 2)) {
            r2 <- (x - knots2[j2]) / (knots2[j2 + ord - 1] - knots2[j2])
            J <- J * exp(params2[j2]) / (r2 * exp(params[j2]))
        }
        # candidate likelihood
        if(curve$name == "hazard") {
            temp <- updatespline(temp, params2)
            candlik <- mklik.spline.haz(temp$spline.par, temp, frailty, regression)
        }
        if(curve$name == "frailty"){
            # for the frailty, adjust the mean to make sure it is 1
            newmean <- exp(params2)%*%temp$spline.basisexp -
                exp(params2[j - 2]) * temp$spline.basisexp[j - 2]
            newinner <- -newmean / temp$spline.basisexp[j - 2]
            if(newinner < 0){
                candlik<- -Inf
            }else{
                params2[j - 2] <- log(newinner)
                temp <- updatespline(temp, params2)
                candlik <- mklik.spline.frail(temp$spline.par, hazard, temp, regression)
            }
        }
        # compute acceptance ratio, and accept-reject
        ratio <- sqrt(2 * pi * curve$spline.priorvar) * abs(J)
        acc <- acceptreject(baselik, candlik, ratio)
            #cat("Lik: ", baselik, candlik, J, ratio, acc, "\n")
        if(any(knots2[(ord + 1):(length(knots2) - ord)] != candknots[occ2 == 1])) browser()
        if(acc) return(temp) else return(curve)
    }

    if(u > pd & u < pd + pb){
        # Birth of a knot
        params <- curve$spline.par
        # choose an unoccupied candidate knot location at random
        birthind <- occind[1]
        while(occ[birthind] > 0) birthind <- floor(runif(1, 1, ncandknots + 1)) + ord
        # get the position and interval of candidate knot
        x <- candknots[birthind]
        i <- findInterval(x, knots)
            #cat("Birth at ", x, " after knot ", i, "\n");
        # Compute new knots and parameters
        knots2 <- c(knots[1:i], x, knots[(i + 1):length(knots)])
        attr(knots2, "boundary") <- attr(knots, "boundary")
        attr(knots2, "order") <- attr(knots, "order")
        attr(knots2, "index") <- 1:(nknots + 1) - ord
        params2 <- c(params[1:(i - ord + 1)], 0, params[(i - ord + 2):length(params)])
        # Even though it is possible to add a knot non-destructively, we need to generate
        # another random number to maintain dimension matching and detailed balance
        for(i2 in (i - ord + 2):i){
            r2 <- (x - knots[i2]) / (knots[i2 + ord - 1] - knots[i2])
            if(i2 == i) r2 <- runif(1)
            params2[i2] <- log(r2 * exp(params[i2]) + (1 - r2) * exp(params[i2 - 1]))
        }
        # Create a temp curve with new knots and params
        temp <- curve
        temp$spline.knots <- knots2
        temp$spline.nknots <- nknots + 1
        occ2 <- occ; occ2[birthind] <- 1
        attr(temp$spline.candknots, "occupied") <- occ2
        temp <- makesplinebasis(temp)
        # Jacobian of the transformation
        J <- (exp(params[i]) - exp(params[i - 1])) / exp(params2[i])
        if(ord > 2) for(i2 in (i - ord + 2):(i - 1)) {
            r2 <- (x - knots[i2]) / (knots[i2 + ord - 1] - knots[i2])
            J <- J * r2 * exp(params[i2]) / exp(params2[i2])
        }
        if(curve$name == "hazard") {
            temp <- updatespline(temp, params2)
            candlik <- mklik.spline.haz(temp$spline.par, temp, frailty, regression)
        }
        if(curve$name == "frailty"){
            # adjust frailty mean
            newmean <- exp(params2)%*%temp$spline.basisexp - exp(params2[i]) *
                temp$spline.basisexp[i]
            newinner <- -newmean / temp$spline.basisexp[i]
            if(newinner < 0){
                candlik<- -Inf
            }else{
                params2[i] <- log(newinner)
                temp <- updatespline(temp, params2)
                candlik <- mklik.spline.frail(temp$spline.par, hazard, temp, regression)
            }
        }
        # acceptance ratio
        ratio <- 1 / sqrt(2 * pi * curve$spline.priorvar) * abs(J)
        acc <- acceptreject(baselik, candlik, ratio)
            #cat("Lik: ", baselik, candlik, J, ratio, acc, "\n")
        if(any(knots2[(ord + 1):(length(knots2) - ord)] != candknots[occ2 == 1])) browser()
        if(acc) return(temp) else return(curve)

    }

    if(u > pd + pb) {
        # Move a knot
        # choose a random knot to move
        moveind <- floor(runif(1, 1, nknots + 1))
        # compute the range of candidate knot indices to which the knot can be moved
        leftknotind <- if(moveind == 1) ord + 1 else occind[moveind - 1] + 1
        rightknotind <- if(moveind == nknots) ncandknots + ord else occind[moveind + 1] - 1
        newknotind <- floor(runif(1, leftknotind, rightknotind + 1))
        oldknotind <- occind[moveind]
        #cat("Moveind: ", moveind, "\n");
        #cat("Old / New: ", candknots[oldknotind], candknots[newknotind], "\n");
        if(newknotind == oldknotind) return(curve)
        # create new knots and parameters
        knots2 <- knots
        params2 <- params
        knots2[moveind + ord] <- candknots[newknotind]
        # update occupied knot indices
        occ2 <- occ
        occ2[newknotind] <- 1
        occ2[occind[moveind]] <- 0
        if(sum(occ == 1) < nknots) browser()
        if(any(diff(knots2) < 0)) browser()
        if(any(knots2[(ord + 1):(length(knots2) - ord)] != candknots[occ2 == 1])) browser()
        temp <- curve
        temp$spline.knots <- knots2
        attr(temp$spline.candknots, "occupied") <- occ2
        # update the spline basis
        temp <- makesplinebasis(temp)
        if(curve$name == "hazard") {
            candlik <- mklik.spline.haz(params2, temp, frailty, regression)
        }
        if(curve$name == "frailty"){
            # adjust frailty mean
            newmean <- exp(params2)%*%temp$spline.basisexp - exp(params2[moveind + ord]) *
                temp$spline.basisexp[moveind + ord]
            newinner <- -newmean / temp$spline.basisexp[moveind + ord]
            if(newinner < 0){
                candlik<- -Inf
            }else{
                params2[moveind + ord] <- log(newinner)
                temp <- updatespline(temp, params2)
                candlik <- mklik.spline.frail(params2, hazard, temp, regression)
            }
        }
        acc <- acceptreject(baselik, candlik, 1)
            #cat("Lik: ", baselik, candlik, acc, "\n")
        if(acc) return(temp) else return(curve)
    }

}
#************ mh.bdm 

#****f* MetropolisHastings/updatepostvar.curve
#  NAME
#    updatepostvar.curve --- update the prior variance for a curve
#  FUNCTION
#    Updates the prior variance for the spline parameters of either the hazard or
#    frailty curve, with an inverse-gamma prior and given hyperparameters.
#  INPUTS
#    curve      an RCurve structure
#  OUTPUTS
#    curve      the updated curve
#  SYNOPSIS
updatepostvar.curve <- function(curve)
#  SOURCE
#
{
    if(curve$hasspline) curve$spline.priorvar <- rinvgamma(1, length(curve$spline.par) /
        2 + curve$spline.hyper[1], scale = curve$spline.penaltyfactor * smoothpen(curve) *
        curve$spline.priorvar + curve$spline.hyper[2])
    if(curve$haspar) curve$param.priorvar <- rinvgamma(1, length(curve$param.par) /
        2 + curve$param.hyper[1], scale = sum(curve$param.par^2) / 2 + curve$param.hyper[2])
    return(curve)
}
#************ updatepostvar.curve 

#****f* MetropolisHastings/updatepostvar.coef
#  NAME
#    updatepostvar.coef --- update prior variance for regression coefficients
#  FUNCTION
#    Updates prior variance of regression coefficients using inverse-gamma prior and
#    given hyperparameters.
#  INPUTS
#    regression     an RRegression structure
#  OUTPUTS
#    regression    updated regression
#  SYNOPSIS
updatepostvar.coef <- function(regression)
#  SOURCE
#
{
    regression$priorvar <- rinvgamma(1, length(regression$coefficients) / 2 +
        regression$hyper[1], scale = sum(regression$coefficients^2) / 2 + regression$hyper[2])
    return(regression)
}
#************ updatepostvar.coef 

}}}

{{{ #C - wrappers
    
#****f* CWrappers/csplinedesign
#  NAME
#    csplinedesign --- wrapper for csplinedesign
#  FUNCTION
#    Computes a design matrix for the spline basis, by evaluating the set
#    of order ord B-splines on the set of knots knots, at values x.
#    Wraps the C implenentation of csplinedesign, which works identically to 
#    splineDesign from the splines package,  but faster.
#  INPUTS
#    knots      vector of knots
#    x          points at which the basis should be evaluated.
#    ord        order of the spline
#  OUTPUTS
#    des        matrix, whose (i,j) entry is spline j evaluated at x[i]
#  SYNOPSIS
csplinedesign <- function(knots, x, ord)
#  SOURCE
#
{
    K <- length(knots) - 2 * ord
    design <- matrix(0, length(x), K + ord)
    out <- .C("csplinedesign",
            des = as.double(design),
            x = as.double(x),
            nx = as.integer(length(x)),
            knots = as.double(knots),
            ord = as.integer(ord),
            K = as.integer(K)
        )
    des <- matrix(out$des, length(x), K + ord)
    return(des)
}
#************ csplinedesign 

#****f* CWrappers/cevalEinte
#  NAME
#    cevalEinte --- wrapper for the C implementation of cevalEinte
#  FUNCTION
#    see evalEinte
#  SYNOPSIS
cevalEinte <- function(knots, ord, N = 1)
#  SOURCE
#
{
    K <- length(knots) - 2 * ord
    einte <- rep(0, K + ord);
    out <- .C("cevalEinte",
            einte = as.double(einte),
            knots = as.double(knots),
            ord = as.integer(ord),
            K = as.integer(K),
            N = as.integer(N)
           )
    einte <- out$einte
    return(einte)
}

#************ cevalEinte 
#****f* CWrappers/cevalBinte
#  NAME
#    cevalBinte --- wrapper for the C implementation of evalBinte
#  FUNCTION
#    see evalBinte
#  SYNOPSIS
cevalBinte <- function(knots, ord)
#  SOURCE
#
{
    K <- length(knots) - 2 * ord;
    binte <- rep(0, K + ord);
    out <- .C("cevalBinte",
            binte = as.double(binte),
            knots = as.double(knots),
            ord = as.integer(ord),
            K = as.integer(K)
          )
    binte <- out$binte
    return(binte)
}
#************ cevalBinte 

#****f* CWrappers/cevalCinte
#  NAME
#    cevalCinte --- wrapper for the C implementation of cevalCinte
#  FUNCTION
#    see evalCinte
#  SYNOPSIS
cevalCinte <- function(knots, ord, obs, Binte)
#  SOURCE
#
{
    K <- length(knots) - 2 * ord;
    cinte <- matrix(0, length(obs), length(Binte))
    out <- .C("cevalCinte",
            cinte = as.double(cinte),
            x = as.double(obs),
            nx = as.integer(length(obs)),
            knots = as.double(knots),
            ord = as.integer(ord),
            K = as.integer(K),
            binte = as.double(Binte)
        )
    cinte <- matrix(out$cinte, length(obs), K + ord)
    return(cinte)
}
#************ cevalCinte 

#****f* CWrappers/cmakePenalty.2deriv
#  NAME
#    cmakePenalty.2deriv --- wrapper for cMakePenalty2diff
#  FUNCTION
#    see makePenalty.2deriv
#  SYNOPSIS
cmakePenalty.2deriv <- function(ord, knots){
#  SOURCE
#
    K <- length(knots) - 2 * ord;
    P <- matrix(0, K + ord, K + ord)
    out <- .C("cMakePenalty2diff",
        P = as.double(P),
        knots = as.double(knots),
        ord = as.integer(ord),
        K = as.integer(K)
    )
    P <- matrix(out$P, K + ord, K + ord)
    return(P)
}
#************ cmakePenalty.2deriv 

#****f* ZZdebug/rmklik.spline.haz
#  NAME
#    rmklik.spline.haz --- R re-implementation of the likelihood function in C
#  FUNCTION 
#    For debugging only, works like cInitLikHazSpline.
#  SYNOPSIS
rmklik.spline.haz <- function(spline.par, status, lp, frailrep, hazParY, 
    hazParYcum, weight, B, C, P, penaltyType, sigma2)
#  SOURCE
#
{
    hazSplineY <- B%*%exp(spline.par)
    hazY <- weight * hazSplineY + (1 - weight) * hazParY
    hazSplineYcum <- C%*%exp(spline.par)
    hazYcum <- weight * hazSplineYcum + (1 - weight) * hazParYcum
    lik <- sum(status * log(hazY) - frailrep * hazYcum * exp(lp))
    lik <- as.numeric(lik)
    if(penaltyType == 1) lik <- lik - t((spline.par))%*%P%*%(spline.par) / (2 * sigma2)
    if(penaltyType == 2) lik <- lik - t(exp(spline.par))%*%P%*%exp(spline.par) / (2 * sigma2)
    return(lik)
}
#************ rmklik.spline.haz 

#****f* CWrappers/cmklik.spline.haz
#  NAME
#    cmklik.spline.haz --- spline hazard likelihood in C wrapper
#  FUNCTION
#    Wrapper for cInitLikHazSpline, likelihood function for use during initialization.
#  INPUTS
#    par        vector of spline parameters whose likelihood should be computed
#    status     vector of event indicators
#    lp         vector of linear predictors, beta%*%Z
#    frailrep   vector of frailties of same length as lp, repeated if necessary
#    hazParY    parametric hazard evaluated at each of the event times
#    hazParYcum parametric cumulative hazards
#    weight     relative weight of parametric and spline component
#    B          spline basis produced by csplinedesign
#    C          cumulative spline basis produced by cevalCinte
#    P          penalty matrix
#    penaltyType    integer, see typePenalty
#    sigma2     prior variance of spline parameters
#    min        minimum allowed value of spline parameters
#  OUTPUTS
#    lik        loglikelihood of par
#  SYNOPSIS
cmklik.spline.haz <- function(par, status, lp, frailrep, hazParY, hazParYcum, weight, 
        B, C, P, penaltyType, sigma2, min)
#  SOURCE
#
{
    lik <- as.double(rep(0, 1))
    out <- .C("cInitLikHazSpline",
            lik = lik, par = par, status = status, lp = lp, frailrep = frailrep,
            hazParY = hazParY, hazParYcum = hazParYcum, weight = weight, B = B, C = C, P = P,
            penaltyType = penaltyType, sigma2 = sigma2, ny = as.integer(length(lp)),
            nj = as.integer(length(par)), DUP = FALSE)
    lik <- out$lik
    lik <- lik - sum(ifelse(par< min, (par - min)^2, 0))    
    return(lik)
}
#************ cmklik.spline.haz 

#****f* ZZdebug/rmkgr.spline.haz
#  NAME
#    rmkgr.spline.haz --- R reimplementation of cInitGrHazSpline
#  FUNCTION
#    For debugging only, works like cInitGrHazSpline
#  SYNOPSIS
rmkgr.spline.haz <- function(spline.par, status, lp, frailrep, hazParY, hazParYcum,
    weight, B, C, P, penaltyType, sigma2)
#  SOURCE
#
{
    B <- matrix(B, length(lp), length(spline.par))
    C <- matrix(C, length(lp), length(spline.par))
    hazSplineY <- B%*%exp(spline.par)
    hazY <- weight * hazSplineY + (1 - weight) * hazParY
    hazSplineYcum <- C%*%exp(spline.par)
    hazYcum <- weight * hazSplineYcum + (1 - weight) * hazParYcum
    gr <- rep(0, length(spline.par))
    status = status / hazY
    lp = exp(lp) * frailrep
    gr <- gr + t(B)%*%status
    gr <- gr - t(C)%*%lp
    gr <- gr * exp(spline.par)
    return(gr)
}
#************ rmkgr.spline.haz 

#****f* CWrappers/cmkgr.spline.haz
#  NAME
#    cmkgr.spline.haz --- wrapper for cInitGrHazSpline
#  FUNCTION
#    Wrapper for cInitGrHazSpline, which computes hazard loglikelihood gradient for
#    initialization only.
#  INPUTS
#    see cmklik.spline.haz for inputs and outputs
#  SYNOPSIS
cmkgr.spline.haz <- function(par, status, lp, frailrep, hazParY, hazParYcum, weight,
    B, C, P, penaltyType, sigma2, min)
#  SOURCE
#
{
    gr <- as.double(rep(0, length(par)))
    out <- .C("cInitGrHazSpline",
            gr = gr, par = par, status = status, lp = lp, frailrep = frailrep,
            hazParY = hazParY, hazParYcum = hazParYcum, weight = weight, B = B, C = C, P = P,
            penaltyType = penaltyType, sigma2 = sigma2, ny = as.integer(length(lp)),
            nj = as.integer(length(par)), DUP = FALSE)
    gr <- out$gr
    gr <- gr - ifelse(par< min, 2 * (par - min), 0)    
    return(gr)
}
#************ cmkgr.spline.haz 

#****f* CWrappers/cmklik.spline.frail
#  NAME
#    cmklik.spline.frail --- wrapper for cInitLikFrailSpline
#  FUNCTION
#    Wraps cInitLikFrailSpline, which computes the frailty spline parameter likelihood
#    during initialization.
#  INPUTS
#    par        vector of spline parameters whose likelihood should be computed
#    fixedind   index of the parameter held fixed for identifiability
#    frailParY  parametric frailty density evaluated at each of the frailties
#    weight     relative weight of parametric and spline component
#    B          spline basis produced by csplinedesign
#    E          vector of 1-means produced by cevalEinte
#    M          large number used to penalize the frailty mean
#    P          penalty matrix
#    penaltyType    integer, see typePenalty
#    sigma2     prior variance of spline parameters
#    min        minimum allowed value of spline parameters
#  OUTPUTS
#  SYNOPSIS
cmklik.spline.frail <- function(par, fixedind, frailParY, weight, B, E, M, P,
            penaltyType, sigma2, min)
#  SOURCE
#
{
    par <- as.double(repairfrailtypar(par, fixedind))
    lik <- as.double(0);
    out <- .C("cInitLikFrailSpline", lik = lik, par = par, frailParY = frailParY,
        weight = weight, B = B, E = E, M = M, P = P, penaltyType = penaltyType, sigma2 = sigma2,
        ny = as.integer(length(frailParY)), nj = as.integer(length(par)), DUP = FALSE)
    lik <- out$lik
    lik <- lik - sum(ifelse(par< min, (par - min)^2, 0))    
    return(lik)
}
#************ cmklik.spline.frail 
}}}

{{{ #S3Methods
##############################################################
# \section{S3Methods} S3 Methods for fitting, printing, summary, etc
##############################################################

################# Methods for fitting

#****f* S3Methods/splinesurv
#  NAME
#    splinesurv --- method dispatch for splinesurv
#  FUNCTION
#    This is the main user-callable function that handles method dispatch
#    to splinesurv.agdata, splinesurv.formula and splinesurv.data.frame.
#    See the package documentation for usage instructions.
#  SYNOPSIS
splinesurv <- function(x, ...)
#  SOURCE
#
{
    UseMethod("splinesurv")
}
#************ splinesurv 

#****f* S3Methods/splinesurv.data.frame
#  NAME
#    splinesurv.data.frame --- splinesurv method for data frames
#  FUNCTION
#    This function takes in a data frame in the format required by splinesurv.agdata
#    except for the sorting order. It simply makes sure that the data frame has the
#    required format and passes it on.
#  INPUTS
#    x  a data frame with columns i, j, time, delta, followed by covariates
#  SYNOPSIS
splinesurv.data.frame <- function(x, ...)
#  SOURCE
#
{
    call <- match.call()
    if(dim(x)[2] > 4 && all(colnames(x)[1:4] == c("i", "j", "time", "delta"))) {
        x <- x[order(x$i), ]
        out <- splinesurv.agdata(x, ...)
        out$call <- call
        return(out)
    } else {
        stop("input data frame needs to have colnames i, j, time, delta")
    }
}
#************ splinesurv.data.frame 

#****f* S3Methods/splinesurv.formula
#  NAME
#    splinesurv.formula --- formula interface for splinesurv
#  FUNCTION
#    This is the main user-facing interface function. It takes a formula and additional
#    parameters, conducts basic input checking and builds a data frame in the format
#    required by splinesurv.agdata. After fitting is done, it constructs the splinesurv
#    output object.
#  INPUTS
#    See package documentation
#  OUTPUTS
#    See package documentation for splinesurv object
#  SYNOPSIS
splinesurv.formula <- function(formula, data = parent.frame(), ...)
#  SOURCE
#
{
    # in part based on coxph function
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, sys.parent()))) m$data <- as.data.frame(data)
    m$...<-NULL
    m[[1]] <- as.name("model.frame")
    special <- "cluster"
    Terms <- if (missing(data)) terms(formula, special) else 
        terms(formula, special, data = data)    
    m$formula <- Terms
    m <- eval(m, sys.parent())
    n <- nrow(m)
    # Check response
    resp <- model.extract(m, "response")
    if (!is.Surv(resp)) stop("model response must be a Surv object")
    if(attr(resp, "type") != "right") stop("right - censored survival data only")
    time <- resp[, "time"]
    delta <- resp[, "status"]
    clusterind <- attr(Terms, "specials")$cluster
    clusterind0 <- NULL
    dropx <- NULL
    # handle cluster special terms
    clusternames <- NULL
    if(length(clusterind) > 0){
        cluster <- m[, clusterind]
        for(j in 1:dim(data)[2]) if(isTRUE(all.equal(data[, j], cluster))) clusterind0 <- j
        if(is.factor(cluster)) clusternames <- levels(cluster)
        if(is.numeric(cluster)) clusternames <- as.character(unique(cluster))
        i <- as.numeric(as.factor(cluster))
        tempc <- untangle.specials(Terms, "cluster", 1:10)
        ord <- attr(Terms, "order")[tempc$terms]
        if (any(ord > 1)) stop("Cluster can not be used in an interaction")
        dropx <- c(dropx, tempc$terms)
    }else{
        i <- rep(1, n)
    }
    Ji <- drop(table(i))
    j <- unlist(as.vector(sapply(Ji, function(x) 1:x)))
    newTerms <- if(length(dropx))  Terms[ - dropx] else Terms
    # construct data frame for splinesurv.agdata
    X <- model.matrix(newTerms, m)
    X <- X[, -1, drop = FALSE]
    agdata <- as.data.frame(cbind(i, j, time, delta, X))
    agdata[, -2] <- agdata[order(agdata$i), -2]
    class(agdata) <- c("agdata", "data.frame")
    fit <- splinesurv.agdata(agdata, ...)
    gcout <- gc()
    # clean up output object
    fit$call <- call
    colnames(fit$history$frailty) <- clusternames
    if(!is.null(fit$posterior.mean)) names(fit$posterior.mean$frailty) <- clusternames
    fit$terms <- newTerms
    attr(fit$terms, "special")$cluster <- clusterind0
    return(fit)
}
#************ splinesurv.formula 

################# Methods for printing

#****f* S3Methods/print.splinesurv
#  NAME
#    print.splinesurv --- print a splinesurv object
#  FUNCTION
#   Prints the posterior mean of splinesurv fit coefficients.
#  SYNOPSIS
print.splinesurv <- function(x, ...)
#  SOURCE
#
{
    coef <- as.matrix(x$posterior.mean$coef, ncol = 1)
    colnames(coef) <- "coef"
    print(coef, ...)
    invisible(x)
}
#************ print.splinesurv 

#****f* S3Methods/summary.splinesurv
#  NAME
#    summary.splinesurv --- creates an object of class summary.splinesurv
#  FUNCTION
#    Summarizes a splinesurv fit. See package documentation for details.
#  SYNOPSIS
summary.splinesurv <- function(object, quantiles = c(.025, .975), ...)
#  SOURCE
#
{
    x <- object
    out <- NULL
    #  Extract components of the fit that should be included in the summary
    out$call <- x$call
    out$coef <- as.matrix(x$posterior.mean$coefficients, ncol = 1)
    colnames(out$coef) <- "mean"
    out$iter <- x$control$iter
    out$burnin <- x$control$burnin
    out$hazard <- x$hazard
    out$frailty <- x$frailty
    # compute frailty variance using two estimators:
    # as the mean of the variances of the posterior densities
    fvar <- post.fvar(x, quantiles)
    # or as the variance of the posterior frailty samples
    fvar2 <- apply(x$history$frailty[(out$burnin + 1):out$iter, ], 1, var)
    out$frailty$spline.fvar <- fvar$mean["spline.fvar", ]
    out$frailty$param.fvar <- fvar$mean["param.fvar", ]
    out$frailty$fvar <- fvar$mean["fvar", ]
    out$frailty$fvar2 <- mean(fvar2)
    out$posterior.mean <- x$posterior.mean
    out$quantiles.coef <- NULL
    if(out$iter < out$burnin) out$burnin <- 0
    # Compute posterior quantiles of regression coefficients
    if(length(quantiles)){
        goodind <- (out$burnin + 1):(out$iter)
        goodcoef <- x$history$coefficients[goodind, ,drop = FALSE]
        for(q in quantiles){
            out$quantiles.coef <- cbind(out$quantiles.coef, 
                apply(goodcoef, 2, function(x) quantile(x, q)))
        }
        # For presentation, make sure the colnames and rownames are nice
        colnames(out$quantiles.coef) <- paste(quantiles * 100, "%", sep = "")
        rownames(out$quantiles.coef) <- rownames(out$coef)
        out$quantiles.fvar <- fvar$quantiles
        out$quantiles.fvar2 <- quantile(fvar2, quantiles)
    }
    # save the dots for printing parameters
    out$dots <- as.list(substitute(list(...)))[ - 1]
    class(out) <- "summary.splinesurv"   
    return(out)
}
#************ summary.splinesurv 

#****f* S3Methods/post.fvar
#  NAME
#    post.fvar --- compute posterior frailty variance
#  FUNCTION
#    Compute the posterior mean frailty variance and its quantiles by weighing posterior
#    spline and parametric component frailty variances together. The former is computed
#    at each iteration as the variance of the density defined by the current set of knots
#    and parameters for the frailty spline.
#  INPUTS
#    x      an object of class splinesurv
#    quantiles  a list of quantiles at which to compute
#  OUTPUTS
#    means  a vector containing the overall frailty variance, the variance of the spline
#           component, and the variance of the parametric component
#    quantiles  a matrix containing the same for each given quantile
#  SYNOPSIS
post.fvar <- function(x, quantiles = c(.025, .975))
#  SOURCE
#
{
    goodind <- (x$control$burnin + 1):(x$control$iter)
    # spline component
    if(hasspline(x$frailty)) spline.fvars <- x$history$frailty.spline.fvar[goodind] else 
        spline.fvars <- NULL
    # parametric component
    if(haspar(x$frailty)) 
        param.fvars <- exp(x$history$frailty.param.par[goodind, ,drop = FALSE]) else 
            param.fvars <- NULL
    # weigh the two components
    if(haspar(x$frailty) & hasspline(x$frailty)){
        weights <- x$history$frailty.weight[goodind, ,drop = F]
        fvars <- weights * spline.fvars + (1 - weights) * param.fvars
    }else{
        if(haspar(x$frailty)) fvars <- param.fvars else fvars <- spline.fvars
    }
    # compute means
    fvar.means <- suppressWarnings(rbind(mean(fvars), mean(spline.fvars), mean(param.fvars)))
    # compute quantiles
    fvar.quantiles <- suppressWarnings(rbind(quantile(fvars, quantiles),
        quantile(spline.fvars, quantiles), quantile(param.fvars, quantiles)))
    rownames(fvar.means) <- rownames(fvar.quantiles) <- c("fvar", "spline.fvar", "param.fvar")
    return(list(mean = fvar.means, quantiles = fvar.quantiles))
}
#************ post.fvar 

#****f* S3Methods/print.summary.splinesurv
#  NAME
#    print.summary.splinesurv --- print summary.splinesurv object
#  FUNCTION
#    Pretty-print the output of summary.splinesurv
#  INPUTS
#    x      an object of class summary.splinesurv
#    ...    additional parameters passed to formatC
#  SYNOPSIS
print.summary.splinesurv <- function(x, ...)   
#  SOURCE
#
{ 
    cat("Call: \n")
    print(x$call)
    # print basic info
    cat("\nIterations: ", x$iter, " (", x$burnin, " discarded as burn - in)\n", sep = "")
    printpars <- paste(names(x$dots), unlist(x$dots), sep = " = ", collapse = ", ")
    if(nchar(printpars)) printpars <- paste(", ", printpars)
    # print regression coefficients and quantiles
    cat("\nRegression parameter posterior:\n")
    eval(parse(text = paste("print(cbind(x$coef, x$quantiles.coef)", printpars, ")")))
    # print frailty variances
    cat("\nFrailty variance:\n")
    fvarout <- rbind(cbind(x$frailty$fvar, x$quantiles.fvar[1, ,drop = F]),
        c(x$frailty$fvar2, x$quantiles.fvar2))
    rownames(fvarout) <- c("fvar", "fvar2"); colnames(fvarout)[1] <- "mean"
    eval(parse(text = paste("print(fvarout", printpars, ")")))
    # print summaries of each of the curves
    cat("\nBaseline hazard:")
    printcurvesummary(x$hazard, x$posterior.mean$hazard.weight,
        x$posterior.mean$hazard.param.par)
    cat("\nFrailty density:")
    printcurvesummary(x$frailty, x$posterior.mean$frailty.weight,
        x$posterior.mean$frailty.param.par)
    invisible(x)
}
#************ print.summary.splinesurv 

#****f* S3Methods/printcurvesummary
#  NAME
#    printcurvesummary --- summarize a curve for summary.splinesurv
#  FUNCTION
#    Routine called by summary.splinesurv to print a nice summary of a curve structure
#    including mean number of knots, boundaries, parametric/spline components and weights, etc.
#  INPUTS
#    curve  an RCurve structure, possibly with portions missing
#    w      weight of splines component
#    param.par  parametric component parameters to use
#  OUTPUTS
#    none, screen output only
#  SYNOPSIS
printcurvesummary <- function(curve, w = NULL, param.par = NULL)
#  SOURCE
#
{
    haspar <- curve$type == "both" || curve$type == "parametric"
    hasspline <- curve$type == "both" || curve$type == "spline"
    # print spline component
    if(hasspline) {
        cat("\n\tSpline")
        if(haspar) cat(" ( weight =", format(w, digits = 3), "):") else cat(":")
        cat("\n\t\tOrder:", curve$spline.ord)
        cat("\n\t\tMean Interior knots:", format(curve$spline.nknots, digits = 2))
        if(curve$name == "frailty") cat("\n\t\tVariance: ", format(curve$spline.fvar,
            digits = 3))
    }
    # print parametric component
    if(haspar){
        cat("\n\tParametric")
        if(hasspline) cat(" ( weight =", format(1 - w, digits = 3), "):") else cat(":")
        dist <- curve$param.dist
        cat("\n\t\tDistribution:", curve$param.dist)
        if(curve$name == "hazard" & dist == "exponential"){
            cat("\n\t\tRate:", format(exp(param.par), digits = 3))
        }
        if(curve$name == "hazard" & dist == "weibull"){
            cat("\n\t\tRate:", format(exp(param.par[1]), digits = 3))
            cat("\n\t\tScale:", format(exp(param.par[2]), digits = 3)) 
        }
        if(curve$name == "frailty" & (dist == "gamma" | dist == "lognormal")){
            cat("\n\t\tVariance:", format(exp(param.par), digits = 3))
        }
        
    }
}
#************ printcurvesummary 


################# Methods for plotting


#****f* S3Methods/plot.splinesurv
#  NAME
#    plot.splinesurv --- plot method for splinesurv objects
#  FUNCTION
#    Plots either a hazard, frailty, survival curve estimate, or coefficient density.
#    There are a large number of options and settings, consult the package documentation
#    for details on the options.
#
#    This function mainly works by calling the predict function to construct the
#    curve fits. Most of the rest is just handling inputs and keeping track of the
#    various parameters needed by different plotting routines.
#  SYNOPSIS
plot.splinesurv <- function(x, which = c("hazard", "survival", "frailty", "coef", "all"),
    newdata = NULL, iter = NULL, fn = mean, marginal = c("none", "mc", "numerical"),
    plotknots = TRUE, npoints = 100, npost = 100,
    alpha=.05, legend = NULL, lty = 1, col = 2, lwd = 2, lty.knots = 1, col.knots = 8,
    lwd.knots = 1, xlab = NULL, ylab = NULL, main = NULL, xlim = NULL, ylim = NULL,
    tk = FALSE, ...)
#  SOURCE
#
{
    
    if(tk) {
        # call the routine that uses tcltk via the tkrplot to plot an interactive
        # curve viewer
        splinesurvtkplot(x, newdata = newdata, fn = fn, marginal = marginal, 
            plotknots = plotknots, npoints = npoints,
            legend = legend, lty = lty, col = col, lwd = lwd, lty.knots = lty.knots,
            col.knots = col.knots, lwd.knots = lwd.knots, xlab = xlab, ylab = ylab,
            main = main, xlim = xlim, ylim = ylim, ...)
        return(invisible());
    }
    # save old parameters
    oldask <- par("ask")
    which <- match.arg(which)
    marginal <- match.arg(marginal)
    if(which == "all") par(ask = TRUE)
    if(!is.null(iter) && iter <= 0) iter <- NULL
    # Hazard and survival curve plotting
    if(which == "hazard" | which == "survival" | which == "all"){
        type <- if(which == "all") "hazard" else which
        # get the set of knots to use. If an iteration is given, use knots for
        # that iteration, otherwise use only boundary knots
        if(is.null(x$hazard$spline.knots)){
            knots <- NULL
        }else{
            if(is.null(iter)) {
                knots <- x$history$hazard.spline.knots[1, ];knots <- knots[knots>-Inf]
                knots <- c(min(knots), max(knots))
            }else{ 
                knots <- x$history$hazard.spline.knots[iter, ];knots <- knots[knots>-Inf]
            }
        }
        if(is.null(knots)) knots <- range(x$data$time)
        if(is.null(xlim)) xlim1 = range(x$data$time) else xlim1 <- xlim
        # set the points at which the curve should be predicted
        times = seq(from = max(xlim1[1], min(knots)), to = min(xlim1[2], max(knots)),
            length = npoints)
        if(is.null(xlab)) xlab1 <- "Time" else xlab1 <- xlab
        # set axis labels
        if(type == "hazard"){
            if(is.null(ylab)) ylab1 <- "Hazard" else ylab1 <- ylab
            if(is.null(main)){
                if(is.null(newdata)) main1 <- "Baseline hazard" 
                else main1 <- "Hazard"
            }else{ main1 <- main }
        }
        if(type == "survival"){
            if(is.null(ylab)) ylab1 <- "Survival" else ylab1 <- ylab
            if(is.null(ylim)) ylim <- c(0,1)
            if(is.null(main)){
                if(is.null(newdata)) main1 <- "Baseline survival" 
                else main1 <- "Survival"
            }else{ main1 <- main }
        }
        # a single curve is plotted if either no newdata is given, or newdata has only one row
        if(is.null(newdata) | (!is.null(newdata) && dim(newdata)[1] < 2)){
            # estimate the hazard curve
            haz <- predict(x, type = type, x = times, marginal = marginal, newdata = newdata, iter = iter, 
                fn = fn, npost = npost, alpha = alpha)
            plot(haz[, 1:2], type = "l", lty = lty, col = col, lwd = lwd, main = main1,
                xlab = xlab1, ylab = ylab1, xlim = xlim1, ylim = ylim, ...)
            # plot knots as vertical lines
            if(plotknots) {
                abline(v = knots, col = col.knots, lty = lty.knots, lwd = lwd.knots, ...)
                lines(haz, lty = lty, col = col, lwd = lwd)
            }
            # plot confidence intervals
            if(!is.null(alpha) & is.null(iter)){
                lines(haz[, c(1, 3)], lty = 2, col = col)
                lines(haz[, c(1, 4)], lty = 2, col = col)
            }
        # if newdata has more than one row, multiple lines have to be predicted and plotted
        }else{
            haz <- times
            # predict a curve for each row of newdata
            for(i in 1:dim(newdata)[1]) haz <- cbind(haz, predict(x, type = type, x = times,
                marginal = marginal, newdata = newdata[i, ,drop = FALSE], iter = iter, fn = fn, npost = npost)[, 2])
            if(length(col) == 1 & length(lty) == 1 & length(lwd) == 1) col = 1:dim(newdata)[1]
            # plot all the new lines
            matplot(haz[, 1], haz[, -1], type = "l", col = col, lwd = lwd, lty = lty,
                main = main1, xlab = xlab1, ylab = ylab1, xlim = xlim1, ylim = ylim, ...)
            if(is.null(legend)) legend <- rownames(newdata)
            # plot knots as vertical lines
            if(plotknots) abline(v = knots, col = col.knots, lty = lty.knots,
                lwd = lwd.knots, ...)
            # add a legend
            if(!is.na(legend)) legend("topright", legend = legend, col = col, lty = lty,
                lwd = lwd)
        }
        if(which == "all") plot(x, which = "survival", newdata = newdata, iter = iter,
            plotknots = plotknots, npoints = npoints, npost = npost, alpha = alpha,
            legend = legend, lty = lty, col = col, lwd = lwd, lty.knots = lty.knots,
            col.knots = col.knots, xlab = xlab, ylab = ylab, main = main, xlim = xlim,
            ylim = ylim, tk = tk, ...)
    }
    # Frailty density plotting
    if(which == "frailty" | which == "all"){
        knots <- x$frailty$spline.knots
        # construct knots to use
        if(is.null(iter)){ 
                knots <- x$history$frailty.spline.knots[1, ];knots <- knots[knots>-Inf]
                knots <- c(min(knots), max(knots))
        } else{ 
            knots <- x$history$frailty.spline.knots[iter, ];knots <- knots[knots>-Inf]
        }
        if(is.null(knots)) knots <- range(x$posterior.mean$frailty)
        # limits on the graph
        if(is.null(xlim)) {
               if(hasspline(x$frailty)) xlim1 = attr(knots, "b") 
               else xlim1 <- range(x$posterior.mean$frailty)
       }else{xlim1 <- xlim}
        # construct the set of values at which the curve shoult be estimated
        Ui = seq(from = max(xlim1[1], min(knots)), to = min(xlim1[2], max(knots)),
            length = npoints)
        if(is.null(xlab)) xlab1 <- "x" else xlab1 <- xlab
        if(is.null(ylab)) ylab1 <- "Density" else ylab1 <- ylab
        if(is.null(main)) main1 <- "Frailty density" else main1 <- main
        # comput the curve
        dens <- predict(x, type = "frailty", x = Ui, iter = iter, fn = fn, npost = npost,
            alpha = alpha)
        # plot the density curve
        plot(dens[, 1:2], type = "l", lty = lty, col = col, lwd = lwd, main = main1,
            xlab = xlab1, ylab = ylab1, xlim = xlim1, ylim = ylim, ...)
        # plot knots as vertical lines
        if(plotknots){
            abline(v = knots, col = col.knots, lty = lty.knots, lwd = lwd.knots, ...)
            lines(dens, type = "l", lty = lty, col = col, lwd = lwd)
        }
        # plot CIs
        if(!is.null(alpha) & is.null(iter)){
            lines(dens[, c(1, 3)], lty = 2, col = col)
            lines(dens[, c(1, 4)], lty = 2, col = col)
        }
    }
    # plot coefficient density
    if(which == "coef" | which == "all"){
        burnin <- x$control$burnin
        if(is.null(burnin)) burnin <- x$call$control$burnin
        if(is.null(burnin)) burnin <- dim(x$history$coefficients)[2]
        if(is.null(xlab)) xlab1 <- "x" else xlab1 <- xlab
        if(is.null(ylab)) ylab1 <- "Posterior density" else ylab1 <- ylab
        coefs <- x$history$coefficients
        coefnames <- colnames(coefs)
        if(length(coefnames) > 1) par(ask = TRUE)
        for(i in 1:dim(coefs)[2]){
            if(is.null(main)) main1 <- paste("Coefficient of", coefnames[i]) else main1 <- main
            betai <- coefs[burnin:(dim(coefs)[1]), i]
            plot(density(betai), lty = lty, col = col, lwd = lwd, main = main1, xlab = xlab1, ylab = ylab1, xlim = xlim, ylim = ylim, ...)
        }     
    }
    par(ask = oldask)
}
#************ plot.splinesurv 

#****f* S3Methods/splinesurvtkplot
#  NAME
#    splinesurvtkplot --- plot the curve using tcltk
#  FUNCTION
#    Uses tcltk and the tkrplot package to plot the curve, with a slider to select
#    the iteration.
#  INPUTS
#    x      a splinesurv object
#    ...    plotting parameters
#  SYNOPSIS
splinesurvtkplot <- function(x, ...)
#  SOURCE
#
{
    library(tkrplot)
	which <- "h"
    # tk canvas
	tt <- tktoplevel()
	tktitle(tt) <- "SplineSurv"
	iter <- 0
    # tcl variables, to hold the selected iteration and plot type
	tcliter <- tclVar(iter)
	tclwhich <- tclVar(which)
	maxiter = x$control$maxiter
	res <- max(1, round(maxiter / 100))
    # plotting canvas
	img <- tkrplot(tt, function() plot(x = x, which = which, iter = iter, ...))	
    
    # an inner function to set the iteration
	setiter <- function(...){
        # get the iteration from the tcliter variable
		thisiter <- round(as.numeric(tclvalue(tcliter)))
		if(iter != thisiter){
			assign("iter", thisiter, inherits = TRUE)
            #replot if the iteration has changed
			tkrreplot(img)
		}	
	}
    
    # an inner function to set the type of plot to "hazard"
	setwhich_h <- function(...) {
		assign("which", "h", inherits = TRUE)
		tkrreplot(img)
	}
    # an inner function to set the type of plot to "survival"
	setwhich_s <- function(...) {
		assign("which", "s", inherits = TRUE)
		tkrreplot(img)
	}
    # an inner function to set the type of plot to "frailty"
	setwhich_f <- function(...){
		assign("which", "f", inherits = TRUE)		
		tkrreplot(img)
	}
    # an inner function to set the current iteration to 0, corresponding to
    # plotting the posterior mean.
	setpost <- function(...){
		assign("iter", 0, inherits = TRUE)
		tclvalue(tcliter) <- 0
		tkrreplot(img)
	}

    # a TK scale control, used to select the iteration to plot
	iter_scale <- tkscale(tt, command = setiter, from = 0, to = maxiter, resolution = res,
        showvalue = T, orient = "horiz", variable = tcliter, length = 400)
    # a frame containing the plot type buttons
	ff <- tkframe(tt, relief = "ridge", borderwidth = 2, width = 150, height = 100)
    # buttons to select the type of plot desired
	which_b_h <- tkradiobutton(ff, command = setwhich_h, text = "hazard",
        variable = tclwhich, value = "h" )
	which_b_s <- tkradiobutton(ff, command = setwhich_s, text = "survival",
        variable = tclwhich, value = "s" )
	which_b_f <- tkradiobutton(ff, command = setwhich_f, text = "frailty",
        variable = tclwhich, value = "f" )
    # button to plot the posterior mean of the curve
	post_b <- tkbutton(tt, command = setpost, text = "Posterior" )
    # layout on the grid
	tkgrid(img, columnspan = 2)
	tkgrid(which_b_h, which_b_s, which_b_f)
	tkgrid(ff, columnspan = 2)
	tkgrid(post_b, iter_scale)
	tkgrid.configure(iter_scale, sticky = "e")
	tkgrid.configure(post_b, sticky = "sw")

}

#************ splinesurvtkplot 

################# predict method

#************ setpost 
#****f* S3Methods/predict.splinesurv
#  NAME
#    predict.splinesurv --- prediction method for splinesurv objects
#  FUNCTION
#    Compute curve estimates, or linear predictor/risk estimates for splinesurv.
# 
#    This routine mainly works recursively. If called with iter=NULL, it calls itself
#    with a subset of iterations multiple times and aggregates the results.
#  INPUTS
#    See package documentation
#  SYNOPSIS
predict.splinesurv <- function(object, type = c("hazard", "survival", "lp", "risk", "frailty"),
    x = NULL, marginal = c("none", "mc", "numerical"), newdata = NULL, iter = NULL, 
    fn = mean, alpha = NULL, npost = 100, ...)
#  SOURCE
#
{
    #browser()
    marginal <- match.arg(marginal)
    type <- match.arg(type)
    fit <- object; haz <- 1
    ntimes <- 100
    # check for valid newdata
    if((type == "hazard" | type == "survival") & !is.null(newdata)) if(dim(newdata)[1] > 1) 
        stop("newdata may only have one row")
        
    # get the set of "x" values at which the curve will be predicted
    if(type == "hazard" | type == "survival") {
        if(is.null(x) | is.character(x))   
            x <- seq(from = min(fit$data$time), to = max(fit$data$time), length = ntimes) 
    }
    if(type == "frailty" & is.null(x)) {
        knots <- fit$history$frailty.spline.knots[1, ]
        knots <- knots[knots>-Inf]
        bounds <- range(knots)
        x <- seq(from = bounds[1], to = bounds[2], length = ntimes) 
    }

    # if iter=NULL and a curve is required, call the routine again with a subset of iters
    if((type == "hazard" | type == "frailty" | type == "survival") & is.null(iter)){
        # get the set of iterations that will be used
        ngooditer <- min(fit$control$maxiter - fit$control$burnin, npost)
        iters <- round(seq(from = fit$control$burnin + 1, to = fit$control$maxiter,
            length = ngooditer))
        # matrix for storage of predictions
        preds <- matrix(0, length(x), ngooditer)
        i <- 1
        # call predict for a single iteration (after this if statement)
        for(iter in iters) {
            thispred <- predict(fit, type, x, marginal=marginal, newdata=newdata, iter=iter, ...)
            preds[, i] <- thispred[, 2]
            i <- i + 1
        }
        # aggregate the set of predictions
        pred <- apply(preds, 1, fn, na.rm = TRUE)
        if(!is.null(alpha)) quantiles <- t(apply(preds, 1, quantile,
            probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)) else quantiles <- NULL
        out <- cbind(x, pred, quantiles)
        colnames(out) <- c(colnames(thispred), colnames(quantiles))
        return(as.data.frame(out))
    }

    # similar to above, but handles "risk" as well
    if(type == "hazard" | type == "survival" | (type == "risk" & !is.null(x))){
        if(is.null(x) | is.character(x))   times <- seq(from = min(fit$data$time),
            to = max(fit$data$time), length = ntimes) else times <- x
        # construct a temp hazard curve
        hazard <- fit$hazard
        hazard$x <- times; hazard$haspar <- haspar(hazard); 
        hazard$hasspline <- hasspline(hazard); hazard$spline.norm <- FALSE
        if(is.null(iter)){
            # if iter=NULL, call the routine for a subset of iters as above
            ngooditer <- min(fit$control$maxiter - fit$control$burnin, npost)
            iters <- round(seq(from = fit$control$burnin, to = fit$control$maxiter,
                length = ngooditer))
            preds <- matrix(0, length(times), ngooditer)
            i <- 1
            for(iter in iters) {
                thispred <- predict(fit, type, times, marginal=marginal, newdata=newdata, iter=iter, ...)
                preds[, i] <- thispred[, 2]
                i <- i + 1
            }
            pred <- apply(preds, 1, fn, na.rm = TRUE)
            if(!is.null(alpha)) quantiles <- apply(preds, 1, quantile, 
                probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE) else quantiles <- NULL
            out <- cbind(times, pred, t(quantiles))
            colnames(out) <- c(colnames(thispred), rownames(quantiles))
            return(as.data.frame(out))
        }else{
            # continue constructing a temporary hazard RCurve
            hazard$spline.par <- fit$history$hazard.spline.par[iter, ]
            hazard$spline.knots <- fit$history$hazard.spline.knots[iter, ]
            hazard$param.par <- fit$history$hazard.param.par[iter, ]
            hazard$spline.par <- hazard$spline.par[hazard$spline.par>-Inf]
            hazard$spline.knots <- hazard$spline.knots[hazard$spline.knots>-Inf]
            hazard$weight <- fit$history$hazard.weight[iter, ]
        }
        # set up the basis and evaluate the parametric and spline components
        hazard <- makesplinebasis(hazard, quick = TRUE)
        hazard <- evalparametric(hazard)
        hazard <- evalspline(hazard, quick = TRUE)
        haz <- hazard$y
        # if newdata is nonzero, compute the frailty and covariate multiplier for each row
        if(!is.null(newdata)){
            risk <- predict(fit, "risk", NULL, marginal=marginal, newdata=newdata, iter=iter)$risk[1]
            haz <- haz * risk
            clusterind <- attr(fit$terms, "special")$cluster
            if(!is.null(clusterind)) cluster <- newdata[[clusterind]] else cluster <- NULL
            if(!is.null(cluster) && !is.na(cluster)){ 
                frail <- fit$history$frailty[iter, cluster] 
                marginal <- FALSE  # if frailty is specified, marginal is not allowed
            }else frail <- 1
            haz <- haz * frail
        }
        if(type == "hazard") return(data.frame(time = times, hazard = haz))
        # compute survival curve by cumulative summing
        if(type == "survival"){
            dx <- c(diff(times), 0)
            surv <- exp(-cumsum(haz * dx))
            if(marginal == "mc")
                surv <- apply(outer(surv, 
                    fit$history$frailty[iter, ], function(x, y) x^y), 1, mean)                
            if(marginal == "numerical") {
                frailpred <- predict(fit, "frailty", iter = iter)
                dx <- frailpred$x[2] - frailpred$x[1]
                frailpred$density <- frailpred$density / sum(frailpred$density * dx)
                survfns <- outer(surv, frailpred$x, function(x, y) x^y)
                surv <- survfns %*% frailpred$density * dx
            }
            return(data.frame(time = times, survival = surv))
        }
    } 
    # for lp and risk, evaluate the model terms and do the prediction
    if(type == "lp" | type == "risk")
    {
        if(is.null(newdata)) {
            vars <- model.matrix(fit$terms, model.frame(fit))[, -1, drop = FALSE]
        }else{
            temp <- cbind(data.frame(i = 0, j = 0, time = 0, delta = 0, newdata))
            vars <- model.matrix(fit$terms, model.frame(fit, data = temp))[, -1, drop = FALSE]
        }
        if(is.null(iter))
            coef <- fit$posterior.mean$coefficients
        else
            coef <- fit$history$coefficients[iter, ]
        # compute lp
        lp <- as.numeric(vars%*%coef)
        if(type == "lp") {
            lp <- data.frame(lp)
            try(rownames(lp) <- if(is.null(newdata)) rownames(fit$data) else 
                rownames(newdata), silent = TRUE)
            return(lp)
        }
    }
    # for risk, given lp and hazard, multiply additionally by the frailty
    if(type == "risk")
    {
        pred <- exp(lp)
        if(is.null(newdata)) {
            Ji <- table(fit$data$i)
            if(is.null(iter))
                frailty <- rep(fit$posterior.mean$frailty, Ji)
            else
                frailty <- rep(fit$history$frailty[iter, ], Ji)
            pred <- frailty * pred
        }
        risk <- outer(haz, pred, "*")
        try(colnames(risk) <- if(is.null(newdata)) rownames(fit$data) else 
            rownames(newdata), silent = TRUE)
        if(any(duplicated(colnames(risk)))) colnames(risk) <- NULL
        if(dim(risk)[1] == 1){
             risk <- as.data.frame(t(risk))
             colnames(risk) <- "risk"
        }else{
            risk <- data.frame(time = times, risk)
        }
        return(risk)
    }
    # for frailty curve
    if(type == "frailty")
    {
        # construct a temp frailty curve that can be evaluated
        frailty <- fit$frailty
        frailty$x <- x; frailty$haspar <- haspar(frailty);
        frailty$hasspline <- hasspline(frailty); frailty$spline.norm <- TRUE
        if(is.null(iter)){
            frailty$spline.par <- fit$posterior.mean$frailty.spline.par
            frailty$param.par <- fit$posterior.mean$frailty.param.par
            frailty$weight <- fit$posterior.mean$frailty.weight
            if(type == "frailty") return(data.frame(x = x, density = 1))
        }else{
            frailty$spline.par <- fit$history$frailty.spline.par[iter, ]
            frailty$spline.knots <- fit$history$frailty.spline.knots[iter, ]
            frailty$param.par <- fit$history$frailty.param.par[iter, ]
            frailty$spline.par <- frailty$spline.par[frailty$spline.par>-Inf]
            frailty$spline.knots <- frailty$spline.knots[frailty$spline.knots>-Inf]
            frailty$weight <- fit$history$frailty.weight[iter, ]
        }
        frailty <- makesplinebasis(frailty, quick = TRUE)
        frailty <- evalparametric(frailty)
        frailty <- evalspline(frailty, quick = TRUE)
        density <- frailty$y        
        return(data.frame(x = x, density = density))
    }
}
#************ predict.splinesurv 

#****f* S3Methods/dic.splinesurv
#  NAME
#    dic -- Deviance Information Criterion for a model
#  FUNCTION
#    Compute the DIC for a fitted model
#  INPUTS
#    object     a splinesurv object
#  SYNOPSIS
dic.splinesurv <- function(object){
#  SOURCE
#
    x <- object
    ll<-x$history$loglik[(x$control$burnin+1):x$control$maxiter]
    lm<-mean(ll[abs(ll-median(ll))<5*abs(median(ll))])
    Dbar <- -2*lm
    ll<-ll[ll>quantile(ll,.025) & ll<quantile(ll,.975)]
    V <- var(-2 *ll)
    Dpost <- -2 * likpost.splinesurv(x)
    pD <- V/2
    dic <- Dbar + pD
    dic2 <- Dpost + V
    list(dic=dic,dic2=dic2,Dbar=Dbar,Dpost=Dpost,pD=pD)
}
#************ dic.splinesurv

#****f* S3Methods/likpost.splinesurv
#  NAME
#    likpost -- posterior log-likelihood
#  FUNCTION
#    Compute the log-likelihood at the posterior mean for a fitted model
#  INPUTS
#    object     a splinesurv object
#  SYNOPSIS
likpost.splinesurv <- function(x){
#  SOURCE
#
    P <- x$post
    D <- x$data

    post.haz <- predict(x, "hazard", x=D$time)$hazard
    post.hazcum <- predict(x, "survival", x=seq(from=min(D$time), to=max(D$time), length=1000))
    post.hazcum <- -log(post.hazcum$surv[findInterval(D$time, post.hazcum$time, all=T)])
    post.frail <- predict(x, "frailty", x=P$frailty)$density
    post.lp <- predict(x, "lp")

    # Evaluate full loglikelihood at posterior
    lik <- 0
    lik <- lik + sum(D$delta * (log(P$frailty[D$i]) + log(post.haz) + post.lp))
    lik <- lik - sum(P$frailty[D$i] * post.hazcum * exp(post.lp))
    lik <- lik + sum(log(post.frail))

    lik <- lik - length(P$coef)/2 * log(P$priorvar$coef) - sum(P$coef^2)/(2*P$priorvar$coef)
    lik <- lik - ((x$control$hyper[1]+1) * log(P$priorvar$coef) +
            x$control$hyper[2]/P$priorvar$coef)
    # spline components
    if(hasspline(x$hazard)){
        lik <- lik - (x$hazard$spline.nknots + 2*x$hazard$spline.ord - 2)/2 * 
            log(P$priorvar$hazard.spline) #TODO: add some proxy for smoothness penalty here
        lik <- lik - (x$hazard$spline.hyper[1]+1) * log(P$priorvar$hazard.spline) +
                x$hazard$spline.hyper[2]/P$priorvar$hazard.spline
    }
    if(hasspline(x$frailty)){
        lik <- lik - (x$frailty$spline.nknots + 2*x$frailty$spline.ord)/2 * 
            log(P$priorvar$frailty.spline) #TODO: add some proxy for smoothness penalty here
        lik <- lik - (x$frailty$spline.hyper[1]+1) * log(P$priorvar$frailty.spline) +
                x$frailty$spline.hyper[2]/P$priorvar$frailty.spline
    }
    # parametric components
    if(haspar(x$hazard))
        lik <- lik - length(P$hazard.param.par)/2 * log(P$priorvar$hazard.param) -
                sum(P$hazard.param.par^2)/(2*P$priorvar$hazard.param) -
                (x$hazard$param.hyper[1]+1) * log(P$priorvar$hazard.param) -
                x$hazard$param.hyper[2]/P$priorvar$hazard.param
    if(haspar(x$frailty))
        lik <- lik - length(P$frailty.param.par)/2 * log(P$priorvar$frailty.param) -
                sum(P$frailty.param.par^2)/(2*P$priorvar$frailty.param) -
                (x$frailty$param.hyper[1]+1) * log(P$priorvar$frailty.param) -
                x$frailty$param.hyper[2]/P$priorvar$frailty.param

    # priors on weights
    if(haspar(x$hazard) & hasspline(x$hazard))
        lik <- lik + (x$hazard$weight.hyper[1]-1) * log(P$hazard.weight) +
                     (x$hazard$weight.hyper[2]-1) * log(1 - P$hazard.weight)
    if(haspar(x$frailty) & hasspline(x$frailty))
        lik <- lik + (x$frailty$weight.hyper[1]-1) * log(P$frailty.weight) +
                     (x$frailty$weight.hyper[2]-1) * log(1 - P$frailty.weight)

    # priors on number/position of knots
    if(hasspline(x$hazard) & x$hazard$spline.adaptive)
        lik <- lik + log(1/choose(x$hazard$spline.maxoccknots,x$hazard$spline.nknots)) +
                log(nknotsPrior(round(x$hazard$spline.nknots),x$hazard))
    if(hasspline(x$frailty) & x$frailty$spline.adaptive)
        lik <- lik + log(1/choose(x$frailty$spline.maxoccknots,x$frailty$spline.nknots)) +
                log(nknotsPrior(round(x$frailty$spline.nknots),x$frailty))
        
    lik
}
#************ likpost.splinesurv

}}}

{{{ #Main
##############################################################
# \section{Main} MAIN FUNCTION
##############################################################


#****f* 00main/splinesurv.agdata
#  NAME
#    splinesurv.agdata --- main estimation function
#  FUNCTION
#    This is the main fitting function for the splinesurv routine. It accepts a data
#    frame in a specified format, conducts the initialization procedure, then
#    either starts the MCMC loop within R or calls compiled C code that does the same
#    thing. 
#
#    See also the call graph linked from 00main.
#   
#    This function should not be called directly, instead the interface provided by
#    splinesurv.formula should be used.
#  INPUTS
#    x      a data frame with the following columns:
#           i       cluster index
#           j       subject index
#           time    event time
#           delta   event indicator
#           ...     covariates
#    See the package documentation for detail on the remaining options
#  OUTPUTS
#    an object of class splinesurv, see package documentation.
#  SYNOPSIS
splinesurv.agdata <- function(x, hazard = NULL, frailty = NULL, regression = NULL,
    control = NULL, coda = FALSE, initial = NULL, verbose = 3, usec = TRUE, ...)
#  SOURCE
#
{
        
    if(verbose >= 1) cat("Initializing...\n")
    
    agdata <- x
    rm(x)
    call <- match.call()
    m <- length(unique(agdata$i))
    if(m == 1) warning("Single cluster: frailty density estimate is meaningless")
    Ji <- table(agdata$i)
    
    if(verbose >= 2) cat("\tSetting initial parameters...\n")
    
    # Parse input (control)
    control.in <- control
    # default control options
    control.default <- list(
        burnin = 500, # Length of the burn - in period
        maxiter = 1000, # Max number of iterations
        thin = 1, # Degree of thinning
        tun.auto = TRUE, # Auto - calibrate tuning parameters
        tun.int = 100 # Interval for calibration of the acceptance rate
    )
    control <- control.default
    controlnames <- names(control)
    innames <- names(control.in)
    if(!is.null(control.in)){
        # replace default control settings by any specified in input
        for(n in innames) eval(parse(text = paste("control$",
            match.arg(n, controlnames), " <- control.in$", n, sep = "")))
    }
    if(control$burnin > control$maxiter) stop("Burnin cannot be greater than maxiter")
     
     # Parse input (frailty)   
    frailty.in <- frailty
    # default settings for frailty RCurve options
    frailty.default <- list(
        type = "spline",
        spline.ord = 4,
        spline.adaptive = TRUE,
        spline.knots = NULL,
        spline.knotspacing = "equal",
        spline.nknots = NULL,
        spline.nknots.prior = NULL,
        spline.nknots.hyper = NULL,
        spline.ncandknots = 100,
        spline.bdmconst=.4,
        spline.maxoccknots = 35,
        spline.maxknot = 5,
        spline.par = NULL,
        spline.min=-100,
        spline.penalty = "none",
        spline.penaltyfactor = 1,
        spline.meanpenalty = 1e10,
        spline.priorvar = 0.1,
        spline.hyper = c(0.01, 0.01),
        spline.tun = 1,
        spline.accept = 0,       
        param.dist = "none",
        param.par = NULL,
        param.priorvar = 0.1,
        param.hyper = c(0.01, 0.01),
        param.tun = 1,
        param.accept = 0,
        weight = 0.5,
        weight.priorvar = 0.1,
        weight.hyper = c(1, 2),
        weight.tun = 0.01,
        weight.accept = 0,
        accept = 0
    )
    frailty <- frailty.default
    frailtynames <- names(frailty)
    if(!is.null(frailty.in)){
        # replace default frailty options by input
        for(n in names(frailty.in)) eval(parse(text = paste("frailty$",
            match.arg(n, frailtynames), " <- frailty.in$", n, sep = "")))
    }
    # set additional RCurve settings
    frailty$type <- match.arg(frailty$type, c("spline", "parametric", "both"))
    frailty$hasspline <- frailty$type == "spline" | frailty$type == "both"
    frailty$haspar <- frailty$type == "parametric" | frailty$type == "both"
    frailty$spline.knotspacing <- match.arg(frailty$spline.knotspacing, c("equal"))
    frailty$spline.penalty <- match.arg(frailty$spline.penalty,
        c("2diff", "2deriv", "log2deriv", "none"))
    frailty$param.dist <- match.arg(frailty$param.dist, c("none", "gamma", "lognormal"))
    if(m == 1 & frailty$haspar) stop("parametric component not allowed for single cluster")
    if(frailty$haspar & frailty$param.dist == "none") {
        warning(paste("no distribution specified for frailty parametric component",
            "-- setting to gamma"))
        frailty$param.dist <- "gamma"
    }
    # match prior settings and set default hyperparameters
    frailty$spline.nknots.prior <- match.arg(frailty$spline.nknots.prior,
        c("poisson", "geometric", "poissonmix", "negbin", "power"))
    if(frailty$spline.nknots.prior == "poisson" & is.null(frailty$spline.nknots.hyper))
        frailty$spline.nknots.hyper <- 10
    if(frailty$spline.nknots.prior == "geometric" & is.null(frailty$spline.nknots.hyper))
        frailty$spline.nknots.hyper <- 0.1
    if(frailty$spline.nknots.prior == "poissonmix" & is.null(frailty$spline.nknots.hyper))
        frailty$spline.nknots.hyper <- c(10, 30)
    if(frailty$spline.nknots.prior == "negbin" & is.null(frailty$spline.nknots.hyper))
        frailty$spline.nknots.hyper <- c(2, .1)
    if(frailty$spline.nknots.prior == "power" & is.null(frailty$spline.nknots.hyper))
        frailty$spline.nknots.hyper <- -1 / 2

    frailty$spline.norm <- TRUE
    frailty$name <- "frailty"
    
     # Parse input (hazard)   
    hazard.in <- hazard
    # default hazard RCurve
    hazard.default <- list(
        type = "spline",
        spline.ord = 4,
        spline.adaptive = TRUE,
        spline.knotspacing = "mixed",
        spline.nknots = NULL,
        spline.nknots.prior = NULL,
        spline.nknots.hyper = NULL,
        spline.ncandknots = 100,
        spline.maxoccknots = 35,
        spline.bdmconst=.4,
        spline.knots = NULL,
        spline.par = NULL,
        spline.min=-100,
        spline.penalty = "none",
        spline.penaltyfactor = 1,
        spline.priorvar = 0.1,
        spline.hyper = c(0.01, 0.01),
        spline.tun = 1,      
        spline.accept = 0, 
        param.dist = "none",
        param.par = NULL,
        param.priorvar = 0.1,
        param.hyper = c(0.01, 0.01),
        param.tun = 1,
        param.accept = 0,
        weight = 0.5,
        weight.priorvar = 0.1,
        weight.hyper = c(1, 2),
        weight.tun = 0.01,
        weight.accept = 0
    )
    hazard <- hazard.default
    haznames <- names(hazard)
    if(!is.null(hazard.in)){
        for(n in names(hazard.in)) eval(parse(text = paste("hazard$",
            match.arg(n, haznames), " <- hazard.in$", n, sep = "")))
    }
    # other hazard settings
    hazard$type <- match.arg(hazard$type, c("spline", "parametric", "both"))
    hazard$hasspline <- hazard$type == "spline" | hazard$type == "both"
    hazard$haspar <- hazard$type == "parametric" | hazard$type == "both"
    hazard$spline.knotspacing <- match.arg(hazard$spline.knotspacing,
        c("quantile", "equal", "mixed"))
    hazard$spline.penalty <- match.arg(hazard$spline.penalty,
        c("2diff", "2deriv", "log2deriv", "none"))
    hazard$param.dist <- match.arg(hazard$param.dist,
        c("none", "exponential", "weibull", "lognormal")) 
    if(hazard$haspar & hazard$param.dist == "none") {
        warning(paste("no distribution specified for hazard parametric component",
            "-- setting to weibull"))
        hazard$param.dist <- "weibull"
    }
    # priors and hyperparameters for the number of knots
    hazard$spline.nknots.prior <- match.arg(hazard$spline.nknots.prior,
        c("poisson", "geometric", "poissonmix", "negbin", "power"))
    if(hazard$spline.nknots.prior == "poisson" & is.null(hazard$spline.nknots.hyper))
        hazard$spline.nknots.hyper <- 10
    if(hazard$spline.nknots.prior == "geometric" & is.null(hazard$spline.nknots.hyper))
        hazard$spline.nknots.hyper <- 0.1
    if(hazard$spline.nknots.prior == "poissonmix" & is.null(hazard$spline.nknots.hyper))
        hazard$spline.nknots.hyper <- c(10, 30)
    if(hazard$spline.nknots.prior == "negbin" & is.null(hazard$spline.nknots.hyper))
        hazard$spline.nknots.hyper <- c(2, .1)
    if(hazard$spline.nknots.prior == "power" & is.null(hazard$spline.nknots.hyper))
        hazard$spline.nknots.hyper <- -1 / 2

    # default weight settings
    if(!hazard$haspar) hazard$weight <- 1
    if(!hazard$hasspline) hazard$weight <- 0
    if(!frailty$haspar) frailty$weight <- 1
    if(!frailty$hasspline) frailty$weight <- 0
    
    hazard$spline.norm <- FALSE
    hazard$name <- "hazard"
    hazard$x <- agdata$time
     
    # Parse input (regression)
    reg.in <- regression
    reg.default <- list(
        priorvar = 0.1,
        hyper = c(0.01, 0.01),
        tun = 1,
        accept = 0,
        loglik = 0
    )
    regression <- reg.default
    regnames <- names(regression)
    if(!is.null(reg.in)) for(n in names(reg.in))
        eval(parse(text = paste("regression$", match.arg(n, regnames),
            " <- reg.in$", n, sep = "")))

     # Default settings for the initial number of knots
    if(is.null(hazard$spline.nknots)) hazard$spline.nknots <- if(hazard$spline.adaptive) 
        min(nknotsPriorMean(hazard), hazard$spline.maxoccknots) else 
        max(min(round(sum(Ji) / 4), 35), 1)
    if(is.null(frailty$spline.nknots)) frailty$spline.nknots <- if(frailty$spline.adaptive)
        min(nknotsPriorMean(frailty), frailty$spline.maxoccknots) else 
        max(1, min(round(m / 4), 35), 1)
    
    
    if(verbose >= 2) cat("\tFitting Cox survival models...\n")
    
    # Cox fit with gamma frailties for initial values of Ui and beta
    varnames <- colnames(agdata)[ - (1:4)]
    qvarnames <- paste("`", varnames, "`", sep = "")
    if(m > 1){
        coxfit <- coxph(as.formula(paste("Surv(time, delta)~",
            paste(qvarnames, collapse = " + "), " + frailty(i)")), data = agdata)
        Ui <- exp(coxfit$frail)
        if(var(Ui) < 1e-5) {
            Ui <- 2 * runif(length(Ui))
            Ui <- 1 + (Ui - mean(Ui))
            Ui[Ui < 0] <- 1
        }
    }else{
        coxfit <- coxph(as.formula(paste("Surv(time, delta)~",
            paste(qvarnames, collapse = " + "))), data = agdata)
        Ui <- 1
    } 
    # set initial frailies and regression structure
    frailty$x <- Ui   
    beta <- coxfit$coef
    regression$m <- m
    regression$Ji <- Ji
    regression$covariates <- as.matrix(agdata[, -(1:4)], sum(Ji), length(beta))
    regression$time <- agdata$time
    regression$status <- agdata$delta
    regression$cluster <- as.integer(agdata$i)
    regression <- updateregression(regression, beta)
        
    rm(coxfit)
    
    # Parametric fits
    if(verbose >= 2 & (frailty$haspar | hazard$haspar)) 
        cat("\tFitting parametric components...\n")
    
    hazard <- fitparametric(hazard, agdata)
    frailty <- fitparametric(frailty, Ui)

    # Spline knots
    if(verbose >= 2 & (frailty$hasspline | hazard$hasspline))
        cat("\tComputing spline knots...\n")
    
    hazbounds <- range(agdata$time) + c(-.1, .1) * diff(range(agdata$time))
    hazbounds[hazbounds < 0] <- 0
    hazard <- makeknots(hazard, agdata$time[agdata$delta == 1], bounds = hazbounds)
    frailty <- makeknots(frailty, Ui, bounds = c(0, max(max(Ui), min(2 * max(Ui),
        frailty$spline.maxknot))))
        
    # Evaluate the splines and integrals
    if(verbose >= 2 & (frailty$hasspline | hazard$hasspline))
        cat("\tConstructing spline basis functions...\n")
    
    hazard <- makesplinebasis(hazard, usec = usec)
    frailty <- makesplinebasis(frailty, usec = usec)
       
    # Penalty matrices
    if(verbose >= 2 & (frailty$hasspline | hazard$hasspline))
        cat("\tInitializing penalty matrices...\n")

    hazard <- makepenalty(hazard, usec = usec)
    frailty <- makepenalty(frailty, usec = usec)    
        

    if(verbose >= 2 & (frailty$hasspline | hazard$hasspline))
        cat("\tObtaining initial values for spline parameters...\n")

    {{{# Initial values for the theta vectors
    
    if(hazard$haspar & hazard$hasspline){
        oldhazweight <- hazard$weight
        hazard$weight <- 1;
        hazard <- weightcurve(hazard)
    }
    if(frailty$haspar & frailty$hasspline){
        oldfrailweight <- frailty$weight
        frailty$weight <- 1;
        frailty <- weightcurve(frailty)
    }

    # Initial values for hazard parameters 
    if(hazard$hasspline){
	    theta.haz <- rep(0, hazard$spline.nknots + hazard$spline.ord)
        if(usec){ # use fast C code to compute likelihoods
            par <- as.double(theta.haz); status <- as.double(regression$status);
            lp <- as.double(regression$lp); frailrep <- as.double(rep(frailty$x, Ji));
            hazParY <- as.double(if(hazard$haspar) hazard$param.y else rep(0, length(lp))); 
            hazParYcum <- as.double(if(hazard$haspar) hazard$param.ycum else
                rep(0, length(lp))); 
            weight <- hazard$weight; B <- as.double(hazard$spline.basis);
            C <- as.double(hazard$spline.basiscum);
            P <- as.double(hazard$spline.penaltyfactor * hazard$spline.penaltymatrix);
            penaltyType <- as.integer(pmatch(hazard$spline.penalty,
                c("none", "2diff", "2deriv", "log2deriv")) - 1);
            splinemin <- as.double(hazard$spline.min)
            sigma2 <- hazard$spline.priorvar
            sigma2 <- 100; sigma2target <- hazard$spline.priorvar
             #compute initial values by slowly decreasing the prior variance for stability
            while(sigma2 > sigma2target){
                sigma2 <- as.double(sigma2 / 10)
                opt.theta.haz <- optim(par, fn = cmklik.spline.haz, gr = cmkgr.spline.haz,
                    status = status, lp = lp, frailrep = frailrep, hazParY = hazParY, 
                    hazParYcum = hazParYcum, weight = weight, B = B, C = C, P = P,
                    penaltyType = penaltyType, sigma2 = sigma2, min = splinemin,
                    method = "BFGS", control = list(fnscale=-1))
                par <- as.double(opt.theta.haz$par)
            }
            opt.theta.haz <- optim(par, fn = cmklik.spline.haz, gr = cmkgr.spline.haz,
                status = status, lp = lp, frailrep = frailrep, hazParY = hazParY,
                hazParYcum = hazParYcum, weight = weight, B = B, C = C, P = P,
                penaltyType = penaltyType, sigma2 = sigma2, min = splinemin,
                method = "BFGS", control = list(fnscale=-1, maxit = 1), hessian = TRUE)
            rm(par, status, lp, hazParY, hazParYcum, weight, B, C, P, penaltyType,
                splinemin, sigma2)
            gcout <- gc()
        }else{ # usec=FALSE
            hazard <- updatespline(hazard, theta.haz)
            gcout <- gc()
            sigma2 <- 100; sigma2target <- hazard$spline.priorvar
             #compute initial values by slowly decreasing the prior variance
            while(sigma2 > sigma2target){
                sigma2 <- sigma2 / 10;
                hazard$spline.priorvar <- sigma2
                opt.theta.haz <- optim(hazard$spline.par,
                        fn = mklik.spline.haz,
                        gr = mkgr.spline.haz,
                        method = "BFGS",
                        control = list(fnscale=-1),
                        hazard = hazard,
                        frailty = frailty,
                        regression = regression,
                        hessian = FALSE)
                hazard <- updatespline(hazard, opt.theta.haz$par)
            }
            opt.theta.haz <- optim(hazard$spline.par,
                    fn = mklik.spline.haz,
                    gr = mkgr.spline.haz,
                    method = "BFGS",
                    control = list(fnscale=-1),
                    hazard = hazard,
                    frailty = frailty,
                    regression = regression,
                    hessian = TRUE)
        }
        gcout <- gc()
        hazard <- updatespline(hazard, opt.theta.haz$par)

    }

        # Initial values for frailty parameters
    if(frailty$hasspline){
        frailty$spline.fixedind <- 1
        # parameter vector not including the fixed index
	    theta.frail <- rep(0, frailty$spline.nknots + frailty$spline.ord - 1)
        {
            opt.theta.frail <- optim(theta.frail,
                    fn = mklik.spline.frail.init,
                    method = "BFGS",
                    control = list(fnscale=-1),
                    hazard = hazard,
                    frailty = frailty,
                    regression = regression,
                    hessian = TRUE)    
        }
        gcout <- gc()
        # add the fixed index back in
        theta.frail <- repairfrailtypar(opt.theta.frail$par, frailty$spline.fixedind)
        badtheta <- TRUE; j <- 0;
        # set the mean of the estimated frailty density to be exactly 1
        while(badtheta){
           j <- j + 1
           thetaj <- suppressWarnings(log(-sum(frailty$spline.basisexp[ - j] * exp(theta.frail[ - j])) / frailty$spline.basisexp[j]))
           if(!is.nan(thetaj)) badtheta <- FALSE
        }
        theta.frail[j] <- thetaj; 
        frailty <- updatespline(frailty, theta.frail)
    }

    gcout <- gc()
    # reweight the inital curves
    if(hazard$hasspline & hazard$haspar){
        hazard$weight <- oldhazweight
        hazard <- weightcurve(hazard);
    }
    if(frailty$hasspline & frailty$haspar){
        frailty$weight <- oldfrailweight
        frailty <- weightcurve(frailty)
    }

    gcout <- gc() 
    # Evaluate variances and hessians for candidate generation
    frailty$tun <- diff(range(Ui))^2 / 6
    hess.coef <- mkhess.coef(regression$coefficients, hazard, frailty, regression)
    Sigma.coef <- solve(-hess.coef)
    regression$candcov <- Sigma.coef
    regression$cholcandcov <- chol(Sigma.coef, pivot = TRUE)
    if(hazard$hasspline){
        hess.haz <- opt.theta.haz$hess
        Sigma.haz <- inverthessian(hess.haz)
        hazard$spline.candcov <- Sigma.haz
        hazard$spline.candsd <- rep(1, hazard$spline.maxoccknots + hazard$spline.ord)
        hazard$spline.cholcandcov <- chol(Sigma.haz, pivot = TRUE)
        rm(hess.haz, Sigma.haz)
    }
    if(frailty$hasspline){
        hess.frail <- opt.theta.frail$hess
        Sigma.frail <- inverthessian(hess.frail)
        frailty$spline.candcov <- Sigma.frail
        frailty$spline.candsd <- rep(1, frailty$spline.maxoccknots + frailty$spline.ord)
        frailty$spline.cholcandcov <- chol(Sigma.frail, pivot = TRUE)
        rm(hess.frail, Sigma.frail)
    }
    # construct a numeric hessian for the parametric parameters
    if(hazard$haspar){ 
        # (I don't remember why it doesn't just use numHess.par like for the frailty, but
        #   there must be a reason for this inelegant setup)
        temphaz <- list(haspar = hazard$haspar, hasspline = hazard$hasspline,
        weight = hazard$weight, spline.y = hazard$spline.y, spline.ycum = hazard$spline.ycum,
        name = hazard$name, param.dist = hazard$param.dist, x = hazard$x, y = hazard$y,
        ycum = hazard$ycum, param.y = hazard$param.y, param.ycum = hazard$param.ycum,
        param.par = hazard$param.par, param.priorvar = hazard$param.priorvar)
        eps <- 1e-5
        par <- temphaz$param.par
        hess <- matrix(0, length(par), length(par))
        for(i in 1:length(par)){
            # this constructs numerical gradients by finite differences
            for(j in 1:length(par)){
            par1 <- par;par2 <- par;par3 <- par;
            if(i == j) {par1[i] <- par[i] + eps;par2 <- par1;par3[i] <- par[i] + 2 * eps}
            else {par1[i] <- par[i] + eps; par2[j] <- par[j] + eps; par3 <- par1;
                par3[j] <- par[j] + eps}
            g1 <- (mklik.param.haz(par1, temphaz, frailty, regression) -
                mklik.param.haz(par, temphaz, frailty, regression)) / eps
            g2 <- (mklik.param.haz(par3, temphaz, frailty, regression) -
                mklik.param.haz(par2, temphaz, frailty, regression)) / eps
            hess[i, j] <- (g2 - g1) / eps
            }
        }
        Sigma.par.haz <- inverthessian(hess)
        hazard$param.candcov <- Sigma.par.haz
        hazard$param.cholcandcov <- chol(Sigma.par.haz, pivot = TRUE)
        rm(Sigma.par.haz, hess, temphaz)
    }
    if(frailty$haspar){
        Sigma.par.frail <- inverthessian(numHess.par(frailty$param.par, mklik.param.frail,
            hazard = hazard, frailty = frailty, regression = regression))
        frailty$param.candcov <- Sigma.par.frail
        frailty$param.cholcandcov <- chol(Sigma.par.frail, pivot = TRUE)
        rm(Sigma.par.frail)
    }
    
    
    }}}

    if(frailty$hasspline) frailty$spline.fvar <- frailtysplinefvar(frailty)

    #browser()
    gcout <- gc()
    # Store initial values in parameter history
    history <- inithistory(hazard, frailty, regression, control)
    avg.tunhist <- NULL
    avg.accepthist <- NULL
    

    # an internal function to prepare the curve structures for calling the C code
    prepforCall <- function(curve){
        if(curve$type == "parametric") return(curve)
        if(curve$spline.nknots < curve$spline.maxoccknots){
            addcols <- curve$spline.maxoccknots - curve$spline.nknots
            curve$spline.knots <- c(curve$spline.knots, rep(-Inf, addcols))
            curve$spline.par <- c(curve$spline.par, rep(-Inf, addcols))
            curve$spline.candsd <- c(curve$spline.candsd, rep(-Inf, addcols))
            curve$spline.basis <- cbind(curve$spline.basis, matrix(-Inf, 
                dim(curve$spline.basis)[1], addcols))
            curve$spline.basisint <- c(curve$spline.basisint, rep(-Inf, addcols))
            if(curve$name == "hazard") curve$spline.basiscum <- cbind(curve$spline.basiscum,
                matrix(-Inf, dim(curve$spline.basiscum)[1], addcols))
            if(curve$name == "frailty") curve$spline.basisexp <- c(curve$spline.basisexp,
                rep(-Inf, addcols))
        }
        curve$spline.candocc <- attr(curve$spline.candknots, "occupied")
        return(curve)
    }
    if(usec){
        hazard <- prepforCall(hazard)
        frailty <- prepforCall(frailty)
    }

     # this does nothing, it just creates a useful marker
     main <- function() {}
        
    if(verbose >= 1) cat("Starting MCMC...\n")
    
    iter <- 1 # Counts the recorded iterations
    iter.el <- 0 # Counts elapsed iterations in each thinning cycle
    
    if(verbose >= 3) cat(iter, " ")

    
    while(iter < control$maxiter)
    {
        if(usec){
            # C version of the main loop
            gcout <- gc()
            nexttunint <- iter - iter%%control$tun.int + control$tun.int
            enditer <- min(nexttunint, control$maxiter)
            out <- .Call("SplineSurvMainLoop", hazard, frailty, regression, history, iter,
                enditer, control$thin, verbose) 
            iter <- enditer
        }else{

            # R version of the main loop

            iter.el <- iter.el + 1

            if(verbose >= 4) cat(iter.el)
                
            # MH update of frailties
            frailty <- mh.frail(hazard, frailty, regression)
            
            # MH update of regression parameters
            regression <- mh.coef(hazard, frailty, regression)

            # MH update of baseline parameters
            hazard <- mh.hazard.spline(hazard, frailty, regression)
            
            # MH update of frailty density parameters
            frailty <- mh.frailty.spline(hazard, frailty, regression)
                    
            # MH update of parametric baseline parameters
            hazard <- mh.hazard.param(hazard, frailty, regression)
            
            # MH update of parametric frailty parameters
            frailty <- mh.frailty.param(hazard, frailty, regression)

            # MH update of weights
            hazard <- mh.weight("hazard", hazard, frailty, regression)
            frailty <- mh.weight("frailty", hazard, frailty, regression)

            # Update of the sigmas / taus
            hazard <- updatepostvar.curve(hazard)
            frailty <- updatepostvar.curve(frailty)
            regression <- updatepostvar.coef(regression)

            # Birth - death - move
            if(hazard$spline.adaptive)
                hazard <- mh.bdm("hazard", hazard, frailty, regression)
            if(frailty$spline.adaptive)
                frailty <- mh.bdm("frailty", hazard, frailty, regression)

            # record history
            if(iter.el == control$thin){        
                iter <- iter + 1
                history <- updatehistory(history, iter, hazard, frailty, regression)
                iter.el <- 0
                if(verbose >= 3) cat(" ", iter, " ", sep = "")
            }
        }

        {{{  # Periodic calibration check
        if(iter%%control$tun.int == 0 & iter < control$maxiter){
  
            if(verbose == 1 | verbose == 2) cat(iter, " ")
            
            calinds <- (iter - control$tun.int + 1):iter   

            # Calibration of the tuning parameters for acceptance rate
            if(control$tun.auto & iter <= control$burnin){
               if(verbose >= 3) cat("\n Calibration ...\n")
                      
                          
                tunnames <- c("regression$tun",
                    "hazard$spline.tun",
                    "frailty$spline.tun",
                    "hazard$param.tun",
                    "frailty$param.tun",
                    "hazard$weight.tun",
                    "frailty$weight.tun",
                    "frailty$tun")             
                alltun <- rep(0, length(tunnames))
                for(i in 1:length(alltun)) eval(parse(text = paste("alltun[", i, "] <- ",
                    tunnames[i])))
                avg.tunhist <- rbind(avg.tunhist, alltun)
                avg.accepthist <- rbind(avg.accepthist,
                                      apply(history$accept[calinds, ], 2, mean))
                for(g in 1:length(alltun)){
                    if(all(avg.accepthist[, g]>.25)) alltun[g] <- alltun[g] * 2
                    if(all(avg.accepthist[, g]<.25)) alltun[g] <- alltun[g] / 2
                    if(any(avg.accepthist[, g]>.25) & any(avg.accepthist[, g]<.25)){
                        # build linear regression model for tuning parameters
                        fit <- lm(y~x, data.frame(x = avg.tunhist[, g], 
                            y = avg.accepthist[, g]))
                        out <- try(max(1e-3, nlm(accrate.predict.lm, alltun[g], m = fit)$est))
                        if(inherits(out, "try-error"))
                            alltun[g] <- avg.tunhist[which.min((avg.accepthist-.25)^2)]
                        else
                            alltun[g] <- out
                    }
                }
                for(i in 1:length(alltun)) eval(parse(text = paste(tunnames[i],
                    " <- alltun[", i, "]")))

                if(verbose >= 4){
                    outmat <- cbind(avg.accepthist, alltun);
                    rownames(outmat) <- tunnames;
                    colnames(outmat) <- c("acceptance", "tuning")
                    print(outmat)
                }
            }
            # Print full iteration info
            if(verbose >= 5){
                cat("Frailties: ", submean(history$frailty, calinds), "\n")
                cat("Coefficients: ", submean(history$coefficients, calinds), "\n")
                cat("Hazard.spline: ", submean(history$hazard.spline.par, calinds), "\n")
                cat("Frailty.spline: ", submean(history$frailty.spline.par, calinds), "\n")
                cat("Hazard.param: ", submean(history$hazard.param.par, calinds), "\n")
                cat("Frailty.param: ", submean(history$frailty.param.par, calinds), "\n")
                cat("Prior variances: ", submean(history$priorvar, calinds), "\n")
            }
            
            #browser()
            
        }
        }}}
        
    }   # end main loop

    gcout <- gc()
    if(verbose > 0) cat("Done!\n")
    hazard <- makeoutputcurve(hazard)
    frailty <- makeoutputcurve(frailty)
   
    gcout <- gc()
   {{{ # Construct output
    if(control$burnin < iter){
        if(hasspline(frailty)){
        frail.par.rs <- rowSums(exp(history$frailty.spline.par))
        history$frailty.spline.par <- exp(history$frailty.spline.par) / frail.par.rs
        }
        sub <- (1:(iter)) > (control$burnin)
        posterior.mean <- list(coefficients = submean(history$coefficients, sub),
                            frailty = submean(history$frailty, sub),
                            hazard.spline.par = submean(history$hazard.spline.par, sub),
                            hazard.param.par = submean(history$hazard.param.par, sub),
                            hazard.weight = submean(history$hazard.weight, sub),
                            frailty.spline.par = submean(history$frailty.spline.par, sub),
                            frailty.param.par = submean(history$frailty.param.par, sub),
                            frailty.weight = submean(history$frailty.weight, sub),
                            priorvar = as.list(submean(history$priorvar, sub)),
                            loglik = submean(history$loglik, sub)
                        )
        # save mean number of knots in spline.nknots for output
        if(hasspline(hazard)) 
            hazard$spline.nknots <- mean(apply(history$hazard.spline.knots[sub, ], 1,
                function(x) sum(x>-Inf))) - 2 * hazard$spline.ord
        if(hasspline(frailty)) {
            frailty$spline.nknots <- mean(apply(history$frailty.spline.knots[sub, ], 1,
                function(x) sum(x>-Inf))) - 2 * frailty$spline.ord
            history$frailty.spline.par <- log(history$frailty.spline.par)
            posterior.mean$frailty.spline.par <- log(posterior.mean$frailty.spline.par)
        }
        colnames(history$coefficients) <- varnames
        names(posterior.mean$coefficients) <- varnames
    }else{
        postmean = NULL
    }
    if(coda){
        # convert history to mcmc objects
        library(coda)
        for(i in 1:length(history)) history[[i]] <- as.mcmc(history[[i]])
    }
    # clear history rownames
    rownames(history$frailty) <- rownames(history$coefficients) <-
        rownames(history$splinepar.haz) <- rownames(history$splinepar.frail) <- NULL
    control$iter <- iter
    control$hyper <- regression$hyper
    }}}

    gcout <- gc()
    # output object
    out <- list(call = call, history = history, posterior.mean = posterior.mean,
        hazard = hazard, frailty = frailty, control = control, data = agdata)
    class(out) <- "splinesurv"
    return(out)
}
#************ splinesurv.agdata 
}}}
