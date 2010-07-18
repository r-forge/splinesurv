/***************************
 * mainloop.c
 * C routines for the main loop of the splinesurv package
 * Emmanuel Sharef
 * ************************/

/****h* /CFitting
 * NAME
 *  CFitting --- model fitting routines in C
 * FUNCTION
 *  These routines are C implementations of the ones described in RFitting, and 
 *  work approximately in the same way, adjusting for C's memory management and
 *  direct BLAS access.
 *  
 *  The main function is SplineSurvMainLoop, which calls several Metropolis-Hastings
 *  functions in the CMetropolisHastings module, which in turn rely on likelihood
 *  functions in the CmakeLikelihood module.
 *
 *  This module is organized into several submodules. The main
 *  MCMC loop is contained in CMetropolisHastings, which relies on likelihood
 *  functions in CmakeLikelihood. In order to evaluate likelihoods, B-spline
 *  curves need to be updated regularly, the functionality for which is provided
 *  by CCurveUpdate, which relies on CsplineUtils.
 * CONTENTS
 *  CcurveUpdate --- update curves
 *  CDefines --- control parameter definitions
 *  CmakeLikelihood --- likelihood functions in C
 *  CMetropolisHastings --- MH steps in C
 *  CmiscUtils --- miscellaneous utilities in C
 *  CsplineUtils --- spline utilities in C
 ********/ 

/****h* CFitting/CmiscUtils
 * NAME
 *  CmiscUtils --- miscellaneous utilities in C
 * FUNCTION
 *  Various useful small utilities, including random number generators.
 * CONTENTS
 *  dcopyWrapper --- wrapper for Fortran call to dcopy
 *  ddotWrapper --- Wrapper for Fortran call to ddot
 *  dfactorial --- compute the factorial of a double
 *  diagmvWrapper --- Wrapper for Fortran call to diagmv
 *  dmax --- maximum of two doubles
 *  dmin --- minimum of two doubles
 *  dnegbin --- negative binomial density
 *  EvalNknotsPrior --- evaluate the prior on the number of knots
 *  getListElement --- get an element of a SEXP list by name
 *  imax --- maximum of two integers
 *  imin --- minimum of two integers
 *  mvrnorm --- multivariate normal random numbers
 *  rinvgamma --- inverse gamma random numbers
 *  UpdateHistory --- update history of parameters
 *********/

/****h* CFitting/CsplineUtils
 * NAME
 *  CsplineUtils --- spline utilities in C
 * FUNCTION
 *  Tools to manage and evaluate B-splines and related integrals
 * CONTENTS
 *  cevalBinte --- compute the integrals of each spline basis function
 *  cevalCinte --- construct matrix of partial integrals of a B-spline basis
 *  cevalCinte2 --- compute partial integrals of B-spline basis functions
 *  cevalCinteOld --- old method for computing cevalCinte
 *  cevalEinte --- compute the N-th moment of the distance from 1 of a B-spline
 *  cnBsmom --- compute the N-h moment of a B-spline basis function
 *  cSplineConvolution --- compute the convolution of two basis functions
 *  cSplineDerivInt --- compute the convolution of the derivatives of two spline basis functions
 *  csplinedesign --- create a B-spline basis
 *  csplineeval --- evaluate a B-spline function
 *********/

/****h* CFitting/CcurveUpdate
 * NAME
 *  CcurveUpdate --- update curves
 * FUNCTION
 *  Functions to re-estimate and update curves in the course of estimation.
 * CONTENTS
 *  EvalCurveAtOnePoint --- evaluate a CCurve at a single point
 *  EvalParamAtOnePoint --- evaluate parametric component at a single point
 *  EvalParametric --- evaluate the parametric component of a curve with new parameters
 *  EvalSpline --- evaluate spline component
 *  EvalSplineAtOnePoint --- evaluate spline at a given value of x
 *  FrailtySplineVar --- compute frailty variance of a spline component
 *  MakeSplineBasis --- make a the spline basis from scratch
 *  PopulateLocalCurve --- populate a CCurve structure from an RCurve
 *  PopulateLocalHistory --- populate a local CHistory curve
 *  PopulateLocalRegression --- populate a local RRegression structure
 *  RemakeSplineBasis --- update spline basis after a birth-death-move operation
 *  ReweightCurve --- weight a curve's parametric and spline components
 *  UpdateCurveX --- change an X value of a CCurve structure
 *  UpdateParamPar --- update parametric component with new parameters
 *  UpdateSplineBasis --- recompute basis if one of the X values changes
 *  UpdateSplinePar --- update the spline parameters
 *********/

/****h* CFitting/CmakeLikelihood
 * NAME
 *  CmakeLikelihood --- likelihood functions in C
 * FUNCTION
 *  Functions to compute the loglikelihood required by the MCMC steps, in C
 * CONTENTS
 *    LikelihoodFrailty --- likelihood of the frailty itself
 *    LikelihoodFrailtyLogSum --- utility to efficiently compute a sum of log-frailties
 *    LikelihoodHazardLogSum --- utility to efficiently compute a sum of log-hazards
 *    LikelihoodParamFrailty --- likelihood of frailty parametric parameter
 *    LikelihoodParamHazard --- likelihood of hazard parametric parameters
 *    LikelihoodRegression --- likelihood of regression coefficients
 *    LikelihoodSplineFrailty --- likelihood of frailty spline parameters
 *    LikelihoodSplineHazard --- likelihood of hazard spline parameters
 *    LikelihoodWeightFrailty --- likelihood of weight for frailty curve
 *    LikelihoodWeightHazard --- likelihood of weight for hazard curve
 *    SmoothnessPenalty --- compute smoothness penalties
 *********/

/****h* CFitting/CMetropolisHastings
 * NAME
 *    CMetropolisHastings --- MH steps in C
 * FUNCTION
 *    Functions to successively update each of the parameter sets by Metropolis-Hastings
 *    or reversible jump M-H. These consitute the MCMC loop.
 * CONTENTS
 *    AcceptReject --- accept-reject step for Metropolis-Hastings
 *    MH_BDM --- birth-death-move steps for spline knots
 *    MH_Frail --- MH step for frailties
 *    MH_ParamFrailty --- MH for frailty parametric component parameters
 *    MH_ParamHazard --- MH for hazard parametric component parameters
 *    MH_Regression --- MH step for regression coefficients
 *    MH_SplineFrailty --- MH for frailty spline parameters
 *    MH_SplineHazard --- MH for hazard spline parameters
 *    MH_Weight --- MH for weight of spline component
 *    UpdatePostvarCurve --- update prior variance for a CCurve
 *    UpdatePostvarRegression --- update prior variance for regression coefficients
 *********/

/****h* CFitting/CDefines
 * NAME
 *    CDefines --- control parameter definitions
 * CONTENTS
 *    DEBUG --- debugging marker
 *    LIK_MOD --- modulus to use for summing likelihoods
 *    MAX_PAR --- hard limit on maximum parameter value allowed
 *    enumPenalty --- an enum for different penalty types
 *    enumDistribution --- an enum for different parametric distributions
 *    enumNknotsPrior --- an enum for different priors on the number of knots
 *********/

#include <R_ext/Linpack.h>    
#include <R_ext/Lapack.h>    
#include <R_ext/Print.h>  
#include <R_ext/Utils.h>  
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <math.h>
#include <stdlib.h>
#include "init.h"

/****d* CDefines/DEBUG
 *  NAME
 *    DEBUG --- debugging marker
 *  FUNCTION
 *    If true, additional debug information is printed to the screen
 *  SOURCE
*/
#define DEBUG
/************ DEBUG */

/****d* CDefines/LIK_MOD
 *  NAME
 *    LIK_MOD --- modulus to use for summing likelihoods
 *  FUNCTION
 *    Modulus used when summing loglikelihoods, see LikelihoodHazardLogSum for its use.
 *  SOURCE
*/
#define LIK_MOD 10
/************ LIK_MOD */

/****d* CDefines/MAX_PAR
 *  NAME
 *    MAX_PAR --- hard limit on maximum parameter value allowed
 *  FUNCTION
 *    For numerical reasons, it's necessary to impose an upper limit on spline
 *    parameter values. It's set to 20, because exp(20) should be a sufficiently
 *    large number.
 *  SOURCE
*/
#define MAX_PAR 20
/************ MAX_PAR */

/****d* CDefines/enumPenalty
 *  NAME
 *    enumPenalty --- an enum for different penalty types
 *  FUNCTION
 *    Index the penalty functions to use. Options are:
 *      0   Gaussian penalty
 *      1   Sum of squared second differences
 *      2   Integrated squared second derivative
 *      3   log of 2
 *  SOURCE
*/
typedef enum {pnone=0, pdiff=1, pderiv=2, plogderiv=3} penalty;
/************/

/****d* CDefines/enumDistribution
 *  NAME
 *    enumDistribution --- an enum for different parametric distributions
 *  FUNCTION
 *    Index parametric component distributions to use for parametric components of
 *    the hazard or density. Options are:
 *      0   None
 *      1   Exponential
 *      2   Weibull
 *      3   Lognormal
 *      4   Gamma
 *  SOURCE
*/
typedef enum {Dnone=0, Dexponential=1, Dweibull=2, Dlognormal=3, Dgamma=4} distribution;
/************/

/****d* CDefines/enumNknotsPrior
 *  NAME
 *    enumNknotsPrior --- an enum for different priors on the number of knots
 *  FUNCTION
 *    Index possible priors on the number of knots for use with adaptive model
 *    selection. Options are:
 *      0   Poisson
 *      1   Geometric
 *      2   Poisson Mixture
 *      3   Negative Binomial
 *      4   Power (x^a)
 *  SOURCE
*/
typedef enum {PrPoisson=0, PrGeometric=1, PrPoissonMix=2,
    PrNegBin=3, PrPower=4} nknotsprior;
/************/

/****d* 01structures/CCurve
 *  NAME
 *    CCurve --- structure to store a curve
 *  FUNCTION
 *    This structure contains all the information about a curve, either the hazard
 *    or the frailty. This includes all parameters for spline components, parametric
 *    components, weights, and knot positions, as well as spline basis functions,
 *    tuning parameters, etc. It is the C analogoue of the RCurve structure.
 *
 *    Note that many of the elements of the structure are pointers. This allows
 *    the memory allocated by R to be used directly, rather than having to re-allocate
 *    memory within C.
 *
 *    The components are described by commments in the source.
 *  SOURCE
*/
typedef struct curve {
    int hasSpline, //has a spline component
         hasPar, //has a parametric component
         isHazard, //whether the curve represents hazard or frailty
         SplineOrd, // order of the spline
         SplineAdaptive, // whether the knots should birth/death/move
         SplineNknots, //spline number of interior knots
         SplineNknotsMax, // max number of knots
         SplineNCandKnots, // number of candidate knots
         nx, //number of observations
         nj, //number of basis functions (spline)
         njmax, // max value of nj
         np, // number of parameters (parametric)
         SplineFixedInd; // fixed index for frailty spline
    double SplineFvar, // frailty variance
           SplineEParSum; // sum of exponentials of spline parameters
    penalty SplinePenaltyType; // (0=none, 1=diff, 2=2deriv, 3=log2der)
    distribution ParDist; // Parametric distribution function
    nknotsprior SplineNknotsPrior;
    double *SplineKnots, //Knots of the spline
           *SplineCandKnots, // Candidate knots for adaptive spline
           *SplineCandOcc, // occupied indices for the spline candidate knots
           *SplineNknotsHyper, // parameter for prior on the number of knots
           *SplineBDMConst, //Tuning parameter for birth-death-move
           *SplineBasis, //Basis of the spline
           *SplineBasisCum, //Cumulative basis
           *SplineBasisInt, //Integral over the basis functions
           *SplineBasisExp, //Expectation of each spline component
           *SplinePar, //Spline parameters (theta)
           *SplineEPar, // exponential of Spline parameters (theta)
           *SplineMin, // minimum recommended parameter value
           *SplinePenaltyMatrix, // Penalty matrix on spline parameters
           *SplinePenaltyFactor, // Weight factor for smoothness penalty
           *SplineMeanPenalty, // penalty on the mean (frailty only)
           *SplinePriorvar, //prior variance
           *SplineHyper, //Spline Hyperparameters
           *SplineCandCov, // Covariance matrix of candidates
           *SplineCandSD, // standard deviations of candidates
           *SplineCholCov, //Covariance matrix Cholesky decomposition
           *SplineTun, // Tuning parameter for spline parameters
           *SplineAccept, // Acceptance indicator for spline parameters
           *ParamPar, //Parametric component parameters
           *ParamPriorvar, // Parametric Prior variance
           *ParamCandCov, // Candidate covariance
           *ParamCholCov, // Candidate covariance Cholesky
           *ParamTun, // Parametric tuning parameter
           *ParamHyper, //Parametric hyperparameters
           *ParamAccept, // Acceptance indicator for parametric parameters
           *Weight, // Weight of spline component
           *WeightPriorvar, //Prior variance of weight
           *WeightTun, //Tuning parameter for weight
           *WeightHyper, //Hyperparameters for weight
           *WeightAccept, // Acceptance for weight
           *Accept, // Acceptance indicator for x (only used by frailty)
           *tun, // Tuning parameter (general)
           *X, // Observations
           *SplineY, // Spline estimates
           *ParamY, // Parametric estimates
           *Y, // estimates
           *SplineYcum, // Spline cumulative
           *ParamYcum, // parametric cumulative
           *Ycum; // Cumulative
} *curveP;
/************/

/****d* 01structures/CRegression
 *  NAME
 *    CRegression --- structure to store regression information
 *  FUNCTION
 *    This structure contains all the information about the regression components
 *    of the problem, including covariates, regression coefficients, status
 *    indicators, etc.
 *
 *    The components are described by commments in the source.
 *  SOURCE
*/
typedef struct regress {
    int m,  // # of clusters
        n,  // # of observations
        p,  // # of regression parameters
        *Ji, // subjects per cluster (length m)
        *Jicum, // cumsum of Ji (length m+1)
        *cluster; // Cluster ID
    double *coefficients, // coefficient estimates
           *covariates, // covariate estimates
           *lp, // linear predictor
           *elp, // exponential of linear predictor
           *status, // event status indicators (note: double)
           *time, // event time
           *frailrep, // repeated frailties
           *frailelp, // frailty * elp
           *CandCov,  // covariance matrix of candidates
           *CholCov, // Cholesky factor of covariance
           *priorvar, // variance of prior
           *hyper,   // hyperparameters of prior
           *Accept, // Acceptance indicator
           *tun;    // tuning parameter
} *regressionP;
/************/

/****d* 01structures/CHistory
 *  NAME
 *    CHistory --- structure to store MCMC history
 *  FUNCTION
 *    This structure contains the state of the chain after every stored iteration.
 *
 *    The components are described by commments in the source.
 *  SOURCE
*/
typedef struct hist {
    int ny; // number of iterations
    double *frailty, // values of the frailties
           *coefficients, // regression coefficients
           *HazardSplinePar, // spline parameters for the hazard
           *HazardSplineKnots, // knots for the hazard spline
           *FrailtySplinePar, // parameters for the frailty spline
           *FrailtySplineKnots, // knots for the frailty spline
           *FrailtySplineFvar, // variance of the frailty spline
           *HazardParamPar, // parametric component parameters for the hazard
           *FrailtyParamPar, // parametric component parameters for the frailty
           *HazardWeight, // weight of the spline component for the hazard
           *FrailtyWeight, // weight of the spline component for the frailty
           *priorvar, // prior error variances
           *accept, // acceptance rates for each of the components
           *loglik; // full log-likelihood history
} *historyP;
/************/

/****f* CmiscUtils/dmin
 *  NAME
 *    dmin --- minimum of two doubles
 *  INPUTS
 *    x, y  doubles
 *  OUTPUTS
 *    min(x,y)
 *  SYNOPSIS
*/
static inline double dmin(double x, double y) 
/*
 * SOURCE
 */
{
    return x<y ? x : y;
}
/************ dmin */

/****f* CmiscUtils/dmax
 *  NAME
 *    dmax --- maximum of two doubles
 *  INPUTS
 *   x,y    doubles
 *  OUTPUTS
 *   max(x,y)
 *  SYNOPSIS
 */
static inline double dmax(double x, double y) 
/*
 *  SOURCE
*/
{
    return x>y ? x : y;
}

/************ dmax */
/****f* CmiscUtils/imin
 *  NAME
 *    imin --- minimum of two integers
 *  INPUTS
 *     x,y  integers
 *  OUTPUTS
 *     min(x,y)
 *  SYNOPSIS
 */
static inline double imin(int x, int y)
/*
 *  SOURCE
*/
{
    return x<y ? x : y;
}

/************ imin */
/****f* CmiscUtils/imax
 *  NAME
 *    imax --- maximum of two integers
 *  INPUTS
 *    x,y   integers
 *  OUTPUTS
 *    max(x,y)
 *  SYNOPSIS
 */
static inline double imax(int x, int y)
/*
 *  SOURCE
*/
{
    return x>y ? x : y;
}
/************ imax */

/****f* CmiscUtils/getListElement
 *  NAME
 *    getListElement --- get an element of a SEXP list by name
 *  FUNCTION
 *    Given an R list structure as a SEXP, get a named list element. This is taken
 *    almost verbatim from "Writing R extensions".
 *  INPUTS
 *    list  a list in R memory
 *    str   the name of the desired element
 *  OUTPUTS
 *    elmt  a pointer to the element of list whose name is str
 *  SYNOPSIS
 */
SEXP getListElement(SEXP list, const char *str)
/*
 *  SOURCE
*/
{
    SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
    int i;
    for (i = 0; i < length(list); i++)
     if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
       elmt = VECTOR_ELT(list, i);
       break;
     }
    return elmt;
}
/************ getListElement */

/****f* CmiscUtils/dcopyWrapper
 *  NAME
 *    dcopyWrapper --- wrapper for Fortran call to dcopy
 *  FUNCTION
 *    Copies the first n elements of vector x to vector y.
 *
 *    Wrapper for BLAS call to dcopy, with fixed incx and incy at 1,
 *    to make it easier to copy vectors.
 *  INPUTS
 *    n     number of elements to copy
 *    x     source vector (unchanged on exit)
 *    y     destination vector
 *  SYNOPSIS
 */
static inline void dcopyWrapper( int n, double *x, double *y)
/*
 *  SOURCE
*/
{
    int c1=1;
    F77_CALL(dcopy)(&n, x, &c1, y, &c1);
}
/************ dcopyWrapper */

/****f* CmiscUtils/ddotWrapper
 *  NAME
 *    ddotWrapper --- Wrapper for Fortran call to ddot
 *  FUNCTION
 *    Compute the dot product of vectors x and y
 *
 *    Wrapper for BLAS call to ddot, with fixed incx and incy at 1.
 *  INPUTS
 *    n     number of elements
 *    x     vector 1
 *    y     vector 2
 *  OUTPUTS
 *    dot product of x and y.
 *  SYNOPSIS
 */
static inline double ddotWrapper( int n, double *x, double *y)
/*
 *  SOURCE
*/
{
    int c1=1;
    return F77_CALL(ddot)(&n, x, &c1, y, &c1);
}
/************ ddotWrapper */

/****f* CmiscUtils/diagmvWrapper
 *  NAME
 *    diagmvWrapper --- Wrapper for Fortran call to diagmv
 *  FUNCTION
 *    Compute the element by element product of vectors v1 and v2.
 *    This is equivalent to computing diag(v1) %*% v2.
 *    
 *    Wrapper for BLAS call to dsbmv for this special case.
 *  INPUTS
 *    n     number of elements to multiply
 *    v1    first vector
 *    v2    second vector
 *    out   storage for output vector
 *  SYNOPSIS
 */
static inline void diagmvWrapper(int n, double *v1, double *v2, double *out)
/*
 *  SOURCE
*/
{
    int c1 = 1;
    int c0 = 0;
    double d1 = 1.0;
    double d0 = 0.0;
    char uplo = 'u';
    F77_CALL(dsbmv)( &uplo, &n, &c0,
        &d1, v1, &c1, 
        v2, &c1,
        &d0, out, &c1);
}
/************ diagmvWrapper */

/****f* CmiscUtils/mvrnorm
 *  NAME
 *    mvrnorm --- multivariate normal random numbers
 *  FUNCTION
 *    Generate multivariate normal random numbers, given a mean and Cholesky-factored
 *    covariance matrix.
 *  INPUTS
 *    n         length of the vector to be generated
 *    out       storage for output of length n
 *    mu        mean vector (length n)
 *    CholSigma Cholesky factorization of the covariance matrix nxn
 *    tun       tuning parameter for the variance
 *  SYNOPSIS
 */
static inline void mvrnorm(int n, double *out, double *mu, double *CholSigma, double tun)
/*
 *  SOURCE
*/
{
    double * temp = (double *) malloc(n * sizeof(double));
    // iid N(0,1) random numbers:
    for(int i=0; i<n; i++) temp[i] = rnorm(0,1);
    int c1=1; 
    char trans='T';
    double c0=0;
    double c1d=1;
    double sqrttun = sqrt(tun);
    // set out = 0*out + t(Ch)%*%temp
    F77_CALL(dgemv)(&trans, &n, &n, &sqrttun, CholSigma, &n, temp, &c1, &c0, out, &c1);
    // out = out + mu
    F77_CALL(daxpy)(&n, &c1d, mu, &c1, out, &c1);
    free(temp);
}
/************ mvrnorm */

/****f* CmiscUtils/dfactorial
 *  NAME
 *    dfactorial --- compute the factorial of a double
 *  FUNCTION
 *    Factorial function, needed by dnegbin to compute negative binomial density.
 *  INPUTS
 *    x     double
 *  OUTPUTS
 *    x!
 *  SYNOPSIS
 */
double dfactorial(double x)
/*
 *  SOURCE
*/
{
    return x>1 ? x*dfactorial(x-1) : x;
}
/************ dfactorial */

/****f* CmiscUtils/rinvgamma
 *  NAME
 *    rinvgamma --- inverse gamma random numbers
 *  FUNCTION
 *    Generate inverse gamma random numbers, using the gamma RNG function in R.
 *  SYNOPSIS
 */
static inline double rinvgamma(double shape, double scale)
/*
 *  SOURCE
*/
{
    double out = 1.0 / rgamma(shape, 1.0/scale);
    return out;
}
/************ rinvgamma */

/****f* CmiscUtils/dnegbin
 *  NAME
 *    dnegbin --- negative binomial density
 *  FUNCTION
 *    Compute the density of a NB(r,p) evaluated at x.
 *  SYNOPSIS
 */
static inline double dnegbin(double x, double r, double p)
/*
 *  SOURCE
*/
{
    double out = gammafn(x+r)/( dfactorial(x))/gammafn(r)*pow(p,r)*pow(1-p,x);
    return out;
}
/************ dnegbin */

/****f* CcurveUpdate/FrailtySplineVar
 *  NAME 
 *    FrailtySplineVar --- compute frailty variance of a spline component
 *  FUNCTION
 *    For the spline component of a frailty density curve, compute the corresponding
 *    frailty variance.
 *  INPUTS
 *    frailty   a frailty CCurve
 *  OUTPUTS
 *    fvar      variance of the frailty density
 *  SYNOPSIS
 */
double FrailtySplineVar(curveP frailty)
/*
 *  SOURCE
*/
{
    double fvar=0;
    double * Moment2 = (double *) malloc( frailty->nj * sizeof(double));
    int ord = frailty->SplineOrd;
    int K = frailty->SplineNknots;
    int N=2;
    // vector of second moments of each basis function
    cevalEinte(Moment2,frailty->SplineKnots,&ord,&K,&N);
    for(int i=0;i<frailty->nj;i++) Moment2[i] = 1-Moment2[i];
    // second moment of the entire density
    fvar = ddotWrapper(frailty->nj, Moment2, frailty->SplineEPar)/frailty->SplineEParSum;
    // subtract mean^2
    fvar = fvar - 1.0;
    return fvar;
}
/************ FrailtySplineVar */

/****f* CcurveUpdate/PopulateLocalCurve
 *  NAME
 *    PopulateLocalCurve --- populate a CCurve structure from an RCurve
 *  FUNCTION
 *    At initialization of the SplineSurvMainLoop, this routine populates a local
 *    curve that is easier to work with in C. Much of this just consists of 
 *    creating a set of pointers that point to the appropriate locations in R memory.
 *  INPUTS
 *    theCurve  a CCurve to be populated
 *    Rcurve    the RCurve that contains the current values
 *  SYNOPSIS
 */
void PopulateLocalCurve( curveP theCurve, SEXP Rcurve)
/*
 *  SOURCE
*/
{
    // has a parametric component?
    theCurve->hasPar = asInteger(getListElement(Rcurve,"haspar"));
    // has a spline component?
    theCurve->hasSpline = asInteger(getListElement(Rcurve,"hasspline"));
    // does it use adaptive knot selection?
    theCurve->SplineAdaptive = asInteger(getListElement(Rcurve,"spline.adaptive"));
    // length of observations (times/frailties)
    theCurve->nx = (int) (length(getListElement(Rcurve,"x")));
    // minimum spline parameter value
    theCurve->SplineMin = REAL(getListElement(Rcurve,"spline.min"));
    // spline prior varianc3
    theCurve->SplinePriorvar =  REAL(getListElement(Rcurve,"spline.priorvar"));
    // spline hyperparameters
    theCurve->SplineHyper =  REAL(getListElement(Rcurve,"spline.hyper"));
    // spline tuning parameters
    theCurve->SplineTun =  REAL(getListElement(Rcurve,"spline.tun"));
    // spline accept/reject of the last step
    theCurve->SplineAccept =  REAL(getListElement(Rcurve,"spline.accept"));
    // parameter prior variance
    theCurve->ParamPriorvar = REAL(getListElement(Rcurve,"param.priorvar"));
    // parametric tuning parameter
    theCurve->ParamTun = REAL(getListElement(Rcurve,"param.tun"));
    // parametric hyperparameters
    theCurve->ParamHyper = REAL(getListElement(Rcurve,"param.hyper"));
    // parametric accept/reject of last step
    theCurve->ParamAccept = REAL(getListElement(Rcurve,"param.accept"));
    // weight of spline component
    theCurve->Weight = REAL(getListElement(Rcurve,"weight"));
    // prior variance of weight
    theCurve->WeightPriorvar = REAL(getListElement(Rcurve,"weight.priorvar"));
    // hyperparameters of weight
    theCurve->WeightHyper = REAL(getListElement(Rcurve,"weight.hyper"));
    // tuning parameter of weight
    theCurve->WeightTun = REAL(getListElement(Rcurve,"weight.tun"));
    // accept/reject for last step of weight
    theCurve->WeightAccept = REAL(getListElement(Rcurve,"weight.accept"));
    // vector of observations
    theCurve->X = REAL(getListElement(Rcurve,"x"));
    // vector of the curve evaluated at X
    theCurve->Y = REAL(getListElement(Rcurve,"y"));
    // hazard curves also need a vector of cumulative hazard integrals
    if(theCurve->isHazard) theCurve->Ycum =REAL(getListElement(Rcurve,"ycum"));
    // frailty curve needs to keep track of acceptance of the frailties themselves
    if(!theCurve->isHazard) theCurve->Accept =REAL(getListElement(Rcurve,"accept"));
    // tuning parameter for the frailties
    if(!theCurve->isHazard) theCurve->tun = REAL(getListElement(Rcurve,"tun"));

    // special terms for spline components
    if(theCurve->hasSpline){  
        // order of the spline
        theCurve->SplineOrd=asInteger(getListElement(Rcurve,"spline.ord"));
        // number of knots of the spline
        theCurve->SplineNknots = asInteger(getListElement(Rcurve,"spline.nknots"));
        // max number of knots
        theCurve->SplineNknotsMax = asInteger(getListElement(Rcurve,"spline.maxoccknots"));
        // hyperparameters for the number of knots
        theCurve->SplineNknotsHyper = REAL(getListElement(Rcurve,"spline.nknots.hyper"));
        // number of candidate knots
        theCurve->SplineNCandKnots = asInteger(getListElement(Rcurve,"spline.ncandknots"));
        // vector of candidate knot positions
        theCurve->SplineCandKnots = REAL(getListElement(Rcurve,"spline.candknots"));
        // vector of occupancies of candidate knots
        theCurve->SplineCandOcc = REAL(getListElement(Rcurve,"spline.candocc"));
        // constant used for birth-death-move step
        theCurve->SplineBDMConst = REAL(getListElement(Rcurve,"spline.bdmconst"));
        // number of spline basis functions
        theCurve->nj = theCurve->SplineNknots + theCurve->SplineOrd;
        // get the type of prior on the number of knots
        const char * charNknotsPrior = 
            CHAR(STRING_ELT(getListElement(Rcurve,"spline.nknots.prior"),0));
        nknotsprior iNknotsPrior;
        if(strcmp(charNknotsPrior,"poisson")==0) iNknotsPrior=PrPoisson;
        if(strcmp(charNknotsPrior,"geometric")==0) iNknotsPrior=PrGeometric;
        if(strcmp(charNknotsPrior,"poissonmix")==0) iNknotsPrior=PrPoissonMix;
        if(strcmp(charNknotsPrior,"negbin")==0) iNknotsPrior=PrNegBin;
        if(strcmp(charNknotsPrior,"power")==0) iNknotsPrior=PrPower;
        theCurve->SplineNknotsPrior = iNknotsPrior;

        // get the type of smoothness penalty
        const char * charPenaltyType = 
            CHAR(STRING_ELT(getListElement(Rcurve,"spline.penalty"),0));
        penalty iPenaltyType;
        if (strcmp(charPenaltyType,"2diff") ==0) { iPenaltyType=pdiff; }
        else if (strcmp(charPenaltyType,"2deriv") ==0) { iPenaltyType=pderiv; }
        else if (strcmp(charPenaltyType,"log2deriv") ==0) { iPenaltyType=plogderiv;}
        else { iPenaltyType=pnone;}
        theCurve->SplinePenaltyType = iPenaltyType;
        if(iPenaltyType != pnone) // if there is a penalty matrix, get it
            theCurve->SplinePenaltyMatrix = 
                REAL(getListElement(Rcurve, "spline.penaltymatrix"));
        // scaling factor for penalty
        theCurve->SplinePenaltyFactor = REAL(getListElement(Rcurve,"spline.penaltyfactor"));
        // vector of spline knots
        theCurve->SplineKnots =  REAL(getListElement(Rcurve,"spline.knots"));
        // spline basis matrix
        theCurve->SplineBasis =  REAL(getListElement(Rcurve,"spline.basis"));
        // vector of integrals of each basis function
        theCurve->SplineBasisInt =  REAL(getListElement(Rcurve,"spline.basisint"));
        // spline parameters
        theCurve->SplinePar =  REAL(getListElement(Rcurve,"spline.par"));
        // exponentials of spline parameters
        theCurve->SplineEPar = (double *) malloc( 
                (theCurve->SplineNknotsMax + theCurve->SplineOrd) * sizeof(double));
        for(int j=0; j< theCurve->nj; j++)
            theCurve->SplineEPar[j] = exp(theCurve->SplinePar[j]);
        // covariance matrix for candidate generation
        theCurve->SplineCandCov =  REAL(getListElement(Rcurve,"spline.candcov"));
        // cholesky factorization of the covariance matrix
        theCurve->SplineCholCov =  REAL(getListElement(Rcurve,"spline.cholcandcov"));
        // standard deviations for candidate generation
        theCurve->SplineCandSD = REAL(getListElement(Rcurve,"spline.candsd"));
        // spline component evaluated at X
        theCurve->SplineY =  REAL(getListElement(Rcurve,"spline.y"));
        if(theCurve->isHazard){
            // hazards additionally need cumulative basis functions and cumulative
            // spline component integrals
            theCurve->SplineBasisCum =  REAL(getListElement(Rcurve,"spline.basiscum"));
            theCurve->SplineYcum =  REAL(getListElement(Rcurve,"spline.ycum"));
        }
        if(!theCurve->isHazard){
            // frailties need integrals of each basis function
            theCurve->SplineBasisInt =  REAL(getListElement(Rcurve,"spline.basisint"));
            // 1-expected value of each basis function
            theCurve->SplineBasisExp =  REAL(getListElement(Rcurve,"spline.basisexp"));
            // sum of exponentials of spline parameters
            theCurve->SplineEParSum = 0;
            for(int j=0; j< theCurve->nj; j++) 
                theCurve->SplineEParSum += theCurve->SplineEPar[j];
            // variance of spline component of frailty
            theCurve->SplineFvar = FrailtySplineVar(theCurve);
            // penalty on the distance of the mean from 1
            theCurve->SplineMeanPenalty = REAL(getListElement(Rcurve,"spline.meanpenalty"));
            // which of the spline parameters is held fixed
            theCurve->SplineFixedInd = 
                asInteger(getListElement(Rcurve, "spline.fixedind")) - 1;
        }
    }

    // parametric component only
    if(theCurve->hasPar){
        // parametric component parameters
        theCurve->ParamPar = REAL(getListElement(Rcurve,"param.par"));
        // parametric component evaluated at X
        theCurve->ParamY = REAL(getListElement(Rcurve,"param.y"));
        // cumulative integral for hazard
        if(theCurve->isHazard)
            theCurve->ParamYcum = REAL(getListElement(Rcurve,"param.ycum"));
        // get parametric distribution
        const char * charParamDist = CHAR(STRING_ELT(getListElement(Rcurve,"param.dist"),0));
        distribution iParamDist;
        if(strcmp(charParamDist,"exponential") == 0) { iParamDist = Dexponential; } 
        else if(strcmp(charParamDist,"weibull") == 0) { iParamDist = Dweibull; } 
        else if(strcmp(charParamDist,"gamma") == 0) { iParamDist = Dgamma; } 
        else if(strcmp(charParamDist,"lognormal") == 0) { iParamDist = Dlognormal; } 
        else { iParamDist = Dnone; } 
        theCurve->ParDist = iParamDist;
        // number of parametric component parameters
        theCurve->np = (int) (length(getListElement(Rcurve,"param.par")));
        // candidate generation for parametric component parameters
        theCurve->ParamCandCov =  REAL(getListElement(Rcurve,"param.candcov"));
        theCurve->ParamCholCov =  REAL(getListElement(Rcurve,"param.cholcandcov"));
    }
}
/************ PopulateLocalCurve */

/****f* CcurveUpdate/PopulateLocalRegression
 *  NAME
 *    PopulateLocalRegression --- populate a local RRegression structure
 *  SYNOPSIS
 *    void PopulateLocalRegression( regressionP theReg, SEXP Rregression){
 *  FUNCTION
 *  INPUTS
 *  OUTPUTS
 *  SOURCE
*/
void PopulateLocalRegression( regressionP theReg, SEXP Rregression){
    theReg->m = asInteger(getListElement(Rregression, "m"));
    theReg->Ji = INTEGER(getListElement(Rregression, "Ji"));
    theReg->Jicum = (int *) malloc( (theReg->m + 1) * sizeof(int));
    theReg->Jicum[0] = 0;
    for(int i=0; i < theReg->m; i++) theReg->Jicum[i+1] = theReg->Jicum[i] + theReg->Ji[i];
    theReg->n = (int) (length(getListElement(Rregression,"cluster")));
    theReg->p = (int) (length(getListElement(Rregression,"coefficients")));
    theReg->cluster = INTEGER(getListElement(Rregression, "cluster"));

    theReg->coefficients = REAL(getListElement(Rregression,"coefficients"));
    theReg->covariates = REAL(getListElement(Rregression,"covariates"));
    theReg->lp = REAL(getListElement(Rregression,"lp"));
    theReg->status = REAL(getListElement(Rregression,"status"));
    theReg->time = REAL(getListElement(Rregression,"time"));
    theReg->frailrep = (double *) malloc( theReg->n * sizeof(double));    
    theReg->frailelp = (double *) malloc( theReg->n * sizeof(double));    
    theReg->CandCov = REAL(getListElement(Rregression,"candcov"));
    theReg->CholCov = REAL(getListElement(Rregression,"cholcandcov"));
    theReg->priorvar = REAL(getListElement(Rregression,"priorvar"));
    theReg->hyper = REAL(getListElement(Rregression,"hyper"));
    theReg->Accept = REAL(getListElement(Rregression,"accept"));
    theReg->tun = REAL(getListElement(Rregression,"tun"));
    theReg->elp = (double *) malloc( theReg->n * sizeof(double));
    for(int i=0; i<theReg->n; i++) theReg->elp[i]=exp(theReg->lp[i]);
}

/************ PopulateLocalRegression */

/****f* CcurveUpdate/PopulateLocalHistory
 *  NAME
 *    PopulateLocalHistory --- populate a local CHistory curve
 *  FUNCTION
 *    Given an RHistory, populate a CHistory structure with pointers to the
 *    memory allocated by R.
 *  INPUTS
 *  OUTPUTS
 *  SYNOPSIS
 */
void PopulateLocalHistory( historyP theHist, SEXP Rhistory)
/*
 *  SOURCE
*/
{
    SEXP elmt;
    // frailties
    elmt = getListElement(Rhistory, "frailty");
    theHist->ny = INTEGER(getAttrib(elmt, R_DimSymbol))[0];
    theHist->frailty = REAL(elmt);
    // regression coefficients
    theHist->coefficients = REAL(getListElement(Rhistory, "coefficients"));
    // Only populate history components if there is a spline/parametric component
    // hazard spline
    elmt = getListElement(Rhistory, "hazard.spline.par");
    if(elmt != R_NilValue ) theHist->HazardSplinePar = REAL(elmt);
    elmt = getListElement(Rhistory, "hazard.spline.knots");
    if(elmt != R_NilValue ) theHist->HazardSplineKnots = REAL(elmt);
    // hazard parametric
    elmt = getListElement(Rhistory, "hazard.param.par");
    if(elmt != R_NilValue ) theHist->HazardParamPar = REAL(elmt);
    // hazard weight
    elmt = getListElement(Rhistory, "hazard.weight");
    if(elmt != R_NilValue ) theHist->HazardWeight = REAL(elmt);
    // frailty spline
    elmt = getListElement(Rhistory, "frailty.spline.par");
    if(elmt != R_NilValue ) theHist->FrailtySplinePar = REAL(elmt);
    elmt = getListElement(Rhistory, "frailty.spline.knots");
    if(elmt != R_NilValue ) theHist->FrailtySplineKnots = REAL(elmt);
    elmt = getListElement(Rhistory, "frailty.spline.fvar");
    if(elmt != R_NilValue ) theHist->FrailtySplineFvar = REAL(elmt);
    // frailty parametric 
    elmt = getListElement(Rhistory, "frailty.param.par");
    if(elmt != R_NilValue ) theHist->FrailtyParamPar = REAL(elmt);
    elmt = getListElement(Rhistory, "frailty.weight");
    if(elmt != R_NilValue ) theHist->FrailtyWeight = REAL(elmt);
    theHist->priorvar = REAL(getListElement(Rhistory, "priorvar"));
    theHist->accept = REAL(getListElement(Rhistory, "accept"));
    theHist->loglik = REAL(getListElement(Rhistory, "loglik"));
}
/************ PopulateLocalHistory */

/****f* CmakeLikelihood/SmoothnessPenalty
 *  NAME
 *    SmoothnessPenalty --- compute smoothness penalties
 *  FUNCTION
 *    Compute the smoothness penalty for a spline CCurve as part of
 *    making likelihoods.
 *
 *    The penalty can be either a penalty on the sum of squared second differences
 *    or the integrated squared second derivative, or a Gaussian penalty. The
 *    penalty matrix and penalty type stored in the CCurve structure should be
 *    of matching type.
 *  INPUTS
 *    theCurve  CCurve structure
 *  OUTPUTS
 *    pen   the smoothness penalty, already divided by 2*(prior variance)
 *          that is, the entire smoothness penalty term.
 *  SYNOPSIS
 */
double SmoothnessPenalty(curveP theCurve)
/*
 *  SOURCE
*/
{
    double pen;
    int c1 = 1;
    // Gaussian penalty
    if(theCurve->SplinePenaltyType == pnone){ 
        pen = pow(F77_CALL(dnrm2)(&(theCurve->nj), theCurve->SplinePar, &c1),2.0);
    }else{
        double * temp = (double *) calloc(theCurve->nj, sizeof(double));
        double * par;
        // 2diff penalty works directly with parameters, 2deriv works with exp(theta)
        if(theCurve->SplinePenaltyType == pdiff) par = theCurve->SplinePar;
        else par = theCurve->SplineEPar;
        double c1d = 1.0;
        char uplo = 'U';
        // compute par %*% P %*% par
        F77_CALL(dsymv)( &uplo, &(theCurve->nj), &c1d, theCurve->SplinePenaltyMatrix,
                &(theCurve->nj), par, &c1, &c1d, temp, &c1);
        pen = F77_CALL(ddot)( &(theCurve->nj), temp, &c1, par, &c1);
        if(theCurve->SplinePenaltyType == plogderiv) pen=log(pen+1);
        free(temp);
    }
    // normalize by prior variance
    pen = pen / (2 * theCurve->SplinePriorvar[0]);
    return pen;
}
/************ SmoothnessPenalty */

/****f* CcurveUpdate/ReweightCurve
 *  NAME
 *    ReweightCurve --- weight a curve's parametric and spline components
 *  FUNCTION
 *    After updating either the parametric or spline component of a curve, or their
 *    relative weight, the total curve must be recomputed.
 *  INPUTS
 *    theCurve  a CCurve structure
 *    i         the index of the observation that should be reweighted, or -1 for all
 *  OUTPUTS
 *    theCurve  the original curve, with Y and Ycum components updated.
 *  SYNOPSIS
 */
void ReweightCurve(curveP theCurve, int i){
/*
 *  SOURCE
*/
    if(theCurve->hasPar & !theCurve->hasSpline) { //parametric only
        // just copy ParamY and ParamYcum to Y and Ycum
        dcopyWrapper(theCurve->nx, theCurve->ParamY, theCurve->Y);
        if(theCurve->isHazard) dcopyWrapper(theCurve->nx, theCurve->ParamYcum, theCurve->Ycum);
        return;
    }
    if(!theCurve->hasPar & theCurve->hasSpline) { //spline only
        // just copy SplineY and SplineYcum to Y and Ycum
        dcopyWrapper(theCurve->nx, theCurve->SplineY, theCurve->Y);
        if(theCurve->isHazard) dcopyWrapper(theCurve->nx, theCurve->SplineYcum, theCurve->Ycum);
        return;
    }
    double w = theCurve->Weight[0];
    if(i>=0){
        // reweight a single observation
        theCurve->Y[i] = w *  theCurve->SplineY[i] + (1-w) * theCurve->ParamY[i];
        if(theCurve->isHazard) theCurve->Ycum[i] = w *  theCurve->SplineY[i] + (1-w) *
            theCurve->ParamY[i];
    }else{
        // reweight the entire curve
        double c0=0;
        int c0i=0;
        int c1=1;
        int n = theCurve->nx;
        double wm1=1-w;
        F77_CALL(dcopy)(&n, &c0, &c0i, theCurve->Y, &c1); //set Y=0
        F77_CALL(daxpy)(&n, &w, theCurve->SplineY, &c1, theCurve->Y, &c1); //add w*splineY
        F77_CALL(daxpy)(&n, &wm1, theCurve->ParamY, &c1, theCurve->Y, &c1); //add (1-w)*paramY
        if(theCurve->isHazard){
            F77_CALL(dcopy)(&n, &c0, &c0i, theCurve->Ycum, &c1); //set Y=0
            F77_CALL(daxpy)(&n, &w, theCurve->SplineYcum, &c1, theCurve->Ycum, &c1); //add w*splineY
            F77_CALL(daxpy)(&n, &wm1, theCurve->ParamYcum, &c1, theCurve->Ycum, &c1); //add (1-w)*paramY
        }
    }
}
/************ ReweightCurve */

/****f* CcurveUpdate/UpdateSplineBasis
 *  NAME
 *    UpdateSplineBasis --- recompute basis if one of the X values changes
 *  FUNCTION
 *    Re-evaluates the spline basis at its X values. This is called indirectly by MH_Frail, as
 *    well as during the Birth-Death-Move steps in MH_BDM.
 *
 *    The function can evaluate only a subset of spline basis functions (as needed when
 *    changing the knot set), or evaluate
 *    at a single observation (as when changing one frailty).
 *  INPUTS
 *    theCurve      CCurve structure to be re-evaluated
 *    i             index of the observation to evaluate, or -1 for all
 *    startj        start index of the basis functions to evaluate
 *    endj          stop index of the basis functions to evaluate
 *  OUTPUTS
 *  SYNOPSIS
 */
void UpdateSplineBasis(curveP theCurve, int i, int startj, int endj)
/*
 *  SOURCE
*/
{
    if(!theCurve->hasSpline) return;
    if(i<0)
        // evaluate every row separately
        for(int ind=0; ind<theCurve->nx; ind++) UpdateSplineBasis(theCurve, ind, startj, endj);
    else{
        for(int j=startj; j<endj; j++){
            // evaluate the basis at entry [i,j]
            theCurve->SplineBasis[i + j * theCurve->nx] = csplineeval(theCurve->X[i], j,
                    theCurve->SplineOrd, theCurve->SplineKnots, theCurve->SplineOrd, theCurve->nj);
            if(!theCurve->isHazard) 
                theCurve->SplineBasis[i+j * theCurve->nx] /= theCurve->SplineBasisInt[j];
        }
        // for the hazard, also evaluate cumulative hazards
        if(theCurve->isHazard) 
            cevalCinte2(theCurve->SplineBasisCum, theCurve->X, theCurve->nx, theCurve->SplineKnots,
                theCurve->nj, theCurve->SplineOrd, theCurve->SplineBasisInt, i, startj, endj);
    }
}
/************ UpdateSplineBasis */

/****f* CcurveUpdate/MakeSplineBasis
 *  NAME
 *    MakeSplineBasis --- make a the spline basis from scratch
 *  FUNCTION
 *    Function to make the spline basis, analogous to makesplinebasis in R. It can
 *    be called when the knots, or X values of the curve have changed.
 *    This turned out to be rather inefficient, and was replaced by RemakeSplineBasis,
 *    but this function remains in the codebase for debugging.
 *  INPUTS
 *    theCurve      a CCurve structure
 *  OUTPUTS
 *    theCurve, with SplineBasis, SplineBasisCum, SplineBasisInt and SplineBasisExp updated
 *  SYNOPSIS
 */
void MakeSplineBasis(curveP theCurve)
/*
 *  SOURCE
*/
{
    if(!theCurve->hasSpline) return;
    double * knots = theCurve->SplineKnots;
    int ord = theCurve->SplineOrd;
    int nknots = theCurve->SplineNknots;
    int nj = theCurve->nj;
    double * x = theCurve->X;
    int c1 = 1;

    // compute integrals over each basis spline
    cevalBinte(theCurve->SplineBasisInt, knots, &ord, &nknots);

    // For frailty, compute the expectations
    if(!theCurve->isHazard) 
        cevalEinte(theCurve->SplineBasisExp, knots, &ord, &nknots, &c1);
    
    // Make basis
    UpdateSplineBasis(theCurve,-1,0,nj); 
    
}
/************ MakeSplineBasis */

/****f* CcurveUpdate/RemakeSplineBasis
 *  NAME
 *    RemakeSplineBasis --- update spline basis after a birth-death-move operation
 *  FUNCTION
 *    This function is called by MH_BDM, to update the spline basis after adding, moving
 *    or deleting a knot. It's a more efficient version of MakeSplineBasis that only
 *    updates the basis functions required by the operation.
 *
 *    For move steps, this means simply updating a subset of basis functions. For birth and move
 *    steps, existing basis functions are moved, and the intermediate ones updated.
 *  INPUTS
 *    theCurve      a CCurve structure
 *    oper          character, either 'm' for move, 'b' for birth or 'd' for death.
 *    j             index of the knot that was moved, added or deleted.
 *  OUTPUTS
 *    theCurve, with SplineBasis, SplineBasisCum, SplineBasisInt and SplineBasisExp updated
 *  SYNOPSIS
 */
void RemakeSplineBasis(curveP theCurve, char oper, int j)
/*
 *  SOURCE
*/
{
    if(!theCurve->hasSpline) return;
    double * knots = theCurve->SplineKnots;
    int ord = theCurve->SplineOrd;
    int nknots = theCurve->SplineNknots;
    int nj = theCurve->nj;
    int nx = theCurve->nx;
    double * x = theCurve->X;
    int c1 = 1;

    // compute integrals over each basis spline
    cevalBinte(theCurve->SplineBasisInt, knots, &ord, &nknots);

    // For frailty, compute the expectations
    if(!theCurve->isHazard) 
        cevalEinte(theCurve->SplineBasisExp, knots, &ord, &nknots, &c1);
        
    if(oper=='m'){
        // Update basis from j to j+ord (j here is the moveind, in internal knot numbering)
        UpdateSplineBasis(theCurve,-1,j,imin(nj,j+ord+1)); 
    }
    if(oper=='d'){
        //move basis j+1 through nj into j through nj-1
        int nmv = ((nj) - (j+1))*nx;
        nmv = nmv>0 ? nmv : 0;
        if(nmv>0){
            F77_CALL(dcopy)(&nmv, theCurve->SplineBasis + (j+1)*nx, &c1,
                    theCurve->SplineBasis + j*nx, &c1);
            if(theCurve->isHazard)
                F77_CALL(dcopy)(&nmv, theCurve->SplineBasisCum + (j+1)*nx,
                    &c1, theCurve->SplineBasisCum + j*nx, &c1);
        }
        // update basis functions j-ord through j
        UpdateSplineBasis(theCurve,-1,imax(0,j-ord),j);
    }
    if(oper=='b'){
        // //move basis j+1 : nj to j+2 : nj+1
        for(int k=nj-2;k>j;k--) F77_CALL(dcopy)(&nx, theCurve->SplineBasis + k*nx, &c1,
                theCurve->SplineBasis + (k+1)*nx, &c1);
        if(theCurve->isHazard) for(int k=nj-2;k>j;k--) F77_CALL(dcopy)(&nx,
                theCurve->SplineBasisCum + k*nx, &c1, theCurve->SplineBasisCum + (k+1)*nx, &c1);
        // update basis from j-ord to j+2
        UpdateSplineBasis(theCurve,-1,imax(0,j-ord),imin(nj,j+2));
    }
}
/************ RemakeSplineBasis */

/****f* CcurveUpdate/EvalSpline
 *  NAME
 *    EvalSpline --- evaluate spline component
 *  FUNCTION
 *    Evaluate the spline component if the spline parameters have been changed. This function
 *    assumes the basis is correctly specified, and only the spline parameters (theta) have changed.
 *
 *    The only components of the curve that are changed are SplineY and SplineYcum; The function
 *    ReweightCurve is called at the end to update Y and Ycum.
 *  INPUTS
 *    theCurve      CCurve structure whose SplinePar and SplineEPar have changed
 *    i             index of observation to update, or -1 for all
 *  OUTPUTS
 *    theCurve      CCurve with Y and Ycum updated
 *  SYNOPSIS
 */
void EvalSpline(curveP theCurve, int i)
/*
 *  SOURCE
*/
{
    if(!theCurve->hasSpline) return;
    if(i>=0){
        // evaluate a single observation
        int c1=1;
        theCurve->SplineY[i] = F77_CALL(ddot)(&(theCurve->nj), theCurve->SplineBasis + i,
                &(theCurve->nx), theCurve->SplineEPar, &c1);
        // normalize frailty density
        if(!theCurve->isHazard) theCurve->SplineY[i] /= theCurve->SplineEParSum;
        // compute cumulative hazard
        if(theCurve->isHazard)
            theCurve->SplineYcum[i] = F77_CALL(ddot)(&(theCurve->nj), theCurve->SplineBasisCum + i,
                &(theCurve->nx), theCurve->SplineEPar, &c1);
    }else{
        // evaluate entire curve
        int c1 = 1;
        double c0 = 0;
        char trans = 'N';
        double scaler = theCurve->isHazard ? 1.0 : 1.0/theCurve->SplineEParSum;
        // set splineY = basis * splinepar + 0*splineY
        F77_CALL(dgemv)(&trans, &(theCurve->nx), &(theCurve->nj), &scaler, theCurve->SplineBasis,
                &(theCurve->nx), theCurve->SplineEPar, &c1, &c0, theCurve->SplineY, &c1);
        if(theCurve->isHazard)
            F77_CALL(dgemv)(&trans, &(theCurve->nx), &(theCurve->nj), &scaler, theCurve->SplineBasisCum,
                    &(theCurve->nx), theCurve->SplineEPar, &c1, &c0, theCurve->SplineYcum, &c1);
    }
    ReweightCurve(theCurve, i);
}
/************ EvalSpline */

/****f* CcurveUpdate/EvalSplineAtOnePoint
 *  NAME
 *    EvalSplineAtOnePoint --- evaluate spline at a given value of x
 *  FUNCTION
 *    Evaluate a spline curve at a given point, even if that point is not stored in its
 *    X component and has not been included in its basis. This differs from EvalSpline in that
 *    the latter must already have a fully computed basis stored. This one does a full curve
 *    evaluation -- in that sense it acts as a wrapper to csplineeval.
 *  INPUTS
 *    theCurve      a CCurve structure
 *  OUTPUTS
 *    x             value at which to evaluate the spline
 *  SYNOPSIS
 */
double EvalSplineAtOnePoint(curveP theCurve, double x)
/*
 *  SOURCE
*/
{
    double splY=0;
    int c1=1;
    if(theCurve->hasSpline){
        // compute the basis at x
        double * tempBasis = (double *) malloc( theCurve->nj * sizeof(double));
        for(int j=0; j<theCurve->nj; j++) tempBasis[j]=csplineeval(x, j, theCurve->SplineOrd,
                theCurve->SplineKnots, theCurve->SplineOrd, theCurve->nj);
        // normalize the frailty basis
        if(!theCurve->isHazard)
            for(int j=0; j<theCurve->nj; j++) tempBasis[j] /= (theCurve->SplineBasisInt[j] *
                    theCurve->SplineEParSum); 
        // compute the spline value at X
        splY = F77_CALL(ddot)( &(theCurve->nj), tempBasis, &c1, theCurve->SplineEPar, &c1);
        free(tempBasis);
    }
    return splY;
}
/************ EvalSplineAtOnePoint */

/****f* CcurveUpdate/UpdateSplinePar
 *  NAME
 *    UpdateSplinePar --- update the spline parameters
 *  FUNCTION
 *    Function called to update the spline parameters in a CCurve. It copies a set
 *    of new parameters into the curve and re-evaluates the curve appropriately.
 *  INPUTS
 *    theCurve      a CCurve structure whose spline parameters should be updated
 *    newpar        new spline parameters
 *    j             index of single parameter to be updated (or -1 for all)
 *  OUTPUTS
 *    theCurve      CCurve with updated Y and Ycum
 *  SYNOPSIS
 */
void UpdateSplinePar(curveP theCurve, double * newpar, int j)
/*
 *  SOURCE
*/
{
    if(!theCurve->hasSpline) return;
    if(j>=0){
        // updateSplinePar can only be called with j>0 for hazard curve
        if(!theCurve->isHazard) Rprintf("Bad call to UpdateSplinePar\n");
        double oldeparj = theCurve->SplineEPar[j];
        double neweparj = exp(newpar[j]);
        double epardiff = neweparj - oldeparj;
        theCurve->SplinePar[j] = newpar[j];
        theCurve->SplineEPar[j] = neweparj;
        theCurve->SplineEParSum += epardiff;
        int c1 = 1;
        // update the curve by changing only the j-th parameter
        F77_CALL(daxpy)(&(theCurve->nx), &epardiff, theCurve->SplineBasis+j*theCurve->nx, &c1,
                theCurve->SplineY, &c1); 
        if(theCurve->isHazard) F77_CALL(daxpy)(&(theCurve->nx), &epardiff,
                theCurve->SplineBasisCum+j*theCurve->nx, &c1, theCurve->SplineYcum, &c1); 
        ReweightCurve(theCurve, -1);
    }else{
        // copy the new parameters into the curve structure
        dcopyWrapper(theCurve->nj, newpar, theCurve->SplinePar);
        // make sure that for the frailty, the parameter at the fixed index = 0
        if(!theCurve->isHazard) 
            for(int i=0;i<theCurve->nj; i++) theCurve->SplinePar[i] -= newpar[theCurve->SplineFixedInd];
        // compute exponentials of parameters
        for(int i=0; i<theCurve->nj; i++) theCurve->SplineEPar[i] = exp(theCurve->SplinePar[i]);
        // normalize parameters for frailty
        if(!theCurve->isHazard){
            theCurve->SplineEParSum = 0;
            for(int j=0; j<theCurve->nj; j++) theCurve->SplineEParSum+=theCurve->SplineEPar[j];
            theCurve->SplineFvar = FrailtySplineVar(theCurve);
        }
        // evaluate the spline at the new parameters
        EvalSpline(theCurve, -1);
    }
}
/************ UpdateSplinePar */

/****f* CcurveUpdate/EvalParamAtOnePoint
 *  NAME
 *    EvalParamAtOnePoint --- evaluate parametric component at a single point
 *  FUNCTION
 *    Called as part of EvalCurveAtOnePoint, evaluates the parametric component of
 *    a CCurve given a point x.
 *  INPUTS
 *    theCurve      CCurve structure
 *    x             double, point at which the curve should be evaluated
 *    cum           integer, if >0, the cumulative integral of the curve is returned
 *  OUTPUTS
 *    y             the value of the parametric component at point x
 *  SYNOPSIS
 */
double EvalParamAtOnePoint(curveP theCurve, double x, int cum)
/*
 *  SOURCE
*/
{
    double parY = 0;
    double parYcum =0;
    if(!theCurve->hasPar) return parY;
    // parametric hazard types
    if(theCurve->isHazard){
        // exponential hazard
        if(theCurve->ParDist == Dexponential){
            double lambda = exp(theCurve->ParamPar[0]);
            parY = lambda;
            parYcum = lambda*x;
        // weibull hazard
        }else if(theCurve->ParDist == Dweibull){
            double lambda = exp(theCurve->ParamPar[0]);
            double alpha = exp(theCurve->ParamPar[1]);
            parY = alpha*lambda*pow(x,alpha-1);
            parYcum = lambda * pow(x,alpha);
        }else{
            Rprintf("distribution not implemented");
        }
    // parametric frailty distributions
    }else{
        // gamma distribution
        if(theCurve->ParDist == Dgamma){
            double alpha = exp(- theCurve->ParamPar[0]);
            parY = dgamma(x, alpha, 1/alpha,0);
        // lognormal distribution
        }else if(theCurve->ParDist == Dlognormal){
            double alpha = exp(theCurve->ParamPar[0]);
            parY = exp(-pow(log(x)+alpha/2,2)/(2*alpha)) / (x*sqrt(2*M_PI*alpha));
        }
    }
    return (cum == 0) ? parY : parYcum;
}
/************ EvalParamAtOnePoint */

/****f* CcurveUpdate/EvalParametric
 *  NAME
 *    EvalParametric --- evaluate the parametric component of a curve with new parameters
 *  FUNCTION
 *    Re-evaluate the parametric component of a curve if the parametric component parameters
 *    have changed.
 *  INPUTS
 *    theCurve      CCurve structure
 *    i             index of the observation that should be evaluated
 *  OUTPUTS
 *    theCurve      CCurve, with Y and Ycum updated
 *  SYNOPSIS
 */
void EvalParametric(curveP theCurve, int i)
/*
 *  SOURCE
*/
{
    if(!theCurve->hasPar) return;
    if(i<0){
        // if all observations should be updated, updated them one by one
        for(int ind=0; ind<theCurve->nx; ind++) EvalParametric(theCurve, ind);
    }else{
        theCurve->ParamY[i] = EvalParamAtOnePoint(theCurve, theCurve->X[i], 0);
        if(theCurve->isHazard)
            theCurve->ParamYcum[i] = EvalParamAtOnePoint(theCurve, theCurve->X[i], 1);
    }
    ReweightCurve(theCurve, i);
}
/************ EvalParametric */

/****f* CcurveUpdate/UpdateParamPar
 *  NAME
 *    UpdateParamPar --- update parametric component with new parameters
 *  FUNCTION
 *    Update the parametric component of a curve, given a new set of parameters. This just
 *    involves copying the new parameters into the CCurve and calling EvalParametric
 *    to re-evaluate the curve.
 *  INPUTS
 *    theCurve  CCurve structure
 *    newPar    new parametric component parameters
 *  OUTPUTS
 *    theCurve  CCurve structure, with newPar copied into ParamPar and Y and Ycum updated.
 *  SYNOPSIS
 */
void UpdateParamPar(curveP theCurve, double * newPar)
/*
 *  SOURCE
*/
{
    if(!theCurve->hasPar) return;
    // copy new parameters into the curve
    for(int i=0; i<theCurve->np; i++)
        theCurve->ParamPar[i] = newPar[i];
    // evaluate the curve
    EvalParametric(theCurve, -1);
}
/************ UpdateParamPar */

/****f* CcurveUpdate/EvalCurveAtOnePoint
 *  NAME
 *    EvalCurveAtOnePoint --- evaluate a CCurve at a single point
 *  FUNCTION
 *    Given a point x, compute the value of the curve at that point, even if it is
 *    not included in the curve's set of observations.
 *  INPUTS
 *    theCurve      CCurve to be evaluated
 *    x             double, point at which to evaluate it.
 *  OUTPUTS
 *    y             CCurve evaluated at x
 *  SYNOPSIS
 */
double EvalCurveAtOnePoint(curveP theCurve, double x)
/*
 *  SOURCE
*/
{
    // this would actually work on frailties, I think, but it never gets called
    // for frailties, so I just put this in as a safeguard.
    if(theCurve->isHazard) Rprintf("Error: using EvalCurveAtOnePoint on hazard");
    // evaluate spline component
    double splY=EvalSplineAtOnePoint(theCurve,x);
    // evaluate parametric component
    double parY=EvalParamAtOnePoint(theCurve,x,0);
    // weight the two components
    if(theCurve->hasPar & theCurve->hasSpline)
        return theCurve->Weight[0] * splY + (1-theCurve->Weight[0])*parY;
    else if(theCurve->hasPar)
        return parY;
    else if(theCurve->hasSpline)
        return splY;
    else return 0;
}
/************ EvalCurveAtOnePoint */

/****f* CcurveUpdate/UpdateCurveX
 *  NAME
 *    UpdateCurveX --- change an X value of a CCurve structure
 *  FUNCTION
 *    Changes theCurve->X[i] to x, and updates Y[i] accordingly
 *  INPUTS
 *    theCurve      a CCurve structure
 *    x             double, to replace theCurve->X[i]
 *    i             index between 0 and theCurve->nx
 *  SYNOPSIS
 */
void UpdateCurveX(curveP theCurve, double x, int i)
/*
 *  SOURCE
*/
{
    if(theCurve->isHazard) Rprintf("Error: using UpdateCurveX on hazard.");
    // copy x into theCurve
    theCurve->X[i] = x;
    // re-evaluate the curve basis and values
    UpdateSplineBasis(theCurve, i, 0, theCurve->nj);
    EvalSpline(theCurve,i);
    EvalParametric(theCurve,i);
}
/************ UpdateCurveX */

/****f* CmiscUtils/EvalNknotsPrior
 *  NAME
 *    EvalNknotsPrior --- evaluate the prior on the number of knots
 *  FUNCTION
 *    Compute the prior probability of a given number of knots for a curve. Used
 *    in MH_BDM.
 *  INPUTS
 *     nknots       the value at which the prior should be evaluated
 *     theCurve     a CCurve structure with SplineNknotsPrior and SplineNknotsHyper components
 *  OUTPUTS
 *     p        prior probability of nknots for a distribution given by SplineNknotsPrior
 *  SYNOPSIS
 */
static inline double EvalNknotsPrior( int nknots, curveP theCurve)
/*
 *  SOURCE
*/
{
    double dnknots = (double) nknots;
    // fetch the prior probability distribution
    nknotsprior prior = theCurve->SplineNknotsPrior;
    // evaluate different types of distribution
    if(prior==PrPoisson)  // Poisson prior
        return dpois( dnknots, theCurve->SplineNknotsHyper[0], 0);
    if(prior==PrGeometric)  // Geometric prior
        return dgeom( dnknots, theCurve->SplineNknotsHyper[0], 0);
    if(prior==PrPoissonMix)  // Poisson mixture
        return 0.5*( dpois( dnknots, theCurve->SplineNknotsHyper[0], 0) +
                dpois( dnknots, theCurve->SplineNknotsHyper[0], 0));
    if(prior==PrNegBin)  // Negative binomial
        return dnegbin( dnknots, theCurve->SplineNknotsHyper[0], theCurve->SplineNknotsHyper[1]);
    if(prior==PrPower)  // Power
        return pow( dnknots, theCurve->SplineNknotsHyper[0]);
}
/************ EvalNknotsPrior */

/****f* CmakeLikelihood/LikelihoodFrailty
 *  NAME
 *    LikelihoodFrailty --- likelihood of the frailty itself
 *  FUNCTION
 *    Compute frailty likelihoods as part of MH_Frail. See also mklik.frail
 *  INPUTS
 *    i             index of the frailty whose likelihood should be evaluated
 *    hazard        CCurve for the hazard
 *    frailty       CCurve for the frailty
 *    regression    CRegression structure
 *  OUTPUTS
 *    lik   loglikelihood of frailty->X[i]
 *  SYNOPSIS
 */
static inline double LikelihoodFrailty(int i, curveP hazard, curveP frailty, regressionP regression)
/*
 *  SOURCE
*/
{
   double Ui = frailty->X[i];
   double lUi = log(Ui);
   double lik = 0;
   int c1 = 1;
   lik += log(frailty->Y[i]);
   int start = regression->Jicum[i];
   int n = regression->Jicum[i+1] - start;
   lik += lUi * F77_CALL(dasum)(&n, regression->status+start, &c1);
   lik -= Ui * ddotWrapper(n, hazard->Ycum+start, regression->elp+start);
   return(lik);
}
/************ LikelihoodFrailty */

/****f* CmakeLikelihood/LikelihoodHazardLogSum
 *  NAME
 *    LikelihoodHazardLogSum --- utility to efficiently compute a sum of log-hazards
 *  FUNCTION
 *    Because the log function is quite slow (especially it seems on the MacIntel architecture),
 *    this function computes the sum of log-hazards by multiplying several together at a time,
 *    then taking logs.
 *  INPUTS
 *    nx        number of observations
 *    Y         vector of length nx
 *    status    vector of length nx
 *  OUTPUTS
 *    lik       sum of status*log(Y)
 *  SYNOPSIS
 */
static inline double LikelihoodHazardLogSum(int nx, double *status, double *Y)
/*
 *  SOURCE
*/
{
    double lik=0;
    double thislik = 1;
    // sum LIK_MOD elements at a time
    for(int i=0; i<nx; i++) {
        thislik *= status[i] > 0 ? Y[i] : 1.0;
        if(i % LIK_MOD == 0){
            lik += log(thislik);
            thislik = 1;
        }
    }  
    lik+=log(thislik);
    return lik;  
}
/************ LikelihoodHazardLogSum */

/****f* CmakeLikelihood/LikelihoodFrailtyLogSum
 *  NAME
 *    LikelihoodFrailtyLogSum --- utility to efficiently compute a sum of log-frailties
 *  FUNCTION
 *    Because the log function is quite slow (especially it seems on the MacIntel architecture),
 *    this function computes the sum of log-frailties by multiplying several together at a time,
 *    then taking logs.
 *  INPUTS
 *    nx        number of observations
 *    Y         vector of length nx
 *  OUTPUTS
 *    lik       sum of log(Y)
 *  SYNOPSIS
 */
static inline double LikelihoodFrailtyLogSum(int nx, double *Y)
/*
 *  SOURCE
*/
{
    double lik=0;
    double thislik = 1;

    for(int i=0; i<nx; i++) {
        thislik *= Y[i];
        if(i % LIK_MOD == 0){
            lik += log(thislik);
            thislik = 1;
        }
    }
    lik += log(thislik);
    return lik;
}
/************ LikelihoodFrailtyLogSum */

/****f* CmakeLikelihood/LikelihoodRegression
 *  NAME
 *    LikelihoodRegression --- likelihood of regression coefficients
 *  FUNCTION
 *    Computes the likelihood of regression coefficients beta, called by MH_Regression. See
 *    also mklik.coef.
 *  INPUTS
 *    hazard        CCurve for the hazard
 *    frailty       CCurve for the frailty
 *    regression    CRegression structure
 *  OUTPUTS
 *    lik    loglikelihood of regression->coefficients
 *  SYNOPSIS
 */
static inline double LikelihoodRegression(curveP hazard,curveP frailty,regressionP regression)
/*
 *  SOURCE
*/
{
    double lik; 
    double * hazYcum = hazard->Ycum;
    double * frailelp = regression->frailelp;
    lik = ddotWrapper(regression->n, regression->status, regression->lp);
    lik -= ddotWrapper(regression->n, frailelp, hazYcum);
    int c1 = 1;
    lik -= pow(F77_CALL(dnrm2)(&(regression->p), regression->coefficients, &c1),2)/
        (2*regression->priorvar[0]);
    return lik;
}
/************ LikelihoodRegression */

/****f* CmakeLikelihood/LikelihoodSplineHazard
 *  NAME
 *    LikelihoodSplineHazard --- likelihood of hazard spline parameters
 *  FUNCTION
 *    Compute loglikelihood of parameters for the spline component of the hazard curve.
 *    See also mklik.spline.haz.
 *  INPUTS
 *    hazard        CCurve for the hazard
 *    frailty       CCurve for the frailty
 *    regression    CRegression structure
 *  OUTPUTS
 *    lik    loglikelihood of hazard->SplinePar
 *  SYNOPSIS
 */
static inline double LikelihoodSplineHazard(curveP hazard, curveP frailty, regressionP regression)
/*
 *  SOURCE
*/
{
    double * frailelp = regression->frailelp;
    double * hazYcum = hazard->Ycum;
    
    // point process likelihood
    double lik = LikelihoodHazardLogSum(hazard->nx, regression->status, hazard->SplineY);
    lik -= ddotWrapper(regression->n, frailelp, hazard->SplineYcum);
    // smoothness penalty
    lik -= hazard->SplinePenaltyFactor[0]*SmoothnessPenalty(hazard); 
    // penalize parameters that are too small or too large
    for(int i=0; i<hazard->nj; i++) lik -= hazard->SplinePar[i]<hazard->SplineMin[0] ?
        pow(hazard->SplinePar[i] - hazard->SplineMin[0],2) : 0.0;
    for(int i=0; i<hazard->nj; i++) lik += hazard->SplinePar[i] > MAX_PAR ? -INFINITY: 0.0; return lik;
}
/************ LikelihoodSplineHazard */

/****f* CmakeLikelihood/LikelihoodSplineFrailty
 *  NAME
 *    LikelihoodSplineFrailty --- likelihood of frailty spline parameters
 *  FUNCTION
 *    Compute loglikelihood of parameters for the spline component of the frailty curve.
 *    See also mklik.spline.frail.
 *  INPUTS
 *    hazard        CCurve for the hazard
 *    frailty       CCurve for the frailty
 *    regression    CRegression structure
 *  OUTPUTS
 *    lik    loglikelihood of frailty->SplinePar
 *  SYNOPSIS
 */
static inline double LikelihoodSplineFrailty(curveP hazard, curveP frailty, regressionP regression)
/*
 *  SOURCE
*/
{
    // sum of log-hazards
    double lik =  LikelihoodFrailtyLogSum(frailty->nx, frailty->SplineY);
    // smoothness penalty
    lik -= frailty->SplinePenaltyFactor[0]*SmoothnessPenalty(frailty); 
    // penalize parameters that are too small or too big
    for(int i=0; i<frailty->nj; i++) lik -= frailty->SplinePar[i]<frailty->SplineMin[0] ?
        pow(frailty->SplinePar[i] - frailty->SplineMin[0],2) : 0.0;
    for(int i=0; i<frailty->nj; i++) lik += frailty->SplinePar[i] > MAX_PAR ? -INFINITY: 0.0; 
    return lik;
}
/************ LikelihoodSplineFrailty */

/****f* CmakeLikelihood/LikelihoodParamHazard
 *  NAME
 *    LikelihoodParamHazard --- likelihood of hazard parametric parameters
 *  FUNCTION
 *    Compute loglikelihood of parameters for the parametric component of the hazard curve.
 *    See also mklik.param.haz.
 *  INPUTS
 *    hazard        CCurve for the hazard
 *    frailty       CCurve for the frailty
 *    regression    CRegression structure
 *  OUTPUTS
 *    lik    loglikelihood of hazard->ParamPar
 *  SYNOPSIS
 */
static inline double LikelihoodParamHazard(curveP hazard, curveP frailty, regressionP regression)
/*
 *  SOURCE
*/
{
    double * frailelp = regression->frailelp;
    double * hazYcum = hazard->Ycum;
   
    // point process likelihood
    double lik = LikelihoodHazardLogSum(hazard->nx, regression->status, hazard->ParamY);
    lik -= ddotWrapper(regression->n, frailelp, hazard->ParamYcum);
    int c1=1;
    // Gaussian prior for parameters
    lik -= pow(F77_CALL(dnrm2)(&(hazard->np), hazard->ParamPar, &c1),2)/(2*hazard->ParamPriorvar[0]);
    return lik;
}
/************ LikelihoodParamHazard */

/****f* CmakeLikelihood/LikelihoodParamFrailty
 *  NAME
 *    LikelihoodParamFrailty --- likelihood of frailty parametric parameter
 *  FUNCTION
 *    Compute loglikelihood of parameters for the parametric component of the frailty curve.
 *    See also mklik.param.frail.
 *  INPUTS
 *    hazard        CCurve for the hazard
 *    frailty       CCurve for the frailty
 *    regression    CRegression structure
 *  OUTPUTS
 *    lik    loglikelihood of frailty->ParamPar
 *  SYNOPSIS
 */
static inline double LikelihoodParamFrailty(curveP hazard, curveP frailty, regressionP regression)
/*
 *  SOURCE
*/
{
    // sum of log-frailties
    double lik = LikelihoodFrailtyLogSum(frailty->nx, frailty->ParamY);

    int c1=1;
    // Gaussian prior for parameters
    lik -= pow(F77_CALL(dnrm2)(&(frailty->np), frailty->ParamPar, &c1),2)/(2*frailty->ParamPriorvar[0]);
    return lik;
}
/************ LikelihoodParamFrailty */

/****f* CmakeLikelihood/LikelihoodWeightHazard
 *  NAME
 *    LikelihoodWeightHazard --- likelihood of weight for hazard curve
 *  FUNCTION
 *    Compute loglikelihood of weight parameter for the hazard curve, if the hazard curve
 *    has both. See also mklik.weight.haz.
 *  INPUTS
 *    hazard        CCurve for the hazard
 *    frailty       CCurve for the frailty
 *    regression    CRegression structure
 *  OUTPUTS
 *    lik   loglikelihood of hazard->Weight
 *  SYNOPSIS
 */
double LikelihoodWeightHazard(curveP hazard, curveP frailty, regressionP regression)
/*
 *  SOURCE
*/
{ 
    double * frailelp = regression->frailelp;
    double * hazYcum = hazard->Ycum;
    
    // point process likelihood
    double lik = LikelihoodHazardLogSum(hazard->nx, regression->status, hazard->Y);
    lik -= ddotWrapper(regression->n, frailelp, hazYcum);
    // prior on the weight
    lik += (hazard->WeightHyper[0] - 1.0) * log(hazard->Weight[0])
          +(hazard->WeightHyper[1] - 1.0) * log(1.0 - hazard->Weight[0]);
    return lik;
}
/************ LikelihoodWeightHazard */

/****f* CmakeLikelihood/LikelihoodWeightFrailty
 *  NAME
 *    LikelihoodWeightFrailty --- likelihood of weight for frailty curve
 *  FUNCTION
 *    Compute loglikelihood of weight parameter for the frailty curve, if the frailty curve
 *    has both. See also mklik.weight.frail.
 *  INPUTS
 *    hazard        CCurve for the hazard
 *    frailty       CCurve for the frailty
 *    regression    CRegression structure
 *  OUTPUTS
 *    lik   loglikelihood of frailty->Weight
 *  SYNOPSIS
 */
double LikelihoodWeightFrailty(curveP hazard, curveP frailty, regressionP regression)
/*
 *  SOURCE
*/
{
    double lik =  LikelihoodFrailtyLogSum(frailty->nx, frailty->Y);

    lik += (frailty->WeightHyper[0] - 1.0) * log(frailty->Weight[0])
          +(frailty->WeightHyper[1] - 1.0) * log(1.0 - frailty->Weight[0]);
    return lik;
}
/************ LikelihoodWeightFrailty */

/****f* CmakeLikelihood/LikelihoodFull
 *  NAME
 *    LikelihoodFull --- Full posterior likelihood
 *  FUNCTION
 *    Compute the full posterior likelihood of parameters at the current iteration. 
 *    Needed for DIC.    
 *  INPUTS
 *    hazard        CCurve for the hazard
 *    frailty       CCurve for the frailty
 *    regression    CRegression structure
 *  OUTPUTS
 *    lik   full loglikelihood
 *  SYNOPSIS
 */
double LikelihoodFull(curveP hazard, curveP frailty, regressionP regression)
/*
 *  SOURCE
*/
{
    double * hazYcum = hazard->Ycum;
    double * frailelp = regression->frailelp;
    double lik =  0.0;
    int c1 = 1;
    // Hazard, frailty and regression (top-level)
    lik += LikelihoodHazardLogSum(regression->n, regression->status, regression->frailrep);
    lik += LikelihoodHazardLogSum(hazard->nx, regression->status, hazard->Y);
    lik += ddotWrapper(regression->n, regression->status, regression->lp);

    lik -= ddotWrapper(regression->n, regression->frailelp, hazard->Ycum);
    lik += LikelihoodFrailtyLogSum(frailty->nx, frailty->Y);

    // priors on regression and spline parameters and their variances
    lik -= (regression->p)/2 * log(regression->priorvar[0]);
    lik -= pow(F77_CALL(dnrm2)(&(regression->p), regression->coefficients, &c1),2)/
        (2.0*regression->priorvar[0]);
    lik -= (regression->hyper[0]+1) * log(regression->priorvar[0])
           + (regression->hyper[1])/(regression->priorvar[0]);
    if(hazard->hasSpline){
        lik -= (hazard->nj)/2 * log(hazard->SplinePriorvar[0])
            + hazard->SplinePenaltyFactor[0]*SmoothnessPenalty(hazard); 
        lik -= (hazard->SplineHyper[0]+1) * log(hazard->SplinePriorvar[0])
            + hazard->SplineHyper[1]/hazard->SplinePriorvar[0];
            }
    if(frailty->hasSpline){
        lik -= (frailty->nj)/2 * log(frailty->SplinePriorvar[0])
            + frailty->SplinePenaltyFactor[0]*SmoothnessPenalty(frailty); 
        lik -= (frailty->SplineHyper[0]+1) * log(frailty->SplinePriorvar[0])
            + frailty->SplineHyper[1]/frailty->SplinePriorvar[0];
    }
    // Priors on parametric components and their variances
    if(hazard->hasPar){
        lik -= (hazard->np)/2 * log(hazard->ParamPriorvar[0]);
        lik -= pow(F77_CALL(dnrm2)(&(hazard->np), hazard->ParamPar, &c1),2)/(2*hazard->ParamPriorvar[0]);
        lik -= (hazard->ParamHyper[0]+1) * log(hazard->ParamPriorvar[0])
            + hazard->ParamHyper[1]/hazard->ParamPriorvar[0];
    }
    if(frailty->hasPar){
        lik -= (frailty->np)/2 * log(frailty->ParamPriorvar[0]);
        lik -= pow(F77_CALL(dnrm2)(&(frailty->np), frailty->ParamPar, &c1),2)/(2*frailty->ParamPriorvar[0]);
        lik -= (frailty->ParamHyper[0]+1) * log(frailty->ParamPriorvar[0])
            + frailty->ParamHyper[1]/frailty->ParamPriorvar[0];
    }

    // Prior on weights
    if(hazard->hasSpline & hazard->hasPar)
        lik += (hazard->WeightHyper[0] - 1.0) * log(hazard->Weight[0])
              +(hazard->WeightHyper[1] - 1.0) * log(1.0 - hazard->Weight[0]);
    if(frailty->hasSpline & frailty->hasPar)
        lik += (frailty->WeightHyper[0] - 1.0) * log(frailty->Weight[0])
              +(frailty->WeightHyper[1] - 1.0) * log(1.0 - frailty->Weight[0]);

    // Prior on number and positions of knots
    if(hazard->hasSpline & hazard->SplineAdaptive){
        lik += log(dfactorial((double) hazard->SplineNknotsMax - hazard->SplineNknots)
                * dfactorial((double) hazard->SplineNknots)
                / dfactorial((double) hazard->SplineNknotsMax));
        lik += log(EvalNknotsPrior(hazard->SplineNknots, hazard));
    }
    if(frailty->hasSpline & frailty->SplineAdaptive){
        lik += log(dfactorial((double) frailty->SplineNknotsMax - frailty->SplineNknots)
                * dfactorial((double) frailty->SplineNknots)
                / dfactorial((double) frailty->SplineNknotsMax));
        lik += log(EvalNknotsPrior(frailty->SplineNknots, frailty));
    }

    return lik;
}
/************ LikelihoodFull */

/****f* CMetropolisHastings/AcceptReject
 *  NAME
 *    AcceptReject --- accept-reject step for Metropolis-Hastings
 *  FUNCTION
 *    Executes the MH accept-reject decision.
 *  INPUTS
 *    baselik   loglikelihood of the base case
 *    candlik   loglikelihood of the candidate
 *    ratio     additional multiplier (e.g. prior ratio)
 *  OUTPUTS
 *    out       integer, if 1, accept, if 0, reject
 *  SYNOPSIS
 */
static inline int AcceptReject(double baselik, double candlik, double ratio)
/*
 *  SOURCE
*/
{
    if(isnan(candlik)) candlik = -DBL_MAX;
    double r = exp(candlik - baselik) * ratio; 
    double pAccept = dmin(1,r);
    int out;
    if(isnan(pAccept)) pAccept = 0;
    if(runif(0,1) < pAccept) out = 1;
    else out = 0;
    return out;
}
/************ AcceptReject */

/****f* CMetropolisHastings/MH_Frail
 *  NAME
 *    MH_Frail --- MH step for frailties
 *  FUNCTION
 *    Update frailty estimates by Metropolis-Hastings. See also mh.frail.
 *  INPUTS
 *    hazard        CCurve for the hazard
 *    frailty       CCurve for the frailty
 *    regression    CRegression structure
 *  SYNOPSIS
 */
void MH_Frail(curveP hazard, curveP frailty, regressionP regression)
/*
 *  SOURCE
*/
{
    double acc = 0;
    double baselik, candlik;
    int j;
    // update the frailties one at a time
    for(int i=0; i < regression->m; i++){
        double u = frailty->X[i];
        double y = frailty->Y[i];
        double v = frailty->tun[0];
        double uj,yj,candj;
        // generate candidate from gamma distribution
        double cand = rgamma( pow(u,2)/v, v/u);
        double pcu = 1; double puc = 1;
        // if the candidate is invalid (e.g. beyond boundary knots), fail
        if(isnan(cand) || cand<1e-5 || ( frailty->hasSpline &&
            (cand > frailty->SplineKnots[frailty->nj + frailty->SplineOrd -1]
             | cand < frailty->SplineKnots[0]) ) ){
            continue;
        }else{
            // find another frailty to move by the same amount
            j = i;
            while(j == i){
                j = (int) floor(runif(0,(double) frailty->nx));
            }
            //Rprintf("i: %d,  j: %d\n",i,j);
            uj = frailty->X[j];
            yj = frailty->Y[j];
            candj = u + uj - cand;
            // if the moved element is invalid, fail
            if( candj<1e-5 || (frailty->hasSpline &&
                (candj > frailty->SplineKnots[frailty->nj + frailty->SplineOrd -1]
                 | candj < frailty->SplineKnots[0]) )) continue;
            // base likelihood
            baselik = LikelihoodFrailty(i, hazard, frailty, regression)
                + LikelihoodFrailty(j, hazard, frailty, regression);
            // update the curve with the candidate, without touching the basis
            frailty->X[i] = cand;
            frailty->Y[i] = EvalCurveAtOnePoint(frailty, cand);
            frailty->X[j] = candj;
            frailty->Y[j] = EvalCurveAtOnePoint(frailty, candj);
            // candidate likelihood
            candlik = LikelihoodFrailty(i, hazard, frailty, regression)
                + LikelihoodFrailty(j, hazard, frailty, regression);
            // transition probabilities
            puc = dgamma(cand, pow(u,2)/v, v/u, 0);
            pcu = dgamma(u, pow(cand,2)/v, v/cand, 0);
        }
        int thisacc = AcceptReject(baselik, candlik, pcu/puc);
        if(thisacc==0) { //Did not accept, so undo the damage
            frailty->X[i] = u;
            frailty->Y[i] = y;
            frailty->X[j] = uj;
            frailty->Y[j] = yj;
        }else{ // accepted, so update the curve fully
            UpdateCurveX(frailty, cand, i);
            UpdateCurveX(frailty, candj, j);
            for(int k=regression->Jicum[i]; k<regression->Jicum[i+1]; k++){
                regression->frailrep[k] = frailty->X[i];
                regression->frailelp[k] = frailty->X[i]*regression->elp[k];
            }
            for(int k=regression->Jicum[j]; k<regression->Jicum[j+1]; k++){
                regression->frailrep[k] = frailty->X[j];
                regression->frailelp[k] = frailty->X[j]*regression->elp[k];
            }
        }
        // update acceptance rate
        acc += (double) thisacc;
    }
    acc /= regression->m;
    frailty->Accept[0] = acc;
}
/************ MH_Frail */

/****f* CMetropolisHastings/MH_Regression
 *  NAME
 *    MH_Regression --- MH step for regression coefficients
 *  FUNCTION
 *    Update the regression coefficients by Metropolis-Hastings. See also mh.coef.
 *  INPUTS
 *    hazard        CCurve for the hazard
 *    frailty       CCurve for the frailty
 *    regression    CRegression structure
 *  SYNOPSIS
 */
void MH_Regression(curveP hazard, curveP frailty, regressionP regression)
/*
 *  SOURCE
*/
{
    double baselik, candlik;
    // base likelihood
    baselik = LikelihoodRegression(hazard, frailty, regression);
    double * cand = (double *) calloc( regression->p, sizeof(double));
    double * oldlp = (double *) malloc( regression->n * sizeof(double));
    double * oldelp = (double *) malloc( regression->n * sizeof(double));
    double * oldfrailelp = (double *) malloc( regression->n * sizeof(double));
    double * oldcoef = (double *) malloc( regression->p * sizeof(double));
    // store the old regression information
    char trans = 'N'; double c0 = 0; int c1 = 1; double c1d = 1;
    dcopyWrapper(regression->p, regression->coefficients, oldcoef);
    dcopyWrapper(regression->n, regression->lp, oldlp);
    dcopyWrapper(regression->n, regression->elp, oldelp);
    dcopyWrapper(regression->n, regression->frailelp, oldfrailelp);
    //generate candidate parameters
    mvrnorm(regression->p, cand, regression->coefficients, regression->CholCov, regression->tun[0]); 
    // Change the regression object with the new lp and elps
    dcopyWrapper(regression->p, cand, regression->coefficients);
    F77_CALL(dgemv)(&trans, &(regression->n), &(regression->p), &c1d, regression->covariates,
            &(regression->n), regression->coefficients, &c1, &c0, regression->lp, &c1);
    for(int i=0; i < regression->n; i++) regression->elp[i] = exp(regression->lp[i]);
    diagmvWrapper(regression->n, regression->frailrep, regression->elp, regression->frailelp);
    candlik = LikelihoodRegression(hazard, frailty, regression);
    int acc = AcceptReject(baselik, candlik, 1);
    if(!acc){
        // not accepted, so restore old regression
        dcopyWrapper(regression->p, oldcoef, regression->coefficients);
        dcopyWrapper(regression->n, oldlp, regression->lp);
        dcopyWrapper(regression->n, oldelp, regression->elp);
        dcopyWrapper(regression->n, oldfrailelp, regression->frailelp);
    }
    regression->Accept[0] = (double) acc; 
    free(cand);
    free(oldlp);
    free(oldelp);
    free(oldfrailelp);
    free(oldcoef);
}
/************ MH_Regression */

/****f* CMetropolisHastings/MH_SplineHazard
 *  NAME
 *    MH_SplineHazard --- MH for hazard spline parameters
 *  FUNCTION
 *    Update the parameters of the hazard spline component by Metropolis-Hastings.
 *    See also mh.hazard.spline.
 *  INPUTS
 *    hazard        CCurve for the hazard
 *    frailty       CCurve for the frailty
 *    regression    CRegression structure
 *  SYNOPSIS
 */
void MH_SplineHazard(curveP hazard, curveP frailty, regressionP regression)
/*
 *  SOURCE
*/
{
    if(!hazard->hasSpline) return;
    double baselik, candlik;
    double sumacc=0;
    // base likelihood
    baselik = LikelihoodSplineHazard(hazard,frailty,regression); 
    double * cand = (double *) calloc( hazard->nj, sizeof(double));
    double * thiscand = (double *) calloc( hazard->nj, sizeof(double));
    // allocate storage for parameters of old spline
    // if the move is not accepted, these will be used to restore it, to make it faster
    double * oldPar = (double *) malloc( hazard->nj * sizeof(double));
    double * oldEPar = (double *) malloc( hazard->nj * sizeof(double));
    double * oldY = (double *) malloc( hazard->nx * sizeof(double));
    double * oldYcum = (double *) malloc( hazard->nx * sizeof(double));
    double * oldSplineY = (double *) malloc( hazard->nx * sizeof(double));
    double * oldSplineYcum = (double *) malloc( hazard->nx * sizeof(double));
    dcopyWrapper(hazard->nj, hazard->SplinePar, oldPar);
    dcopyWrapper(hazard->nj, hazard->SplineEPar, oldEPar);
    dcopyWrapper(hazard->nj, hazard->SplinePar, thiscand);
    dcopyWrapper(hazard->nx, hazard->Y, oldY);
    dcopyWrapper(hazard->nx, hazard->SplineY, oldSplineY);
    dcopyWrapper(hazard->nx, hazard->Ycum, oldYcum);
    dcopyWrapper(hazard->nx, hazard->SplineYcum, oldSplineYcum);
    // create candidate parameters
    for(int j=0;j<hazard->nj;j++)
        cand[j] = hazard->SplinePar[j]+hazard->SplineTun[0]*
            rnorm(0,hazard->SplineCandSD[j]);
    // update spline parameters one at a time
    for(int j=0; j < hazard->nj; j++){
        thiscand[j] = cand[j];
        UpdateSplinePar(hazard,thiscand,j); 
        // candidate likelihod
        candlik = LikelihoodSplineHazard(hazard,frailty,regression);
        int acc = AcceptReject(baselik, candlik, 1);
        if(acc){ // accepted, so store the new information
            baselik = candlik;
            sumacc++;
            dcopyWrapper(hazard->nj, hazard->SplineEPar, oldEPar);
            dcopyWrapper(hazard->nx, hazard->Y, oldY);
            dcopyWrapper(hazard->nx, hazard->SplineY, oldSplineY);
            dcopyWrapper(hazard->nx, hazard->Ycum, oldYcum);
            dcopyWrapper(hazard->nx, hazard->SplineYcum, oldSplineYcum);
        }else{ // rejected, so restore the old information
            thiscand[j]=oldPar[j];
            dcopyWrapper(hazard->nj, thiscand, hazard->SplinePar);
            dcopyWrapper(hazard->nj, oldEPar, hazard->SplineEPar);
            dcopyWrapper(hazard->nx, oldY, hazard->Y);
            dcopyWrapper(hazard->nx, oldSplineY, hazard->SplineY);
            dcopyWrapper(hazard->nx, oldYcum, hazard->Ycum);
            dcopyWrapper(hazard->nx, oldSplineYcum, hazard->SplineYcum);
        }
    }
    // acceptance rate
    hazard->SplineAccept[0] =  sumacc / ((double) hazard->nj);
    // free memory
    free(cand);
    free(thiscand);
    free(oldPar);
    free(oldEPar);
    free(oldY);
    free(oldYcum);
    free(oldSplineY);
    free(oldSplineYcum);
}
/************ MH_SplineHazard */

/****f* CMetropolisHastings/MH_SplineFrailty
 *  NAME
 *    MH_SplineFrailty --- MH for frailty spline parameters
 *  FUNCTION
 *    Update frailty spline parameters by Metropolis-Hastings. See also mh.frailty.spline.
 *  INPUTS
 *    hazard        CCurve for the hazard
 *    frailty       CCurve for the frailty
 *    regression    CRegression structure
 *  SYNOPSIS
 */
void MH_SplineFrailty(curveP hazard, curveP frailty, regressionP regression)
/*
 *  SOURCE
*/
{
    if(!frailty->hasSpline) return;
    double baselik, candlik;
    double sumacc=0;
    // base likelihood
    baselik = LikelihoodSplineFrailty(hazard,frailty,regression); 
    double * cand = (double *) calloc( frailty->nj, sizeof(double));
    // allocate storage for parameters of old spline
    double * oldPar = (double *) malloc( (frailty->nj) * sizeof(double));
    double * oldEPar = (double *) malloc( frailty->nj * sizeof(double));
    double * oldY = (double *) malloc( frailty->nx * sizeof(double));
    double * oldSplineY = (double *) malloc( frailty->nx * sizeof(double));
    dcopyWrapper(frailty->nj, frailty->SplinePar, oldPar);
    dcopyWrapper(frailty->nj, frailty->SplineEPar, oldEPar);
    dcopyWrapper(frailty->nx, frailty->Y, oldY);
    dcopyWrapper(frailty->nx, frailty->SplineY, oldSplineY);
    double oldSplineEParSum = frailty->SplineEParSum;
    double oldSplineFvar = frailty->SplineFvar;
    int ord2 = frailty->SplineOrd / 2;

    // update the frailty spline parameters one by one
    for(int j=0; j<frailty->nj; j++){
        dcopyWrapper(frailty->nj, oldPar, cand);
        // choose which parameter will compensate for j
        // to ensure the frailty mean remains 1
        int k = j;
        while(j == k | k<ord2-1 | k>frailty->nj-ord2-1) 
            k = (int) floor(runif(0,(double) frailty->nj));
        
        // Generate candidate parameter at j
        cand[j] = dmax(frailty->SplinePar[j],frailty->SplineMin[0])+frailty->SplineTun[0]*
            rnorm(0,frailty->SplineCandSD[j]);
        if(cand[j]<frailty->SplineMin[0]) cand[j] = 1e-10;
        // Try to compute value at k to compensate for change at j
        double newmean = frailty->SplineBasisExp[j] * (exp(cand[j])-exp(oldPar[j]));
        double candk = log(oldEPar[k] - newmean/frailty->SplineBasisExp[k]);
        if(isnan(candk)) continue;
        cand[k] = candk;
        // Compute candidate likelihood
        UpdateSplinePar(frailty,cand,-1); 
        candlik = LikelihoodSplineFrailty(hazard,frailty,regression);
        int acc = AcceptReject(baselik, candlik, 1);
        if(acc){ // accepted, save and continue
            baselik = candlik;
            sumacc++;
            dcopyWrapper(frailty->nj, frailty->SplinePar, oldPar);
            dcopyWrapper(frailty->nj, frailty->SplineEPar, oldEPar);
            dcopyWrapper(frailty->nx, frailty->Y, oldY);
            dcopyWrapper(frailty->nx, frailty->SplineY, oldSplineY);
            oldSplineEParSum = frailty->SplineEParSum;
            oldSplineFvar = frailty->SplineFvar;
        }else{ // rejected, restore old values and continue
            dcopyWrapper(frailty->nj, oldPar, frailty->SplinePar);
            dcopyWrapper(frailty->nj, oldEPar, frailty->SplineEPar);
            dcopyWrapper(frailty->nx, oldY, frailty->Y);
            dcopyWrapper(frailty->nx, oldSplineY, frailty->SplineY);
            frailty->SplineEParSum = oldSplineEParSum ;
            frailty->SplineFvar = oldSplineFvar ;
        }
    }
    // acceptance rate
    frailty->SplineAccept[0] = sumacc / ((double) frailty->nj);
    // free memory
    free(cand);
    free(oldPar);
    free(oldEPar);
    free(oldY);
    free(oldSplineY);
}
/************ MH_SplineFrailty */


/****f* CMetropolisHastings/MH_ParamHazard
 *  NAME
 *    MH_ParamHazard --- MH for hazard parametric component parameters
 *  FUNCTION
 *    Update the parameters indexing the parametric component by Metropolis-Hastings,
 *    see also mh.hazard.param
 *  INPUTS
 *    hazard        CCurve for the hazard
 *    frailty       CCurve for the frailty
 *    regression    CRegression structure
 *  SYNOPSIS
 */
void MH_ParamHazard(curveP hazard, curveP frailty, regressionP regression)
/*
 *  SOURCE
*/
{
    if(!hazard->hasPar) return;
    double baselik, candlik;
    // base likelihoood
    baselik = LikelihoodParamHazard(hazard,frailty,regression);
    // storage for old settings
    double * cand = (double *) calloc( hazard->np, sizeof(double));
    double * oldPar = (double *) malloc( hazard->np * sizeof(double));
    dcopyWrapper(hazard->np, hazard->ParamPar, oldPar);
    mvrnorm(hazard->np, cand, hazard->ParamPar, hazard->ParamCholCov, hazard->ParamTun[0]);
    // update curve at candidate
    UpdateParamPar(hazard,cand);
    // candidate likelihodo
    candlik = LikelihoodParamHazard(hazard,frailty,regression);
    int acc = AcceptReject(baselik, candlik, 1);
    if(!acc) // rejected, undo damage
        UpdateParamPar(hazard,oldPar);
    hazard->ParamAccept[0] = (double) acc;
    free(cand);
    free(oldPar);
}
/************ MH_ParamHazard */

/****f* CMetropolisHastings/MH_ParamFrailty
 *  NAME
 *    MH_ParamFrailty --- MH for frailty parametric component parameters
 *  FUNCTION
 *    Update the parameters indexing the parametric component by Metropolis-Hastings,
 *    see also mh.frailty.param
 *  INPUTS
 *    hazard        CCurve for the hazard
 *    frailty       CCurve for the frailty
 *    regression    CRegression structure
 *  SYNOPSIS
 */
void MH_ParamFrailty(curveP hazard, curveP frailty, regressionP regression)
/*
 *  SOURCE
*/
{ 
    if(!frailty->hasPar) return;
    double baselik, candlik;
    // base likelihood
    baselik = LikelihoodParamFrailty(hazard,frailty,regression);
    double * cand = (double *) calloc( frailty->np, sizeof(double));
    double * oldPar = (double *) malloc( frailty->np * sizeof(double));
    dcopyWrapper(frailty->np, frailty->ParamPar, oldPar);
    // generate candidate
    mvrnorm(frailty->np, cand, frailty->ParamPar, frailty->ParamCholCov, frailty->ParamTun[0]);
    // update curve with candidate parameter
    UpdateParamPar(frailty,cand);
    // candidate likelihood
    candlik = LikelihoodParamFrailty(hazard,frailty,regression);
    int acc = AcceptReject(baselik, candlik, 1);
    if(!acc) // rejected, undo damage
        UpdateParamPar(frailty,oldPar);
    frailty->ParamAccept[0] = (double) acc;
    free(cand);
    free(oldPar);
}
/************ MH_ParamFrailty */

/****f* CMetropolisHastings/MH_Weight
 *  NAME
 *    MH_Weight --- MH for weight of spline component
 *  FUNCTION
 *    Update the relative weight of the spline component in a curve with both spline
 *    and parametric components by Metropolis-Hastings. See also mh.weight.
 *  INPUTS
 *    theCurve      CCurve to be updated, can be either hazard or frailty
 *    hazard        CCurve for the hazard
 *    frailty       CCurve for the frailty
 *    regression    CRegression structure
 *  SYNOPSIS
 */
void MH_Weight(curveP theCurve, curveP hazard, curveP frailty, regressionP regression)
/*
 *  SOURCE
*/
{
    if(!theCurve->hasPar | !theCurve->hasSpline) return;
    // get the likelihood function, depending on whether the curve is the hazard or frailty
    double ( *likfun )(curveP, curveP, regressionP);
    if(theCurve->isHazard) likfun = &LikelihoodWeightHazard;
    else likfun = &LikelihoodWeightFrailty;  
    double oldW = theCurve->Weight[0];
    double w = dmin(dmax(oldW, .01), 0.99);
    double v = theCurve->WeightTun[0];
    double alpha = w*(w*(1-w)/v-1);
    double beta = (1-w)/w*alpha;
    // generate candidate as beta
    double cand = rbeta(alpha, beta);
    if(isnan(cand)){
        theCurve->WeightAccept[0]=0;
        return;
    }
    double alphac = cand*(cand*(1-cand)/v-1);
    double betac = (1-cand)/cand*alphac;
    //base likelihood
    double baselik = likfun(hazard, frailty, regression);
    theCurve->Weight[0] = cand;
    ReweightCurve(theCurve, -1);
    // candidate likelihood
    double candlik = likfun(hazard, frailty, regression);
    // transition rates
    double puc = dbeta(cand, alpha, beta, 0);
    double pcu = dbeta(w, alphac, betac, 0);
    int acc = AcceptReject(baselik, candlik, pcu/puc);
    if(!acc){
        theCurve->Weight[0] = oldW;
        ReweightCurve(theCurve, -1);
    }
    theCurve->WeightAccept[0] = (double) acc;
}
/************ MH_Weight */

/****f* CMetropolisHastings/MH_BDM
 *  NAME
 *    MH_BDM --- birth-death-move steps for spline knots
 *  FUNCTION
 *    Reversible-Jump MCMC steps for adding, deleting, or moving knots as part of the
 *    adaptive knot selection procedure. Can work with either the hazard or frailty
 *    curve, see also mh.bdm.
 *  INPUTS
 *    which         character, can be either 'h' or 'f', for hazard or frailty respectively
 *    hazard        CCurve for the hazard
 *    frailty       CCurve for the frailty
 *    regression    CRegression structure
 *  SYNOPSIS
 */
void MH_BDM(char which, curveP hazard, curveP frailty, regressionP regression)
/*
 *  SOURCE
*/
{
    curveP theCurve;
    if(which == 'h') theCurve = hazard;
    if(which == 'f') theCurve = frailty;
    if(!theCurve->hasSpline) return;
    // Pointers to useful components
    int ord = theCurve->SplineOrd;
    int nknots = theCurve->SplineNknots;
    double * candknots = theCurve->SplineCandKnots;
    int ncandknots = theCurve->SplineNCandKnots;
    double * occ = theCurve->SplineCandOcc; // occupied candidate knots
    int * occind = malloc(nknots * sizeof(int));
    // allocate temporary storage
    double * oldKnots = calloc(nknots+2*ord, sizeof(double));
    double * oldPar = calloc(nknots+ord, sizeof(double));
    double * knots2 = calloc(nknots+2*ord+1, sizeof(double));
    double * params2 = calloc(nknots+ord+1, sizeof(double));
    double * knots = theCurve->SplineKnots;
    double * params = theCurve->SplinePar;
    dcopyWrapper(nknots+2*ord,theCurve->SplineKnots,oldKnots);
    dcopyWrapper(nknots+ord,theCurve->SplinePar,oldPar);
    dcopyWrapper(nknots+2*ord,theCurve->SplineKnots,knots2);
    dcopyWrapper(nknots+ord,theCurve->SplinePar,params2);

    int i=0;
    // compute probability of birth-death and move steps
    for(int j=ord;j<ncandknots+2*ord; j++) if(occ[j]==1) occind[i++]=j;
    double pk = EvalNknotsPrior( nknots, theCurve );
    double pkp1 = (nknots < theCurve->SplineNknotsMax) ? EvalNknotsPrior(nknots+1, theCurve) : 0;
    double pkm1 = (nknots > 1) ? EvalNknotsPrior(nknots-1, theCurve) : 0;
    double pb = theCurve->SplineBDMConst[0] * dmin(1.0,pkp1/pk); // P(birth)
    double pd = theCurve->SplineBDMConst[0] * dmin(1.0,pkm1/pk); // P(death)
    double pm = dmax(0.0,1.0-pb-pd); // P(move)
    double u = runif(0,1);
    double baselik,candlik;
    // base likelihood
    if(theCurve->isHazard) baselik = LikelihoodSplineHazard(theCurve,frailty,regression);
    if(!theCurve->isHazard) baselik = LikelihoodSplineFrailty(hazard,theCurve,regression);
    if(u<pd){ // death step
        int j  = (int) floor(runif(0,(double) nknots))+ord; // index of dying knot
        double x = knots[j];  // value of dying knot
            //Rprintf("Remove knot %d at %f (occ: %d)\n",j,knots[j],(int) occ[occind[j-ord]]);
        for(int i=0;i<nknots+ord;i++) if(i>=j-1) params2[i] = params2[i+1]; //remove params2[j];
        for(int i=0;i<nknots+ord;i++) if(i>=j) knots2[i] = knots2[i+1]; //remove knots2[j2];
        // update spline parameters for knot deletion
        if(ord>2) for(int j2=j-ord+1; j2<j-1; j2++){
            double r2 = (x-knots2[j2])/(knots2[j2+ord-1]-knots2[j2]);
            double inner = 1/r2 * exp(params2[j2])-(1-r2)/r2*exp(params2[j2-1]);
            if(inner>0) params2[j2]=log(inner); else params2[j2] = theCurve->SplineMin[0];
        }
        theCurve->SplineNknots--; theCurve->nj--;
        // Jacobian for the transformation
        double J = exp(params2[j-1]) / (exp(params[j-1]) - exp(params[j-2]));
        if(ord>2) for(int j2 = j-ord+2; j2<j-1; j2++){
            double r2 = (x-knots2[j2])/(knots2[j2+ord-1]-knots2[j2]);
            J = J * exp(params2[j2])/(r2*exp(params[j2]));
        }
        // update the curve with the new knot set
        dcopyWrapper(nknots+2*ord, knots2,theCurve->SplineKnots);
        RemakeSplineBasis(theCurve,'d',j);
        // update curve with new parameters
        if(theCurve->isHazard) {
            UpdateSplinePar(theCurve,params2,-1);
            candlik = LikelihoodSplineHazard(theCurve,frailty,regression);
        }
        if(!theCurve->isHazard){
            // for the frailty curve, make sure that the mean is 1
            double newmean = -exp(params2[j-2])*theCurve->SplineBasisExp[j-2];
            for(int i=0;i<nknots+ord-1;i++) newmean += exp(params2[i])*theCurve->SplineBasisExp[i];
            double newinner = -newmean/theCurve->SplineBasisExp[j-2];
            if(newinner<0) candlik = -INFINITY;
            else{
                params2[j-2] = log(newinner);
                UpdateSplinePar(theCurve,params2,-1);
                candlik = LikelihoodSplineFrailty(hazard,theCurve,regression);
            }
        }
        // prior ratio * Jacobian
        double ratio = sqrt(2*M_PI*theCurve->SplinePriorvar[0])*fabs(J);
        int acc = AcceptReject(baselik, candlik, ratio);
        if(acc){
            // Update occupied index and candsd
            occ[occind[j-ord]] = 0;
        }else{ //undo damage
            theCurve->SplineNknots++; theCurve->nj++;
            dcopyWrapper(nknots+2*ord, oldKnots, theCurve->SplineKnots);
            RemakeSplineBasis(theCurve,'b',j);
            UpdateSplinePar(theCurve,oldPar,-1);
        }
    }
    if(u>pd & u<pd+pb){
        // Birth
        int birthind = occind[0];
        // choose a random unoccupied location for the knot to be born
        while(occ[birthind]) birthind = (int) floor(runif(0,(double) ncandknots))+ord;
        double x = candknots[birthind]; // value of the new knot
        int j = 0; while(knots[j]<x) j++; j--; // find the interval in which x lies
            //Rprintf("Birth at %f after knot %d\n",x,j);
        // update knot set and spline parameters
        for(int i=nknots+2*ord;i>j;i--) knots2[i]=knots2[i-1];
        knots2[j+1]=x;
        for(int i=nknots+ord;i>j-ord+1;i--) params2[i]=params2[i-1];
        for(int j2 = j-ord+2; j2<j+1; j2++){
            double r2 = (x-knots[j2])/(knots[j2+ord-1]-knots[j2]);
            if(j2==j) r2 = runif(0,1);
            params2[j2] = log(r2*exp(params[j2])+(1-r2)*exp(params[j2-1]));
        }
            //Rprintf("Birthind,x,j: %d %f %d\n",birthind,x,j);
        theCurve->SplineNknots++; theCurve->nj++;
        // Jacobian of the transformation
        double J = (exp(params[j])-exp(params[j-1]))/exp(params2[j]);
        if(ord>2) for(int j2 = j-ord+2; j2<j; j2++){
            double r2 = (x-knots[j2])/(knots[j2+ord-1]-knots[j2]);
            J = J * r2 * exp(params[j2])/exp(params2[j2]);
        }
        // update curve with new knots and parameters
        dcopyWrapper(nknots+2*ord+1, knots2,theCurve->SplineKnots);
        RemakeSplineBasis(theCurve,'b',j);
        if(theCurve->isHazard){
            UpdateSplinePar(theCurve,params2,-1);
            candlik = LikelihoodSplineHazard(theCurve,frailty,regression);
        }
        if(!theCurve->isHazard){
            // for frailty, make sure the mean is 1
            double newmean = -exp(params2[j])*theCurve->SplineBasisExp[j];
            for(int i=0;i<nknots+ord+1;i++) newmean += exp(params2[i])*theCurve->SplineBasisExp[i];
            double newinner = -newmean/theCurve->SplineBasisExp[j];
            if(newinner<0) candlik = -INFINITY;
            else{
                params2[j] = log(newinner);
                UpdateSplinePar(theCurve,params2,-1);
                candlik = LikelihoodSplineFrailty(hazard,theCurve,regression);
            }
        }
        // prior ratio * Jacobian
        double ratio = 1/sqrt(2*M_PI*theCurve->SplinePriorvar[0]) * fabs(J);
        int acc = AcceptReject(baselik, candlik, ratio);
            //Rprintf("Lik: %f %f %f %f %d\n",baselik,candlik, J, ratio,acc);
        if(acc){
            // Update set of occupied candidate knots
            occ[birthind] = 1;
        }else{ //undo damage
            theCurve->SplineNknots--; theCurve->nj--;
            dcopyWrapper(nknots+2*ord, oldKnots, theCurve->SplineKnots);
            RemakeSplineBasis(theCurve,'d',j+1);
            UpdateSplinePar(theCurve,oldPar,-1);
        }
    }
    if(u>pb+pd) { // Move a knot
        // choose a random knot to move
        int moveind=0;
        moveind  = (int) floor(runif(0,(double) nknots));
        // find range of movement
        int leftknotind = (moveind == 0) ? ord : occind[moveind-1]+1;
        int rightknotind = (moveind == nknots-1) ? ncandknots + ord - 1 : occind[moveind+1]-1;
        // choose a candidate knot to place the new knot in
        int newknotind = (int) floor(runif((double) leftknotind,(double) rightknotind+1));
        // if the knot doesn't stay the same, update the curve
        if(candknots[newknotind] != oldKnots[moveind+ord]){ 
            theCurve->SplineKnots[moveind+ord] = candknots[newknotind];
            // update spline basis with new knot position
            RemakeSplineBasis(theCurve,'m',moveind);
            if(theCurve->isHazard){ //hazard
                EvalSpline(theCurve,-1);
                candlik = LikelihoodSplineHazard(theCurve,frailty,regression);
            }
            if(!theCurve->isHazard){ //frailty
                // frailty needs to ensure the mean is 1
                double newmean = -exp(params2[moveind+ord])*theCurve->SplineBasisExp[moveind+ord];
                for(int i=0;i<nknots+ord;i++) newmean += exp(params2[i])*theCurve->SplineBasisExp[i];
                double newinner = -newmean/theCurve->SplineBasisExp[moveind+ord];
                if(newinner<0) candlik = -INFINITY;
                else{
                    params2[moveind+ord] = log(newinner);
                    UpdateSplinePar(theCurve,params2,-1);
                    candlik = LikelihoodSplineFrailty(hazard,theCurve,regression);
                }
            }
            int acc = AcceptReject(baselik,candlik,1.0);
            //Rprintf("Lik: %f %f %d\n",baselik,candlik,acc);
            if(!acc){
                // not accepted, undo the damage
                theCurve->SplineKnots[moveind+ord] = oldKnots[moveind+ord];
                RemakeSplineBasis(theCurve,'m',moveind);
                if(theCurve->isHazard){ // hazard
                    EvalSpline(theCurve,-1);
                }else{      // frailty
                    UpdateSplinePar(theCurve,oldPar,-1);
                }
            }else{ // update occupied index
                occ[occind[moveind]]=0;
                occ[newknotind] = 1;
            }
        }
    }

    free(oldKnots);
    free(oldPar);
    free(knots2);
    free(params2);
    free(occind);
}
/************ MH_BDM */

/****f* CMetropolisHastings/UpdatePostvarCurve
 *  NAME
 *    UpdatePostvarCurve --- update prior variance for a CCurve
 *  FUNCTION
 *    Update the prior variance for the spline and parametric parameters of a hazard or frailty curve,
 *    with an inverse-gamma prior on the hyperparameters. See also updatepostvar.curve.
 *  INPUTS
 *    theCurve      a CCurve structure
 *  OUTPUTS
 *    theCurve      the input, with SplinePriorVar component updated
 *  SYNOPSIS
 */
void UpdatePostvarCurve(curveP theCurve)
/*
 *  SOURCE
*/
{
    int c1 = 1;
    // spline prior variance
    if(theCurve->hasSpline) theCurve->SplinePriorvar[0] = rinvgamma(
            (double) theCurve->nj / 2 + theCurve->SplineHyper[0],
            theCurve->SplinePenaltyFactor[0] * SmoothnessPenalty(theCurve) * theCurve->SplinePriorvar[0] 
                + theCurve->SplineHyper[1] );

    // parametric component prior variance
    if(theCurve->hasPar) theCurve->ParamPriorvar[0] = rinvgamma(
            (double) theCurve->np / 2 + theCurve->ParamHyper[0],
            pow(F77_CALL(dnrm2)(&(theCurve->np), theCurve->ParamPar, &c1),2)/2 + 
                theCurve->ParamHyper[1] );
}
/************ UpdatePostvarCurve */

/****f* CMetropolisHastings/UpdatePostvarRegression
 *  NAME
 *    UpdatePostvarRegression --- update prior variance for regression coefficients
 *  FUNCTION
 *    Update the prior varaince of regression coefficients using inverse-gamma prior and
 *    the given hyperparameters.
 *  INPUTS
 *    theReg        a CRegression structure
 *  OUTPUTS
 *    theReg        the input, with theReg->priorvar updated
 *  SYNOPSIS
 */
void UpdatePostvarRegression(regressionP theReg)
/*
 *  SOURCE
*/
{
    int c1 = 1;
    theReg->priorvar[0] = rinvgamma(
            (double) theReg->p / 2 + theReg->hyper[0],
            pow(F77_CALL(dnrm2)(&(theReg->p), theReg->coefficients, &c1),2)/2+ theReg->hyper[1] );
}
/************ UpdatePostvarRegression */

/****f* CmiscUtils/UpdateHistory
 *  NAME
 *    UpdateHistory --- update history of parameters
 *  FUNCTION
 *    Store the state of the chain at the current iteration in the CHistory structure
 *  INPUTS
 *    hazard        CCurve for the hazard
 *    frailty       CCurve for the frailty
 *    regression    CRegression structure
 *    history       a CHistory structure
 *    iter          the current iteration counter
 *  SYNOPSIS
 */
void UpdateHistory(curveP hazard, curveP frailty, regressionP regression, historyP history, int iter)
/*
 *  SOURCE
*/
{
    int c1=1;
    int ny = history->ny;
    // store frailties
    F77_CALL(dcopy)(&(frailty->nx), frailty->X, &c1, history->frailty + iter-1, &(history->ny));
    // store regression coefficients
    F77_CALL(dcopy)(&(regression->p), regression->coefficients, &c1, history->coefficients + iter-1,
            &(history->ny));

    if(frailty->hasSpline){
        // store frailty spline knots and parameters
        int lknots = frailty->SplineNknots + 2*frailty->SplineOrd;
        F77_CALL(dcopy)(&(frailty->nj), frailty->SplinePar, &c1, history->FrailtySplinePar + iter-1,
                &(history->ny));
        F77_CALL(dcopy)(&(lknots), frailty->SplineKnots, &c1, history->FrailtySplineKnots + iter-1,
                &(history->ny));
        history->FrailtySplineFvar[iter-1] = frailty->SplineFvar;
    }
    if(hazard->hasSpline){
        // store hazard spline knots and parameters
        int lknots = hazard->SplineNknots + 2*hazard->SplineOrd;
        F77_CALL(dcopy)(&(hazard->nj), hazard->SplinePar, &c1, history->HazardSplinePar + iter-1,
                &(history->ny));
        F77_CALL(dcopy)(&(lknots), hazard->SplineKnots, &c1, history->HazardSplineKnots + iter-1,
                &(history->ny));
    }
    // store parametric component parameters
    if(frailty->hasPar)
        F77_CALL(dcopy)(&(frailty->np), frailty->ParamPar, &c1, history->FrailtyParamPar + iter-1,
                &(history->ny));
    if(hazard->hasPar)
        F77_CALL(dcopy)(&(hazard->np), hazard->ParamPar, &c1, history->HazardParamPar + iter-1,
                &(history->ny));
    // store weights
    if(hazard->hasPar & hazard->hasSpline)
        history->HazardWeight[iter-1] = hazard->Weight[0];
    if(frailty->hasPar & frailty->hasSpline)
        history->FrailtyWeight[iter-1] = frailty->Weight[0];

    // store prior variances
    history->priorvar[iter-1 + ny*0] = regression->priorvar[0];
    history->priorvar[iter-1 + ny*1] = hazard->SplinePriorvar[0];
    history->priorvar[iter-1 + ny*2] = frailty->SplinePriorvar[0];
    history->priorvar[iter-1 + ny*3] = hazard->ParamPriorvar[0];
    history->priorvar[iter-1 + ny*4] = frailty->ParamPriorvar[0];
    history->priorvar[iter-1 + ny*5] = hazard->WeightPriorvar[0];
    history->priorvar[iter-1 + ny*6] = frailty->WeightPriorvar[0];

    // store acceptance rates
    history->accept[iter-1 + ny*0] = regression->Accept[0];
    history->accept[iter-1 + ny*1] = hazard->SplineAccept[0];
    history->accept[iter-1 + ny*2] = frailty->SplineAccept[0];
    history->accept[iter-1 + ny*3] = hazard->ParamAccept[0];
    history->accept[iter-1 + ny*4] = frailty->ParamAccept[0];
    history->accept[iter-1 + ny*5] = hazard->WeightAccept[0];
    history->accept[iter-1 + ny*6] = frailty->WeightAccept[0];
    history->accept[iter-1 + ny*7] = frailty->Accept[0];

    // store full loglikelihood
    history->loglik[iter-1] = LikelihoodFull(hazard, frailty, regression);
}
/************ UpdateHistory */


/****f* 00main/SplineSurvMainLoop
 *  NAME
 *    SplineSurvMainLoop --- main loop in C
 *  FUNCTION
 *    This is the main fitting loop in C, which is called by splinesurv.agdata if it
 *    is called with usec=TRUE (default). Everything in the set of CFitting routines is
 *    implemented to match the RFitting routines, with the caveat of C's memory management.
 *
 *    See also the call graph linked from 00main.
 *  INPUTS
 *    Rhazard       RCurve structure for the hazard
 *    Rfrailty      RCurve structure for the frailty
 *    Rregression   RRegression structure for the regression information
 *    Rhistory      RHistory structure for the history of the parameters
 *    Rstartiter    integer, for the starting iteration from which the loop should begin
 *    Renditer      integer, for the iteration until which the loop should run
 *    Rthin         integer, for the number of iterations that should be discarded each time
 *    Rverbose      integer, containing the amount of output that should be printed to the
 *                  screen. This is not as well-behaved as the verbose option within R.
 *  SYNOPSIS
 */
SEXP SplineSurvMainLoop( SEXP Rhazard, SEXP Rfrailty, SEXP Rregression, 
        SEXP Rhistory, SEXP Rstartiter, SEXP Renditer, SEXP Rthin, 
        SEXP Rverbose)
/*
 *  SOURCE
*/
{
    
    // initialize random number generator.
    GetRNGstate();
    // Create local curve and regression structures, as CCurve, CRegression 
    // and CHistory structures. These consist primarily of pointers to the memory
    // allocated by R, but they make the addressing easier.
    curveP hazard = (struct curve *) R_alloc(1,sizeof(struct curve));
    curveP frailty = (struct curve *) R_alloc(1,sizeof(struct curve));
    regressionP regression = (struct regress *) R_alloc(1,sizeof(struct regress));
    historyP history = (struct hist *) R_alloc(1,sizeof(struct hist));
    hazard->isHazard=1;
    frailty->isHazard=0;
    PopulateLocalCurve(hazard, Rhazard);
    PopulateLocalCurve(frailty, Rfrailty);
    PopulateLocalRegression(regression, Rregression);
    PopulateLocalHistory(history,Rhistory);
    for(int j=0; j<regression->m; j++) 
        for(int i=regression->Jicum[j]; i<regression->Jicum[j+1]; i++)
            regression->frailrep[i] = frailty->X[j];
    // compute the value of the frailty*lp vector
    diagmvWrapper(regression->n, regression->frailrep, regression->elp, regression->frailelp);

    // Begin main loop
    int iter = asInteger(Rstartiter); // iteration counter
    int enditer = asInteger(Renditer);
    int thin = asInteger(Rthin);
    int verbose = asInteger(Rverbose);
    int iterEl = 0; // counts the elapsed iterations between history updates
    while(iter < enditer)
    {
        R_CheckUserInterrupt();
        iterEl++;
        if(verbose >= 4) Rprintf("%d", iterEl);
        
        // Metropolis-Hastings steps
        MH_Frail(hazard,frailty,regression);

        MH_Regression(hazard,frailty,regression);
        
        MH_SplineHazard(hazard,frailty,regression);

        MH_SplineFrailty(hazard,frailty,regression);
        
        MH_ParamHazard(hazard,frailty,regression);

        MH_ParamFrailty(hazard,frailty,regression);

        MH_Weight(hazard, hazard, frailty, regression);
        MH_Weight(frailty, hazard, frailty, regression);

        UpdatePostvarCurve(hazard);
        UpdatePostvarCurve(frailty);
        UpdatePostvarRegression(regression);

        // RJMCMC steps for adaptive knot selection
        if(hazard->SplineAdaptive) MH_BDM('h',hazard,frailty,regression);
        if(frailty->SplineAdaptive) MH_BDM('f',hazard,frailty,regression);
        
        // store current state in CHistory
        if(iterEl == thin){
            iter++;
            UpdateHistory(hazard, frailty, regression, history, iter);
            iterEl = 0;
            if(verbose>=3) Rprintf(" %d ",iter);
        }
        
    }
    REAL(getListElement(Rhazard,"spline.nknots"))[0] = (double) hazard->SplineNknots;
    REAL(getListElement(Rfrailty,"spline.nknots"))[0] = (double) frailty->SplineNknots;
    SEXP returnval = R_NilValue;
    PutRNGstate();
    return returnval;
}

/************ SplineSurvMainLoop */
