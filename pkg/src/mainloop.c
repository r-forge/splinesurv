
#include <R_ext/Linpack.h>    
#include <R_ext/Lapack.h>    
#include <R_ext/Print.h>  
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <math.h>
#include <stdlib.h>
#include "init.h"

#define DEBUG

typedef enum {pnone=0, pdiff=1, pderiv=2, plogderiv=3} penalty;
typedef enum {Dnone=0, Dexponential=1, Dweibull=2, Dlognormal=3, Dgamma=4} distribution;

typedef struct curve {
    int hasSpline, //has a spline component
         hasPar, //has a parametric component
         isHazard, //whether the curve represents hazard or frailty
         SplineOrd, // order of the spline
         SplineNknots, //spline number of knots
         nx, //number of observations
         nj, //number of basis functions (spline)
         np, // number of parameters (parametric)
         SplineFixedInd; // fixed index for frailty spline
    double SplineEParSum;
    penalty SplinePenaltyType; // (0=none, 1=diff, 2=2deriv, 3=log2der)
    distribution ParDist; // Parametric distribution function
    double *SplineKnots, //Knots of the spline
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

typedef struct regress {
    int m,  // # of clusters
        n,  // # of observations
        p,  // # of regression parameters
        *Ji, // subjects per cluster
        *cluster; // Cluster ID
    double *coefficients, // coefficient estimates
           *covariates, // covariate estimates
           *lp, // linear predictor
           *elp, // exponential of linear predictor
           *status, 
           *time,
           *CandCov,  // covariance matrix of candidates
           *CholCov, // Cholesky factor of covariance
           *priorvar, // variance of prior
           *hyper,   // hyperparameters of prior
           *Accept, // Acceptance indicator
           *tun;    // tuning parameter
} *regressionP;

typedef struct hist {
    int ny; // number of iterations
    double *frailty,
           *coefficients,
           *HazardSplinePar,
           *FrailtySplinePar,
           *HazardParamPar,
           *FrailtyParamPar,
           *HazardWeight,
           *FrailtyWeight,
           *priorvar,
           *accept;
} *historyP;

inline double dmin(double x, double y) {
    return x<y ? x : y;
}

inline double dmax(double x, double y) {
    return x>y ? x : y;
}

SEXP getListElement(SEXP list, const char *str)
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

void dcopyWrapper( int n, double *x, double *y)
{
    int c1=1;
    F77_CALL(dcopy)(&n, x, &c1, y, &c1);
}

double ddotWrapper( int n, double *x, double *y)
{
    int c1=1;
    return F77_CALL(ddot)(&n, x, &c1, y, &c1);
}

void mvrnorm(int n, double *out, double *mu, double *CholSigma, double tun)
{
    double * temp = (double *) malloc(n * sizeof(double));
    GetRNGstate();
    for(int i=0; i<n; i++) temp[i] = rnorm(0,1);
    PutRNGstate();
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

double rinvgamma(double shape, double scale)
{
    GetRNGstate();
    double out = 1.0 / rgamma(shape, 1.0/scale);
    PutRNGstate();
    return out;
}

// Populate a local curve structure based on an R expresion
void PopulateLocalCurve( curveP theCurve, SEXP Rcurve)
{
    theCurve->hasPar = asInteger(getListElement(Rcurve,"haspar"));
    theCurve->hasSpline = asInteger(getListElement(Rcurve,"hasspline"));
    theCurve->nx = (int) (length(getListElement(Rcurve,"x")));
    theCurve->SplineMin = REAL(getListElement(Rcurve,"spline.min"));
    theCurve->SplinePriorvar =  REAL(getListElement(Rcurve,"spline.priorvar"));
    theCurve->SplineHyper =  REAL(getListElement(Rcurve,"spline.hyper"));
    theCurve->SplineTun =  REAL(getListElement(Rcurve,"spline.tun"));
    theCurve->SplineAccept =  REAL(getListElement(Rcurve,"spline.accept"));
    theCurve->ParamPriorvar = REAL(getListElement(Rcurve,"param.priorvar"));
    theCurve->ParamTun = REAL(getListElement(Rcurve,"param.tun"));
    theCurve->ParamHyper = REAL(getListElement(Rcurve,"param.hyper"));
    theCurve->ParamAccept = REAL(getListElement(Rcurve,"param.accept"));
    theCurve->Weight = REAL(getListElement(Rcurve,"weight"));
    theCurve->WeightPriorvar = REAL(getListElement(Rcurve,"weight.priorvar"));
    theCurve->WeightHyper = REAL(getListElement(Rcurve,"weight.hyper"));
    theCurve->WeightTun = REAL(getListElement(Rcurve,"weight.tun"));
    theCurve->WeightAccept = REAL(getListElement(Rcurve,"weight.accept"));
    theCurve->X = REAL(getListElement(Rcurve,"x"));
    theCurve->Y = REAL(getListElement(Rcurve,"y"));
    if(theCurve->isHazard) theCurve->Ycum =REAL(getListElement(Rcurve,"ycum"));
    if(!theCurve->isHazard) theCurve->Accept =REAL(getListElement(Rcurve,"accept"));
    if(!theCurve->isHazard) theCurve->tun = REAL(getListElement(Rcurve,"tun"));
    if(theCurve->hasSpline){
        theCurve->SplineOrd=asInteger(getListElement(Rcurve,"spline.ord"));
        theCurve->SplineNknots = (int) (length(getListElement(Rcurve,"spline.knots")));
        theCurve->nj = theCurve->SplineNknots - theCurve->SplineOrd;
        const char * charPenaltyType = CHAR(STRING_ELT(getListElement(Rcurve,"spline.penalty"),0));
        penalty iPenaltyType;
        if (strcmp(charPenaltyType,"2diff") ==0) { iPenaltyType=pdiff; }
        else if (strcmp(charPenaltyType,"2deriv") ==0) { iPenaltyType=pderiv; }
        else if (strcmp(charPenaltyType,"log2deriv") ==0) { iPenaltyType=plogderiv;}
        else { iPenaltyType=pnone;}
        theCurve->SplinePenaltyType = iPenaltyType;
        if(iPenaltyType != pnone)
            theCurve->SplinePenaltyMatrix = REAL(getListElement(Rcurve, "spline.penaltymatrix"));
        theCurve->SplinePenaltyFactor =  REAL(getListElement(Rcurve,"spline.penaltyfactor"));
        theCurve->SplineKnots =  REAL(getListElement(Rcurve,"spline.knots"));
        theCurve->SplineBasis =  REAL(getListElement(Rcurve,"spline.basis"));
        theCurve->SplineBasisInt =  REAL(getListElement(Rcurve,"spline.basisint"));
        theCurve->SplinePar =  REAL(getListElement(Rcurve,"spline.par"));
        theCurve->SplineEPar = (double *) malloc( theCurve->nj * sizeof(double));
        theCurve->SplineCandCov =  REAL(getListElement(Rcurve,"spline.candcov"));
        theCurve->SplineCholCov =  REAL(getListElement(Rcurve,"spline.cholcandcov"));
        for(int j=0; j< theCurve->nj; j++) theCurve->SplineEPar[j] = exp(theCurve->SplinePar[j]);
        theCurve->SplineY =  REAL(getListElement(Rcurve,"spline.y"));
        if(theCurve->isHazard){
            theCurve->SplineBasisCum =  REAL(getListElement(Rcurve,"spline.basiscum"));
            theCurve->SplineYcum =  REAL(getListElement(Rcurve,"spline.ycum"));
        }
        if(!theCurve->isHazard){
            theCurve->SplineBasisInt =  REAL(getListElement(Rcurve,"spline.basisint"));
            theCurve->SplineBasisExp =  REAL(getListElement(Rcurve,"spline.basisexp"));
            theCurve->SplineEParSum = 0;
            for(int j=0; j<theCurve->nj; j++) theCurve->SplineEParSum+=theCurve->SplineEPar[j];
            theCurve->SplineMeanPenalty = REAL(getListElement(Rcurve,"spline.meanpenalty"));
            theCurve->SplineFixedInd = asInteger(getListElement(Rcurve, "spline.fixedind")) - 1;
        }
    }
    if(theCurve->hasPar){
        theCurve->ParamPar = REAL(getListElement(Rcurve,"param.par"));
        theCurve->ParamY = REAL(getListElement(Rcurve,"param.y"));
        if(theCurve->isHazard)
            theCurve->ParamYcum = REAL(getListElement(Rcurve,"param.ycum"));
        const char * charParamDist = CHAR(STRING_ELT(getListElement(Rcurve,"param.dist"),0));
        distribution iParamDist;
        if(strcmp(charParamDist,"exponential") == 0) { iParamDist = Dexponential; } 
        else if(strcmp(charParamDist,"weibull") == 0) { iParamDist = Dweibull; } 
        else if(strcmp(charParamDist,"gamma") == 0) { iParamDist = Dgamma; } 
        else if(strcmp(charParamDist,"lognormal") == 0) { iParamDist = Dlognormal; } 
        else { iParamDist = Dnone; } 
        theCurve->ParDist = iParamDist;
        theCurve->np = (int) (length(getListElement(Rcurve,"param.par")));
        theCurve->ParamCandCov =  REAL(getListElement(Rcurve,"param.candcov"));
        theCurve->ParamCholCov =  REAL(getListElement(Rcurve,"param.cholcandcov"));
    }
}

void FreeCurveMem(curveP theCurve)
{
    if(theCurve->hasSpline)
        free(theCurve->SplineEPar);
}

void PopulateLocalRegression( regressionP theReg, SEXP Rregression){
    theReg->m = asInteger(getListElement(Rregression, "m"));
    theReg->Ji = INTEGER(getListElement(Rregression, "Ji"));
    theReg->n = (int) (length(getListElement(Rregression,"cluster")));
    theReg->p = (int) (length(getListElement(Rregression,"coefficients")));
    theReg->cluster = INTEGER(getListElement(Rregression, "cluster"));

    theReg->coefficients = REAL(getListElement(Rregression,"coefficients"));
    theReg->covariates = REAL(getListElement(Rregression,"covariates"));
    theReg->lp = REAL(getListElement(Rregression,"lp"));
    theReg->status = REAL(getListElement(Rregression,"status"));
    theReg->time = REAL(getListElement(Rregression,"time"));
    theReg->CandCov = REAL(getListElement(Rregression,"candcov"));
    theReg->CholCov = REAL(getListElement(Rregression,"cholcandcov"));
    theReg->priorvar = REAL(getListElement(Rregression,"priorvar"));
    theReg->hyper = REAL(getListElement(Rregression,"hyper"));
    theReg->Accept = REAL(getListElement(Rregression,"accept"));
    theReg->tun = REAL(getListElement(Rregression,"tun"));
    theReg->elp = (double *) malloc( theReg->n * sizeof(double));
    for(int i=0; i<theReg->n; i++) theReg->elp[i]=exp(theReg->lp[i]);
}

void FreeRegMem(regressionP theReg)
{
    free(theReg->elp);
}


void PopulateLocalHistory( historyP theHist, SEXP Rhistory)
{
    SEXP elmt;
    elmt = getListElement(Rhistory, "frailty");
    theHist->ny = INTEGER(getAttrib(elmt, R_DimSymbol))[0];
    theHist->frailty = REAL(elmt);
    theHist->coefficients = REAL(getListElement(Rhistory, "coefficients"));
    // Only populate history components if there is a spline/parametric component
    elmt = getListElement(Rhistory, "hazard.spline.par");
    if(elmt != R_NilValue ) theHist->HazardSplinePar = REAL(elmt);
    elmt = getListElement(Rhistory, "hazard.param.par");
    if(elmt != R_NilValue ) theHist->HazardParamPar = REAL(elmt);
    elmt = getListElement(Rhistory, "hazard.weight");
    if(elmt != R_NilValue ) theHist->HazardWeight = REAL(elmt);
    elmt = getListElement(Rhistory, "frailty.spline.par");
    if(elmt != R_NilValue ) theHist->FrailtySplinePar = REAL(elmt);
    elmt = getListElement(Rhistory, "frailty.param.par");
    if(elmt != R_NilValue ) theHist->FrailtyParamPar = REAL(elmt);
    elmt = getListElement(Rhistory, "frailty.weight");
    if(elmt != R_NilValue ) theHist->FrailtyWeight = REAL(elmt);
    theHist->priorvar = REAL(getListElement(Rhistory, "priorvar"));
    theHist->accept = REAL(getListElement(Rhistory, "accept"));
}

double SmoothnessPenalty(curveP theCurve)
{
    if(theCurve->SplinePenaltyType == pnone) return 0;
    double * temp = (double *) calloc(theCurve->nj, sizeof(double));
    double * par;
    double pen;
    if(theCurve->SplinePenaltyType == pdiff) par = theCurve->SplinePar;
    else par = theCurve->SplineEPar;
    int c1 = 1;
    double c1d = 1.0;
    char uplo = 'U';
    F77_CALL(dsymv)( &uplo, &(theCurve->nj), &c1d, theCurve->SplinePenaltyMatrix, &(theCurve->nj), par, &c1, &c1d, temp, &c1);
    pen = F77_CALL(ddot)( &(theCurve->nj), temp, &c1, par, &c1);
    if(theCurve->SplinePenaltyType == plogderiv) pen=log(pen+1);
    pen = pen / (2 * theCurve->SplinePriorvar[0]);
    free(temp);
    return pen;
}

void ReweightCurve(curveP theCurve, int i){
    if(theCurve->hasPar & !theCurve->hasSpline) { //parametric only
        dcopyWrapper(theCurve->nx, theCurve->ParamY, theCurve->Y);
        if(theCurve->isHazard) dcopyWrapper(theCurve->nx, theCurve->ParamYcum, theCurve->Ycum);
        return;
    }
    if(!theCurve->hasPar & theCurve->hasSpline) { //spline only
        dcopyWrapper(theCurve->nx, theCurve->SplineY, theCurve->Y);
        if(theCurve->isHazard) dcopyWrapper(theCurve->nx, theCurve->SplineYcum, theCurve->Ycum);
        return;
    }
    double w = theCurve->Weight[0];
    if(i>=0){
        theCurve->Y[i] = w *  theCurve->SplineY[i] + (1-w) * theCurve->ParamY[i];
        if(theCurve->isHazard) theCurve->Ycum[i] = w *  theCurve->SplineY[i] + (1-w) * theCurve->ParamY[i];
    }else{
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

void UpdateSplineBasis(curveP theCurve, int i)
{
    if(!theCurve->hasSpline) return;
    if(theCurve->isHazard) Rprintf("Cannot update the hazard basis");
    if(i<0)
        for(int ind=0; ind<theCurve->nx; ind++) UpdateSplineBasis(theCurve, ind);
    else{
        for(int j=0; j<theCurve->nj; j++)
            theCurve->SplineBasis[i + j * theCurve->nx] = csplineeval(theCurve->X[i], j, theCurve->SplineOrd, theCurve->SplineKnots, theCurve->SplineOrd, theCurve->nj);
        if(!theCurve->isHazard)
            for(int j=0; j<theCurve->nj; j++) theCurve->SplineBasis[i+j * theCurve->nx] /= theCurve->SplineBasisInt[j];
    }
}


void EvalSpline(curveP theCurve, int i)
{
    if(!theCurve->hasSpline) return;
    if(i>=0){
        int c1=1;
        theCurve->SplineY[i] = F77_CALL(ddot)(&(theCurve->nj), theCurve->SplineBasis + i, &(theCurve->nx), theCurve->SplineEPar, &c1);
        if(!theCurve->isHazard) theCurve->SplineY[i] /= theCurve->SplineEParSum;
        if(theCurve->isHazard)
            theCurve->SplineYcum[i] = F77_CALL(ddot)(&(theCurve->nj), theCurve->SplineBasisCum + i, &(theCurve->nx), theCurve->SplineEPar, &c1);
    }else{
        int c1 = 1;
        double c0 = 0;
        char trans = 'N';
        double scaler = theCurve->isHazard ? 1.0  : 1.0/ theCurve->SplineEParSum;
        // set splineY = basis * splinepar + 0*splineY
        F77_CALL(dgemv)(&trans, &(theCurve->nx), &(theCurve->nj), &scaler, theCurve->SplineBasis, &(theCurve->nx), theCurve->SplineEPar, &c1, &c0, theCurve->SplineY, &c1);
        if(theCurve->isHazard)
            F77_CALL(dgemv)(&trans, &(theCurve->nx), &(theCurve->nj), &scaler, theCurve->SplineBasisCum, &(theCurve->nx), theCurve->SplineEPar, &c1, &c0, theCurve->SplineYcum, &c1);
    }
    ReweightCurve(theCurve, i);
}

double EvalSplineAtOnePoint(curveP theCurve, double x)
{
    double splY=0;
    int c1=1;
    if(theCurve->hasSpline){
        double * tempBasis = (double *) malloc( theCurve->nj * sizeof(double));
        for(int j=0; j<theCurve->nj; j++) tempBasis[j]=csplineeval(x, j, theCurve->SplineOrd,
                theCurve->SplineKnots, theCurve->SplineOrd, theCurve->nj);
        if(!theCurve->isHazard)
            for(int j=0; j<theCurve->nj; j++) tempBasis[j]/= (theCurve->SplineBasisInt[j] * theCurve->SplineEParSum);
        splY = F77_CALL(ddot)( &(theCurve->nj), tempBasis, &c1, theCurve->SplineEPar, &c1);
        free(tempBasis);
    }
    return splY;
}

void UpdateSplinePar(curveP theCurve, double * newpar)
{
    if(!theCurve->hasSpline) return;
    dcopyWrapper(theCurve->nj, newpar, theCurve->SplinePar);
    for(int i=0; i<theCurve->nj; i++) theCurve->SplineEPar[i] = exp(theCurve->SplinePar[i]);
    EvalSpline(theCurve, -1);
}

double EvalParamAtOnePoint(curveP theCurve, double x, int cum)
{
    double parY = 0;
    double parYcum =0;
    if(!theCurve->hasPar) return parY;
    if(theCurve->isHazard){
        if(theCurve->ParDist == Dexponential){
            double lambda = exp(theCurve->ParamPar[0]);
            parY = lambda;
            parYcum = lambda*x;
        }else if(theCurve->ParDist == Dweibull){
            double lambda = exp(theCurve->ParamPar[0]);
            double alpha = exp(theCurve->ParamPar[1]);
            parY = alpha*lambda*pow(x,alpha-1);
            parYcum = lambda * pow(x,alpha);
        }else{
            Rprintf("distribution not implemented");
        }
    }else{
        if(theCurve->ParDist == Dgamma){
            double alpha = exp(- theCurve->ParamPar[0]);
            parY = dgamma(x, alpha, 1/alpha,0);
        }else if(theCurve->ParDist == Dlognormal){
            double alpha = exp(theCurve->ParamPar[0]);
            parY = exp(-pow(log(x)+alpha/2,2)/(2*alpha)) / (x*sqrt(2*M_PI*alpha));
        }
    }
    return (cum == 0) ? parY : parYcum;
}

void EvalParametric(curveP theCurve, int i)
{
    if(!theCurve->hasPar) return;
    if(i<0){
        for(int ind=0; ind<theCurve->nx; ind++) EvalParametric(theCurve, ind);
    }else{
        theCurve->ParamY[i] = EvalParamAtOnePoint(theCurve, theCurve->X[i], 0);
        if(theCurve->isHazard)
            theCurve->ParamYcum[i] = EvalParamAtOnePoint(theCurve, theCurve->X[i], 1);
    }
    ReweightCurve(theCurve, i);
}

void UpdateParamPar(curveP theCurve, double * newPar)
{
    if(!theCurve->hasPar) return;
    for(int i=0; i<theCurve->np; i++)
        theCurve->ParamPar[i] = newPar[i];
    EvalParametric(theCurve, -1);
}

double EvalCurveAtOnePoint(curveP theCurve, double x)
{
    if(theCurve->isHazard) Rprintf("Error: using EvalCurveAtOnePoint on hazard");
    double splY=EvalSplineAtOnePoint(theCurve,x);
    double parY=EvalParamAtOnePoint(theCurve,x,0);
    if(theCurve->hasPar & theCurve->hasSpline)
        return theCurve->Weight[0] * splY + (1-theCurve->Weight[0])*parY;
    else if(theCurve->hasPar)
        return parY;
    else if(theCurve->hasSpline)
        return splY;
    else return 0;
}

void UpdateCurveX(curveP theCurve, double x, int i)
{
    if(theCurve->isHazard) Rprintf("Error: using UpdateCurveX on hazard.");
    theCurve->X[i] = x;
    UpdateSplineBasis(theCurve, i);
    EvalSpline(theCurve,i);
    EvalParametric(theCurve,i);
}

double LikelihoodFrailty(int i, curveP hazard, curveP frailty, regressionP regression)
{
   double Ui = frailty->X[i];
   double lUi = log(Ui);
   int clust = i+1;
   double lik = 0;
   lik += log(frailty->Y[i]);
   for(int ind = 0; ind < hazard->nx; ind++){
        if(regression->cluster[ind] < clust) continue;
        if(regression->cluster[ind] > clust) break;
        lik += (regression->status[ind] * lUi);
        lik -= Ui * hazard->Ycum[ind] * regression->elp[ind];
   }
   return(lik);
}

double LikelihoodRegression(curveP hazard,curveP frailty,regressionP regression)
{
   double lik; 
   double * elp = regression->elp;
   int * cluster = regression->cluster;
   double * frailX = frailty->X;
   double * hazYcum = hazard->Ycum;

   lik = ddotWrapper(regression->n, regression->status, regression->lp);
   for(int i=0; i< regression->n; i++) lik-=frailX[cluster[i]-1] * hazYcum[i] * elp[i];
   int c1=1;
   lik -= pow(F77_CALL(dnrm2)(&(regression->p), regression->coefficients, &c1),2)/(2*regression->priorvar[0]);
   return lik;
}

double LikelihoodSplineHazard(curveP hazard, curveP frailty, regressionP regression)
{
    double lik=0;
    double * elp = regression->elp;
    int * cluster = regression->cluster;
    double * frailX = frailty->X;
    double * hazYcum = hazard->Ycum;
    //unMOD double * hazYcum = hazard->SplineYcum;

    for(int i=0; i<hazard->nx; i++) lik+=regression->status[i]*log(hazard->Y[i]);
    //unMOD for(int i=0; i<hazard->nx; i++) lik+=regression->status[i]*log(hazard->SplineY[i]);
    for(int i=0; i< regression->n; i++) lik-=frailX[cluster[i]-1] * hazYcum[i] * elp[i];
    lik -= hazard->SplinePenaltyFactor[0]*SmoothnessPenalty(hazard); 
    for(int i=0; i<hazard->nj; i++) lik -= hazard->SplinePar[i]<hazard->SplineMin[0] ? pow(hazard->SplinePar[i] - hazard->SplineMin[i],2) : 0.0;
    return lik;
}

double LikelihoodSplineFrailty(curveP hazard, curveP frailty, regressionP regression)
{
    double lik=0;

    for(int i=0; i<frailty->nx; i++) lik+=log(frailty->Y[i]);
    //unMOD for(int i=0; i<frailty->nx; i++) lik+=log(frailty->SplineY[i]);
    lik -= frailty->SplinePenaltyFactor[0]*SmoothnessPenalty(frailty); 
    lik -= frailty->SplineMeanPenalty[0] * pow(ddotWrapper(frailty->nj, frailty->SplineBasisExp, frailty->SplineEPar),2);
    for(int i=0; i<frailty->nj; i++) lik -= frailty->SplinePar[i]<frailty->SplineMin[0] ? pow(frailty->SplinePar[i] - frailty->SplineMin[i],2) : 0.0;
    return lik;
}

double LikelihoodParamHazard(curveP hazard, curveP frailty, regressionP regression)
{
    double lik=0;
    double * elp = regression->elp;
    int * cluster = regression->cluster;
    double * frailX = frailty->X;
    double * hazYcum = hazard->Ycum;
    //unMOD double * hazYcum = hazard->ParamYcum;
    for(int i=0; i<hazard->nx; i++) lik+=regression->status[i]*log(hazard->Y[i]);
    //unMOD for(int i=0; i<hazard->nx; i++) lik+=regression->status[i]*log(hazard->ParamY[i]);
    for(int i=0; i< regression->n; i++) lik-=frailX[cluster[i]-1] * hazYcum[i] * elp[i];
    int c1=1;
    lik -= pow(F77_CALL(dnrm2)(&(hazard->np), hazard->ParamPar, &c1),2)/(2*hazard->ParamPriorvar[0]);
    return lik;
}

double LikelihoodParamFrailty(curveP hazard, curveP frailty, regressionP regression)
{
    double lik=0;
    for(int i=0; i<frailty->nx; i++) lik+=log(frailty->Y[i]);
    //unMOD for(int i=0; i<frailty->nx; i++) lik+=log(frailty->ParamY[i]);
    int c1=1;
    lik -= pow(F77_CALL(dnrm2)(&(frailty->np), frailty->ParamPar, &c1),2)/(2*frailty->ParamPriorvar[0]);
    return lik;
}

double LikelihoodWeightHazard(curveP hazard, curveP frailty, regressionP regression)
{ 
    double lik=0;
    double * elp = regression->elp;
    int * cluster = regression->cluster;
    double * frailX = frailty->X;
    double * hazYcum = hazard->Ycum;
    for(int i=0; i<hazard->nx; i++) lik+=regression->status[i]*log(hazard->Y[i]);
    for(int i=0; i< regression->n; i++) lik-=frailX[cluster[i]-1] * hazYcum[i] * elp[i];
    //lik -= pow(hazard->Weight[0],2)/(2*hazard->WeightPriorvar[0]);
    lik += (hazard->WeightHyper[0] - 1.0) * log(hazard->Weight[0])
          +(hazard->WeightHyper[1] - 1.0) * log(1.0 - hazard->Weight[0]);
    return lik;
}

double LikelihoodWeightFrailty(curveP hazard, curveP frailty, regressionP regression)
{
    double lik=0;
    for(int i=0; i<frailty->nx; i++) lik+=log(frailty->Y[i]);
    //lik -= pow(frailty->Weight[0],2)/(2*frailty->WeightPriorvar[0]);
    lik += (frailty->WeightHyper[0] - 1.0) * log(frailty->Weight[0])
          +(frailty->WeightHyper[1] - 1.0) * log(1.0 - frailty->Weight[0]);
    return lik;
}

int AcceptReject(double baselik, double candlik, double ratio)
{
    if(isnan(candlik)) candlik = -DBL_MAX;
    double r = exp(candlik - baselik) * ratio; 
    double pAccept = dmin(1,r);
    int out;
    if(isnan(pAccept)) pAccept = 0;
    GetRNGstate();
    if(runif(0,1) < pAccept) out = 1;
    else out = 0;
    PutRNGstate();
    return out;
}

void MH_Frail(curveP hazard, curveP frailty, regressionP regression)
{
    double acc = 0;
    double baselik, candlik;
    for(int i=0; i < regression->m; i++){
        baselik = LikelihoodFrailty(i, hazard, frailty, regression);
        double u = frailty->X[i];
        double y = frailty->Y[i];
        double v = frailty->tun[0];
        GetRNGstate();
        double cand = rgamma( pow(u,2)/v, v/u);
        PutRNGstate();
        double pcu = 1; double puc = 1;
        if(isnan(cand) || ( frailty->hasSpline &&
            (cand > frailty->SplineKnots[frailty->SplineNknots-1]
             | cand < frailty->SplineKnots[0]) ) ){
            candlik = -DBL_MAX;
        }else{
            frailty->X[i] = cand;
            frailty->Y[i] = EvalCurveAtOnePoint(frailty, cand);
            candlik = LikelihoodFrailty(i, hazard, frailty, regression);
            puc = dgamma(cand, pow(u,2)/v, v/u, 0);
            pcu = dgamma(u, pow(cand,2)/v, v/cand, 0);
        }
        int thisacc = AcceptReject(baselik, candlik, pcu/puc);
        if(thisacc==0) { //Did not accept, so undo the damage
            frailty->X[i] = u;
            frailty->Y[i] = y;
        }else{
            UpdateCurveX(frailty, cand, i);
        }
        acc += (double) thisacc;
    }
    acc /= regression->m;
    frailty->Accept[0] = acc;
}

void MH_Regression(curveP hazard, curveP frailty, regressionP regression)
{
    double baselik, candlik;
    baselik = LikelihoodRegression(hazard, frailty, regression); 
    double * cand = (double *) calloc( regression->p, sizeof(double));
    double * oldlp = (double *) malloc( regression->n * sizeof(double));
    double * oldelp = (double *) malloc( regression->n * sizeof(double));
double * oldcoef = (double *) malloc( regression->p * sizeof(double));
    // store the old regression information
    char trans = 'N'; double c0 = 0; int c1 = 1; double c1d = 1;
    dcopyWrapper(regression->p, regression->coefficients, oldcoef);
    dcopyWrapper(regression->n, regression->lp, oldlp);
    dcopyWrapper(regression->n, regression->elp, oldelp);
    //generate candidate parameters
    mvrnorm(regression->p, cand, regression->coefficients, regression->CholCov, regression->tun[0]); 
    // Change the regression object with the new lp and elps
    dcopyWrapper(regression->p, cand, regression->coefficients);
    F77_CALL(dgemv)(&trans, &(regression->n), &(regression->p), &c1d, regression->covariates, &(regression->n), regression->coefficients, &c1, &c0, regression->lp, &c1);
    for(int i=0; i < regression->n; i++) regression->elp[i] = exp(regression->lp[i]);
    candlik = LikelihoodRegression(hazard, frailty, regression);
    int acc = AcceptReject(baselik, candlik, 1);
    if(!acc){
        // not accepted, so restore old regression
        dcopyWrapper(regression->p, oldcoef, regression->coefficients);
        dcopyWrapper(regression->n, oldlp, regression->lp);
        dcopyWrapper(regression->n, oldelp, regression->elp);
    }
    regression->Accept[0] = (double) acc; 
    free(cand);
    free(oldlp);
    free(oldelp);
    free(oldcoef);
}

void MH_SplineHazard(curveP hazard, curveP frailty, regressionP regression)
{
    if(!hazard->hasSpline) return;
    double baselik, candlik;
    baselik = LikelihoodSplineHazard(hazard,frailty,regression); 
    double * cand = (double *) calloc( hazard->nj, sizeof(double));
    // allocate storage for parameters of old spline
    double * oldPar = (double *) malloc( hazard->nj * sizeof(double));
    dcopyWrapper(hazard->nj, hazard->SplinePar, oldPar);
    // create candidate parameters
    mvrnorm(hazard->nj, cand, hazard->SplinePar,hazard->SplineCholCov,hazard->SplineTun[0]);
    UpdateSplinePar(hazard,cand); 
    candlik = LikelihoodSplineHazard(hazard,frailty,regression);
    int acc = AcceptReject(baselik, candlik, 1);
    if(!acc)
        UpdateSplinePar(hazard,oldPar);
    hazard->SplineAccept[0] = (double) acc;
    free(cand);
    free(oldPar);
}

void MH_SplineFrailty(curveP hazard, curveP frailty, regressionP regression)
{
    if(!frailty->hasSpline) return;
    double baselik, candlik;
    baselik = LikelihoodSplineFrailty(hazard,frailty,regression); 
    double * cand = (double *) calloc( frailty->nj, sizeof(double));
    double * candShort = (double *) calloc( (frailty->nj-1), sizeof(double));
    // allocate storage for parameters of old spline
    double * oldPar = (double *) malloc( (frailty->nj) * sizeof(double));
    dcopyWrapper(frailty->nj, frailty->SplinePar, oldPar);
    double * candMean = (double *) malloc ((frailty->nj-1) * sizeof(double));
    // The candidate must include every entry except the fixed index
    int ind=0;
    for(int i = 0; i<frailty->nj; i++){
        if(i != frailty->SplineFixedInd){
            candMean[ind]=frailty->SplinePar[i];
            ind++;
        }
    }
    // create candidate parameters
    mvrnorm(frailty->nj-1, candShort, candMean,frailty->SplineCholCov,frailty->SplineTun[0]);
    // Fix up the candidate by setting the fixed entry to 0
    ind=0;
    for(int i=0; i<frailty->nj; i++){
        if(i != frailty->SplineFixedInd){
            cand[i] = candShort[ind];
            ind++;
        }else{
            cand[i] = 0;
        }
    }
    UpdateSplinePar(frailty,cand); 
    candlik = LikelihoodSplineFrailty(hazard,frailty,regression);
    int acc = AcceptReject(baselik, candlik, 1);
    if(!acc)
        UpdateSplinePar(frailty,oldPar);
    frailty->SplineAccept[0] = (double) acc;
    free(cand);
    free(candShort);
    free(candMean);
    free(oldPar);
}

void MH_ParamHazard(curveP hazard, curveP frailty, regressionP regression)
{
    if(!hazard->hasPar) return;
    double baselik, candlik;
    baselik = LikelihoodParamHazard(hazard,frailty,regression);
    double * cand = (double *) calloc( hazard->np, sizeof(double));
    double * oldPar = (double *) malloc( hazard->np * sizeof(double));
    dcopyWrapper(hazard->np, hazard->ParamPar, oldPar);
    mvrnorm(hazard->np, cand, hazard->ParamPar, hazard->ParamCholCov, hazard->ParamTun[0]);
    UpdateParamPar(hazard,cand);
    candlik = LikelihoodParamHazard(hazard,frailty,regression);
    int acc = AcceptReject(baselik, candlik, 1);
    if(!acc)
        UpdateParamPar(hazard,oldPar);
    hazard->ParamAccept[0] = (double) acc;
    free(cand);
    free(oldPar);
}

void MH_ParamFrailty(curveP hazard, curveP frailty, regressionP regression)
{ 
    if(!frailty->hasPar) return;
    double baselik, candlik;
    baselik = LikelihoodParamFrailty(hazard,frailty,regression);
    double * cand = (double *) calloc( frailty->np, sizeof(double));
    double * oldPar = (double *) malloc( frailty->np * sizeof(double));
    dcopyWrapper(frailty->np, frailty->ParamPar, oldPar);
    mvrnorm(frailty->np, cand, frailty->ParamPar, frailty->ParamCholCov, frailty->ParamTun[0]);
    UpdateParamPar(frailty,cand);
    candlik = LikelihoodParamFrailty(hazard,frailty,regression);
    int acc = AcceptReject(baselik, candlik, 1);
    if(!acc)
        UpdateParamPar(frailty,oldPar);
    frailty->ParamAccept[0] = (double) acc;
    free(cand);
    free(oldPar);
}

void MH_Weight(curveP theCurve, curveP hazard, curveP frailty, regressionP regression)
{
    if(!theCurve->hasPar | !theCurve->hasSpline) return;
    double (*likfun)(curveP, curveP, regressionP);
    if(theCurve->isHazard) likfun = &LikelihoodWeightHazard;
    else likfun = &LikelihoodWeightFrailty;  
    double oldW = theCurve->Weight[0];
    double w = dmin(dmax(oldW, .01), 0.99);
    double v = theCurve->WeightTun[0];
    double alpha = w*(w*(1-w)/v-1);
    double beta = (1-w)/w*alpha;
    GetRNGstate();
    double cand = rbeta(alpha, beta);
    PutRNGstate();
    if(isnan(cand)){
        theCurve->WeightAccept[0]=0;
        return;
    }
    double alphac = cand*(cand*(1-cand)/v-1);
    double betac = (1-cand)/cand*alphac;
    double baselik = likfun(hazard, frailty, regression);
    theCurve->Weight[0] = cand;
    ReweightCurve(theCurve, -1);
    double candlik = likfun(hazard, frailty, regression);
    double puc = dbeta(cand, alpha, beta, 0);
    double pcu = dbeta(w, alphac, betac, 0);
    int acc = AcceptReject(baselik, candlik, pcu/puc);
    if(!acc){
        theCurve->Weight[0] = oldW;
        ReweightCurve(theCurve, -1);
    }
    theCurve->WeightAccept[0] = (double) acc;
}

void UpdatePostvarCurve(curveP theCurve)
{
    int c1 = 1;
    if(theCurve->hasSpline) theCurve->SplinePriorvar[0] = rinvgamma(
            (double) theCurve->nj / 2 + theCurve->SplineHyper[0],
            theCurve->SplinePenaltyFactor[0] * SmoothnessPenalty(theCurve) * theCurve->SplinePriorvar[0] + theCurve->SplineHyper[1] );

    if(theCurve->hasPar) theCurve->ParamPriorvar[0] = rinvgamma(
            (double) theCurve->np / 2 + theCurve->ParamHyper[0],
            pow(F77_CALL(dnrm2)(&(theCurve->np), theCurve->ParamPar, &c1),2)/2+ theCurve->ParamHyper[1] );

   // if(theCurve->hasPar & theCurve->hasSpline) theCurve->WeightPriorvar[0] = rinvgamma(
   //         0.5 + theCurve->WeightHyper[0],
   //         pow(theCurve->Weight[0],2)/2 + theCurve->WeightHyper[1] );
}

void UpdatePostvarRegression(regressionP theReg)
{
    int c1 = 1;
    theReg->priorvar[0] = rinvgamma(
            (double) theReg->p / 2 + theReg->hyper[0],
            pow(F77_CALL(dnrm2)(&(theReg->p), theReg->coefficients, &c1),2)/2+ theReg->hyper[1] );
}

void UpdateHistory(curveP hazard, curveP frailty, regressionP regression, historyP history, int iter)
{
    int c1=1;
    int ny = history->ny;
    F77_CALL(dcopy)(&(frailty->nx), frailty->X, &c1, history->frailty + iter-1, &(history->ny));
    F77_CALL(dcopy)(&(regression->p), regression->coefficients, &c1, history->coefficients + iter-1, &(history->ny));
    if(frailty->hasSpline)
        F77_CALL(dcopy)(&(frailty->nj), frailty->SplinePar, &c1, history->FrailtySplinePar + iter-1, &(history->ny));
    if(hazard->hasSpline)
        F77_CALL(dcopy)(&(hazard->nj), hazard->SplinePar, &c1, history->HazardSplinePar + iter-1, &(history->ny));
    if(frailty->hasPar)
        F77_CALL(dcopy)(&(frailty->np), frailty->ParamPar, &c1, history->FrailtyParamPar + iter-1, &(history->ny));
    if(hazard->hasPar)
        F77_CALL(dcopy)(&(hazard->np), hazard->ParamPar, &c1, history->HazardParamPar + iter-1, &(history->ny));
    if(hazard->hasPar & hazard->hasSpline)
        history->HazardWeight[iter-1] = hazard->Weight[0];
    if(frailty->hasPar & frailty->hasSpline)
        history->FrailtyWeight[iter-1] = frailty->Weight[0];

    history->priorvar[iter-1 + ny*0] = regression->priorvar[0];
    history->priorvar[iter-1 + ny*1] = hazard->SplinePriorvar[0];
    history->priorvar[iter-1 + ny*2] = frailty->SplinePriorvar[0];
    history->priorvar[iter-1 + ny*3] = hazard->ParamPriorvar[0];
    history->priorvar[iter-1 + ny*4] = frailty->ParamPriorvar[0];
    history->priorvar[iter-1 + ny*5] = hazard->WeightPriorvar[0];
    history->priorvar[iter-1 + ny*6] = frailty->WeightPriorvar[0];

    history->accept[iter-1 + ny*0] = regression->Accept[0];
    history->accept[iter-1 + ny*1] = hazard->SplineAccept[0];
    history->accept[iter-1 + ny*2] = frailty->SplineAccept[0];
    history->accept[iter-1 + ny*3] = hazard->ParamAccept[0];
    history->accept[iter-1 + ny*4] = frailty->ParamAccept[0];
    history->accept[iter-1 + ny*5] = hazard->WeightAccept[0];
    history->accept[iter-1 + ny*6] = frailty->WeightAccept[0];
    history->accept[iter-1 + ny*7] = frailty->Accept[0];
}


SEXP SplineSurvMainLoop( SEXP Rhazard, SEXP Rfrailty, SEXP Rregression, 
        SEXP Rhistory, SEXP Rstartiter, SEXP Renditer, SEXP Rverbose)
{
    //Create locally stored curve structures
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
/*#ifdef DEBUG
    Rprintf("HAZARD:\n");
    if(hazard->hasSpline){
        Rprintf("HasSpline: %i\n",hazard->hasSpline);
        Rprintf("Nknots: %i\n",hazard->SplineNknots);
        Rprintf("Ord: %i\n",hazard->SplineOrd);
        Rprintf("Knots: "); for(int i=0;i<hazard->SplineNknots;i++) Rprintf("%f ",hazard->SplineKnots[i]);Rprintf("\n");
        Rprintf("PenaltyType:");  Rprintf("%d\n",hazard->SplinePenaltyType);
        Rprintf("CandCov:"); for(int i=0; i<41; i++) Rprintf("%f ",hazard->SplineCandCov[i]); Rprintf("\n");
        Rprintf("CandChol:"); for(int i=0; i<41; i++) Rprintf("%f ",hazard->SplineCholCov[i]); Rprintf("\n");
    }
    if(hazard->hasPar){
        Rprintf("ParamPar: %f\n", hazard->ParamPar[0]);
    }
    Rprintf("\nX: "); for(int i=0;i<10;i++) Rprintf("%f ",hazard->X[i]);
    Rprintf("\nY: "); for(int i=0;i<10;i++) Rprintf("%f ",hazard->Y[i]);
#endif */

    // Begin main loop
    int iter = asInteger(Rstartiter);
    int enditer = asInteger(Renditer);
    int verbose = asInteger(Rverbose);
    while(iter < enditer)
    {
        iter++;
        if(verbose >= 3) Rprintf("%d ", iter);
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

        UpdateHistory(hazard, frailty, regression, history, iter);
        
    }
    SEXP returnval = R_NilValue;
    FreeCurveMem(hazard);
    FreeCurveMem(frailty);
    FreeRegMem(regression);
    return returnval;
}

SEXP mychol( SEXP x){
    double * xmat = REAL(x);
    int n = (int) sqrt(length(x));
    int job=1;
    int * pivot = (int *) malloc( n * sizeof(int));
    int info=0;
    int c1=1;
    char uplo='u';
    double * work = (double *) malloc( n * sizeof(double));
    F77_CALL(dchdc)(xmat, &n, &n, work, pivot, &job, &info);
//    F77_CALL(dpotrf)(&uplo, &n, xmat, &n, &info);
    for(int i=1; i< n; i++)
        for(int j=0; j<i; j++)
            xmat[i + n * j]=0;
    free(work);
    free(pivot);
    SEXP returnval = R_NilValue;
    return returnval;
}
