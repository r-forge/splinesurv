/****h* initRoutine/CinitRoutine
 *  NAME
 *    CinitRoutine --- additional C functions used for initialization only
 *  FUNCTION
 *    Set of functions, mainly likelihood and gradient functions called during 
 *    initialization of the splinesurv.agdata function.
 *  CONTENTS
 *    addInitSmoothnessPenaltyGr --- adds gradient of smoothness penalty to an input vector
 *    cInitGrHazSpline --- gradient of likelihood of spline parameters for hazard
 *    cInitLikFrailSpline --- loglikelihood of parameters for spline component of frailty curve
 *    cInitLikHazSpline --- spline hazard likelihood for initialization
 *    cInitMultAndWeight --- multiply a spline component and add it to a weighted parametric component
 *    cMakePenalty2diff --- compute a penalty matrix on second derivatives
 *    InitSmoothnessPenalty --- compute the smoothness penalty during initialization
 *********/

#include <R_ext/BLAS.h>    
#include <R_ext/Print.h>    
#include <math.h>
#include <stdlib.h>
#include "init.h"

/****f* CsplineUtils/cnBsmom
 *  NAME
 *    cnBsmom --- compute the N-h moment of a B-spline basis function
 *  FUNCTION
 *   Computes the N-th moment of the j-th B-spline basis function of order ord defined
 *   on the given set of knots. See also nBsmom.
 *  INPUTS
 *    N       moment to compute
 *    ord     order of the spline
 *    j       which basis function to compute the moment for
 *    knots   set of knots to use
 *    knotord order for which the knots were constructed
 *  OUTPUTS
 *  SYNOPSIS
 */
double cnBsmom(const int N, const int ord, const int j, const double *knots, const int knotord)
/*
 *  SOURCE
*/
{
    double lknot=knots[j-ord+knotord];
    double rknot=knots[j+knotord];
    if(lknot==rknot){
        return pow(rknot,N);
    }
    double out;
    // magic recursive formula
    out = ord/(rknot-lknot)/(N+1);
    out = out*(-cnBsmom(N+1,ord-1,j-1,knots,knotord)+cnBsmom(N+1,ord-1,j,knots,knotord));
    return out;
}
/************ cnBsmom */

/****f* CsplineUtils/cevalEinte
 *  NAME
 *    cevalEinte --- compute the N-th moment of the distance from 1 of a B-spline
 *  FUNCTION
 *    C implementation of evalEinte, see its description for details. The difference is
 *    that this function can compute general moments, not just the first moment. The second
 *    moment is needed for computing the frailty variance.
 *  INPUTS
 *    einte     vector of length K+ord to store the output
 *    knots     vector of knots, of length K+2*ord
 *    ord       order of the spline
 *    K         number of interior knots
 *    N         moment to compute
 *  SYNOPSIS
 */
void cevalEinte(double *einte, double *knots, int *ord, int *K, int *N)
/*
 *  SOURCE
*/
{
    int j;
    for(j=0; j < *K+*ord; j++)
        einte[j]=1.0-cnBsmom(*N,*ord,j,knots,*ord);
}
/************ cevalEinte */

/****f* CsplineUtils/cevalBinte
 *  NAME
 *    cevalBinte --- compute the integrals of each spline basis function
 *  FUNCTION
 *    C re-implmentation of evalBinte function. See its description for detail.
 *  INPUTS
 *    binte     vector of length K+ord to store output
 *    knots     vector of knots of length K+2*ord
 *    ord       order of the spline
 *    K         number of interior knots
 *  SYNOPSIS
 */
void cevalBinte(double *binte, double *knots, int *ord, int *K)
/*
 *  SOURCE
*/
{
    for(int j=0; j<*K+*ord; j++)
        binte[j]=(knots[j+*ord]-knots[j]) / *ord;
}
/************ cevalBinte */


/****f* CsplineUtils/csplineeval
 *  NAME
 *    csplineeval --- evaluate a B-spline function
 *  FUNCTION
 *    Evaluates a B-spline function at a given value. This is the lowest-level
 *    B-spline evaluator in this code. It uses the simple recursive formulation of
 *    a B-spline basis function.
 *
 *    Given a set of knots, this routine evaluates at x the basis function of order ord
 *    that has support on knots [j-ord,j]. Note that because negative array indices are not
 *    allowed, knot -ord is stored in knots[0] here.
 *  INPUTS
 *    x     the value at which it should be evaluated
 *    j     which of the basis functions should be evaluated
 *    ord   order of the spline
 *    knots knots on which the spline is defined
 *    splord order of the spline, but remains unchanged throughout evaluation
 *    nj    number of basis functions defined on the knots (length(knots)-splord)
 *  OUTPUTS
 *    The j-th basis function evaluated at x: B_j(x)
 *  SYNOPSIS
 */
double csplineeval( double x, int j, int ord, double *knots, int splord, int nj)
/*
 *  SOURCE
*/
{
    double knotj=knots[j+splord];
    double knotjm1=knots[j-1+splord];
    // base case for ord=1
    if(ord==1) {
        if((x >= knotjm1) & (x < knotj)) return 1.0;
        if((x==knotj) & (j+splord==nj)) return 1.0;
        return 0.0;
    }
    // recursive formula otherwise
    double out = 0;
    double knotjmn=knots[j-ord+splord];
    double denom1 = (knotjm1 - knotjmn);
    if(denom1>0){
        double numer1 = (x-knotjmn);
        double rec1 = csplineeval(x,j-1,ord-1,knots,splord,nj);
        out += numer1/denom1*rec1;
    }
    double knotjmn1=knots[j-ord+1+splord];
    double denom2 = (knotj - knotjmn1);
    if(denom2>0){
        double numer2 = (knotj-x);
        double rec2 = csplineeval(x,j,ord-1,knots,splord,nj);
        out += numer2/denom2*rec2;
    }
    return out;
}
/************ csplineeval */

/****f* CsplineUtils/csplinedesign
 *  NAME
 *    csplinedesign --- create a B-spline basis
 *  FUNCTION
 *    Compute a full B-spline basis given a set of knots and a set of points at which to 
 *    evaluate the basis functions.
 *    Plays the role of spline_basis in the splines package.
 *  INPUTS
 *    des       matrix of size nx x (K+ord) to hold the output
 *    x         set of points at which to evaluate the basis functions
 *    nx        length of x
 *    knots     vector of knot points, with appropriate boundary knots
 *    ord       order of the spline
 *    K         number of interior knots
 *  OUTPUTS
 *    des[i,j] contains the j-th B-spline basis function evaluated at x[i]
 *  SYNOPSIS
 */
void csplinedesign(double *des, double *x, int *nx, double *knots, int *ord, int *K)
/*
 *  SOURCE
*/
{
    int nj= *K + *ord;
    for(int i=0;i<*nx;i++){
        for(int j=0; j< nj; j++)
            des[i + j* *nx] = csplineeval( x[i], j, *ord, knots, *ord, nj);
    }
}
/************ csplinedesign */

/****f* CsplineUtils/cevalCinteOld
 *  NAME
 *    cevalCinteOld --- old method for computing cevalCinte
 *  FUNCTION
 *    See evalCinte. This function works analogously to the R implementation, but
 *    was too slow, so it was replaced by cevalCinte2.
 *  SYNOPSIS
 */
void cevalCinteOld(double *cinte, double *x, int *nx, double *knots, int *ord, int *K, double *binte)
/*
 *  SOURCE
*/
{
    int nj = *K + *ord;
    int nk = nj + *ord;
    double * knots2 = (double *) malloc((nk +2) * sizeof(double));
    knots2[0]=knots[0];
    for(int j=0; j<nk; j++)
        knots2[j+1]=knots[j];
    knots2[nk+1]=knots[nk-1];
#ifdef DEBUGCINTE
    Rprintf("knots2:\n ");
    for(int j=0; j<nk+2; j++) Rprintf("%f ",knots2[j]);
    Rprintf("\n");
#endif
    int ord2 = *ord+1;
    double * bs2 = (double *) malloc((nj+1)* *nx * sizeof(double));
    csplinedesign(bs2,x,nx,knots2,&ord2,K);
#ifdef  DEBUGCINTE
    Rprintf("bs2:\n");
    for(int j=0; j<nj+1; j++) Rprintf("%f ",bs2[j * *nx]);
    Rprintf("\n");
    for(int j=0; j<nj+1; j++) Rprintf("%f ",bs2[1+ j * *nx]);
    Rprintf("\n");
#endif
    for(int i=0; i< *nx; i++){
        for(int j=0; j<nj; j++){
            cinte[i+j* *nx]=0;
            if(x[i]>=knots[j+ *ord]) cinte[i + j* *nx] = binte[j];
            if((x[i]<knots[j+ *ord]) & (x[i]>=knots[j]))
            {
                for(int k=j+1;k<nj+1;k++) 
                    cinte[i + j* *nx]+=binte[j]*bs2[i + k* *nx];
            }
        }
    }
    free(bs2);
    free(knots2);
}
/************ cevalCinteOld */


/****f* CsplineUtils/cevalCinte2
 *  NAME
 *    cevalCinte2 --- compute partial integrals of B-spline basis functions
 *  FUNCTION
 *    This is needed to compute the cumulative baseline hazard. This function updates the i-th
 *    row of the input matrix cinte, from jstart to jstop, so that cinte[i,j] contains the
 *    integral of the j-th B-spline basis function of order ord, defined on the given knots,
 *    from 0 to x[i].
 *  INPUTS
 *    cinte     output matrix of dimension nx x nj
 *    x         vector of observations at which to evaluate
 *    nx        length of x
 *    knots     vector of knot positions, length nj+ord
 *    ord       order of the spline
 *    binte     vector of integrals of each B-spline basis function, see cevalBinte
 *    i         row of cinte to update
 *    jstart    starting column of cinte to update
 *    jstop     ending column of cinte to update
 *  SYNOPSIS
 */
void cevalCinte2(double *cinte, double *x, int nx, double *knots, int nj, int ord, double *binte,
        int i, int jstart, int jstop)
/*
 *  SOURCE
*/
{
    int nk = nj + ord;
    double * knots2 = (double *) malloc((nk +2) * sizeof(double));
    knots2[0]=knots[0];
    for(int j=0; j<nk; j++) knots2[j+1]=knots[j];
    knots2[nk+1]=knots[nk-1];
    int ord2 = ord+1;
    double * bs2 = (double *) malloc((nj+1)*  sizeof(double));
    // compute a basis of order ord+1
    for(int j=jstart;j<nj+1;j++) bs2[j] = csplineeval(x[i],j,ord2,knots2,ord2,nj);
    // compute as in evalCinte for the row and set of columns given
    for(int j=jstart; j<jstop; j++){
        cinte[i+j* nx]=0;
        if(x[i]>=knots[j+ord]) cinte[i + j*nx] = binte[j];
        if((x[i]<knots[j+ord]) & (x[i]>=knots[j]))
        {
            for(int k=j+1;k<nj+1;k++) 
                cinte[i + j*nx]+=binte[j]*bs2[k];
        }
    }
    free(bs2);
    free(knots2);

}
/************ cevalCinte2 */

/****f* CsplineUtils/cevalCinte
 *  NAME
 *    cevalCinte --- construct matrix of partial integrals of a B-spline basis
 *  FUNCTION
 *    This is needed to compute the cumulative baseline hazard. This function updates the input matrix 
 *    cinte, so that cinte[i,j] contains the
 *    integral of the j-th B-spline basis function of order ord, defined on the given knots,
 *    from 0 to x[i]. It's just a wrapper for cevalCinte2 that evaluates the entire matrix.
 *  INPUTS
 *    cinte     output matrix of dimension nx x nj
 *    x         vector of observations at which to evaluate
 *    nx        length of x
 *    knots     vector of knot positions, length nj+ord
 *    ord       order of the spline
 *    binte     vector of integrals of each B-spline basis function, see cevalBinte
 *  SYNOPSIS
 */
void cevalCinte(double *cinte, double *x, int *nx, double *knots, int *ord, int *K, double *binte)
/*
 *  SOURCE
*/
{
    int nj = *K + *ord;
    for(int i=0; i<*nx; i++) cevalCinte2(cinte, x, *nx, knots, nj, *ord, binte, i, 0, nj);
}
/************ cevalCinte */

/****f* CsplineUtils/cSplineConvolution
 *  NAME
 *    cSplineConvolution --- compute the convolution of two basis functions
 *  FUNCTION
 *    This function computes the convolution of two spline basis functions, that is,
 *       int_0^infty ( x^k B_1(x) * B_2(x) ) dx
 *    where the splines may be of different orders, but are defined on the same set of knots.
 *    See also splineconv.
 *  INPUTS
 *    k      the power of x used in the convolution
 *    n1     order of the first spline
 *    j1     largest knot of the first spline
 *    n2     order of the second spline
 *    j2     largest knot of the second spline
 *    knots  set of knots on which the splines are defined
 *    splord order of the splines
 *  OUTPUTS
 *    The k-th order convolution of the splines defined by the input parameters
 *  SYNOPSIS
 */
double cSplineConvolution(int k, int n1,int j1, int n2, int j2, double *knots, int splord)
/*
 *  SOURCE
*/
{
    double out=0;
    // if the splines don't overlap, the convolution is 0
    if((j1 -n1>=j2) | (j2 -n2>=j1)) return out;
    // if both splines are first-order, the convolution is trivial
    if((n1==1) & (n2==1)){
        out= 1.0/(k+1.0)*(pow(knots[j1+splord],k+1.0)-pow(knots[j1-1+splord],k+1.0));
        return out;
    }
    // assume that n1>n2 wlog, if this is not true, switch the indices
    if(n2>n1){
        int n3=n1; n1=n2; n2=n3;
        int j3=j1; j1=j2; j2=j3;
    }
    // use a magic recursive formula!
    out=0;
    double denom1 = knots[j1-1+splord]-knots[j1-n1+splord];
    double denom2 = knots[j1+splord] - knots[j1-n1+1+splord];
    if(denom1>0){
        out += 1.0/denom1 * cSplineConvolution(k+1,n1-1,j1-1,n2,j2,knots,splord);
        out -= knots[j1-n1+splord]/denom1 * cSplineConvolution(k,n1-1,j1-1,n2,j2,knots,splord);
    }
    if(denom2>0){
        out += knots[j1+splord]/denom2* cSplineConvolution(k,n1-1,j1,n2,j2,knots,splord);
        out -= 1.0/denom2 * cSplineConvolution(k+1,n1-1,j1,n2,j2,knots, splord);
    }
    return out;

}
/************ cSplineConvolution */

/****f* CsplineUtils/cSplineDerivInt
 *  NAME
 *    cSplineDerivInt --- compute the convolution of the derivatives of two spline basis functions
 *  FUNCTION
 *    Used to compute the penalty matrix for the integrated squared second derivative. This
 *    routine computes the integral from 0 to infinity of the l1 derivative and the l2
 *    derivative of the j1 and j2-th splines of order n1 and n2 defined on a set of knots.
 *
 *    See also splinederivint for the R implementation.
 *  INPUTS
 *    l1     derivative of the first spline
 *    n1     order of the first spline
 *    j1     index of the first spline
 *    l2     derivative of the second spline
 *    n2     order of the second spline
 *    j2     index of the second spline
 *    knots  set of knots on which the splines are defined
 *    splord order of the splines
 *  OUTPUTS
 *  SYNOPSIS
 */
double cSplineDerivInt(int l1, int n1, int j1, int l2, int n2, int j2, double *knots, int splord)
/*
 *  SOURCE
*/
{
    if((j1-n1>=j2) | (j2-n2>=j1)) return 0;
    if((l1==0) & (l2==0)) return cSplineConvolution(0,n2,j1,n2,j2,knots,splord);
    if(l2>l1){
        int l3=l1; l1=l2; l2=l3;
        int n3=n1; n1=n2; n2=n3;
        int j3=j1; j1=j2; j2=j3;
    }
    double out=0;
    double denom1 = knots[j1-1+splord]-knots[j1-n1+splord];
    double denom2 = knots[j1+splord] - knots[j1-n1+1+splord];
    if(denom1>0)
        out += (n1-1.0)/denom1 * cSplineDerivInt(l1-1,n1-1,j1-1,l2,n2,j2,knots,splord);
    if(denom2>0)
        out -= (n1-1.0)/denom2 * cSplineDerivInt(l1-1,n1-1,j1,l2,n2,j2,knots,splord);
    return out;
}
/************ cSplineDerivInt */

/****f* CinitRoutine/cMakePenalty2diff
 *  NAME
 *    cMakePenalty2diff --- compute a penalty matrix on second derivatives
 *  FUNCTION
 *    Generates a penalty matrix that penalizes the integral of the squared second derivative
 *    of a curve. See also makePenalty.2deriv.
 *    
 *    This function should be called cMakePenalty2deriv, but I'm scared of changing its name.
 *  INPUTS
 *    P     a matrix of dimension (K+ord) x (K+ord)
 *    knots set of knots on which the spline is defined
 *    ord   order of the spline
 *    K     number of interior knots
 *  OUTPUTS
 *  SYNOPSIS
 */
void cMakePenalty2diff(double *P, double *knots, int *ord, int *K)
/*
 *  SOURCE
*/
{
    // P should be (K+ord)^2
    int nj = *K + *ord;
    for(int j1=0; j1<nj; j1++)
    {
        for(int j2=j1; j2<nj; j2++){
            P[j1 + j2*nj]=cSplineDerivInt(2, *ord, j1, 2, *ord, j2, knots, *ord);
            P[j2 + j1*nj]=P[j1 + j2*nj];
        }
    }
}
/************ cMakePenalty2diff */

/****f* CinitRoutine/cInitMultAndWeight
 *  NAME
 *    cInitMultAndWeight --- multiply a spline component and add it to a weighted parametric component
 *  FUNCTION
 *    Given a B-spline basis and a set of parameters, this curve multiplies the basis by the parameters
 *    to produce an estimate of the spline component, then weights it and adds it to an existing
 *    parametric component. It's inelegant, but fast for initialization.
 *  INPUTS
 *    y         parametric curve component, on output holds the total curve
 *    B         B-spline basis
 *    par       weights of each basis element
 *    weight    weight of the spline component
 *    ny        length of y
 *    nj        number of basis functions
 *  SYNOPSIS
 */
void cInitMultAndWeight( double *y,  double *B,  double *par,  double *weight,  int  *ny,  int *nj)
/*
 *  SOURCE
*/
{
    // compute w*B%*%exp(par) + (1-w)y
    const int c1 = 1;
    const char trans = 'N';
    const double onemw=1-*weight;
    F77_CALL(dgemv)(&trans,ny,nj,weight,B,ny,par,&c1,&onemw,y,&c1);
}
/************ cInitMultAndWeight */

/****f* CinitRoutine/InitSmoothnessPenalty
 *  NAME
 *    InitSmoothnessPenalty --- compute the smoothness penalty during initialization
 *  FUNCTION
 *    Special function to compute the smoothness penalty during initialization. This differs from
 *    SmoothnessPenalty in that this routine does not take a CCurve as input, rather only the
 *    relevant components, and can thus be called directly from R.
 *  INPUTS
 *    pen       output storage
 *    par       parameters of the spline
 *    P         penalty matrix
 *    penaltyType   integer indexing the type of penalty (0=gaussian, 1=2diff, 2=2deriv, 3=log2deriv)
 *    sigma2    prior variance
 *    nj        number of B-spline basis functions
 *  SYNOPSIS
 */
void InitSmoothnessPenalty(double *pen, double *par, double *P, int *penaltyType, double *sigma2, int *nj)
/*
 *  SOURCE
*/
{
    if(*penaltyType==0){  // Gaussian
        int c1 = 1;
        *pen = pow(F77_CALL(dnrm2)(nj, par, &c1),2.0);
        *pen = *pen / (2 * *sigma2);
    }
    // second differences or second derivative
    if(*penaltyType==1 || *penaltyType==2 || *penaltyType==3){   
        double * temp = (double *) calloc((*nj) , sizeof(double));
        int c1 = 1;
        double c1d = 1.0;
        char uplo = 'u';
        // compute par %*% P %*% par
        F77_CALL(dsymv)( &uplo, nj, &c1d, P, nj, par, &c1, &c1d, temp, &c1);
        *pen = F77_CALL(ddot)( nj, temp, &c1, par, &c1);
        if(*penaltyType==3) *pen=log(*pen+1);
        *pen = *pen / (2 * *sigma2);
        free( temp );
    }
}
/************ InitSmoothnessPenalty */

/****f* CinitRoutine/addInitSmoothnessPenaltyGr
 *  NAME
 *    addInitSmoothnessPenaltyGr --- adds gradient of smoothness penalty to an input vector
 *  FUNCTION
 *    Computes the gradient of the smoothness penaly over the spline parameters and adds it to
 *    an input vector.
 *  INPUTS
 *    gr        gradient vector
 *    penpar    parameters of the spline
 *    P         penalty matrix
 *    penaltyType   integer indexing the type of penalty (0=gaussian, 1=2diff, 2=2deriv, 3=log2deriv)
 *    sigma2    prior variance
 *    nj        number of B-spline basis functions
 *  OUTPUTS
 *  SYNOPSIS
 */
void addInitSmoothnessPenaltyGr(double *gr, double * penpar, double *P, int *penaltyType,
        double * sigma2, int *nj)
/*
 *  SOURCE
*/
{
    int c1 = 1;
    if(*penaltyType==0){
        double scalefactor = -1 / *sigma2;
        F77_CALL(daxpy)(nj, &scalefactor, penpar, &c1, gr, &c1);
    }
    // second differences or second derivative
    if(*penaltyType==1 || *penaltyType==2 || *penaltyType==3){   
        double * temp = (double *) calloc((*nj) , sizeof(double));
        double c1d = 1.0;
        char uplo = 'u';
        // compute smoothness penaly gradient.
        F77_CALL(dsymv)( &uplo, nj, &c1d, P, nj, penpar, &c1, &c1d, temp, &c1);
        for(int i=0; i<*nj; i++) temp[i]=temp[i]*penpar[i] / *sigma2;
        double scalefactor=-1;
        if(*penaltyType==3) {
            double pen=0;
            int c2 = 2;
            InitSmoothnessPenalty( &pen, penpar, P, &c2, sigma2, nj);
            scalefactor= -1/((2* *sigma2) * pen + 1);
        }
        // add to gradient
        F77_CALL(daxpy)(nj, &scalefactor, temp, &c1, gr, &c1);
        free( temp );
    }
}
/************ addInitSmoothnessPenaltyGr */

/****f* CinitRoutine/cInitLikHazSpline
 *  NAME
 *    cInitLikHazSpline --- spline hazard likelihood for initialization
 *  FUNCTION
 *    Compute the loglikelihood of parameters of the hazard spline components. Called by
 *    cmklik.spline.haz during initialization.
 *  INPUTS
 *    see cmklik.spline.haz for inputs and outputs
 *  SYNOPSIS
 */
void cInitLikHazSpline(double *lik, double *par, double *status, double *lp, double *frailrep,
    double *hazParY, double *hazParYcum, double *weight, double *B, double *C,
    double *P, int *penaltyType, double *sigma2, int *ny, int *nj)
/*
 *  SOURCE
*/
{
    double * epar = (double *) malloc((*nj) * sizeof(double));
    double * hazY = (double *) malloc((*ny) * sizeof(double));
    double * hazYcum = (double *) malloc((*ny) * sizeof(double));
    const int c1 = 1;
    F77_CALL(dcopy)(ny, hazParY, &c1, hazY, &c1);
    F77_CALL(dcopy)(ny, hazParYcum, &c1, hazYcum, &c1);
    for(int i=0; i< *nj; i++) epar[i]=exp(par[i]);
    cInitMultAndWeight( hazY, B, epar, weight, ny, nj);
    cInitMultAndWeight( hazYcum, C, epar, weight, ny, nj);
    double out=0;
    for(int i=0; i< *ny; i++) hazY[i]=log(hazY[i]);
    out += F77_CALL(ddot)(ny,status,&c1,hazY,&c1);
    for(int i=0; i< *ny; i++)
        out -= frailrep[i]*hazYcum[i]*exp(lp[i]);
    double pen = 0;
    double *penpar = ((*penaltyType==2) | (*penaltyType==3))  ? epar : par;
    InitSmoothnessPenalty(&pen, penpar, P, penaltyType, sigma2, nj);
    out-=pen;
    *lik = out;
    free(epar);
    free(hazY);
    free(hazYcum);
}
/************ cInitLikHazSpline */

/****f* CinitRoutine/cInitGrHazSpline
 *  NAME
 *    cInitGrHazSpline --- gradient of likelihood of spline parameters for hazard
 *  FUNCTION
 *    Compute the gradient of the loglikelihood of parameters of the hazard spline component,
 *    called by cmkgr.spline.haz during initialization only.
 *  INPUTS
 *   see cmkgr.spline.haz for inputs and outputs.
 *  SYNOPSIS
 */
void cInitGrHazSpline(double *gr, double *par, double *status, double *lp, double *frailrep,
    double *hazParY, double *hazParYcum, double *weight, double *B, double *C,
    double *P, int *penaltyType, double *sigma2, int *ny, int *nj)
/*
 *  SOURCE
*/
{
    double * epar = (double *) malloc((*nj) * sizeof(double));
    double * hazY = (double *) malloc((*ny) * sizeof(double));
    double * temp = (double *) malloc((*ny) * sizeof(double));
    const int c1 = 1;
    F77_CALL(dcopy)(ny, hazParY, &c1, hazY, &c1);
    for(int i=0; i< *nj; i++) epar[i]=exp(par[i]);
    const double c1d=1;
    const double cm1d=-1;
    cInitMultAndWeight( hazY, B, epar, weight, ny, nj);
    const char trans = 'T';
    for(int i=0; i<*ny; i++) hazY[i]=status[i]/hazY[i]; 
    for(int i=0; i<*ny; i++) temp[i]=exp(lp[i])*frailrep[i]; 
    F77_CALL(dgemv)(&trans,ny,nj,&c1d,B,ny,hazY,&c1,&c1d,gr,&c1);
    F77_CALL(dgemv)(&trans,ny,nj,&cm1d,C,ny,temp,&c1,&c1d,gr,&c1);
    for(int i=0; i<*nj; i++) gr[i]=gr[i]*epar[i] * *weight;
    double *penpar = ((*penaltyType==2) | (*penaltyType==3))  ? epar : par;
    addInitSmoothnessPenaltyGr(gr, penpar, P, penaltyType, sigma2, nj);
    free(epar);
    free(hazY);
    free(temp);
}
/************ cInitGrHazSpline */

/****f* CinitRoutine/cInitLikFrailSpline
 *  NAME
 *    cInitLikFrailSpline --- loglikelihood of parameters for spline component of frailty curve
 *  FUNCTION
 *    Compute the loglikelihood of parameters of the frailty spline component,
 *    called by cmklik.spline.frail during initialization only.
 *  INPUTS
 *    see cmklik.spline.frail
 *  SYNOPSIS
 */
void cInitLikFrailSpline(double *lik, double *par, double *frailParY, double *weight, double *B,
        double *E, double *M, double *P, int * penaltyType, double *sigma2, int *ny, int *nj)
/*
 *  SOURCE
*/
{
    double out=0;
    int c1=1;
    double * epar = (double *) malloc((*nj) * sizeof(double));
    double * eparnorm = (double *) malloc((*nj) * sizeof(double));
    double * frailY = (double *) malloc((*ny) * sizeof(double));
    F77_CALL(dcopy)(ny, frailParY, &c1, frailY, &c1);
    for(int i=0; i< *nj; i++) epar[i]=exp(par[i]);
    double  eparsum = 0;
    for(int i=0; i< *nj; i++) eparsum+=epar[i];
    for(int i=0; i< *nj; i++) eparnorm[i]=epar[i]/eparsum;
    cInitMultAndWeight( frailY, B, eparnorm, weight, ny, nj);
    for(int i=0; i<*ny; i++) out+=log(frailY[i]); 
    double pen = 0;
    double *penpar = ((*penaltyType==2) | (*penaltyType==3))  ? epar : par;
    InitSmoothnessPenalty(&pen, penpar, P, penaltyType, sigma2, nj);
    out-=pen;
    out -= *M * pow((F77_CALL(ddot)(nj, E, &c1, epar, &c1)),2);
    *lik = out;
    free(epar);
    free(eparnorm);
    free(frailY);
}
/************ cInitLikFrailSpline */
