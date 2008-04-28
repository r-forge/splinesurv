/* ************
 * Initialization functions
 * The functions in this file should only be called directly from R 
 * during the initialization phase.
 ***************/
#include <R_ext/BLAS.h>    
#include <R_ext/Print.h>    
#include <math.h>
#include <stdlib.h>
#include "init.h"

double cnBsmom(const int N, const int ord, const int j, const double *knots, const int knotord)
{
    double lknot=knots[j-ord+knotord];
    double rknot=knots[j+knotord];
    if(lknot==rknot){
        return pow(rknot,N);
    }
    double out;
    out = ord/(rknot-lknot)/(N+1);
    out = out*(-cnBsmom(N+1,ord-1,j-1,knots,knotord)+cnBsmom(N+1,ord-1,j,knots,knotord));
    return out;
}

void cevalEinte(double *einte, double *knots, int *ord, int *K, int *N)
{
    int j;
    for(j=0; j < *K+*ord; j++)
        einte[j]=1.0-cnBsmom(*N,*ord,j,knots,*ord);
}

void cevalBinte(double *binte, double *knots, int *ord, int *K)
{
    for(int j=0; j<*K+*ord; j++)
        binte[j]=(knots[j+*ord]-knots[j]) / *ord;
}


double csplineeval( double x, int j, int ord, double *knots, int splord, int nj)
{
    double knotj=knots[j+splord];
    double knotjm1=knots[j-1+splord];
    if(ord==1) {
        if((x >= knotjm1) & (x < knotj)) return 1.0;
        if((x==knotj) & (j+splord==nj)) return 1.0;
        return 0.0;
    }
    double knotjmn=knots[j-ord+splord];
    double knotjmn1=knots[j-ord+1+splord];
    double numer1 = (x-knotjmn);
    double denom1 = (knotjm1 - knotjmn);
    double numer2 = (knotj-x);
    double denom2 = (knotj - knotjmn1);
    double out = 0;
    if(denom1>0){
        double rec1 = csplineeval(x,j-1,ord-1,knots,splord,nj);
        out += numer1/denom1*rec1;
    }
    if(denom2>0){
        double rec2 = csplineeval(x,j,ord-1,knots,splord,nj);
        out += numer2/denom2*rec2;
    }
    return out;
}

void csplinedesign(double *des, double *x, int *nx, double *knots, int *ord, int *K)
{
    int nj= *K + *ord;
    for(int i=0;i<*nx;i++){
        for(int j=0; j< nj; j++)
            des[i + j* *nx] = csplineeval( x[i], j, *ord, knots, *ord, nj);
    }
}

void cevalCinte(double *cinte, double *x, int *nx, double *knots, int *ord, int *K, double *binte)
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


double cSplineConvolution(int k, int n1,int j1, int n2, int j2, double *knots, int splord)
{
    double out=0;
    if((j1 -n1>=j2) | (j2 -n2>=j1)) return out;
    if((n1==1) & (n2==1)){
        out= 1.0/(k+1.0)*(pow(knots[j1+splord],k+1.0)-pow(knots[j1-1+splord],k+1.0));
        return out;
    }
    if(n2>n1){
        int n3=n1; n1=n2; n2=n3;
        int j3=j1; j1=j2; j2=j3;
    }
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

double cSplineDerivInt(int l1, int n1, int j1, int l2, int n2, int j2, double *knots, int splord)
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

void cMakePenalty2diff(double *P, double *knots, int *ord, int *K)
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

void cInitMultAndWeight( double *y,  double *B,  double *par,  double *weight,  int  *ny,  int *nj)
{
    const int c1 = 1;
    const char trans = 'N';
    const double onemw=1-*weight;
    // compute w*B%*%exp(par) + (1-w)y
    F77_CALL(dgemv)(&trans,ny,nj,weight,B,ny,par,&c1,&onemw,y,&c1);
}

void InitSmoothnessPenalty(double *pen, double *par, double *P, int *penaltyType, double *sigma2, int *nj)
{
    if(*penaltyType==0) *pen=0; //No penalty
    if(*penaltyType==1 || *penaltyType==2 || *penaltyType==3){   // second differences or second derivative
        double * temp = (double *) calloc((*nj) , sizeof(double));
        int c1 = 1;
        double c1d = 1.0;
        char uplo = 'u';
        F77_CALL(dsymv)( &uplo, nj, &c1d, P, nj, par, &c1, &c1d, temp, &c1);
        *pen = F77_CALL(ddot)( nj, temp, &c1, par, &c1);
        if(*penaltyType==3) *pen=log(*pen+1);
        *pen = *pen / (2 * *sigma2);
        free( temp );
    }
}

void addInitSmoothnessPenaltyGr(double *gr, double * penpar, double *P, int *penaltyType, double * sigma2, int *nj)
{
    if(*penaltyType==1 || *penaltyType==2 || *penaltyType==3){   // second differences or second derivative
        double * temp = (double *) calloc((*nj) , sizeof(double));
        int c1 = 1;
        double c1d = 1.0;
        char uplo = 'u';
        F77_CALL(dsymv)( &uplo, nj, &c1d, P, nj, penpar, &c1, &c1d, temp, &c1);
        for(int i=0; i<*nj; i++) temp[i]=temp[i]*penpar[i] / *sigma2;
        double scalefactor=-1;
        if(*penaltyType==3) {
            double pen=0;
            int c2 = 2;
            InitSmoothnessPenalty( &pen, penpar, P, &c2, sigma2, nj);
            scalefactor= -1/((2* *sigma2) * pen + 1);
        }
        F77_CALL(daxpy)(nj, &scalefactor, temp, &c1, gr, &c1);
        free( temp );
    }
}

void cInitLikHazSpline(double *lik, double *par, double *status, double *lp, double *frailrep,
    double *hazParY, double *hazParYcum, double *weight, double *B, double *C,
    double *P, int *penaltyType, double *sigma2, int *ny, int *nj)
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

void cInitGrHazSpline(double *gr, double *par, double *status, double *lp, double *frailrep,
    double *hazParY, double *hazParYcum, double *weight, double *B, double *C,
    double *P, int *penaltyType, double *sigma2, int *ny, int *nj)
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

void cInitLikFrailSpline(double *lik, double *par, double *frailParY, double *weight, double *B, double *E, double *M, double *P, int * penaltyType, double *sigma2, int *ny, int *nj)
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
