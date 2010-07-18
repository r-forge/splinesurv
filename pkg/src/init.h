#ifndef SSINIT_H
#define SSINIT_H
double cnBsmom(const int N, const int ord, const int j, const double *knots, const int knotord);
void cevalEinte(double *einte, double *knots, int *ord, int *K, int *N);
void cevalBinte(double *binte, double *knots, int *ord, int *K);
void cevalCinte(double *cinte, double *x, int *nx, double *knots, int *ord, int *K, double *binte);
void cevalCinte2(double *cinte, double *x, int nx, double *knots, int nj, int ord, double *binte, int i, int jstart, int jstop);
double csplineeval( double x, int j, int ord, double *knots, int splord , int n);
//double csplinecumeval(double x, int j, int ord, int nj, double *knots, double * binte);
void csplinedesign(double *des, double *x, int *nx, double *knots, int *ord, int *K);
#endif
