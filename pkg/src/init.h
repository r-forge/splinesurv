#ifndef SSINIT_H
#define SSINIT_H
double cnBsmom(const int N, const int ord, const int j, const double *knots, const int knotord);
void cevalEinte(double *einte, double *knots, int *ord, int *K);
void cevalBinte(double *binte, double *knots, int *ord, int *K);
double csplineeval( double x, int j, int ord, double *knots, int splord , int n);
void csplinedesign(double *des, double *x, int *nx, double *knots, int *ord, int *K);
#endif
