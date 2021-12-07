#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

void nrerror(char error_text[]);
double *vector(long nl, long nh);
void free_vector(double *v, long nl, long nh);
double **matrix(long nrl, long nrh, long ncl, long nch);
void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);
void mmid(double y[], double dydx[], int nvar, double xs, double htot, int nstep,
    double yout[], void (*derivs)(double, double[], double[]));
void pzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv);
void bsstep(double y[], double dydx[], int nv, double *xx, double htry, double eps,
    double yscal[], double *hdid, double *hnext,
    void (*derivs)(double, double [], double []));
void rkck(double y[], double dydx[], int n, double x, double h, double yout[],
            double yerr[], void (*derivs)(double, double [], double []));
void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
    double yscal[], double *hdid, double *hnext,
    void (*derivs)(double, double [], double []));
void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
    double hmin, int *nok, int *nbad,
    void (*derivs)(double, double [], double []),
    void (*rkqs)(double [], double [], int, double *, double, double, double [],
    double *, double *, void (*)(double, double [], double [])));
