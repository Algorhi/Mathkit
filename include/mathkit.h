
#ifndef __MATHKIT_H__
#define __MATHKIT_H__

#ifdef __cpluplus
extern "C" {
#endif

#include"mkmatrix.h"
#include"mktrigono.h"
typedef void (*mkODE)(double t, double *y, double *dy, double *user);
int rk4fixstep(double t, double h, double *y, double *user, int n ,mkODE fun);
int RungeKutta(double t, double h, double *y, double *user, int n ,mkODE fun);


double spline2d(double x[], double y[], int n, const double t);
void interp1(double x[], double y[], int n, double t[], int m,double z[], char* method);
void interp1linear(double x[], double y[], int n, double t[], int m,double z[]);
void interp1cubic(double x[], double y[], int n, double t[], int m,double z[]);
void interp1spline(double x[], double y[], int n, double t[], int m,double z[]);
double spline3d(double *x, double *y, double *z, int n, int m, const double u, const double v);
void linspace(double d1, double d2, int n, double* y);
void mappingsort(double* X, double* Y, int* Num);
double norm(double* vector,int p);


#ifdef __cpluplus
}
#endif
#endif
