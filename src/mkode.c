#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>

typedef void (*ODE_mk)(double t,  double *y, double *dy, double *user);
const double DEFAULTEPS = 1.0e-10;
int rk4fixstep(double t, double h, double *y, double *user, int n, ODE_mk fun){
  //double *K1 = new double[n];
  //double *K2 = new double[n];
  //double *K3 = new double[n];
  //double *K4 = new double[n];
  //double *tempY = new double[n];
  double *K1 = (double*)malloc(sizeof(double)*n);
  double *K2 = (double*)malloc(sizeof(double)*n);
  double *K3 = (double*)malloc(sizeof(double)*n);
  double *K4 = (double*)malloc(sizeof(double)*n);
  double *tempY= (double*)malloc(sizeof(double)*n);
  int i;

  if(!K1 || !K2 || !K3 || !K4 || !tempY)
    return 0;

  fun(t,y,K1,user);

  for (i=0;i<n;i++)
    tempY[i] = y[i] + K1[i]*h/2.0;
  fun(t+h/2.0,tempY,K2,user);

  for (i=0;i<n;i++)
    tempY[i] = y[i] + K2[i]*h/2.0;
  fun(t+h/2.0,tempY,K3,user);

  for (i=0;i<n;i++)
    tempY[i] = y[i] + K3[i]*h;
  fun(t+h,tempY,K4,user);

  for (i=0;i<n;i++)
    tempY[i] = y[i] + K1[i]*h/2.0;
  fun(t+h,tempY,K4,user);

  for(i=0;i<n;i++)
    y[i] = y[i]+h/6.0*(K1[i]+2.0*K2[i]+2.0*K3[i]+K4[i]);
  
  free(K1);
  free(K2);
  free(K3);
  free(K4);
  free(tempY);
  return 1;

}

int RungeKutta(double t, double h, double *y, double *user, int n,ODE_mk Fun)
{
  int m,i,j,k;
  double hh,p,dt,x,tt,q,a[4];
  double *g=(double*)malloc(sizeof(int)*n);
  double *b=(double*)malloc(sizeof(int)*n);
  double *c=(double*)malloc(sizeof(int)*n);
  double *d=(double*)malloc(sizeof(int)*n);
  double *e=(double*)malloc(sizeof(int)*n);
  hh=h;
  m=1;
  p=1.0+DEFAULTEPS;
  x=t;
  for (i=0;i<=n-1;i++) c[i]=y[i];
  while (p>=DEFAULTEPS)
  {
    a[0]=hh/2.0;
    a[1]=a[0];
    a[2]=hh;
    a[3]=hh;
    for (i=0;i<=n-1;i++)
    {
      g[i]=y[i];
      y[i]=c[i];
    }
    dt=h/m;
    t=x;
    for (j=0;j<=m-1;j++)
    {
      Fun(t,y,d,user);
      for (i=0; i<=n-1; i++)
      {
        b[i]=y[i];
        e[i]=y[i];
      } 
      for (k=0; k<=2; k++)
      {
        for (i=0; i<=n-1; i++)
        {
          y[i]=e[i]+a[k]*d[i];
          b[i]=b[i]+a[k+1]*d[i]/3.0;
        }
        tt=t+a[k];
        Fun(tt,y,d,user);
      }
      for (i=0; i<=n-1; i++) y[i]=b[i]+hh*d[i]/6.0;
      t=t+dt;
      
    }
    p=0.0;
    for (i=0; i<=n-1; i++)
    {
      q=fabs(y[i]-g[i]);
      if (q>p) p=q;
    }
    hh=hh/2.0; 
    m=m+m;
  }
  free(g); 
  free(b); 
  free(c); 
  free(d);
  free(e);
  return 1;
}
