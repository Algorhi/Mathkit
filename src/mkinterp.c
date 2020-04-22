#include<math.h>
#include<malloc.h>
#include"mathkit.h"

static void getindex(double x[], int s, double point, int* bef, int* aft);

double spline2d(double *x, double *y, int N, const double p){
  double yy,t;
  int j,k;
  yy = 0.0;
  for(j=0;j<N;j++)
  {
    t=1.0;
    for(k=0;k<N;k++)
      if(k!=j)
        t=t*(p-x[k])/(x[j]-x[k]);
    yy = yy + t*y[j];
  }
  return yy;
}
void interp1(double* x, double* y, int n, double* t, int m, double* z, char* method){
  if(strcmp(method,"linear")==0)
    interp1linear(x, y, n, t, m, z);
  if(strcmp(method,"spline")==0)
    interp1spline(x, y, n, t, m, z);
}

void interp1cubic(double* x, double* y, int n, double* t, int m, double* z){
  double txx;
}

void interp1linear(double* x, double* y, int n, double* t, int m, double* z){
  int bef, aft;
  for(int i=0;i<m;i++)
  {
    getindex(x,n,*(t+i),&bef,&aft);
    *(z+i) = y[bef]+(t[i]-x[bef])*((y[bef]-y[aft])/(x[bef]-x[aft]));
  }
  return;
}

void interp1spline(double* x, double* y, int n, double* t, int m, double* z){
  // 第二种边界条件下的三次样条插值、微商和积分 
  //
  //
  double* dy  =(double*)malloc(sizeof(double)*n);
  int i, j;
  double h0,h1,alpha,beta,g,*s;
  double* ddy  =(double*)malloc(sizeof(double)*n);
  //ddy[0] = -0.279319;
  //ddy[n-1] = 0.0111560;
  ddy[0] = 0.0;
  ddy[n-1] = 0.0;
  //
  s = (double*)malloc(sizeof(double)*n);
  dy[0]=-0.5;
  h0=x[1]-x[0];
  s[0]=3.0*(y[1]-y[0])/(2.0*h0)-ddy[0]*h0/4.0;
  for(j=1;j<=n-2;j++)
  {
    h1=x[j+1]-x[j];
    alpha=h0/(h0+h1);
    beta=(1.0-alpha)*(y[j]-y[j-1])/h0;
    beta=3.0*(beta+alpha*(y[j+1]-y[j])/h1);
    dy[j]=-alpha/(2.0+(1.0-alpha)*dy[j-1]);
    s[j]=(beta-(1.0-alpha)*s[j-1]);
    s[j]=s[j]/(2.0+(1.0-alpha)*dy[j-1]);
    h0=h1;
  }
  dy[n-1]=(3.0*(y[n-1]-y[n-2])/h1+ddy[n-1]*h1/2.0-s[n-2])/(2.0+dy[n-2]);
  for(j=n-2;j>=0;j--)
    dy[j]=dy[j]*dy[j+1]+s[j];

  for(j=0;j<=n-2;j++)
    s[j]=x[j+1]-x[j];

  for(j=0;j<=n-2;j++)
  {
    h1=s[j]*s[j];
    ddy[j]=6.0*(y[j+1]-y[j])/h1-2.0*(2.0*dy[j]+dy[j+1])/s[j];
  }

  h1 = s[n-2]*s[n-2];
  ddy[n-1]=6.0*(y[n-2]-y[n-1])/h1+2.0*(2.0*dy[n-1]+dy[n-2])/s[n-2];
  g=0.0;

  for(i=0;i<=n-2;i++)
  {
    h1=0.5*s[i]*(y[i]+y[i+1]);
    h1=h1-s[i]*s[i]*s[i]*(ddy[i]+ddy[i+1])/24.0;
    g=g+h1;
  }

  for(j=0;j<=m-1;j++)
  {
    if(t[j]>=x[n-1])
      i=n-2;
    else
    {
      i=0;
      while(t[j]>x[i+1])
        i=i+1;
    }
    h1=(x[i+1]-t[j])/s[i];
    h0=h1*h1;
    z[j]=(3.0*h0-2.0*h0*h1)*y[i];
    z[j]=z[j]+s[i]*(h0-h0*h1)*dy[i];
    h1=(t[j]-x[i])/s[i];
    h0=h1*h1;
    z[j]=z[j]+(3.0*h0-2.0*h0*h1)*y[i+1];
    z[j]=z[j]-s[i]*(h0-h0*h1)*dy[i+1];
  }
  free(s);
  free(dy);
  free(ddy);
}
double spline3d(double *x, double *y, double *z, int n, int m, const double u, const double v){
  int ip,ipp,i,j,l,iq,iqq,k;
  double h,w,b[10];
  if (u<=x[0]) { ip=1; ipp=4;}
  else if (u>=x[n-1]) { ip=n-3; ipp=n;}
  //else     //X方向取u前后4个坐标
  else 
  { 
    i=1; j=n;
    while (((i-j)!=1)&&((i-j)!=-1))
    { 
      l=(i+j)/2;
      if (u<x[l-1]) j=l;
      else i=l;
    }
    ip=i-3; ipp=i+4;
  }
  if (ip<1) ip=1;
  if (ipp>n) ipp=n;
  if (v<=y[0]) { iq=1; iqq=4;}
  else if (v>=y[m-1]) { iq=m-3; iqq=m;}
  else    //Y方向取v前后4个坐标
  { 
    i=1; j=m;
    while (((i-j)!=1)&&((i-j)!=-1))
    { 
      l=(i+j)/2;
      if (v<y[l-1]) j=l;
      else i=l;
    }
    iq=i-3; iqq=i+4;
  }
  if (iq<1) iq=1;
  if (iqq>m) iqq=m;
  for (i=ip-1;i<=ipp-1;i++)
  { 
    b[i-ip+1]=0.0;
    for (j=iq-1;j<=iqq-1;j++)
    { 
      h=z[m*i+j];
      for (k=iq-1;k<=iqq-1;k++)
      if (k!=j) h=h*(v-y[k])/(y[j]-y[k]);
      b[i-ip+1]=b[i-ip+1]+h;
    }
  }
  w=0.0;
  for (i=ip-1;i<=ipp-1;i++)
  { 
    h=b[i-ip+1];
    for (j=ip-1;j<=ipp-1;j++)
    if (j!=i) h=h*(u-x[j])/(x[i]-x[j]);
    w=w+h;
  }
  return(w);
}

void linspace(double d1, double d2, int n, double* y){
  //LINSPACE linearing spaced vector.
  for (int i=0;i<n;i++)
  {
    *(y+i) = d1 + (d2-d1)/(n-1)*i;
  }
}

void mappingsort(double* X, double* Y, int* Num){
  /* Make the sequence Y change the sort way as X, 
     where X is arranged in the order of arrival from small to small. */
  int n = *Num;
  double Xtemp;
  double Ytemp;
  for(int i=0;i<n-1;i++)
  {
    for(int j=i+1;j<n;j++)
    {
      if(X[i]>X[j])
      {
        Xtemp = X[i];
        Ytemp = Y[i];
        X[i]  = X[j];
        Y[i]  = Y[j];
        X[j]  = Xtemp;
        Y[j]  = Ytemp;
      }
    }
  }
}

void getindex(double x[], int s, double t, int* m, int* n){
  if(s<=2) printf("s must >= 2!!\n");
  if(t<=MinDouble(x,s))
  {
    *m = 0;
    *n = 1;
    printf("extern point!\n");
  }
  if(t>=MaxDouble(x,s))
  {
    *m = s-2;
    *n = s-1;
    printf("extern point!\n");
  }
  for(int i=0;i<s-1;i++)
  {
    if(t>x[i] && t<= x[i+1])
    {
      *m = i;
      *n = i+1;
    }
  }
  return;
}
