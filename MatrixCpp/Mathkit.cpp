
#include <math.h>
#include <stdlib.h> 
#include "Mathkit.h"
#include <stdio.h>
 
//定步长改进欧拉方法解一阶常微分方程
int mathkit::Fixstep_Euler(double t, double h, double *y, double *user, int n, ODE Fun)
{
	double *K1=new double[n];                              
	double *K2=new double[n];                              
	double *tempY=new double[n];                           
	int i;

	if(!K1 || !K2 || !tempY) 
		return 0;
	
	Fun(t,y,K1,user);                                        

	for(i=0;i<n;i++)
		tempY[i]=y[i]+K1[i]*h;

	Fun(t+h,tempY,K2,user);                              

	for(i=0;i<n;i++)
		y[i]=y[i]+h/2.0*(K1[i]+K2[i]);

	delete [] K1;
	delete [] K2;
	delete [] tempY;

	return 1;
}


//定步长4阶龙格－库塔方法解一阶常微分方程
int mathkit::Fixstep_RungeKutta(double t, double h, double *y, double *user, int n, ODE Fun)
{
	double *K1=new double[n];                              
	double *K2=new double[n];                              
	double *K3=new double[n];                              
	double *K4=new double[n];                              
	double *tempY=new double[n];                           
	int i;

	if(!K1 || !K2 || !K3 || !K4 || !tempY) 
		return 0;
	
	Fun(t,y,K1,user);                                        

	for(i=0;i<n;i++)
		tempY[i]=y[i]+K1[i]*h/2.0;

	Fun(t+h/2.0,tempY,K2,user);                              

	for(i=0;i<n;i++)
		tempY[i]=y[i]+K2[i]*h/2.0;

	Fun(t+h/2.0,tempY,K3,user);                                

	for(i=0;i<n;i++)
		tempY[i]=y[i]+K3[i]*h;

	Fun(t+h,tempY,K4,user);                                  

	for(i=0;i<n;i++)
		y[i]=y[i]+h/6.0*(K1[i]+2.0*K2[i]+2.0*K3[i]+K4[i]);

	delete [] K1;
	delete [] K2;
	delete [] K3;
	delete [] K4;
	delete [] tempY;

	return 1;
}

int mathkit::RungeKutta(double t, double h, double *y, double *user, int n, ODE Fun)
{
    int m,i,j,k;
    double hh,p,dt,x,tt,q,a[4];
    double *g=new double[n];
    double *b=new double[n];
    double *c=new double[n];
    double *d=new double[n];
    double *e=new double[n];
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
    delete[] g; 
    delete[] b; 
    delete[] c; 
    delete[] d;
    delete[] e;
    return 1;
    
}


//圆周率
double mathkit::Pi()
{
	return   3.1415926535897932384626433832795;
}


//弧度转化为角度
double mathkit::Rad2Deg(double Alpha)
{
	return Alpha/Pi()*180.;
}


//角度转化为弧度
double mathkit::Deg2Rad(double Alpha)
{
	return Alpha/180.*Pi();
}

//符号函数
int mathkit::sign(double x)
{
	if(x >= 0.)
		return +1;
	else
		return -1;
}



void mathkit::Cross_Multip(double X[], double Y[], double Z[])
{
	Z[0] = -X[2]*Y[1] + X[1]*Y[2];
	Z[1] =  X[2]*Y[0] - X[0]*Y[2];
	Z[2] = -X[1]*Y[0] + X[0]*Y[1];	
}

double mathkit::spline(double x[],double y[], int n, const double t)
{
    printf("--------\n");
    return(0);
    
    
}

double mathkit::spline(double *x, double *y, double *z, int n, int m, const double u, const double v)
{
	  int ip,ipp,i,j,l,iq,iqq,k;
      double h,w,b[10];
      if (u<=x[0]) { ip=1; ipp=4;}
      else if (u>=x[n-1]) { ip=n-3; ipp=n;}
      else     //X方向取u前后4个坐标
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





