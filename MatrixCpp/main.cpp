/******************************************************************************

FSD

*******************************************************************************/
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include "Mathkit.h"
#include <cmath>
#include "RunFlight.h"
using namespace std;
using namespace mathkit;
#include"Matrix.h"
int main ()
{
  double ys[]={4,9,16,25}; 
  CMatrix M(2,2,ys);
  CMatrix C = M.Pow(0.5);
  CMatrix TD = M.Transpose();//GetTransposed  Transpose
  CMatrix TD2=TD.GetTransposed();
  int rank = C.Rank();
  printf("%d\n",rank);
  //CMatrix *pTD=&TD;
  
  
  //CMatrix TD2(*pTD);
  for (int i = 0; i < 4; i++)
  cout << "  y(" << i << ") = " << setw (8) << *(TD.GetData()+i);
  printf("\n%f\n",*(C.GetData()+1));
  /*
  const NormalEarth Earth;
  double at = Earth.Get_Azimuth (1, 1, 1, 0);
  double gcd = Earth.Get_Great_Circle_Distance(0.1,0.2,1.8,1.6);
  printf ("at=%f\n", at);
  printf ("gcd=%f\n", gcd);
  double r0 = 0;
  double theta0 = 0;
  double psi0 = 0;
  double V0 = 0;
  double gamma0 = 0;
  double phi0 = 0;
  double Y0[] = { r0, theta0, psi0, gamma0, phi0 };
  EVER EV (Y0);
  double a = Pi ();
  //printf("a=%f\n",a);


  double t, h, y[3];
  double user[3] = { 10, 8 / 3, 13 };
  y[0] = 1.0;
  y[1] = 1.0;
  y[2] = 1.0;
  t = 0.0;
  h = 0.1;
  int n = 3;
  for (int j = 1; j <= 100; j++)
    {
      //Fixstep_RungeKutta(t,h,y,user,n,Lorenz);
      RungeKutta (t, h, y, user, n, Lorenz);
      t = t + h;
      cout << "t=" << setw (4) << t;
      for (int i = 0; i < n; i++)
	cout << "  y(" << i << ") = " << setw (8) << y[i];
      cout << endl;
    }

/*
	  int i,j,n,m;
      double u,v,w;
      n=10;
      m=11;
      double x[n],y[m],z[n][m];
      for (i=0;i<n;i++)
      { x[i]=0.1*i; y[i]=x[i];}
      for (i=0;i<n;i++)
      for (j=0;j<m;j++)
        z[i][j]=exp(-(x[i]-y[j]));
      u=0.35; v=0.65;
      w=spline(x,y,*z,11,11,u,v);
	  cout <<"x = " <<u <<",   y = " <<v <<"     z(x,y) = " <<w <<endl;
      u=0.45; v=0.55;
      w=spline(x,y,*z,11,11,u,v);
	  cout <<"x = " <<u <<",   y = " <<v <<"     z(x,y) = " <<w <<endl;
      w=spline(y,*(z+1),m,u);
	  cout <<"x = " <<u <<",   y = " <<v <<"     z(x,y) = " <<w <<endl;	  
	  double ta=1.0;
	  double *pta=&ta;
	  
	  printf("%f\n",pta[1]);
	  printf("%f\n",*z[1]); 
	  
*/

  return 0;
}
