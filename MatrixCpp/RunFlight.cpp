#include"Vehicle.h"
#include"NormalEarth.h"
#include<cmath>
using namespace std;


void Dynamic(double t, double *y, double *dy, double *user)
{
    
    double r=y[0];
    double theta=y[1];
    double phi=y[2];
    double V=y[3];
    double gamma=y[4];
    double psi=y[5];
    
    double D =0;
    
    dy[0]=V*sin(gamma);
    dy[1]=V*cos(gamma)*sin(psi)/(r*cos(phi));
    dy[2]=V*cos(gamma)*cos(psi)/r;
    dy[3]=-D-(sin(gamma)/(r*r));
    
    
}


//test Functions
void Func(double t, double *y, double *dy, double *user)
{
    
    *dy = *y + t;
}
void rktf(double t, double y[], double dy[],double *user)
{ 
    dy[0]=y[1]; 
    dy[1]=-y[0]+1;
}
void Lorenz(double t, double y[], double dy[],double *user)
{
    double sigma= user[0];
    double rou  = user[1];
    double beta = user[2];
    dy[0]=sigma*(y[1]-y[2]);
    dy[1]=y[0]*(rou-y[2])-y[1];
    dy[2]=y[0]*y[1]-beta*y[2];
}