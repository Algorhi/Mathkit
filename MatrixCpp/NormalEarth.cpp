#include "NormalEarth.h"
#include <cmath>
using namespace std;
//#include "Matrix.h"
//#include "Mathkit.h"
//using namespace mathkit;
#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
NormalEarth::NormalEarth ()
{

}
NormalEarth::NormalEarth (double _g,double _R,double _ome):g0(_g),R0(_R),Omega(_ome)
{

}

NormalEarth::~NormalEarth ()
{

}


/// h

double
NormalEarth::Get_Core_Lat (double B)	// 
{
  /*double phi=atan( tan(B)*pow(1-alphae(),2) );

     return phi; */
  double phi = atan (tan (B) * be () * be () / (ae () * ae ()));

  return phi;

}

//h
double
NormalEarth::Get_Core_Distance (double phi)
{
  /*double cs2 = 1- sin(phi)*sin(phi);

     double R = (1.-alphae())*ae()/sqrt(1.+alphae()*(alphae()-2.)*cs2);

     return R; */
  double sn2 = sin (phi) * sin (phi);
  double cs2 = 1 - sn2;

  double R = ae () * be () / sqrt (ae () * ae () * sn2 + be () * be () * cs2);

  return R;

}

//h.!g
void
NormalEarth::Get_R0 (double A0, double B0, double h0,
		      double &R0, double &R0x, double &R0y, double &R0z)
{
  double phi0 = Get_Core_Lat (B0);

  double R0_temp = Get_Core_Distance (phi0);
  R0x = -R0_temp * sin ((B0 - phi0)) * cos (A0);
  R0y = R0_temp * cos ((B0 - phi0)) + h0;
  R0z = R0_temp * sin ((B0 - phi0)) * sin (A0);
  R0 = sqrt (R0x * R0x + R0y * R0y + R0z * R0z);

}


//f1
double
NormalEarth::Get_Height (double r, double phi)
{
  //
  double R = Get_Core_Distance (phi);

  //
  double h = r - R;

  return h;
}

double
NormalEarth::Get_Great_Circle_Distance (const double phi0,
					 const double theta0,
					 const double phif,
					 const double thetaf) const 
{
    double xtemp;
    xtemp = cos(phi0)*cos(phif)*cos(theta0-thetaf)+sin(phi0)*sin(phif);
  return R0*acos(xtemp);
}

double
NormalEarth::Get_Azimuth (const double phi0, const double theta0,
			   const double phif, const double thetaf) const 
{
  double x = sin (thetaf - theta0);
  double y = cos (phi0) * tan (phif) - sin (phi0) * cos (thetaf - theta0);
  return atan2 (x, y);
}
