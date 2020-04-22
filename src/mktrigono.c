#include<math.h>
#include<stdio.h>
#include"mktrigono.h"
//static double abs(double x){
//  if (x>0)
//    return x;
//  else
//    return -x;
//}

double Pi(){
  return 3.1415926535897932384626433832795;
}
double Rad2Deg(double Alpha){
  return Alpha/Pi()*180;
}

double Deg2Rad(double Alpha){
  return Alpha/180.*Pi();
}

void Sphere(double theta1, double phi1, double theta2, double phi2, double* arc, double* LOS){

  double Dlon = theta2 - theta1;
  double Dlat = phi2 - phi1;
  int flag = 0;
  if (Dlon == 0 && Dlat == 0)
  {
    flag = 1;
    *arc = 0.0;
    *LOS = 0.0;
  }
  else if (fabs(Dlon) < 0.000000000003 && fabs(Dlat) <0.000000000003)
  {
    flag = 1;
    double Dlongi = Dlon*cos(phi1);
    *arc = sqrt(Dlongi*Dlongi + Dlat*Dlat);
    if (Dlat > 0.0)
      *LOS = atan(Dlongi/Dlat);
    else if (Dlat <= 0.0)
      *LOS = atan(Dlongi/Dlat) + Pi();
    else if (Dlat == 0.0 && Dlon  > 0.0)
      *LOS = 1.5707963;
    else 
      *LOS = -1.5707963;
  }

  else 
  {
    double sn1 = sin(phi1);
    double sn2 = sin(phi2);
    double cn1 = cos(phi1);
    double cn2 = cos(phi2);
    double arc0 = sn2*sn1 + cn2*cn1*cos(Dlon);
    double arc1 = acos(arc0);
    *arc = arc1;
    double sn = sin(Dlon)*cn2/sin(arc1);
    double cs = (sn2*cn1-cn2*sn1*cos(Dlon)/sin(arc1));
    if(cs >= 0.0 && sn >= 0.0)
    {
      *LOS = asin(sn);
      flag = 1;
    }
    else if (cs<=0.0 && sn >= 0.0)
    {
      *LOS = acos(cs);
      flag = 1;
    }
    else if (cs >= 0.0 && sn <= 0.0)
    {
      *LOS = asin(sn);
      flag = 1;
    }
    else if (cs <= 0.0 && sn <= 0.0)
    {
      *LOS = -acos(cs);
      flag = 1;
    }
  }

  if (flag==0)
  {
    printf("Warning! arc or azimuth is not returned.\n");
    printf("theta1=%6.2f, phi1=%6.2f, theta2=%6.2f, phi2=%6.2f.\n",theta1,phi1,theta2,phi2);
    printf("Dlon=%6.2f,Dlat=%6.2f.\n",Dlon,Dlat);
  }

}

void Bearing(double* PointA, double* PointB, double* BearingA, double* BearingB){
  double longiA = *(PointA);
  double latiA = *(PointA+1);
  double longiB = *(PointB);
  double latiB = *(PointB+1);
  *BearingA = atan2(sin(longiB-longiA)*cos(latiB),cos(latiA)*sin(latiB)-sin(latiA)*cos(latiB)*cos(longiB-longiA));
  double longiT = longiA;
  double latiT = latiA;
  longiA = longiB;
  latiA = latiB;
  longiB = longiT;
  latiB = latiT;
  *BearingB = atan2( sin(longiB-longiA)*cos(latiB), cos(latiA)*sin(latiB)-sin(latiA)*cos(latiB)*cos(longiB-longiA) );
  *BearingB = fmod(*BearingB+Pi(),2*Pi());
}

double GreatCircleDistance(double* statpoint, double* despoint){
  double longiA = *(statpoint);
  double latiA = *(statpoint+1);
  double longiB = *(despoint);
  double latiB = *(despoint+1);
  double a;
  a = pow((sin((latiB-latiA)/2)),2) + cos(latiA)*cos(latiB)*pow(sin((longiB-longiA)/2),2);
  double c;
  c = 2*atan2(sqrt(a),sqrt(1-a));
  return c;
}
