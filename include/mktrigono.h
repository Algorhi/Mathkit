#ifndef __MAKTRIGONO_H__
#define __MAKTRIGONO_H__
#ifdef __cplusplus
extern "C" {
#endif

double Pi();
double Rad2Deg(double Alpha);
double Deg2Rad(double Alpha);
void Sphere(double theta1, double phi1, double theta2,double phi2, double* arc, double* LOS);
void Bearing(double* PointA, double* PointB, double* BearingA, double* BearingB);
double GreatCircleDistance(double* statpoint, double* despoint);

#ifdef __cplusplus
}
#endif
#endif
