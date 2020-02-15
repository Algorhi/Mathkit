#include"Vehicle.h"
#include<cmath>
using namespace std;
EVER::EVER(double *X){
    r    = X[0];
    theta= X[1];
    psi  = X[2];
    V    = X[3];
    gamma= X[4];
    phi  = X[5];
};
EVER::~EVER(){}