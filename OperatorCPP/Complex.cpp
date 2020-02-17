#include<Complex.h>
#include<iostream>
#include<math.h>
using namespace Algorhi;
Complex::Complex(double re,double im):real(re),imag(im){
  std::cout<<"Constructer called" << std::endl;
};

Complex operator+(Complex& a, Complex& b){
  return Complex(a.real + b.real, a.imag + b.imag);
}
