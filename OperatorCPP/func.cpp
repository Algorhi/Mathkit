#include"func.h"
#include<iostream>
using namespace std;
using namespace Algorhi;

Func::Func(int a):swap(a){
  cout<<"first construction function"<<endl;
}

Func::Func(int a,int b):m(a),n(b){
  cout<<"second construction function"<<endl;
}

Func::Func(const Func& other){
  this->m = other.m;
  this->n = other.n;
  cout<<"copied construction function"<<endl;
}

Func::~Func(){
  cout<<"Disconstruction function"<<endl;
}

void Func::get(){
  Func result(*this);
  cout<<m<<n<<endl;
}
