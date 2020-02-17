#include<iostream>
using namespace std;
void func_a(int&);
void func_b(int*&);
int main(){
  int a = 0;
  //int b[4] = {1,2,3,4};
  int* b;
  *b = 0;
  func_a(a);
  func_b(b);
  cout<<"a="<<a<<endl;
  cout<<"b="<<*b<<endl;
}

void func_a(int& a){
  a = 2; 
};
void func_b(int*& b){
  *b = 2; 
};
