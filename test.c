#include<stdio.h>
#include<math.h>
#include"mathkit.h"

int main(void){
double h = 10.0;
double t = 0;
int n = 2;
double s[2]={2.0,3.0};
double user[2] ={2.0,3.0};

// mappingsort function test.
//double X[]={1,2,4,3,6,5};
//double Y[]={1,3,5,2,4,6};
//int Num = 6;
//for(int i=0;i<Num;i++) printf("X[%d]=%f  ",i,X[i]);
//printf("\n");
//for(int i=0;i<Num;i++) printf("Y[%d]=%f  ",i,Y[i]);
//printf("\n");
//mappingsort(X,Y,&Num);
//for(int i=0;i<Num;i++) printf("X[%d]=%f  ",i,X[i]);
//printf("\n");
//for(int i=0;i<Num;i++) printf("Y[%d]=%f  ",i,Y[i]);
//printf("\n");

double valA[]={1,2,3,4,5,6};
double valB[]={1,3,5,2,4,6};
//double val[]={10,0,0,0,-2,0,3,9,0,0,0,3,0,7,8,7,0,0,3,0,8,7,5,0,0,8,0,9,9,13,0,4,0,0,2,-1};
double val[]={10,0,0,0,-2,0, 3,9,0,0,0,3, 0,7,0,7,0,0, 3,0,0,7,5,0, 0,8,0,9,9,13, 0,4,0,0,2,-1};
//double val[]={1,0,0,0,2,2.1,3,0,3.2};
Matrix* mat = MatrixInit(6,6);
MatrixSetData(mat,val);
MatrixPrint(mat);
MatrixBlockPrint(mat,1,5,2,4);
double Value[17];
int Row_ind[17];
int Col_ptr[7];
//printf("number of Value and Row_ind is %d\n\n",MatrixNumNonzero(mat));
Matrix2CCS(mat,Value,Row_ind,17,Col_ptr,7);
//Matrix* MATT=MatrixZeros(6,6);
//Matrix2CCS(MATT,Value,Row_ind,17,Col_ptr,7);
//MatrixPrint(MATT);

for(int i=0;i<19;i++) printf("Value[%d]=%f  ",i,Value[i]);
printf("\n");
for(int i=0;i<19;i++) printf("Row_ind[%d]=%d  ",i,Row_ind[i]);
printf("\n");
for(int i=0;i<7;i++) printf("Col_ptr[%d]=%d  ",i,Col_ptr[i]);
printf("\n");
double Gpr[]={-1,-1,-1,-1,-1,-1,-1/sqrt(2),-1/sqrt(2),-1/sqrt(2),-1/sqrt(2),-1};
int Gir[] = {0,1,2,3,4,5,6,7,6,7,8};
int Gjc[] = {0,1,2,3,4,5,6,8,10,11};
int n_Gir = 11;
int n_Gjc = 10;
int mat_row = MaxInt(Gir,n_Gir)+1;
int mat_col = n_Gjc-1;
printf("mat_row=%d,mat_col=%d\n",mat_row,mat_col);
Matrix* G=MatrixZeros(mat_row,mat_col);
MatrixFromCCS(G,Gpr,Gir,n_Gir,Gjc,n_Gjc);
MatrixPrint(G);
//Matrix* G_TMP=MatrixZeros(mat_row,mat_col);
//MatrixSetElement(G_TMP,0,6,100);
//MatrixFromCCS(G_TMP,Gpr,Gir,n_Gir,Gjc,n_Gjc);
//MatrixPrint(G_TMP);

//Test interp1spline function
//
double x[]={0.52,8,17.95,28.65,50.65,104.6,156.6,260.7,364.4,468,507,520};
double y[]={5.28794,13.84,20.2,24.9,31.1,36.5,36.6,31,20.9,7.8,1.5,0.2};
//double x[]={0.52,8,17.95};
//double y[]={5.28794,13.84,20.2};
double ti[8]={4,14,30,60,130,230,450,555};
double z[8];
interp1(x,y,12,ti,8,z,"spline");
//interp1spline(x,y,12,ti,8,z);
//interp1linear(x,y,3,ti,4,z);
for(int i=0;i<8;i++)
  printf("z[%d]=%f\n",i,z[i]);

////Matrix* mat1=MatrixInit(mat1,3,3);
//double row = 3;
//double col = 2;
//Matrix* matA=MatrixInit(row,col);
////Matrix* matA=MatrixInit(3,2);
//MatrixSetData(matA,valA);
//Matrix* matB=MatrixInit(2,3);
//MatrixSetData(matB,valB);
//Matrix* matC=MatrixInit(3,3);
//matC = MatrixMulti(matA,matB);
//
//Matrix* matadd = MatrixAdd(matC,matC);
//MatrixPrint(matA);
//printf("--------\n");
//MatrixPrint(matB);
//printf("--------\n");
//MatrixPrint(matC);
//printf("--------\n");
//MatrixPrint(matadd);
//printf("--------\n");
//Matrix* matAT=MatrixTrans(matA);
//printf("--------A Trans------\n");
//MatrixPrint(matAT);

//MatrixInit(&mat1,2,2);
//MatrixInit(pmat,2,2);
////MatrixTranspose(&mat1,&mat2);
//printf("%f\n",pmat->data[1][1]);
//printf("%d\n",pmat->init);

void func(double t, double *s, double *ds, double *user);
int flag;
flag = rk4fixstep(t,h,s,user,n,func);
//flag = RungeKutta(t,h,s,user,n,func);
printf("flag = %d\n",flag);
printf("y[this] = %f\n",s[0]);
printf("y[next] = %f\n",s[1]);
double ttin = 23.445;
double tt;
tt = Deg2Rad(ttin);
printf("t[next] = %f\n",tt);
}
void func(double t, double *s, double *ds, double *user){
  double theta;
  double V;
  theta = user[0];
  V = user[1];
  ds[0] = V*cos(theta);
  ds[1] = V*sin(theta);
}

