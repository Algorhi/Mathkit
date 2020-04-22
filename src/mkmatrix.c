#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<assert.h>
#include"mkmatrix.h"


//#define MAT_INIT_FLAG 0x5C
//
//typedef struct {
//  int row, col;
//  double *pData;
//  unsigned char init;
//}Matrix;

int MaxInt(int* array, int n){
  int max = array[0];
  for(int i=0;i<n;i++)
  {
    if (max < array[i])
    max = array[i];
  }
  return max;
}

double MaxDouble(double* array, int n){
  double max = array[0];
  for(int i=0;i<n;i++)
  {
    if (max < array[i])
    max = array[i];
  }
  return max;
}

int MinInt(int* array, int n){
  int min = array[0];
  for(int i=0;i<n;i++)
  {
    if (min > array[i])
    min = array[i];
  }
  return min;
}

double MinDouble(double* array, int n){
  double min = array[0];
  for(int i=0;i<n;i++)
  {
    if (min > array[i])
    min = array[i];
  }
  return min;
}

static void swapab(int *a, int *b){
  int m;
  m = *a;
  *a = *b;
  *b = m;
}
//
//static void perm(int list[], int k, int m, int* p, Matrix* mat, float* det){
//  int i;
//  if(k > m)
//  {
//    float res = mat->pData[0][list[0]];
//    
//    for(i = 1; i < mat->row ; i++){
//      res *= mat->data[i][list[i]];
//    }
//    if(*p%2){
//      //odd is negative
//      *det -= res;
//    }else{
//      //even is positive
//      *det += res;
//    }
//  }
//  else
//  {
//    perm(list, k + 1, m, p, mat, det);
//    for(i = k+1; i <= m; i++)
//    {
//      swap(&list[k], &list[i]); 
//      *p += 1;
//      perm(list, k + 1, m, p, mat, det);
//      swap(&list[k], &list[i]);
//      *p -= 1; 
//    }
//  }
//}


//Matrix* MatrixInit(Matrix* mat, int row, int col){
Matrix* MatrixInit(int row, int col){
  Matrix* mat; 

#ifdef MAT_LEGAL_CHECKING
  if(mat->init == MAT_INIT_FLAG){
    if(mat->row == row && mat->col ==col)
      return mat;
    else
      MatrixFree(mat);
  }
#endif

  if (row>0 && col>0){
    mat = (Matrix*)malloc(sizeof(Matrix));
    mat->row = row;
    mat->col = col;
    mat->init = MAT_INIT_FLAG;
    mat->pData=(double*)malloc(sizeof(double)*row*col);
    return mat;
  }
  else
    return NULL;
  //if (mat->pData == NULL){
  //  return NULL;
  //}
}

Matrix* MatrixEye(int n){
  Matrix* mat;
  if (n>1)
  {
    mat = (Matrix*)malloc(sizeof(Matrix));
    mat->row = n;
    mat->col = n;
    mat->init = MAT_INIT_FLAG;
    mat->pData=(double*)malloc(sizeof(double)*n*n);
    for(int i=0;i<n;i++)
    {
      for(int j=0;j<n;j++)
      {
        if(i==j)
          MatrixSetElement(mat,i,j,1.);
        else
          MatrixSetElement(mat,i,j,0.);
        //printf("%.4f\n",MatrixGetElement(mat,i,j));
      }
    }
    return mat;
  }
  else
    return NULL;
}

Matrix* MatrixZeros(int nRow, int nCol){
  Matrix* mat;
  if (nRow>0 && nCol>0)
  {
    mat = (Matrix*)malloc(sizeof(Matrix));
    mat->row = nRow;
    mat->col = nCol;
    mat->init = MAT_INIT_FLAG;
    mat->pData=(double*)malloc(sizeof(double)*nRow*nCol);
    for(int i=0;i<nRow;i++)
    {
      for(int j=0;j<nCol;j++)
      {
        MatrixSetElement(mat,i,j,0);
      }
    }
    return mat;
  }
  else
    return NULL;
}

void MatrixCopy(Matrix* src, Matrix* dst){
#ifdef MAT_LEGAL_CHECKING
  if (mat->init != MAT_INIT_FLAG){
    return;
  }
#endif
  dst->row = src->row;
  dst->col = src->col;
  int matsize = src->row*src->col;
  memcpy(dst->pData, src->pData, matsize*sizeof(double));
}

void MatrixFree(Matrix* mat){
#ifdef MAT_LEGAL_CHECKING
  if (mat->init != MAT_INIT_FLAG){
    return;
  }
#endif
  free(mat->pData);
  if (mat->pData == NULL)
    printf("Matrix Free!\n");
  else
    printf("Matrix Free defeat!\n");
}

void MatrixSetData(Matrix* mat, double* data){
  if (mat->pData != NULL){
    memcpy(mat->pData, data, mat->row*mat->col*sizeof(double*));
  }
}

void MatrixSetElement(Matrix* mat, int nRow, int mCol, double value){
  assert((nRow>=0) && (nRow<mat->row) && (mCol>=0) && (mCol<mat->col) );
  //array bounds error 
  assert(mat->pData != NULL);
  mat->pData[mCol+nRow*mat->col] = value;
  return;
}

//void MatrixSetBlock(Matrix* matA, int* RowIndex, int* ColIndex, Matrix* matB){
void MatrixSetBlock(Matrix* matA, int RowBegin, int RowEnd, int ColBegin, int ColEnd, Matrix* matB){
  int row0 = RowBegin;
  int rowf = RowEnd;
  int col0 = ColBegin;
  int colf = ColEnd;
  //assert();
  double value_temp;
  for(int i=0;i<matB->row;i++)
  {
    for(int j=0;j<matB->col;j++)
    {
      value_temp = MatrixGetElement(matB,i,j);
      MatrixSetElement(matA,row0+i,col0+j,value_temp);
    }
  }
  return;
}

double MatrixGetElement(Matrix* mat, int nRow, int mCol){
  assert((nRow>=0) && (nRow<mat->row) && (mCol>=0) && (mCol<mat->col) );
  //array bounds error 
  assert(mat->pData != NULL);
  double array;
  memcpy(&array, mat->pData+mCol+nRow*mat->col, sizeof(double));
  return array;
  //return mat->pData[mCol+nRow*mat->col];
}

//void MatrixGetBlock(Matrix* matA, int* RowIndex, int* ColIndex, Matrix* matB){
void MatrixGetBlock(Matrix* matA, int RowBegin, int RowEnd, int ColBegin, int ColEnd, Matrix* matB){
  int row0 = RowBegin;
  int rowf = RowEnd;
  int col0 = ColBegin;
  int colf = ColEnd;
  //assert();
  matB->row = rowf-row0+1;
  matB->col = colf-col0+1;
  double value_temp;
  for(int i=0;i<matB->row;i++)
  {
    for(int j=0;j<matB->col;j++)
    {
      value_temp = MatrixGetElement(matA,row0+i,col0+j);
      MatrixSetElement(matB,i,j,value_temp);
    }
  }
  return;
}

int MatrixGetRowVector(Matrix* mat, int nRow, double* pVector){
  //pVector=(double*)malloc(sizeof(double)*mat->col);
  assert(pVector != NULL);
  for(int j=0; j<mat->col;j++)
    pVector[j]=MatrixGetElement(mat,nRow,j);
  return mat->col;
}

int MatrixGetColVector(Matrix* mat, int mCol, double* pVector){
  //pVector=(double*)malloc(sizeof(double)*mat->row);
  assert(pVector != NULL);
  for(int i=0; i<mat->row;i++)
    pVector[i]=MatrixGetElement(mat,i,mCol);
  return mat->row;
}

int MatrixSetRowVector(Matrix* mat, int nRow, double* pVector){
  assert(pVector != NULL);
  for(int j=0; j<mat->col;j++)
    MatrixSetElement(mat,nRow,j,pVector[j]);
  return mat->col;
}

int MatrixSetColVector(Matrix* mat, int mCol, double* pVector){
  assert(pVector != NULL);
  for(int i=0; i<mat->row;i++)
    MatrixSetElement(mat,i,mCol,pVector[i]);
  return mat->col;
}

void MatrixPrint(const Matrix* mat){
  int matsizes=mat->row*mat->col;
#ifdef MAT_LEGAL_CHECKING
  if(mat == NULL){
    return;
  }
#endif
  printf(" |");
  for(int i=0;i<mat->col;i++) printf("    %d\t    ",i);
  printf("\n");
  //for(int i=0;i<mat->row;i++)
  //{
  //  for(int j=0;j<mat->col;j++)
  //  {
  //    printf("%lf\t",mat-pData[])
  //  }
  //}
  for(int i=0;i<matsizes;i++){
    printf("%lf\t",mat->pData[i]);
    if((i+1)%mat->col == 0)
    {
      printf("\n");
    }
  }
}

void MatrixBlockPrint(const Matrix* mat, int RowBegin, int RowEnd, int ColBegin, int ColEnd){
#ifdef MAT_LEGAL_CHECKING
  if(mat == NULL){
    return;
  }
#endif
  printf(" |");
  for(int p=ColBegin;p<=ColEnd;p++) printf("    %d\t    ",p);
  printf("\n");
  for(int i=RowBegin;i<=RowEnd;i++)
  {
    for(int j=ColBegin;j<=ColEnd;j++)
    {
      double val = MatrixGetElement(mat,i,j);
      printf("%lf\t",val);
    }
    printf("\n");
  }
}


Matrix* MatrixMulti(Matrix* matA, Matrix* matB){
  if(matA->col == matB->row)
  {
    Matrix* matC = MatrixInit(matA->row, matB->col);
    double value;
    for (int i=0;i<matC->row;++i)
    {
      for (int j=0;j<matB->col;++j)
      {
        value = 0.;
        for (int k=0;k<matA->col;++k)
        {
          //matC->pData[i*matC->col+j] += matA->pData[i*matA->col+k] * matB->pData[k*matB->row + j];
          value += MatrixGetElement(matA,i,k)*MatrixGetElement(matB,k,j);
        }
        MatrixSetElement(matC,i,j,value);
      }
    }
    return matC;
  }
  else
  {
    printf("oprate defeat!\n");
    return NULL;
  }
}

Matrix* MatrixNumMulti(Matrix* matA, double value){
  Matrix* dst = MatrixInit(matA->row,matA->col);
  for(int i=0;i<matA->row;++i)
  {
    for(int j=0;j<matA->col;++j)
      MatrixSetElement(dst,i,j,MatrixGetElement(matA,i,j)*value);
  }
  return dst;
}
Matrix* MatrixAdd(Matrix* matA, Matrix* matB){
  assert((matA->row==matB->row) && (matA->col==matB->col));
  Matrix* matC = MatrixInit(matA->row,matA->col);
  for (int i=0;i<matA->row;++i)
    for(int j=0;j<matA->col;++j)
    {
      MatrixSetElement(matC,i,j,MatrixGetElement(matA,i,j)+MatrixGetElement(matB,i,j));
    }

  return matC;
}

Matrix* MatrixTrans(Matrix* src){
  Matrix* dst = MatrixInit(src->col, src->row);
//#ifdef MAT_LEGAL_CHECKING
//  if (src->row != dst->col || src-> col != dst->row ){
//    printf("err check, unmatch matrix for MatrixTranspose.\n");
//    MatrixFree(src);
//    MatrixFree(dst);
//    return NULL;
//  }
//#endif
  for (int i = 0; i < src->row; i++)
    for (int j = 0; j < src->col; j++)
      MatrixSetElement(dst,j,i,MatrixGetElement(src,i,j));

  return dst;
}

/*
 * Compressed Column Storage
 */
 int MatrixNumNonzero(Matrix* mat){
   int Nnonzero = 0;
   for(int i=0;i<mat->row;i++){
     for(int j=0;j<mat->col;j++){
       if (MatrixGetElement(mat,i,j)!=0) Nnonzero++;
     }
   }
   return Nnonzero;

 }

 void Matrix2CCS(Matrix* mat, double* pValue, int* pRow_ind, int m, int* pCol_ptr, int n){
   int Row = mat->row;
   int Col = mat->col;
   *pCol_ptr = 0;
   int vid = 0;
   int pid = 1;
   int ptr = 0;
   double val_ij;
   for (int j=0;j<Col;j++)
   {
     for(int i=0;i<Row;i++)
     {
       val_ij = MatrixGetElement(mat,i,j);
       if(val_ij != 0)
       {
         pValue[vid] = val_ij;
         pRow_ind[vid] = i;
         vid++;
         ptr++;
       }
     }
     pCol_ptr[pid] = ptr;
     pid++;
   }
   if (vid != m) printf("Size of pValue is wrong, that must be %d.\n",vid);
   if (pid != n) printf("Size of pCol_ptr is wrong, that must be %d.\n",pid);
   //double Gpr[11] = { -1,-1,-1,-1,-1,-1,-1 / sqrt(2),-1 / sqrt(2),-1 / sqrt(2),1 / sqrt(2),-1 };
   //for(int i=0;i<11;i++)  pValue[i] = Gpr[i];
 }


void MatrixFromCCS(Matrix* mat, double* pValue, int* pRow_ind, int m, int* pCol_ptr, int n){
  int Row = mat->row;
  int Col = mat->col;
  int CCS_Row = MaxInt(pRow_ind,m)+1;
  int CCS_Col = n-1;
  //printf("=====%d,%d=========\n",CCS_Row,CCS_Col);
  if (CCS_Row!=Row) printf("Row is error, that must be %d.\n",CCS_Row);
  if (CCS_Col!=Col) printf("Col is error, that must be %d.\n",CCS_Col);
  //MatrixZeros(mat,CCS_Row,CCS_Col);
  int ptra,ptrb;
  double setval;
  int seti,setj;
  for(int j=0;j<n-1;j++)
  {
    ptra = *(pCol_ptr+j);
    ptrb = *(pCol_ptr+j+1);
    setj   = j;
    for(int i=ptra;i<ptrb;i++)
    {
      setval = *(pValue+i);
      seti   = *(pRow_ind+i);
      //printf("ptra=%d,ptrb=%d\n",ptra,ptrb);
      //printf("setval=%f,seti=%d,setj=%d\n",setval,seti,setj);
      MatrixSetElement(mat,seti,setj,setval);
    }
  }
}
