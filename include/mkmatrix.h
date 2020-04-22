/*
 * mathkit matrix: C code implementation for basic matrix operation
 * 2020 by Algorhi.com
 *
 */
#ifndef __MKMATRIX__H__
#define __MKMATRIX__H__

#ifdef __cplusplus
extern "C" {
#endif

#define MAT_INIT_FLAG 0x5C

typedef struct{
  int row,col;
  double *pData;
  unsigned char init;
}Matrix;
//Matrix* MatrixInit(Matrix* mat, int row, int col);
Matrix* MatrixInit(int nRow, int mCol);
Matrix* MatrixEye(int n);
Matrix* MatrixZeros(int nRow, int mCol);
void MatrixSetData(Matrix* mat, double *array);
void MatrixSetElement(Matrix* mat, int nRow, int mCol, double value);
//void MatrixSetBlock(Matrix* matA, int* RowIndex, int* ColIndex, Matrix* matB);
void MatrixSetBlock(Matrix* matA, int RowBegin, int RowEnd, int ColBegin, int ColEnd, Matrix* matB);
int MatrixSetRowVector(Matrix* mat, int nRow, double* pVector);
int MatrixSetColVector(Matrix* mat, int nCol, double* pVector);
double MatrixGetElement(Matrix* mat, int nRow, int mCol);
//void MatrixGetBlock(Matrix* matA, int* RowIndex, int* ColIndex, Matrix* matB);
void MatrixGetBlock(Matrix* matA, int RowBegin, int RowEnd, int ColBegin, int ColEnd, Matrix* matB);
int MatrixGetRowVector(Matrix* mat, int nRow, double* pVector);
int MatrixGetColVector(Matrix* mat, int nCol, double* pVector);
void MatrixFree(Matrix* mat);
void MatrixCopy(Matrix* src,Matrix* dst);
void MatrixPrint(const Matrix* mat);
void MatrixBlockPrint(const Matrix* mat, int RowBegin, int RowEnd, int ColBegin, int ColEnd);
Matrix* MatrixMulti(Matrix* matA, Matrix* matB);
Matrix* MatrixNumMulti(Matrix* matA, double value);
Matrix* MatrixAdd(Matrix* matA, Matrix* matB);
Matrix* MatrixTrans(Matrix* src);

int MaxInt(int* array, int n);
double MaxDouble(double* array, int n);
int MinInt(int* array, int n);
double MinDouble(double* array, int n);
/*
 * Matrix Compressed Storage
 */
int MatrixNumNonzero(Matrix* mat);
void Matrix2CCS(Matrix* mat, double* pValue, int* pRow_ind, int m, int* pCol_ptr, int n);
void MatrixFromCCS(Matrix* mat, double* pValue, int* pRow_ind, int m, int* pCol_ptr, int n);
#ifdef __cplusplus
}
#endif
#endif
