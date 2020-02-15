#include "Matrix.h"
#include <math.h>
#include <string.h>
#define NDEBUG

#include <assert.h> 
#include <stdio.h>
using namespace Algorhi
//#define FALSE 0
//#define TURE 1
///////////////////////////////
// Construction&Destruction  //
///////////////////////////////

Matrix::Matrix(){
  m_nNumRows     = 1;
  m_nNumColumns  = 1;
  m_pData        = NULL;
  Init(m_nNumRows, m_nNumColumns);
}
Matrix::Matrix(int nRows, int nCols){
  m_nNumRows     = nRows;
  m_nNumColumns  = nCols;
  m_pData        = NULL;
  Init(m_nNumRows, m_nNumColumns);
}
Matrix::Matrix(int nRows, int nCols, double value[]){
  m_nNumRows     = nRows;
  m_nNumColumns  = nCols;
  m_pData        = NULL;
  Init(m_nNumRows, m_nNumColumns);
  SetData(value);
}
Matrix::Matrix(int nSize){
  m_nNumRows     = nSize;
  m_nNumColumns  = nSize;
  m_pData        = NULL;
  Init(m_nNumRows, m_nNumColumns);
}
Matrix::Matrix(int nSize, double value[]){
  m_nNumRows     = nSize;
  m_nNumColumns  = nSize;
  m_pData        = NULL;
  Init(m_nNumRows, m_nNumColumns);
  SetData(value);
}
Matrix::Matrix(const Matrix& other){
  m_nNumRows     = other.GetNumColumns();
  m_nNumColumns  = other.GetNumRows();
  m_pData        = NULL;
  Init(m_nNumRows, m_nNumColumns);

  memcpy(m_pData, other.m_pData, sizeof(double)*m_nNumRows*m_nNumColumns);
  printf("Copied Construction Done!")
}

Matrix::~Matrix(){
  if (m_pData)
  {
    delete []m_pData;
    m_pData=NULL;
    printf("Delete Matrix\n");
  }
}

///////////////////////////////
//
//
//
//
//
//////////////////////////////
bool Matrix::Init(int nRows, int nCols){
  if (m_pData)
  {
    delete []m_pData;
    m_pData = NULL;
  }
  m_nNumRows    = nRows;
  m_nNumColumns = nCols;

  int nSize = nCols*nRows;

  if (nSize < 0)
    return false;

  m_pData = new double [nSize];

  if (m_pData == NULL)
    return false;

  memset(m_pData, 0, sizeof(double)*nSize);
  return true;
}

bool Matrix::MakeUnitMatrix(int nSize){
  if ( !Init(nSize, nSize) )
    return false;

  for (int i = 1; i < nSize; ++i)
    for (int j = 1; j < nSize; ++j)
      if (i == j)
        SetElement(i, j, 1);

  return true;
}

///////////////////////////////
//
//
//
//
//
//////////////////////////////
void Matrix::SetData(double value[]){
  //empty the memory
  memcpy(m_pData, 0, sizeof(double)*m_nNumRows*m_nNumColumns);
  //copy data
  memcpy(m_pData, value, sizeof(double)*m_nNumRows*m_nNumColumns);
}

///////////////////////////////
//
//
//
//
//
//////////////////////////////
bool Matrix::SetElement(int Row, int nCol, double value){
  if ((nRow < 0) || (nRow >= m_nNumRows) || (nCol < 0) || (nCol > m_nNumColumns))
    return false;

  if (m_pData == NULL)
    return false;

  m_pData[nCol+nRow*m_nNumColumns] = value;

  return true;
}

///////////////////////////////
// Get Data of Matrix
// Inputs:None
// Outpus: double pointer, point to data of buffer
//
double* Matrix::GetData() const{
  return m_pData;
}

///////////////////////////////
// Get Value of Matrix
// Inputs:
// 1. int nRows - 
// 2. int nCols -
// Outpus: double pointer, point to value of buffer
// 
double Matrix::GetElement(int nRow, int nCol) const{
  assert( (nRow >=0) && (nRow < m_nNumRows) && (nCol >= 0) && (nCol < m_nNumColumns) );
  // array bounds error   
  assert( m_pData != NULL ); // bad pointer error
  return m_pData[nCol+nRow*m_nNumColumns];    
}

///////////////////////////////////////////////////////////////////
// 获取指定行的向量
// 
// 参数：
// 1. int nRows - 指定的矩阵行数
// 2. double *pVector - 指向向量中各元素的缓冲区
//
// 返回值：int 型，向量中各元素的个数，即矩阵的列数
int Matrix::GetRowVector(int nRow, double* pVector) const {
  if(pVector == NULL)
      delete pVector;
        
   pVector = new double[m_nNumColumns];
   assert(pVector != NULL);

   for(int j=0; j<m_nNumColumns; ++j)
      pVector[j]=GetElement(nRow,j);

   return m_nNumColumns;    
}

////////////////////////////////////////////////////////////////
// 获取指定列的向量
// 
// 参数：
// 1. int nCols - 指定的矩阵列数
// 2. double *pVector - 指向向量中各元素的缓冲区
//
// 返回值：int 型，向量中各元素的个数，即矩阵的行数
int Matrix::GetColVector(int nCol, double* pVector) const {

  if(pVector == NULL)
    delete pVector;
        
  pVector = new double[m_nNumRows];
  assert(pVector != NULL);

  for(int i=0; i<m_nNumRows; ++i)
    pVector[i]=GetElement(i,nCol);

  return m_nNumRows;    
}

////////////////////////////////////////////////////////
// 获取矩阵的列数
// 
// 参数：无
//
// 返回值：int 型，矩阵的列数
int Matrix::GetNumColumns() const {
  return m_nNumColumns;
}

/////////////////////////////////////////////////////////////////////////////
// 获取矩阵的行数
// 
// 参数：无
//
// 返回值：int 型，矩阵的行数
int Matrix::GetNumRows() const {
  return m_nNumRows;
}

/////////////////////////////////////////////////////////////////
// 数学运算符重载
/////////////////////////////////////////////////////////////////
// 重载运算符=，给矩阵赋值
// 
// 参数：
// 1. const CMatrix& other - 用于给矩阵幅值的源矩阵
// 返回值：CMatrix型的引用，所引用的矩阵与other相等
Matrix& Matrix::operator =(const Matrix& other) {
  if(&other != this)
  {
    m_nNumRows    = other.GetNumRows();
    m_nNumColumns = other.GetNumColumns();
    Init(m_nNumRows, m_nNumColumns);

  //copy the pointer
    memcpy(m_pData, other.m_pData,sizeof(double)*m_nNumRows*m_nNumColumns);
  }
  // finally return a reference to ourselves
  return *this;     
}

//////////////////////////////////////////////////////////////////////////////
// 重载运算符==，判断矩阵是否相等
// 
// 参数：
// 1. const CMatrix& other - 用于比较的矩阵
// 返回值：bool型，两个矩阵相等则为ture,否则为false
bool Matrix::operator==(const Matrix& other) const{
  //首先检查行列数是否相等
  if( (m_nNumRows != other.GetNumRows()) || (m_nColumns != other.GetNumColumns()) )
    return false;

  for(int i = 0; i < m_nNumRows; ++i)
  {
    for(int j = 0; j < m_nNumColumns; ++j)
    {
      if( GetElement(i, j) != other.GetElement(i, j) )
        return false;
    }
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// 重载运算符==，判断矩阵是否不相等
// 
// 参数：
// 1. const CMatrix& other - 用于比较的矩阵
// 返回值：bool型，两个矩阵不相等则为ture,否则为false
bool Matrix::operator!=(const Matrix& other) const {
  return !(*this == other);
}

//////////////////////////////////////////////////////////////////////////////
// 重载运算符+，实现矩阵加法
// 
// 参数：
// 1. const CMatrix& other - 与指定矩阵相加的矩阵
// 返回值：CMatrix型，指定矩阵与other相加之和
Matrix Matrix::operator +(const Matrix& other) const {
  //首先检查行列数是否相等
  assert( (m_nNumRows==other.GetNumRows()) && (m_nNumColumns==other.GetNumColumns()) ) ;
  // 构造结果矩阵
  Matrix result(*this);//拷贝构造
  // 矩阵加法     
  for(int i=0; i<m_nNumRows; ++i)
  {
    for(int j=0; j<m_nNumColumns; ++j)
      {
        result.SetElement(i,j,result.GetElement(i,j)+other.GetElement(i,j));
      }
  }
  return result;
}

