#include"Matrix.h"
#include<math.h> 
#include<string.h>
#define NDEBUG
#include <assert.h> /* 此时断言功能关闭 */
#include<stdio.h>
//#define FALSE 0
//#define TRUE 1
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CMatrix::CMatrix()
{
	m_nNumRows    = 1;
	m_nNumColumns = 1;
	m_pData       = NULL;
    Init(m_nNumRows, m_nNumColumns);
}
// 指定行列构造函数
//
// 参数：
// 1. int nRows - 指定的矩阵行数
// 2. int nCols - 指定的矩阵列数
CMatrix::CMatrix(int nRows, int nCols)
{
	m_nNumRows    = nRows;
	m_nNumColumns = nCols;
	m_pData       = NULL;
    Init(m_nNumRows, m_nNumColumns);
}
// 指定值构造函数
//
// 参数：
// 1. int nRows - 指定的矩阵行数
// 2. int nCols - 指定的矩阵列数
// 3. double value[] - 一维数组，长度为nRows*nCols,存储矩阵各元素的值
CMatrix::CMatrix(int nRows, int nCols, double value[])
{
	m_nNumRows    = nRows;
	m_nNumColumns = nCols;
	m_pData       = NULL;
    Init(m_nNumRows, m_nNumColumns);

	SetData(value);
}
// 方阵构造函数
//
// 参数：
// 1. int nSize - 指定方阵列数
CMatrix::CMatrix(int nSize)
{
	m_nNumRows    = nSize;
	m_nNumColumns = nSize;
	m_pData       = NULL;
    Init(m_nNumRows, m_nNumColumns);
}
// 方阵构造函数
//
// 参数：
// 1. int nSize - 指定方阵列数
// 2. double value[] - 一维数组，长度为nSize*nSize,存储方阵各元素的值
CMatrix::CMatrix(int nSize, double value[])
{
	m_nNumRows    = nSize;
	m_nNumColumns = nSize;
	m_pData       = NULL;
    Init(m_nNumRows, m_nNumColumns);

	SetData(value);
}

// 拷贝构造函数
//
// 参数：
// 1. const CMatrix& other - 源矩阵
CMatrix::CMatrix(const CMatrix& other)
{
    m_nNumColumns = other.GetNumColumns();
    m_nNumRows    = other.GetNumRows();
	m_pData       = NULL;
    Init(m_nNumRows, m_nNumColumns);

	memcpy(m_pData, other.m_pData,sizeof(double)*m_nNumRows*m_nNumColumns);
	printf("copied Construction done!");
}
// 析构函数
CMatrix::~CMatrix()
{
	if(m_pData)
	{
		delete []m_pData;
		m_pData=NULL;
		printf("delete Matrix\n");
	}
}
///////////////////////////////////////////////////////////////////////////////
// 初始化函数
// 
// 参数：
// 1. int nRows - 指定的矩阵行数
// 2. int nCols - 指定的矩阵列数
//
// 返回值：bool 型，初始化是否成功
bool CMatrix::Init(int nRows, int nCols)
{
	if(m_pData)
	{
		delete []m_pData;
		m_pData = NULL;
	}

	m_nNumRows    = nRows;
	m_nNumColumns = nCols;
	
	int nSize = nCols*nRows;

	if(nSize < 0)
		return false;

	m_pData = new double [nSize];  //分配内存

	if(m_pData == NULL)
		return false;

	//if(IsBadReadPtr(m_pData, sizeof(double)*nSize))
	//	return FALSE;

	memset(m_pData, 0, sizeof(double)*nSize);

	return true;
} 

bool CMatrix::MakeUnitMatrix(int nSize)
{
	if( !Init(nSize,nSize) )
		return false;

	for(int i=0; i<nSize; ++i)
		for(int j=0; j<nSize; ++j)
			if(i == j)
				SetElement(i,j,1);

	return true;
}
////////////////////////////////////////////////////////////////////////////////
// 设置矩阵各元素的值
// 
// 参数：
// 1. doulbe value[] - 一维数组，长度为m_nNumColumns*m_nNumRows,存储
//                     矩阵各元素的值
// 2. int nCols - 指定的矩阵列数
//
// 返回值：无
void CMatrix::SetData(double value[])
{
    //empty  the memory
    memset(m_pData,0,sizeof(double)*m_nNumColumns*m_nNumRows);
    //copy data
    memcpy(m_pData,value,sizeof(double)*m_nNumColumns*m_nNumRows);
}
////////////////////////////////////////////////////////////////////////////////
// 设置指定元素的值
// 
// 参数：
// 1. int nRows - 指定的矩阵行数
// 2. int nCols - 指定的矩阵列数
// 3. double value - 指定元素的值
//
// 返回值：bool 型，说明设置是否成功
bool CMatrix::SetElement(int nRow, int nCol, double value)
{
	if( (nRow < 0) || (nRow >= m_nNumRows) || (nCol < 0) || (nCol >= m_nNumColumns) )
		return false;   //array bounds error
	
	if(m_pData == NULL)
		return false;   // bad pointer error

	m_pData[nCol+nRow*m_nNumColumns] = value;

	return true;
}
////////////////////////////////////////////////////////////////////////////////
// 获取矩阵的数据
// 
// 参数：无
//
// 返回值：double型指针，指向矩阵各元素的数据缓冲区
double* CMatrix::GetData()const
{
	return m_pData;
}
///////////////////////////////////////////////////////////////////////////////
// 获取指定元素的值
// 
// 参数：
// 1. int nRows - 指定的矩阵行数
// 2. int nCols - 指定的矩阵列数
//
// 返回值：double 型，指定元素的值
double CMatrix::GetElement(int nRow, int nCol) const
{
	assert( (nRow >=0) && (nRow < m_nNumRows) && (nCol >= 0) && (nCol < m_nNumColumns) );
	// array bounds error   
	assert( m_pData != NULL ); // bad pointer error

	return m_pData[nCol+nRow*m_nNumColumns];    
}
///////////////////////////////////////////////////////////////////////////////
// 获取指定行的向量
// 
// 参数：
// 1. int nRows - 指定的矩阵行数
// 2. double *pVector - 指向向量中各元素的缓冲区
//
// 返回值：int 型，向量中各元素的个数，即矩阵的列数
int CMatrix::GetRowVector(int nRow, double* pVector) const
{
	if(pVector == NULL)
		delete pVector;
	
	pVector = new double[m_nNumColumns];
	assert(pVector != NULL);

	for(int j=0; j<m_nNumColumns; ++j)
		pVector[j]=GetElement(nRow,j);

	return m_nNumColumns;    
}
///////////////////////////////////////////////////////////////////////////////
// 获取指定列的向量
// 
// 参数：
// 1. int nCols - 指定的矩阵列数
// 2. double *pVector - 指向向量中各元素的缓冲区
//
// 返回值：int 型，向量中各元素的个数，即矩阵的行数
int CMatrix::GetColVector(int nCol, double* pVector) const
{
	if(pVector == NULL)
		delete pVector;
	
	pVector = new double[m_nNumRows];
	assert(pVector != NULL);

	for(int i=0; i<m_nNumRows; ++i)
		pVector[i]=GetElement(i,nCol);

	return m_nNumRows;    
}
///////////////////////////////////////////////////////////////////////////////
// 获取矩阵的列数
// 
// 参数：无
//
// 返回值：int 型，矩阵的列数
int CMatrix::GetNumColumns() const
{
    return m_nNumColumns;
}
///////////////////////////////////////////////////////////////////////////////
// 获取矩阵的行数
// 
// 参数：无
//
// 返回值：int 型，矩阵的行数
int CMatrix::GetNumRows() const
{
    return m_nNumRows;
}

////////////////////////////////////////////////////////////////////////////////
/// 数学运算符重载
///////////////////////////////////////////////////////////////////////////////
// 重载运算符=，给矩阵赋值
// 
// 参数：
// 1. const CMatrix& other - 用于给矩阵幅值的源矩阵
// 返回值：CMatrix型的引用，所引用的矩阵与other相等
CMatrix& CMatrix::operator =(const CMatrix& other)
{
	if(&other != this)
	{
		m_nNumRows    = other.GetNumRows();
		m_nNumColumns = other.GetNumColumns();

		Init(m_nNumRows, m_nNumColumns);
       // copy the pointer
		memcpy(m_pData, other.m_pData,sizeof(double)*m_nNumRows*m_nNumColumns);
	}
      // finally return a reference to ourselves
	return *this;     
}
////////////////////////////////////////////////////////////////////////////////
// 重载运算符==，判断矩阵是否相等
// 
// 参数：
// 1. const CMatrix& other - 用于比较的矩阵
// 返回值：bool型，两个矩阵相等则为ture,否则为false
bool CMatrix::operator==(const CMatrix& other) const
{
    // 首先检查行列数是否相等
	if( (m_nNumRows != other.GetNumRows()) || (m_nNumColumns != other.GetNumColumns()) )
		return false;

	for(int i=0; i<m_nNumRows; ++i)
	{
		for(int j=0; j<m_nNumColumns; ++j)
		{
			if( GetElement(i,j) != other.GetElement(i,j) )
				return false;
		}
	}

	return true;    
}
////////////////////////////////////////////////////////////////////////////////
// 重载运算符==，判断矩阵是否不相等
// 
// 参数：
// 1. const CMatrix& other - 用于比较的矩阵
// 返回值：bool型，两个矩阵不相等则为ture,否则为false
bool CMatrix::operator!=(const CMatrix& other) const
{
    return !(*this == other);
}

////////////////////////////////////////////////////////////////////////////////
// 重载运算符+，实现矩阵加法
// 
// 参数：
// 1. const CMatrix& other - 与指定矩阵相加的矩阵
// 返回值：CMatrix型，指定矩阵与other相加之和
CMatrix CMatrix::operator +(const CMatrix& other) const
{
    // 首先检查行列数是否相等
	assert( (m_nNumRows==other.GetNumRows()) && (m_nNumColumns==other.GetNumColumns()) ) ;
    // 构造结果矩阵
	CMatrix result(*this);//拷贝构造
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

////////////////////////////////////////////////////////////////////////////////
// 重载运算符-，实现矩阵减法
// 
// 参数：
// 1. const CMatrix& other - 与指定矩阵相减的矩阵
// 返回值：CMatrix型，指定矩阵与other相加之差
CMatrix  CMatrix::operator -(const CMatrix& other) const
{
    // 首先检查行列数是否相等
	assert( (m_nNumRows==other.GetNumRows()) && (m_nNumColumns==other.GetNumColumns()) ) ;
    // 构造结果矩阵
	CMatrix result(*this);
    // 矩阵加法 
	for(int i=0; i<m_nNumRows; ++i)
	{
		for(int j=0; j<m_nNumColumns; ++j)
		{
			result.SetElement(i,j,result.GetElement(i,j)-other.GetElement(i,j));
		}
	}

	return result;    
}

////////////////////////////////////////////////////////////////////////////////
// 重载运算符*，实现矩阵数乘
// 
// 参数：
// 1. double value - 与指定矩阵相乘的实数
// 返回值：CMatrix型，指定矩阵与value相乘之积
CMatrix CMatrix::operator *(double value) const
{
    // 构造目标矩阵
	CMatrix result(*this);
    // 进行数乘
	for(int i=0; i<m_nNumRows; ++i)
	{
		for(int j=0; j<m_nNumColumns; ++j)
		{
			result.SetElement(i,j,result.GetElement(i,j)*value);
		}
	}

	return result;    
}
////////////////////////////////////////////////////////////////////////////////
// 重载运算符*，实现矩阵乘法
// 
// 参数：
// 1. const CMatrix& other - 与指定矩阵相乘的矩阵
// 返回值：CMatrix型，指定矩阵与other相乘之积
CMatrix CMatrix::operator *(const CMatrix& other) const
{
    // 首先检查行列数是否符合要求
    assert(m_nNumColumns == other.GetNumRows());
    // construct the object we are going to return
	CMatrix result(m_nNumRows,other.GetNumColumns());
	//
	//
	//
	//
	//
	double value;

	for(int i=0; i<result.GetNumRows(); ++i)
	{
		for(int j=0; j<other.GetNumColumns(); ++j)
		{
			value =0.;

			for(int k=0; k<m_nNumColumns; ++k)
			{
				value += GetElement(i,k)*other.GetElement(k,j);
			} 

			result.SetElement(i,j,value);
		}
	} 
	
	return result;
}
////////////////////////////////////////////////////////////////////////////////
// 重载运算符+=，实现矩阵+=
// 
// 参数：
// 1. const CMatrix& other - 与指定矩阵+=的矩阵
// 返回值：CMatrix型，指定矩阵与other相加之和
void CMatrix::operator+=(const CMatrix& other)
{
	assert( (m_nNumRows==other.GetNumRows()) && (m_nNumColumns==other.GetNumColumns()) ) ;
	
	for(int i=0; i<m_nNumRows; ++i)
	{
		for(int j=0; j<m_nNumColumns; ++j)
		{
			SetElement(i,j, GetElement(i,j)+other.GetElement(i,j));
		}
	}    
}
////////////////////////////////////////////////////////////////////////////////
// 实现矩阵开方
// 
// 参数：无
// 返回值：CMatrix型，矩阵开方结果
CMatrix CMatrix::Sqrt() const
{
	CMatrix Sq(m_nNumColumns, m_nNumRows);

	for(int i=0; i<m_nNumRows; ++i)
	{
		for(int j=0; j<m_nNumColumns; ++j)
		{
			double old   = GetElement(i,j);
			double value = ::sqrt(old);//::全局作用的函数
			Sq.SetElement(i, j, value);
		}
	}

	return Sq;    
}
////////////////////////////////////////////////////////////////////////////////
// 实现矩阵指数运算
// 
// 参数：double p
// 返回值：CMatrix型，矩阵pow结果
CMatrix CMatrix::Pow(double p) const
{
	CMatrix Pw(m_nNumColumns, m_nNumRows);

	for(int i=0; i<m_nNumRows; ++i)
	{
		for(int j=0; j<m_nNumColumns; ++j)
		{
			double old   = GetElement(i,j);
			double value = std::pow(old,p);//::全局作用的函数
			Pw.SetElement(i, j, value);
		}
	}

	return Pw;    
}

////////////////////////////////////////////////////////////////////////////////
// 实现矩阵的转置
// 
// 参数：无
// 返回值：CMatrix型，矩阵转置结果
CMatrix CMatrix::Transpose() const
{
	CMatrix Trans(m_nNumColumns, m_nNumRows);
    //CMatrix Trans(CMatrix(m_nNumColumns, m_nNumRows));
	for(int i=0; i<m_nNumRows; ++i)
	{
		for(int j=0; j<m_nNumColumns; ++j)
		{
			Trans.SetElement(j,i,GetElement(i,j));
		}
	}

	return Trans;   
// 	return *this; // 调用拷贝构造函数了！
}

///////////////////////////////////////////////////////////////////////////////
CMatrix CMatrix::GetTransposed() const
{
	CMatrix	transposed(*this) ;		// make a copy of ourselves 调用复制构造函数！

	//transposed.Transpose();

	return transposed.Transpose();    
}

////////////////////////////////////////////////////////////////////////////////
// 求矩阵的迹
// 
// 参数：无
// 返回值：double型，矩阵的迹
double CMatrix::Trace() const
{
	assert(m_nNumColumns == m_nNumRows);

    double tr = 0.;
    for(int i = 0; i < m_nNumRows; i++)
		tr += GetElement(i,i);

	return tr;    
}

////////////////////////////////////////////////////////////////////////////////
// 实现矩阵的逆
// 
// 参数：无
// 返回值：CMatrix型，矩阵转置结果
CMatrix CMatrix::GetInverted() const
{
	// matrix inversion will only work on square matrices
	if (m_nNumColumns != m_nNumRows)
		throw "Matrix must be square." ;
	// return this matrix inverted
	
	CMatrix	copy(*this) ;   // make a copy of ourselves 调用复制构造函数！
	copy.Invert() ;

	return copy ;    
}
////////////////////////////////////////////////////////////////////////////////
// 实现矩阵的逆
// 
// 参数：无
// 返回值：bool型，矩阵转置是否成功
bool CMatrix::Invert()
{
	int *pnRow, *pnCol, i,j,k,l,u,v;
	double d=0, p=0;

	pnRow = new int[m_nNumColumns];
	pnCol = new int[m_nNumRows];
	
	if(pnRow == NULL || pnCol ==NULL)
		return false;

    for(k=0; k<=m_nNumColumns-1; k++)
	{
		d=0.;

		for(i=k; i<=m_nNumColumns-1; i++)
		{
			for(j=k; j<=m_nNumColumns-1; j++)
			{
				l=i*m_nNumColumns+j;
				p=fabs(m_pData[l]);

				if(p>d)
				{
					d=p;
					pnRow[k]=i;
					pnCol[k]=j;
				}
			}
		}

		if(d == 0.)
		{
			delete[] pnRow;
			delete[] pnCol;
			return false;
		}

		if(pnRow[k] != k)
		{
			for(j=0; j<=m_nNumColumns-1; j++)
			{
				u=k*m_nNumColumns+j;
				v=pnRow[k]*m_nNumColumns+j;
				p=m_pData[u];
				m_pData[u]=m_pData[v];
				m_pData[v]=p;
			}
		}

		if(pnCol[k] != k)
		{
			for(i=0; i<=m_nNumColumns-1; i++)
			{
				u=i*m_nNumColumns+k;
				v=i*m_nNumColumns+pnCol[k];
				p=m_pData[u];
				m_pData[u]=m_pData[v];
				m_pData[v]=p;
			}
		}

		l=k*m_nNumColumns+k;
		m_pData[l]=1.0/m_pData[l];
		
		for(j=0; j<=m_nNumColumns-1; j++)
		{
			if(j != k)
			{
				u=k*m_nNumColumns+j;
				m_pData[u]=m_pData[u]*m_pData[l];
			}
		}

		for(i=0; i<=m_nNumColumns-1; i++)
		{
			if(i != k)
			{
				for(j=0; j<=m_nNumColumns-1; j++)
				{
					if(j != k)
					{
				       u=i*m_nNumColumns+j;
 				       m_pData[u]=m_pData[u]-m_pData[i*m_nNumColumns+k]*m_pData[k*m_nNumColumns+j];
					}
				}
			}
		}

		for(i=0; i<=m_nNumColumns-1; i++)
		{
			if(i != k)
			{
				u=i*m_nNumColumns+k;
				m_pData[u]=-m_pData[u]*m_pData[l];
			}
		}
	
	}

	for(k=m_nNumColumns-1; k>=0; k--)
	{
		if(pnCol[k] != k)
		{
			for(j=0; j<=m_nNumColumns-1; j++)
			{
				u=k*m_nNumColumns+j;
				v=pnCol[k]*m_nNumColumns+j;
				p=m_pData[u];
				m_pData[u]=m_pData[v];
				m_pData[v]=p;
			}
		}

		if(pnRow[k] != k)
		{
			for(i=0; i<=m_nNumColumns-1; i++)
			{
				u=i*m_nNumColumns+k;
				v=i*m_nNumColumns+pnRow[k];
				p=m_pData[u];
				m_pData[u]=m_pData[v];
				m_pData[v]=p;
			}
		}
	}

	delete[] pnRow;
	delete[] pnCol;
	
	return true;
}

////////////////////////////////////////////////////////////////////////////////
// 求矩阵的秩
// 
// 参数：无
// 返回值：int型，矩阵的秩
int CMatrix::Rank() const
{
    double *a = m_pData;
    int m = m_nNumRows;
    int n = m_nNumColumns;
    int i,j,k,nn,is,js,l,ll,u,v;
    double q,d;
    nn=m;
    if (m>=n) nn=n;
    k=0;
    for (l=0; l<=nn-1; l++)
      { q=0.0;
        for (i=l; i<=m-1; i++)
        for (j=l; j<=n-1; j++)
          { ll=i*n+j; d=fabs(a[ll]);
	    if (d>q) { q=d; is=i; js=j;}
          }
        if (q+1.0==1.0) return(k);
        k=k+1;
        if (is!=l)
          { for (j=l; j<=n-1; j++)
              { u=l*n+j; v=is*n+j;
                d=a[u]; a[u]=a[v]; a[v]=d;
              }
          }
        if (js!=l)
          { for (i=l; i<=m-1; i++)
              { u=i*n+js; v=i*n+l;
                d=a[u]; a[u]=a[v]; a[v]=d;
              }
          }
        ll=l*n+l;
        for (i=l+1; i<=n-1; i++)
          { d=a[i*n+l]/a[ll];
            for (j=l+1; j<=n-1; j++)
              { u=i*n+j;
                a[u]=a[u]-d*a[l*n+j];
              }
          }
      }
    return(k);

}

