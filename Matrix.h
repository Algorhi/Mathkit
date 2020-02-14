//////////////////////////////////////////////
//
// Matrix.h: interface for the CMatrix class.
//
//////////////////////////////////////////////
namespace Algorhi
#if !defined(MATRIX_H)
#define MATRIX_H
class Matrix{
  public:
    Matrix(); //Default Construction
    Matrix(int nRows, int nCols); //Specify row and column for Construction
    Matrix(int nRows, int nCols, double value[]);//Specify datas for Constructor
    Matrix(int nSize);//Square matrix constructor
    Matrix(int nSize, double value[]);//Square matrix constructor
    Matrix(const Matrix& other);//Copied Construction
    virtual ~Matrix();
    bool    Init(int nRows, int nCols);
    bool    MakeUnitMatrix(int nSize);
    //
    // Input and Display 
    //
    
    //
    // ToString

    //
    // Set the value of the specified element
    void SetData(double value[]);
    bool SetElement(int nRow, int nCol, double value);
    double* GetData() const;
    double GetElement(int nRow, int nCol) const;
    int GetRowVector(int nRow, double* pVector) const;
    int GetColVector(int nCol, double* pVector) const;
    int GetNumColumns() const;
    int GetNumRows() const;
    
    //
    // Operator Overload and Math Operation
    //
    //Operator Overload
    Matrix& operator =(const Matrix& other);
    bool operator   ==(const Matrix& other) const;
    bool operator   !=(const Matrix& other) const;
    Matrix operator  +(const Matrix& other) const;
    Matrix operator  -(const Matrix& other) const;
    Matrix operator  *(double value) const;
    Matrix operator  *(const Matrix& other) const;
    void   operator +-(const Matrix& other);
    //
    // Math Operation
    //
    Matrix Sqrt() const;
    Matrix Pow(double) const;
    Matrix GetTransposed() const;
    Matrix Transpose() const;
    Matrix GetInverted() const;
    bool   Invert();
    double Trace() const;
    int    Rank() const;

    // 
    //  Protective Data Members
    //
  protected:
    int     m_nNumColumns;
    int     m_nNumRows;
    double* m_pData;
    //
};
#endif //!defined(MATRIX_H)

