/////////////////////////////////////////////////////////////////
//
//   Matrix.h: interface for the CMatrix class.
//
////////////////////////////////////////////////////////////////
                                                             //
#if !defined(MATRIX_H)                                      //
#define MATRIX_H                                           //
class CMatrix                                             //
{                                                        //
    public:                                             //
    CMatrix();                                         // Default Construction
    CMatrix(int nRows,int nCols);                     // Specify row and column for Construction
    CMatrix(int nRows, int nCols, double value[]);   // Specify datas for constructor
	CMatrix(int nSize);                             // Square matrix constructor
	CMatrix(int nSize, double value[]);            // Specify row and column for Square Matrix Construction
	CMatrix(const CMatrix& other);                // Copied Construction
	virtual ~CMatrix();                          // Destruction
	bool    Init(int nRows, int nCols);         // Initialization Matrix
	bool    MakeUnitMatrix(int nSize);         // Initialize a square matrix into a unit matrix
	//                                                     //
	//  Input and Display                                  //
	//                                                     //  
	//ToString
	
	
	//                                                     //
	// 	Set the value of the specified element             //
	void    SetData(double value[]);                       //
    bool    SetElement(int nRow, int nCol, double value);  //
	double* GetData()const;                                //       
	double  GetElement(int nRow, int nCol) const;
    int     GetRowVector(int nRow, double* pVector)const;
    int     GetColVector(int nCol, double* pVector)const;
	int     GetNumColumns()const;
	int     GetNumRows()const;
	
	//                                                    //
	//  Operator Overload and Math Operation              //
	//                                                    //
	//  Operator Overload
    CMatrix& operator =(const CMatrix& other);
	bool     operator==(const CMatrix& other)const;
	bool     operator!=(const CMatrix& other)const;
	CMatrix  operator +(const CMatrix& other)const;
	CMatrix  operator -(const CMatrix& other)const;
	CMatrix  operator *(double value)const;
	CMatrix  operator *(const CMatrix& other)const;
    void     operator+=(const CMatrix& other);	
	//  Math Operation
	CMatrix Sqrt()const;
	CMatrix Pow(double)const;
    CMatrix GetTransposed()const;
	CMatrix Transpose()const;
    CMatrix GetInverted()const;
	bool    Invert();
    double  Trace()const;
    int     Rank()const;
    
    //
    // Protective Data Members
    //
    
    protected:
	int     m_nNumColumns;
	int     m_nNumRows;
	double* m_pData;
	
	//
	//
	//
	
	
};
// class CVector
// {
    
// };

#endif // !defined(MATRIX_H)