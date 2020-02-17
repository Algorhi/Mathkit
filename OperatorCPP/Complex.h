namespace Algorhi;
class Complex{
  public:
    Complex( double re = 0,double im = 0);
    Complex(const Complex&);
    double* get();
    virtual ~Complex();
    
  private:
    double real,imag;
};
