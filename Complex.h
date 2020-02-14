// Complex 
class Complex{
  public:
    Complex();
    Complex(double dblx,double dblY);
    Complex(const Complex& other);
    virtual ~Complex() {};

    void SetReal(double dblX);
    void SetImag(double dblY);
    void GetReal();
    void GetImag();
    String ToString() const;
    void FromString(String s, const String& sDelir = "");

    //
    //
    //

    


}
