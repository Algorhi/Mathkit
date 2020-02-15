#if !defined(NORMALEARTH_H)
#define NORMALEARTH_H

class NormalEarth  
{
public:
	NormalEarth();
	NormalEarth(double,double,double);
	virtual ~NormalEarth();
	
public:
   
    double Get_Core_Lat(double B);
    void Get_R0(double A0, double B0, double h0, 
		        double &R0,double &R0x,double &R0y,double &R0z);

    double Get_Core_Distance(double phi);
    double Get_Height(double r, double phi);
    
    double Get_Great_Circle_Distance(const double,const double,const double,const double) const;
    double Get_Azimuth(const double, const double,const double,const double) const;
    
	

public:
	
	double GM()      { return 3.986004418e14;  };
	double J2()      { return 1.08263e-3;   };
	double J4()      { return -2.37091e-6;   };
	double ae()      { return 6378137.;     };
	double be()      { return 6356752.3141; };
	
	
	double Get_Omega() const { return Omega; };
	double Get_Gravity() const { return g0; };
	double Get_Radius() const { return R0; };
	
	void Set_Omega(double _omega) {Omega=_omega;};
	void Set_Gravity(double _g) { g0 =_g; };
	void Set_Radius(double _R) {R0 = _R; }; 

	private:
	
	double Omega = 7.292115e-5;
	double R0 = 6378135;
	double g0 = 9.81;

};

#endif // !defined(NORMALEARTH_H)
