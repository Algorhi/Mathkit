// MathKit.h: interface for the MATHUTILITY
//
//////////////////////////////////////////////////////////////////////
#if !defined(MATHKIT_H)
#define MATHKIT_H

namespace mathkit{
	//Default Precision Control
	const double DEFAULTEPS = 1.e-10;

	//The Right Function Pointer of Differential Equations for Solving Ordinary Differential Equations by Runge-Kutta Method
	typedef void (*ODE)(double t, double *y, double *dy, double *user);  //定义数据类型！！！！

	
	//Generate uniformly distributed random numbers between [0,1]
	double ProduceRandNumber(double *seed);
	
	//Generating a sequence of uniformly distributed random numbers between [0,1]
	void ProduceRandNumbers(double *seed, double r[], int n);
	
	//Generating Normal Distribution Random Numbers with Mean u and Mean Variance g
	double ProduceRandnNumber(double *seed, double u, double g);
	
	//生成均值为u, 均方差为g的正态分布随机序列
	void ProduceRandnNumbers(double *seed, double u, double g, double r[], int n);
	
	//用混洗法产生[0, 1]的均匀分布随机数(性能很好) 
	double Even_rand(int& idum);
	
	//用变换法产生均值为u，均方差为g的正态分布随机数(性能很好)
	double Gauss_rand(int& idum, double u, double g);
	
	//Fixed Step Modified Euler Method for Solving First Order Ordinary Differential Equations
	int Fixstep_Euler(double t, double h, double *y, double *user, int n, ODE Fun);
	
	//Fourth-order Runge-Kutta method with fixed step size for solving first-order ordinary differential equations
	int Fixstep_RungeKutta(double t, double h, double *y, double *user, int n, ODE Fun);
	//Fourth-order Runge-Kutta method with variable step size for solving first-order ordinary differential equations
	int RungeKutta(double t, double h, double *y, double *user, int n, ODE Fun);
	
	//Pi
	double Pi();
	
	//求反三角函数时的限制范围函数
	double limit(double x);
	
	//角度的规范化:
	// fg=0: [0,2*Pi]
	// else: [-Pi,Pi]
	void Norm_Angle_RAD(double& Alpha,int fg);
	
	//角度的规范化:
	// fg=0: [0,360]
	// else: [-180,180]
	void Norm_Angle_DEG(double& Alpha,int fg);
	
	//弧度转化为角度
	double Rad2Deg(double Alpha);
	
	//角度转化为弧度
	double Deg2Rad(double Alpha);
	
	//符号函数
	int sign(double x);

	void Cross_Multip(double X[], double Y[], double Z[]);    
	
	//Interpolation function
	double spline(double x[],double y[], int n, const double t);
	//x[n],y[m],z[n][m]
	double spline(double *x,double *y, double *z,int n, int m, const double u, const double v);	
    
    
}
namespace mathmat{
    
}



#endif // !defined(MATHKIT_H)
