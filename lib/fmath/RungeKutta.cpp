#include "RungeKutta.h"
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_integration.h>

namespace {
	void inline incremento(int entero, int left, int right, double &A1, double &A2, double &A3, double R)
	{
			double floaT = entero;
			double aux = fmod(floaT,2.0);
			
			if(entero==left || entero==right) //si el entero coincide con los extremos
				{A1 += R ;}
			else if(aux==0)				//si el entero es par
				{A2 += R ;}
			else 
				{A3 += R ;}    //si el entero es impar

	}
}

RungeKuttaOpt DefOpt_RungeKutta{ 50, 50 };
double RungeKutta(double a, double b, fun1 c, fun1 d, fun2 f, const RungeKuttaOpt& opt)    //paso como argumento los cuatro limites
//de las integrales y la funcion a integrar																														
{

	int n_x = opt.samples_x;  //recordar que el numero de puntos en los que interpolo 
	//no puede ser mayor al numero de puntos de la funcion

	int n_y = opt.samples_y;

	double x_int = pow((b / a), (1.0 / n_x));

	double J1(0.0), J2(0.0), J3(0.0);

	if (a < b)
	{
		double x = a;
        //#pragma omp parallel for reduction(+:J1,+:J2,+:J3)
		for (int i = 0; i < n_x; ++i)     //le saco el n para que se multiplique n veces y no n+1
		{
			double dx = x*(x_int - 1);

			double sup = d(x);
			double inf = c(x);

			if (inf < sup)
				//		{ if (E > sup)
			{

				double y_int = pow((sup / inf), (1.0 / n_y));

				double K1(0.0), K2(0.0), K3(0.0);

				double y = inf;

				for (int j = 0; j < n_y; ++j)
				{
					double dy = y*(y_int - 1);

					double L1 = f(x, y)*dy;

					if (L1 > 0.0) { incremento(j, 0, n_y - 1, K1, K2, K3, L1); }

					y = y*y_int;
				}

				double L2 = (K1 + 2 * K2 + 4 * K3)*dx / 3;

				if (L2 > 0.0) { incremento(i, 0, n_x - 1, J1, J2, J3, L2); }
				//		}
			}

			x = x*x_int;

		}

		return (J1 + 2.0 * J2 + 4.0 * J3) / 3.0;
	}
	else
	{
		return 0.0;
	}


}

RungeKuttaOpt DefOpt_RungeKuttaSimple{ 50, -1 };

double RungeKuttaStep(fun2 f, double x, double y, double dx, const RungeKuttaOpt& opt)
{
	double k1 = dx*f(x,y);
	double k2 = dx*f(x+0.5*dx,y+0.5*k1);
	double k3 = dx*f(x+0.5*dx,y+0.5*k2);
	double k4 = dx*f(x+dx,y+k3);
	return (1.0/6.0)*(k1+2*k2+2*k3+k4);
}

double integSimpson(double a, double b, fun1 f, size_t n, const RungeKuttaOpt& opt)
{
	a += 0.0001*abs(a);
	b -= 0.0001*abs(b);
	double dx = (b-a)/n;
	double sum_even = 0.0;
	for (size_t j=1; j <= n/2-1; j++)
		sum_even += f(a+2*j*dx);
	double sum_odd = 0.0;
	for (size_t j=1; j <= n/2; j++)
		sum_odd += f(a+(2*j-1)*dx);
	return dx/3.0 * (f(a)+f(b)+2.0*sum_even+4.0*sum_odd);
}

double integSimpsonLog(double a, double b, fun1 f, size_t n, const RungeKuttaOpt& opt)
{
	if (a > 0.0 && b > a) {
		return integSimpson(log(a),log(b),[&](double logx)
				{
					double x = exp(logx);
					return x*f(x);
				},n);
	} else
		return 0.0;
}

double RungeKuttaSimple(double a, double b, fun1 f, const RungeKuttaOpt& opt)
{
	int RK_N = opt.samples_x;
	int n = RK_N;

	double x_int = pow((b / a), (1.0 / n));

	double J1(0.0), J2(0.0), J3(0.0);

	if (a < b)
	{
		double x = a;

		for (int i = 0; i < n; ++i)
		{
			double dx = x*(x_int - 1.0);

			double L1 = f(x)*dx;

			if (L1 > 0.0) { incremento(i, 0, n - 1, J1, J2, J3, L1); }

			x = x*x_int;

		}
	}

	return (J1 + 2.0 * J2 + 4.0 * J3) / 3.0;
}


double intSimple(double a, double b, fun1 f, const RungeKuttaOpt& opt)
{
	int RK_N = opt.samples_x;
	int n = RK_N;

	double x_int = pow((b / a), (1.0 / n));

	double L1 = 0.0;

	if (a < b)
	{
		double x = a;

		for (int i = 0; i < n; ++i)
		{
			double dx = x*(x_int - 1.0);

			L1 += f(x)*dx;

			x = x*x_int;

		}
	}

	return L1;
}



//===========================================================================
double qRomberg(double a, double b, fun1 Func, double eps)
//---------------------------------------------------------------------------
// Integrates function Func on interval [a,b] with relative precision eps
// using the adaptive Romberg method
//---------------------------------------------------------------------------
{
	const int kmax = 30;			// max. no. of step halving iterations
	double r1[kmax+1];				// two consecutive lines
	double r2[kmax+1];				// from the method table
	double f, h, sum;
	long i, n;
	int j, k;
	h = b-a; n = 1;
	r1[0] = 0.5*h*(Func(a) + Func(b));					// initial approximation
	for (k=1; k<=kmax; k++) {							// step halving loop
		sum = 0e0;
		for (i=1; i<=n; i++) sum += Func(a+(i-0.5)*h);
		r2[0] = 0.5*(r1[0] + h*sum);					// trapezoid formula
		f = 1e0;
		for (j=1; j<=k; j++) {							// increase quadrature order
			f *= 4;
			r2[j] = (f*r2[j-1] - r1[j-1])/(f-1);		// new approximation
		}
		if (k > 1) {									// convergence check
			if (fabs(r2[k]-r1[k-1]) <= eps*fabs(r2[k])) break;
			if (fabs(r2[k]) <= eps && fabs(r2[k]) <= fabs(r2[k]-r1[k-1])) break;
		}
		h *= 0.5; n *= 2;								// halve integration step
		for (j=0; j<=k; j++) r1[j] = r2[j];				// shift table lines
	}
	if (k > kmax) {
		printf("qRomberg: max. no. of iterations exceeded !\n");
		k--;
	}
	return r2[k];
}


//===========================================================================
double qImprop2(double a, double b, fun1 Func, double eps)
//---------------------------------------------------------------------------
// Integrates function Func on interval [a,b] with a and/or b singular
// integrable points.
// Calls: qRomberg
//---------------------------------------------------------------------------
{
	double h, h0, s, s1, x;
	h0 = 0.1e0 * (b-a);						// extent of vicinities of singular boundaries
	s = 0e0;
	if (fabs(Func(a)) > 9e99) {				// a is singular point?
		h = h0;
		s1 = 1e0;
		while(fabs(s1) > eps*fabs(s)) {
			h *= 0.5;							// halve interval
			x = a + h;							// left boundary of [x,x+h]
			s1 = qRomberg(x,x+h,Func,eps);		// partial integral on [x,x+h]
			s += s1;							// add contribution to total integral
		}
		a += h0;								// new left boundary of core interval
	}
	if (fabs(Func(b)) > 9e99) {					// b is singular point?
		h = h0;
		s1 = 1e0;
		while(fabs(s1) > eps*fabs(s)) {
			h *= 0.5;							// halve interval
			x = b - h;							// right boundary of [x-h,x]
			s1 = qRomberg(x-h,x,Func,eps);		// partial integral on [x-h,x]
			s += s1;							// add contribution to total integral
		}
		b -= h0;								// new right boundary of core interval
	}
	s += qRomberg(a,b,Func,eps);				// add integral on core interval
	return s;
}

double qImpropLog(double a, double b, fun1 Func, double eps)
{
	if (a > 0.0 && b > a) {
		return qImprop2(log(a),log(b),[&Func](double logx)
				{
					double x = exp(logx);
					return x*Func(x);
				},eps);
	} else
		return 0.0;
}

//===========================================================================
double qMidPoint(double a, double b, fun1 Func, double eps)
//---------------------------------------------------------------------------
// Integrates function Func on interval (a,b) with relative precision eps
// using the adaptive midpoint rule
//---------------------------------------------------------------------------
{
	const int kmax = 19;
	double h, s, s0, sum;
	double f1p6 = 1e0/6e0, f5p6 = 5e0/6e0;
	long i, n;
	int k;
	h = b-a; n = 1;
	s0 = h * Func(a+0.5*h);					// max. no. of divisions
											// initial approximation
	for (k=1; k<=kmax; k++) {				// step division loop
		sum = 0e0;
		for (i=1; i<=n; i++) sum += Func(a+(i-f5p6)*h) + Func(a+(i-f1p6)*h);
		s = (s0 + h*sum)/3;								// new approximation
		if (fabs(s - s0) <= eps*fabs(s)) break;			// convergence check
		h /= 3; n *= 3;									// reduce step
		s0 = s;
	}
	if (k > kmax) printf("qMidPoint: max. no. of iterations exceeded !\n");
	return s;
}

double qMidPointLog(double a, double b, fun1 Func, double eps)
{
	if (a > 0.0 && b > a) {
		return qMidPoint(log(a),log(b),[&Func](double logx)
				{
					double x = exp(logx);
					return x*Func(x);
				},eps);
	} else
		return 0.0;
}


//===========================================================================
double Legendre(int n, double x, double &d)
//---------------------------------------------------------------------------
// Evaluates the n-th order Legendre polynomial and its derivative d in x
// using the recurrence relation
//---------------------------------------------------------------------------
{
	double f, fm1, fm2;
	int i;
	if (n == 0) {
		f = 1e0; d = 0e0;
	} else {
		f = x; fm1 = 1e0;
		for (i=2; i<=n; i++) {
			fm2 = fm1; fm1 = f;
			f = ((2*i-1)*x*fm1 - (i-1)*fm2)/i;
		}
		d = (x*x-1e0) ? n*(x*f-fm1)/(x*x-1e0) : 0.5*n*(n+1)*f/x;
	}
	return f;
}

//===========================================================================
void xGaussLeg(double a, double b, double x[], double w[], int n)
//---------------------------------------------------------------------------
// Calculates abscissas x[] and weights w[] for the n-point Gauss-Legendre
// quadrature on interval [a,b]
// Calls: Legendre (from specfunc.h)
//---------------------------------------------------------------------------
{
	const double eps = 1e-14;							// relative precision of zeros
	double d, f, xc, xi;
	int i, n2;
	n2 = n/2;
	for (i=1; i<=n2; i++) {
		xi = cos(pi*(i-0.25e0)/(n+0.5e0));				// initial approximation for zeros
		f = 9e99;
		while (fabs(f) > eps*fabs(xi)) {				// Newton-Raphson refinement
			f = Legendre(n,xi,d) / d;
			xi -= f;
		}
		x[i] = -xi; x[n-i+1] = xi;						// symmetrical zeros
		w[i] = w[n-i+1] = 2e0/((1e0-xi*xi)*d*d);		// equal weights
	}
	if (n % 2 == 1) {
		Legendre(n,0e0,d);
		x[n2+1] = 0e0;
		w[n2+1] = 2e0/(d*d);
	}
	f = 0.5e0*(b-a); xc = 0.5e0*(b+a);
	for (i=1; i<=n; i++) {
		x[i] = f*x[i] + xc;
		w[i] = f*w[i];
	}													// odd no. of mesh points
														// scaling to interval [a,b]Integration of Functions
}

//===========================================================================
double *dVector(int imin, int imax)
//---------------------------------------------------------------------------
// Allocates a double vector with indices in the range [imin,imax]
//---------------------------------------------------------------------------
{
	double *p;				// assign block start to array pointer
	p = (double*) malloc((size_t) ((imax-imin+1)*sizeof(double)));
	if (!p) {
		printf("Vector: allocation error !\n");
		exit(1);
	}
	return p - imin;		// adjust for offset
}

//===========================================================================
void FreeVector(double *p, int imin)
//---------------------------------------------------------------------------
// Deallocates a double vector allocated with Vector, with offset imin
//---------------------------------------------------------------------------
{
	free((void*) (p+imin));		// compensate for offset
}


//===========================================================================
double qGaussLeg(double a, double b, fun1 Func, int n)
//---------------------------------------------------------------------------
// Integrates function Func on interval [a,b] using n-point Gauss-Legendre
// quadratures
// Calls: xGaussLeg
//---------------------------------------------------------------------------
{
	double s, *x, *w;
	int i;
	x = dVector(1,n);
	w = dVector(1,n);
	xGaussLeg(a,b,x,w,n);
	s = 0e0;
	for (i=1; i<=n; i++) s += w[i] * Func(x[i]);
	FreeVector(x,1);
	FreeVector(w,1);
	return s;
}

double qGaussLegLog(double a, double b, fun1 Func, int n)
{
	if (a > 0.0 && b > a) {
		return qGaussLeg(log(a),log(b),[&Func](double logx)
				{
					double x = exp(logx);
					return x*Func(x);
				},n);
	} else
		return 0.0;
}