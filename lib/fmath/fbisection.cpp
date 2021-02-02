#include "fbisection.h"
#define NTRY 100
#define FACTOR 1.6


int zbrac2(fun1 func,double& x1,double& x2)
    /* Given a function func and an initial guessed range x1 to x2, the routine expands the range
    geometrically until a root is bracketed by the returned values x1 and x2 (in which case zbrac
    returns 1) or until the range becomes unacceptably large (in which case zbrac returns 0). */
{
    void nrerror(char error_text[]);
    int j;
    float f1,f2;
    if (x1 == x2) nrerror("Bad initial range in zbrac");
    f1=func(x1);
    f2=func(x2);
    for (j=1;j<=NTRY;j++) {
        if (f1*f2 < 0.0) return 1;
        if (fabs(f1) < fabs(f2))
            f1=func(x1 += FACTOR*(x1-x2));
        else
            f2=func(x2 += FACTOR*(x2-x1));
    }
    return 0;
}

#define JMAX 40

//===========================================================================
int Bisect(fun1 Func, double a, double b, double &x)
//---------------------------------------------------------------------------
// Determines a real root x of function Func isolated in interval [a,b] by
// the bisection method
// Error code: 0 - normal execution
{
	const double eps = 1e-10;		// relative precision of root
	const int itmax = 100;			// max. no. of iterations
	double fa, fb, fx;
	int it;
	x = a; fa = Func(x);
	if (fabs(fa) == 0e0) return 0;	// is a the root?
	x = b; fb = Func(x);
	if (fabs(fb) == 0e0) return 0;	// is b the root?
	if (fa*fb > 0) return 1;		// [a,b] does not contain a root
									// or contains several roots
	for (it=1; it<=itmax; it++) {
		x = 0.5e0 * (a + b);		// new approximation
		fx = Func(x);
		if (fa*fx > 0) a = x; else b = x;	// choose new bounding interval
		if (((b-a) <= eps*fabs(x)) || (fabs(fx) <= eps)) return 0;
	}
	cout << "Bisect: max. no. of iterations exceeded !" << endl;
	return 2;
}