#include "RungeKutta.h"


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