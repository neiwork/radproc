#include "lossesAnisotropicIC.h"

//#include "crossSectionInel.h"
//#include "dataLosses.h"
//#include <fmath\RungeKutta.h>
#include "modelParameters.h"
#include <fmath\physics.h>



void incremento(int entero, int left, int right, double &A1, double &A2, double &A3, double R);

double intTriple(double E, double eps_min, double eps_max, double r);



double b_theta(double theta, double w0, double E)
{
	//corregido es 1.0 - cos(theta') escrito en funcion de theta sin primar
	double Dstar = 1.0 / (Gamma - 1.0);
	double corregido = (1.0 - cos(theta));// *Dlorentz*Dstar;
	double b = 2.0 * (corregido)*w0*E / P2(electronMass*cLight2);
	return b;
}


double cAniIC(double w0)   //limite inferior
{
	return w0;
}

double dAniIC(double w0, double E, double theta)   //limite superior     
{

	double s = b_theta(theta, w0, E);  //reemplazo el Gamma = 4*u*E / P2(mass*cLight2), para incuir la dep theta
	return s*E/(1+s);
}



double difN(double theta, double w, double w0, double E, double r)   //funcion a integrar
{
	// eps -> variable de afuera
	// eps_0 -> variable integracion
	
	double b = b_theta(theta, w0, E);
	double z = w / E;

	//defino F(z)
	double F = 1.0 + P2(z) / (2.0*(1.0 - z)) - 2.0*z / (b*(1.0 - z)) + 2.0*P2(z) / P2(b*(1.0 - z));

	double nph_0 = blackBody(w, r); 
	double invariant = nph_0 / w; 

	double result = (3.0*thomson / (16.0*pi)) * P2(electronMass*cLight2/E) * invariant * F;

	return result;

//	double Nterm = tpf(u); //state.targetPhotonField(Ep);
//	std::cout << logE << "\t" << Nterm << std::endl;  //<< log10(E*6.25e11) << "\t"  

}


double lossesAnisotropicIC(double E, Particle& particle, double r)
{
	double mass = particle.mass;


	double starT = 1.0e5; //VER
	double a  = 1.0e-3*boltzmann*starT;      //energia minima de los fotones en erg
	double b  = 1.0e3*boltzmann*starT;    //energia maxima de los fotones en erg

	double integral = intTriple(E, a, b, r);// double eps_min, double eps_max, double r)
	
	//double integral = RungeKutta(a,b,&cIC,
	//	[mass, E](double u){return dAniIC(u, E, mass);},
	//	[E, mass, tpf](double u, double t){ return fIC(u, t, E, mass, tpf); });    //le asigno a la variable integral el resultado de la integracion   

			
	return integral*cLight;

	}


//para losses(Ee)  -> intTriple(E, photon_Emin, photon_Emax, cAnimin, dAnimax)
double intTriple(double E, double eps_min, double eps_max, double r)  
{
	int n = 30;

	double x_int = pow((eps_max / eps_min), (1.0 / n));

	double X1(0.0), X2(0.0), X3(0.0);

	if (eps_min < eps_max)
	{
		double x = eps_min;

		for (int i_x = 0; i_x < n; ++i_x)     //le saco el n para que se multiplique n veces y no n+1
		{
			double dx = x*(x_int - 1);

			double t_min =  pi / 2;//   //t_min=0.0
			double t_max =  pi;

			//	if (t_min < t_max)
			//	{

			double t_int = pow((t_max / t_min), (1.0 / n));

			double T1(0.0), T2(0.0), T3(0.0);

			double t = t_min;

			for (int i_t = 0; i_t < n; ++i_t)
			{
				double dt = t*(t_int - 1);

				double inf = cAniIC(x); //limite inferior de eps_1
				double sup = dAniIC(x, E, t); //limite superior (double u, double E, double theta) 

				if (inf < sup)
				{

					double y_int = pow((sup / inf), (1.0 / n));

					double Y1(0.0), Y2(0.0), Y3(0.0);

					double y = inf;

					for (int i_y = 0; i_y < n; ++i_y)
					{
						double dy = y*(y_int - 1);

						////////////// calculo la funcion a integrar

						// eps -> variable de afuera = x
						// eps_0 -> variable integracion =y
						//difN(double theta, double eps, double eps_0, double E, double r);
						double Q = (y - x)*difN(t, x, y, E, r);//
						//double Q = 1.0;
							 
						////////////////////////////

						double L1 = Q*dy;

						//if (L1 > 0.0) { 
						incremento(i_y, 0, n - 1, Y1, Y2, Y3, L1);// }

						y = y*y_int;
					}

					double L2 = (Y1 + 2.0 * Y2 + 4.0 * Y3) *(0.5*sin(t))*dt / 3.0; //agregue int sobre ang solido
					//double L2 = (Y1 + 2.0 * Y2 + 4.0 * Y3) *(E*t)*dt / 3.0; //prube de usar b variable int

//					//if (L2 > 0.0) { 
						incremento(i_t, 0, n - 1, T1, T2, T3, L2); //}			
				}
		

				t = t*t_int;

			}

//			//}
//
			double L3 = (T1 + 2.0 * T2 + 4.0 * T3)*dx / 3.0; 

			//if (L3 > 0.0) { 
				incremento(i_x, 0, n - 1, X1, X2, X3, L3); //}

			x = x*x_int;
		}
		
		double L4 = (X1 + 2.0 * X2 + 4.0 * X3) / 3.0;
		return L4;

	}//cierra el primer if
	else
	{
		return 0.0;
	}

}


void incremento(int entero, int left, int right, double &A1, double &A2, double &A3, double R)
{
	double floaT = entero;
	double aux = fmod(floaT, 2.0);

	if (entero == left || entero == right) //si el entero coincide con los extremos
	{
		A1 += R;
	}
	else if (aux == 0)				//si el entero es par
	{
		A2 += R;
	}
	else
	{
		A3 += R;
	}    //si el entero es impar

}