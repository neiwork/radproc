#include "distribution.h"

#include "injection.h"
#include "losses.h"
#include "timeDistribution.h"
#include <fmath\interpolation.h>
#include <fmath\elimiGaussiana.h>
#include <fmath\matrixInit.h>
#include <fmath\physics.h>


void distribution(Particle& p, State& st)
{

	for (int z_i = 0; z_i < p.ps[1].values.size(); z_i++) {
		
		double z = p.ps[1][z_i]; //valor d ela variable
		int z_position = z_i; // posicion en la dimension z 

	//p.ps.iterate([&p, &st](const SpaceIterator& j){  //reemplace este iterate con el for de arriba
		//double z = j.par.R; //valor de la variable que recorro en este iterate
		//int z_position = j.coord[1];

		//bloque 1

		double Eeff = 0.0;

	//	timeDistribution(p, st, z_position, Eeff);  //descomentar

		//bloque 2 

		ParamSpaceValues Nprevt(p.ps), Naux(p.ps);
		//genero dos psv auxiliares para el siguiente iterate
		p.ps.iterate([&p, &Nprevt, &Naux](const SpaceIterator& i){
						
//			if (i.its[2].canPeek(-1))     //descomentar el bloque
//			{
//				SpaceCoord prevt = i.moved({ 0, 0, -1 });
//				Nprevt.set(i, p.distribution.get(prevt));  //N evaluado en t-1
//			}
			Naux.set(i, p.distribution.get(i)); //copia del N  //aca se rompe
		
		}, { -1, z_position, -1 }); //fijo z //copio el iterador que sigue asi los aux se completan bien


		p.ps.iterate([&p, &st, &z, &Eeff, &Nprevt, &Naux](const SpaceIterator& i){

			double E = i.par.E;
			double ratioLosses = losses(Eeff, p, st) / losses(E, p, st);

			double dist;
			if (i.its[2].canPeek(-1))
			{
				//dist = N1(Eeff, t[i - 1]);  
				dist = Nprevt.dimInterpolator(0)(Eeff); //interpola en la dimension 0, y evalua en Eeff
			}
			else{
				dist = 0.0; //ver que pasa aca, sería en t[0]
			}

			double  Ntot = Naux.get(i) + dist*ratioLosses;
			/*Naux.get(i)      = N1
			  dist*ratioLosses = N2 == N2(E[l], t[i]) = N1(Eeff, t[i - 1])*losses(Eeff, p, st) / losses(E[l], p, st);
			  distribution(E[l], z[k], t[i]) = N1(E[l], z[k], t[i]) + N2(E[l], z[k], t[i]);*/

			p.distribution.set(i, Ntot);

		}, { -1, z_position, -1 }); //fijo z

	}
//}, { 0, -1, 0 }); //afuera itero solo sobre z reemplace este iterate con un for sobre r(i)


	//bloque 3

	p.ps.iterate([&p, &st](const SpaceIterator& j){

		int E_position = j.coord[0];
		int t_position = j.coord[2];

		if (j.its[2].canPeek(1)){

			double delta_t = (j.its[2].peek(1) - j.par.T) / Gamma;//    delta_t = (t[i] - t[i - 1]) / Gamma; en realidad t(i+1)-t(i)
			//std::cout << " (last index: " << i.coord.dims[0] - 1 << ", value: " << i.its[0].peek(-1) << ")";

			p.ps.iterate([&p, &st, &delta_t](const SpaceIterator& i){

				double dist;
				// example: gets value at the previous position of the iterator of dimension 0
				// (also checks if we can -- i.e. current index > 0)
				if (i.its[1].canPeek(1))
				{

					double delta_xk = i.its[1].peek(1) - i.par.R; //delta_xk = z[k] - z[k - 1];

					if (i.its[1].canPeek(2)){

						double delta_xk_1 = i.its[1].peek(2) - i.its[1].peek(1);// z[k - 1] - z[k - 2];

						dist = 0.0; //Ver como hago lo de abajo
						
						SpaceCoord prevZ = i.moved({ 0, -1 }); //ver el -1, cambarlo de la mismaa forma que los deltas !!!!!!!!!!!!!!!!!!!!!
						p.distribution.get(prevZ);
						
						//distribution(E[l], z[k], t[i]) = distribution(E[l], z[k - 1], t[i])*(delta_t*cLight / delta_xk - 1)
						//							- distribution(E[l], z[k], t[i])*(delta_t*cLight / delta_xk);;
					}
					else{
						dist = 0.0; //Ver como hago lo de abajo
						//distribution(E[l], z[k], t[i]) = distribution(E[l], z[k - 1], t[i])*(delta_t*cLight / delta_xk - 1)
						//				aca esta linea no			- distribution(E[l], z[k], t[i])*(delta_t*cLight / delta_xk);;
					}
				}
				else {
					dist = 0.0; //este sería para el zmax, ver que pasa	
				}

				p.distribution.set(i, dist);

			}, { E_position, -1, t_position });
		}
		else {

			p.ps.iterate([&p, &st](const SpaceIterator& i){
				//aca estoy en el caso de tmax, ver que pasa
				double dist = 0.0;

				p.distribution.set(i, dist);

			}, { E_position, -1, t_position });

		}

	}, { -1, 0, -1 }); //VER necesito fijar alguna coordenada??? , { 0, 0, 0 }); //afuera itero solo sobre t
}

/*	p.ps.iterate([&p, &st](const SpaceIterator& i){
		
		i.par.magneticField;
		double total = 0.0, E = i.par.E, z = i.par.R, t = i.par.T;

	// agarra injection de injection (por ahora);
	Vector& Q = p.injection.values;
	//Vector Q = particle.injection;

	Vector& distribution = p.distribution.values;

	// { asume una unica dimension 0 que es la energia }
	Vector& E = p.ps.dimensions[0]->values;
	//Vector E = particle.energyPoints;

	Vector& z = p.ps.dimensions[1]->values;  //particle.zPoints;
	Vector& t = p.ps.dimensions[2]->values;  //particle.zPoints;

	int ePoints = E.size() - 1;    //el -1 es para que sea el nEnergyPoints
	int zPoints = z.size() - 1;
	int tPoints = t.size() - 1;*/


//hacer una funcion para encontrar el incide en el vector, dadas las posiciones E[l], z[k], t[i]

//void distribution(Particle& p, State& state)
//{
//
//	// agarra injection de injection (por ahora);
//	Vector& injection = p.injection.values;
//	Vector& distribution = p.distribution.values;
//
//	// { asume una unica dimension 0 que es la energia }
//	Vector& Ee = p.ps.dimensions[0]->values;
//
//	int ne = Ee.size() - 1;    //el -1 es para que sea el nEnergyPoints y 
//
//	Matrix a;
//
//	matrixInit(a, Ee.size(), Ee.size() + 1, 0);
//
//	for (int i = 0; i <= ne; ++i)	{     //ver como acomodo los limites
//
//		if (i == ne)	{
//			a[i][i] = 1;    //esta es la condicion inicial: N(Emax)=0
//			a[i][ne + 1] = 0;
//		}
//		else	{
//			double h = Ee[i + 1] - Ee[i];
//
//			for (int j = 0; j <= ne + 1; ++j){
//				if (j == i)	{                               //escapeTime esta en [s]
//					a[i][j] = losses(Ee[i], p, state) + 0.5*h / escapeTime(Ee[i], p);
//				}
//				else if (j == i + 1)	{
//					a[i][j] = -losses(Ee[i + 1], p, state) + 0.5*h / escapeTime(Ee[i], p);
//				}
//				else if (j == ne + 1)	{
//
//					a[i][j] = 0.5*h*(injection[i] + injection[i + 1]);
//				}
//				else	{
//					a[i][j] = 0;
//				}
//			}
//		}
//	}
//
//	elimiGaussiana(ne + 1, a, distribution);
//
//	if (p.type == PT_proton){
//		distribution[0] = distribution[1] = distribution[2]; //esto lo hago porque el primer punto me da enorme
//	}                                                      //y arrastro el error para piones y muones
//}
//
//
//


