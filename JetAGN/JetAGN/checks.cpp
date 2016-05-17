#include "checks.h"

#include "injection.h"
#include "distribution.h"

#include <fparameters/Dimension.h>
#include <fparameters/SpaceIterator.h>

#include <fmath/interpolation.h>

#include <iostream>
#include <iomanip>
void check_vec(const Vector& a, const Vector& b)
{
	//size_t N = a.size();
	//bool error = false;
	//for (int i = 0; i < N; ++i) {
	//	if (std::abs((a[i] - b[i]) / b[i]) >= 1.0e-15) {
	//		std::cout << (a[i] - b[i]) / std::max(std::abs(a[i]), std::abs(b[i])) << " diff at " << i << std::endl;
	//		error = true;
	//	}
	//}
	//if (error) {
	//	throw;
	//}
	//std::cout << "OK." << std::endl;
}

void check_inj(const Particle& p)
{
	//check_vec(p.injection, p.injection.values);
}

void check_refactor(const Particle& p)
{
	check_inj(p);
	check_dist(p);
}

void check_dist(const Particle& p)
{
	//check_vec(p.distribution, p.distribution.values);
}

void zeroToN(Vector& v) {
	for (int i = 0; i < (int)v.size(); ++i) {
		v[i] = i;
	}
}


void testFillAndDump() {

	ParamSpace ps;
	ps.add(new Dimension(2, zeroToN));
	ps.add(new Dimension(3, zeroToN));
	ps.add(new Dimension(2, zeroToN));

	// create the PSVs
	ParamSpaceValues m(ps); // some PSV

	// iterate the dimensions, and for each combination of values; 
	// compute a value and store it appropriate location in PSV
	int c = 0;
	m.fill([&c](const SpaceIterator& i){
		return c++;
	});

	ps.iterate([&m, &c](const SpaceIterator& i){
		for (auto coord : i.coord.dims) {
			std::cout << coord;
		}
		std::cout << ":" << m.get(i) << std::endl;
	});

}

void test3DPrint() {

	ParamSpace ps;
	ps.add(new Dimension(4, zeroToN));
	ps.add(new Dimension(4, zeroToN));
	ps.add(new Dimension(4, zeroToN));

	// create the PSV
	int c = 0;
	ParamSpaceValues m(ps, [&c](const SpaceIterator& i){
		return c++;
	});

	ps.iterate([&m](const SpaceIterator& i){
		int plane = i.coord.dims[2];
		std::cout << plane << std::endl;
		std::cout << "  ";
		i.ps.iterate([&m](const SpaceIterator& j){
			if (j.coord.dims[0] == 0) {
				std::cout << std::endl << "  ";
			}
			std::cout << std::setw(3) << m.get(j);
		}, { -1, -1, plane });
		std::cout << std::endl;
	}, { 0, 0 });
}

void testIteratePlane() {

	ParamSpace ps;
	ps.add(new Dimension(4, zeroToN));
	ps.add(new Dimension(4, zeroToN));

	// create the PSV
	int c = 0;
	ParamSpaceValues m(ps);

	/*
		El tercer parametro opcional de iterate() es una lista de 
		dimensiones que seran fijadas al iterar. Se iteraran las
		dimensiones que no fueron fijadas. Se especifican los indices
		en los cuales se fijara cada dimension en orden. Por ejemplo,
		fijar la primera dimension en el primero de sus valores se
		escribe:

		  { 0 }

		Fijar la segunda dimension en el segundo de sus valores, ademas
		de la primera ya fijada:

		  { 0, 1 }

		Si se necesita fijar la segunda dimension pero no la primera,
		se saltea con -1:

		  { -1, 1 }

		Se puede fija cualquier cantidad de dimensiones, incluso todas,
		siendo la iteracion de una unica celda en el hipercubo PSV.	
	*/
	
	ps.iterate([&m, &c](const SpaceIterator& i){
		m.set(i, 1);
	}, {-1,1});

	ps.iterate([&m, &c](const SpaceIterator& i){
		m.set(i, 2);
	}, { 2 });

	ps.iterate([&m, &c](const SpaceIterator& i){
		m.set(i, 3);
	}, { 3,3 });

	m.debug_dump();

}

void testCanPeek() {

	ParamSpace ps;
	ps.add(new Dimension(4, zeroToN));
	ps.add(new Dimension(4, zeroToN));

	// create the PSV
	int c = 0;
	ParamSpaceValues m(ps,[&c](const SpaceIterator& i){
		return c++;
	});

	m.debug_dump();

	std::cout << std::endl;

	ps.iterate([&m, &c](const SpaceIterator& i){
		SpaceCoord coord = i.moved({ -1, -1 });

		if (i.canPeek(coord)) {
			m.set(i, m.get(coord));
		}
	});

	m.debug_dump();

}

void testSpaceIterator2() {
	//setParameters();
	//factor_qrel   = 3.0; 

	double electronLogEmin = 6.0;
	double electronLogEmax = 15.0;
	Particle electron(PT_electron, electronMass, electronLogEmin, electronLogEmax, 10);

	electron.ps.add(new Dimension(5, zeroToN));
	electron.ps.add(new Dimension(4, zeroToN));

	// create the PSVs
	//electron.initialize();

	ParamSpaceValues m(electron.ps); // some PSV

	// iterate the dimensions, and for each combination of values; 
	// compute a value and store it appropriate location in PSV
	electron.ps.iterate([&m](const SpaceIterator& i){
		m.set(i, (i.val(DIM_E) + parameters.magneticField) * 2);
	});


	// DEBUG >>>>>>>>>>>>>
	electron.ps.iterate([](const SpaceIterator& i){
		// read values from the parameters object (shows how they iterate);
		std::cout << i.val(DIM_E) << ", " << parameters.magneticField;

		// example: gets value at the previous position of the iterator of dimension 0
		// (also checks if we can -- i.e. current index > 0)
		if (i.its[0].canPeek(-1)) {
			std::cout << " (last index: " << i.coord.dims[0] - 1 << ", value: " << i.its[0].peek(-1) << ")";
		}

		std::cout << std::endl;
	});

	m.debug_dump();
	//electron.distribution = m;
	//std::cout << electron.dist(1.5) << std::endl;
	// <<<<<<<<<<<<<<<<<<<

}

void testBinarySearch(){

	int size = 10;
	Vector Q(size, 0.0);
	Vector E(size, 0.0);

	for (size_t i = 0; i < Q.size(); ++i){
		Q[i]=P2(i);
		E[i] = i;
	}

	double Eint = 3.5;
	double Qint1 = interpol(Eint, E, Q, Q.size());
	double Qint2 = interpolNew(Eint, E, Q, Q.size());

	std::cout << Qint1 << 't' << Qint2 << std::endl;


}

void testMultiDimensionalInterpolate(){

	ParamSpace ps;
	ps.add(new Dimension(2, zeroToN));
	ps.add(new Dimension(3, zeroToN));
	ps.add(new Dimension(2, zeroToN));

	// create the PSVs
	ParamSpaceValues m(ps); // some PSV

	// iterate the dimensions, and for each combination of values; 
	// compute a value and store it appropriate location in PSV
	int c = 0;
	m.fill([&c](const SpaceIterator& i){
		return c++;
	});

	ps.iterate([&m, &c](const SpaceIterator& i){
		for (auto coord : i.coord.dims) {
			std::cout << coord;
		}
		std::cout << ":" << m.get(i) << std::endl;
	});

	std::cout << std::endl;
	

	SpaceCoord p({ 0, 0, 0 });
	
	m.interpolate({ { 0, 0.1 }, { 1, 0.5 }},&p);
	m.interpolate({ { 0, 0.1 }, { 1, 0.5 }, { 2, 0 } });
}