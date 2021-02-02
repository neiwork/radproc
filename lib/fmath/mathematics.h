#pragma once

#include <cmath>
//#include <math>
#include <vector>
#include <functional>

// c++ math functions plus added math functions

//namespace physics {

	inline double P2(double x) { return x*x; }
	inline double P3(double x) { return x*x*x; }

	typedef std::vector< std::vector< double > > Matrix;
	typedef std::vector<double> Vector;
	typedef std::vector<float> fVector;

	typedef std::function<double(double)> fun1;
	typedef std::function<double(double, double)> fun2;
	typedef std::function<double(double, double, double)> fun3;

class Fun {
public:
	virtual double call(double v)=0;
};

struct zero_params {};
struct one_d_params {double p;};
struct two_d_params {double p1,p2;};
struct three_d_params {double p1,p2,p3;};
struct four_d_params {double p1,p2,p3,p4;};

//}

//using namespace physics;




//}

//using namespace physics;


