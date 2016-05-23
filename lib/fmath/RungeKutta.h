#pragma once
#include "physics.h"

class RungeKuttaOpt {
public:
	int samples_x;
	int samples_y;
};

/** Defino una funcion que resuelve una integral doble usando el método
    de Runge-Kutta de orden cuatro
*/
extern RungeKuttaOpt DefOpt_RungeKutta;
double RungeKutta(double a, double b, fun1 c, fun1 d, fun2 f, const RungeKuttaOpt& opt = DefOpt_RungeKutta);

/** Esta función calcula una integral simple*/ 
extern RungeKuttaOpt DefOpt_RungeKuttaSimple;
double RungeKuttaSimple(double a, double b, fun1 f, const RungeKuttaOpt& opt = DefOpt_RungeKuttaSimple);
