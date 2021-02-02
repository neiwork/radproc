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
double RungeKuttaStep(fun2 f, double x, double y, double dx, const RungeKuttaOpt& opt = DefOpt_RungeKutta);
double integSimpson(double a, double b, fun1 f, size_t n, const RungeKuttaOpt& opt = DefOpt_RungeKutta);
double integSimpsonLog(double a, double b, fun1 f, size_t n, const RungeKuttaOpt& opt = DefOpt_RungeKutta);
double intSimple(double a, double b, fun1 f, const RungeKuttaOpt& opt = DefOpt_RungeKuttaSimple);
double qImprop2(double a, double b, fun1 f, double eps);
double qImpropLog(double a, double b, fun1 f, double eps);
double qMidPoint(double a, double b, fun1 f, double eps);
double qMidPointLog(double a, double b, fun1 f, double eps);
double qGaussLeg(double a, double b, fun1 f, int n);
double qGaussLegLog(double a, double b, fun1 f, int n);