#pragma once
#include "physics.h"

/** Defino una funcion que resuelve una integral doble usando el método
    de Runge-Kutta de orden cuatro
*/
double RungeKutta(double a, double b, fun1 c, fun1 d, fun2 f);

/** Esta función calcula una integral simple*/ 
double RungeKuttaSimple(double a, double b, fun1 f);
