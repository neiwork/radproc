#pragma once


double bessj0(double x);
// Returns the Bessel function J_0(x) for any real x.

double bessy0(double x);
// Returns the Bessel function Y_0(x) for positive x.

double bessj1(double x);
// Returns the Bessel function J_1(x) for any real x.

double bessy1(double x);
// Returns the Bessel function Y_1(x) for positive x.

double bessy(int n, double x);
// Returns the Bessel function Y_n(x) for positive x and n≥2.

double bessj(int n, double x);
// Returns the Bessel function J_n(x) for any real x and n≥2.

double bessi0(double x);
// Returns the modified Bessel function I_0(x) for any real x.

double bessk0(double x);
// Returns the modified Bessel function K_0(x) for positive real x.

double bessi1(double x);
// Returns the modified Bessel function I_1(x) for any real x.

double bessk1(double x);
// Returns the modified Bessel function K_1(x) for positive real x.	

double bessk(int n, double x);
// Returns the modified Bessel function K_n(x) for positive x and n≥2.

double bessi(int n, double x);
// Returns the modified Bessel function I_n(x) for any real x and n≥2.