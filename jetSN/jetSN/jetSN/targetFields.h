#pragma once


//double wph(double z);

double nph_ICani(double Ep, double z, double r, double gamma, const char* id);

double nph_ICani2(double E, double gamma, double eta, const char* id);

double starBlackBody(double Ep, double r, double gamma);

double starIR(double Ep, double z, double gamma);

double wphIR(double z, const char* id);

