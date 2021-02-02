#pragma once

#include <cstddef>
#include "mathematics.h"
#include <iomanip>

using namespace std;

void matrixInit(Matrix& m, size_t width, size_t height, double initValue);

/*  ver las siguientes funciones
void matrixInit4(Matrix& m1, Matrix& m2, Matrix& m3, Matrix& m4,
				size_t width, size_t height, double initValue);
void matrixInitCopy(Matrix& m, size_t width, size_t height, Matrix copy);
void matrixInitSum3(Matrix& m, size_t width, size_t height, Matrix sum1, Matrix sum2, Matrix sum3);
void matrixRead(const string name, Matrix& m, size_t rows, size_t columns);
void vectorRead(const string name, Vector& v, size_t size);
void matrixWrite(const string name, Matrix m, size_t rows, size_t columns);
void vectorWrite(const string name, Vector v, size_t size);
void readBin(const char *name, float *v, unsigned size);
void matrixInitTwoVec(Matrix& m, size_t height, Vector x, Vector y);*/