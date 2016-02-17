#include "matrixInit.h"

void matrixInit(Matrix& m, int height, int width, double initValue ) {

		Vector row(width,initValue);
		Matrix aux(height,row);
		m = aux;

	}