#include "matrixInit.h"

void matrixInit(Matrix& m, size_t height, size_t width, double initValue )
{
		Vector row(width,initValue);
		Matrix aux(height,row);
		m = aux;
}

void matrixInitCopy(Matrix& m, size_t height, size_t width, Matrix copy)
{
		Vector row(width,0.0);
		Matrix aux(height,row);
		m = aux;
		
		for (std::size_t i=0;i<height;i++) {
			for (std::size_t j=0;j<width;j++) {
				m[i][j] = copy[i][j];
			}
		}
}

void matrixInitTwoVec(Matrix& m, size_t height, Vector x, Vector y)
{
	for (std::size_t i=0;i<height;i++) {
		m[i][0] = x[i];
		m[i][1] = y[1];
	}
}

void matrixInitSum3(Matrix& m, size_t height, size_t width, Matrix sum1, Matrix sum2, Matrix sum3)
{
		Vector row(width,0.0);
		Matrix aux(height,row);
		m = aux;
		
		for (std::size_t i=0;i<height;i++) {
			for (std::size_t j=0;j<width;j++) {
					m[i][j] = sum1[i][j]+sum2[i][j]+sum3[i][j];
			}
		}
}

void matrixInit4(Matrix& m1, Matrix& m2, Matrix& m3, Matrix& m4,
				size_t height, size_t width, double initValue)
{
		Vector row(width,initValue);
		Matrix aux(height,row);
		m1 = aux;  m2 = aux;  m3 = aux;  m4 = aux;
}

using namespace std;
#include <fstream>

void matrixRead(const string filename, Matrix& m, size_t rows, size_t columns) {
    
    ifstream f(filename);
    
    double aux;
    
    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < columns; j++)
            f >> m[i][j];
	f.close();
}

void matrixWrite(const string filename, Matrix m, size_t rows, size_t columns) {
	
	ofstream f(filename);
	
	for (size_t i=0;i<rows;i++) {
		for (size_t j=0;j<columns; j++)
			f << m[i][j] << "\t";
		f << endl;
	}
	f.close();
}

void vectorRead(const string filename, Vector& v, size_t size) 
{    
    ifstream f(filename);
    for (size_t i = 0; i < size; i++)
            f >> v[i];
	f.close();
}

void vectorWrite(const string filename, Vector v, size_t size)
{
	ofstream f(filename);
	for (size_t i=0;i<size;i++)
		f << v[i] << endl;
	
	f.close();
}

/*void readBin(const string filename, float *v, unsigned size)
{
	ifstream file;
	file.open(filename,ios::binary);
	file.read();
	fclose(ptr);
}*/