#include "interpolation.h"
#include <iostream>


void locate(Vector xx, size_t n, double x, size_t& j)
{
	// Given an array xx[0..n-1] , and given a value x , returns a value j such that x is 
	//between xx[j] and xx[j+1] . xx must be monotonic, either increasing or decreasing. 
	// j=-1 or j=n-1 is returned to indicate that x is out of range.
	
	size_t ju,jm,jl;
	int ascnd;
	
	jl=-1;             // Initialize lower
	ju=n;           // and upper limits.

	ascnd = (xx[n-1] >= xx[0]);
	while (ju-jl > 1) {            // If we are not yet done,
		jm = (ju+jl) >> 1;         // compute a midpoint,
		if (x >= xx[jm] == ascnd)
			jl=jm;                 // and replace either the lower limit
		else
			ju=jm;                 // or the upper limit, as appropriate.
	}                              // Repeat until the test condition is satisfied.
	if (x == xx[0]) j=0;           // Then set the output
	else if(x == xx[n-1]) j=n-2;
	else j=jl;
}

double interpolMod(double& E, const Vector& ener, const Vector& lum, const int last)
{
	int i = 1;   //esto esta modificado porque las dist de particulas en la posicion 0 dan cero
	
    while(E > ener[i] && i < last)    //&&=and
		i=i+1;

	if ( (E > ener[last]) | (E < ener[0]) ){
		return 0.0;}
	else if(i == 1){               //en estod dos uno va 0
		return lum[1];}
	else if(i >= last){
		return lum[last];}
	else{
		return (E-ener[i-1])*(lum[i]-lum[i-1])/(ener[i]-ener[i-1])+lum[i-1];}

}


int binarySearch(const Vector& array, int lowerbound, int upperbound, double key)
{                               
	/*array es el arreglo de la dimension en la que interpolo
	key es el valor en el que interpolo
	*/

	int position = 0;
	int comparisonCount = 1;    //count the number of comparisons (optional)

	// To start, find the subscript of the middle position.
	
//  i=0
//	while (E > key[i] && i < (size - 1))
	//	i = i + 1;

	while (!(array[position] < key && array[position + 1] > key))
	{
		position = (lowerbound + upperbound) / 2;
		if (array[position] > key){
			upperbound = position;
		}
		else{
			lowerbound = position;
		}
	}

	return position;

	/*if (array[position] < key) return position;
	else if (array[position] > key) return position-1; //me quedo con el anterior al valor de la key
	else std::cout << "error";*/ 

}


int fbinarySearch(float *array, int lowerbound, int upperbound, double key)
{                               
	/*array es el arreglo de la dimension en la que interpolo
	key es el valor en el que interpolo
	*/

	int position = 0;
	int comparisonCount = 1;    //count the number of comparisons (optional)

	// To start, find the subscript of the middle position.
	
//  i=0
//	while (E > key[i] && i < (size - 1))
	//	i = i + 1;

	if (position < upperbound-1) {
		while (!(array[position] < key && array[position + 1] > key))
		{
			cout << position << "\t" << array[position] << "\t" << array[position+1] << "\t" << key << endl;
			cout << lowerbound << "\t" << upperbound << endl;	
			position = (lowerbound + upperbound) / 2;
			if (array[position] > key){
				upperbound = position;
			}
			else{
				lowerbound = position;
			}
		}
		return position;
	} else {
		cout << "ERROR RETURN -1" << endl;
		return -1;
	}
	
	/*if (array[position] < key) return position;
	else if (array[position] > key) return position-1; //me quedo con el anterior al valor de la key
	else std::cout << "error";*/ 
}



//double interpol(double& E, const Vector& ener, const Vector& lum, const int last, const int first=0);
double interpol(double& E, const Vector& key, const Vector& val, const int size, const int base) //const int first=0)
{                                                                             //el igual indica q puedo llamar sin poner
	                                                                          //el ultimo arg
	int i = 0;
	int last = base+size-1;

	//binarySearch(const Vector& array, int lowerbound, int upperbound, int key)
	//int s = binarySearch(key, 0, last, E); //VER porque no da bien
	
		// TODO binary search here.
	while (E > key[i] && i < (size-1))
		i=i+1;

	if ( (E > key[size-1]) || (E < key[0]) ){
		return 0.0;
	} else if(i == 0) {
		return val[base];
	} else if(i >= last) {
		return val[last];
	} else {
		return (E-key[i-1])*(val[base+i]-val[base+i-1])/(key[i]-key[i-1])+val[base+i-1];
	}
}

double interpolNew(double& E, const Vector& key, const Vector& val, const int size, const int base) //const int first=0)
{                                                                             //el igual indica q puedo llamar sin poner
	//el ultimo arg
	//int i = 0;
	int last = base + size - 1;

	//binarySearch(const Vector& array, int lowerbound, int upperbound, int key)
	int i = binarySearch(key, 0, last, E); //VER  no da bien

	// TODO binary search here.
	//while (E > key[i] && i < (size - 1))
	//	i = i + 1;

	if ((E > key[size - 1]) || (E < key[0])){
		return 0.0;
	}
	else if (i == 0) {
		return val[base];
	}
	else if (i >= last) {
		return val[last];
	}
	else {
		return (E - key[i - 1])*(val[base + i] - val[base + i - 1]) / (key[i] - key[i - 1]) + val[base + i - 1];
	}
}


double interpolDoble(double E, double t, const Vector& ener, const Vector& time, const Vector& dist)
{

	Vector value(4, 0.0);
	Vector ratio(4, 0.0);
	
	int ne = ener.size()-1;
	int nt = time.size()-1;
	
	int i = 0;  
	int j = 0;
	
    while(E >= ener[i] && i < ne)    //&&=and
		i=i+1;

	while(t >= time[j] && j < nt)    //&&=and
		j=j+1;

	if((E < ener[0]) | (E > ener[ne]) | (t < time[0]) | (t > time[nt]))
	{
		return 0.0;
	}

	value[0] = dist[(j-1)*(ne+1)+(i-1)];  //(grid[energyIndexBelow][zIndexBelow];
	value[1] = dist[(j-1)*(ne+1)+(i)];    //grid[energyIndexAbove][zIndexBelow];
	value[2] = dist[(j)*(ne+1)+(i-1)];    //grid[energyIndexBelow][zIndexAbove];
	value[3] = dist[(j)*(ne+1)+(i)];      //grid[energyIndexAbove][zIndexAbove];
		
/*	ratio[0] = (energyPoints[energyIndexAbove] - particleEnergy)/
			   (energyPoints[energyIndexAbove] - 
			    energyPoints[energyIndexBelow]);*/ 
	ratio[0] = (ener[i]-E)/(ener[i]-ener[i-1]);

/*	ratio[1] = (particleEnergy - energyPoints[energyIndexBelow])/
			   (energyPoints[energyIndexAbove] - 
			    energyPoints[energyIndexBelow]);*/ 
	ratio[1] = (E-ener[i-1])/(ener[i]-ener[i-1]);

/*	ratio[2] = (zPoints[zIndexAbove] - z)/
			   (zPoints[zIndexAbove] - zPoints[zIndexBelow]);*/
	ratio[2] = (time[j]-t)/(time[j]-time[j-1]);

/*	ratio[3] = (z - zPoints[zIndexBelow])/
			   (zPoints[zIndexAbove] - zPoints[zIndexBelow]);*/
	ratio[3] = (t-time[j-1])/(time[j]-time[j-1]);

	double valueEnergyZBelow = ratio[0]*value[0] + ratio[1]*value[1];
	double valueEnergyZAbove = ratio[0]*value[2] + ratio[1]*value[3];

	double valueEnergyZ = ratio[2]*valueEnergyZBelow + ratio[3]*valueEnergyZAbove; 
	
	return valueEnergyZ;

}









/*	if(i == 0 && j == 0){    
		return dist[0];}
	else if(i >= lastE && j >= lastt){
		return dist[lastt*(lastE+1)+lastE];}   //este es el ultimo elemento
	else if(j == 0 && i != 0){    
		return (E-ener[i-1])*(dist[i]-dist[i-1])/(ener[i]-ener[i-1])+dist[i-1];}
    else if(i == 0 && j != 0){    
		return (t-time[j-1])*(dist[j*(lastE+1)]-dist[(j-1)*(lastE+1)])/(time[j]-time[j-1])+dist[(j-1)*(lastE+1)];}
	else if(j >= lastt && i < lastE){    
		return (E-ener[i-1])*(dist[i]-dist[i-1])/(ener[i]-ener[i-1])+dist[i-1];}
    else if(i >= lastE && j >= lastt){    
		return (t-time[j-1])*(dist[j*lastE]-dist[(j-1)*lastE])/(time[j]-time[j-1])+dist[(j-1)*lastE];}*/ 