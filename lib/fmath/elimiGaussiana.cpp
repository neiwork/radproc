/*  GAUSSIAN ELIMINATION WITH BACKWARD SUBSTITUTION ALGORITHM 
     TO SOLVE THE N BY N LINEAR SYSTEM

   E1:  A(1,1) X(1) + A(1,2) X(2) + ... + A(1,N) X(N) = A(1,N+1)
   E2:  A(2,1) X(1) + A(2,2) X(2) + ... + A(2,N) X(N) = A(2,N+1)
   :
   .
   EN:  A(N,1) X(1) + A(N,2) X(2) + ... + A(N,N) X(N) = A(N,N+1)

     INPUT:   NUMBER OF UNKNOWNS AND EQUATIONS N; AUGMENTED
              MATRIX A = (A(I,J)) WHERE 1<=I<=N AND 1<=J<=N+1.

     OUTPUT:  SOLUTION X(1),X(2),...,X(N) OR A MESSAGE THAT THE LINEAR
              SYSTEM HAS NO UNIQUE SOLUTION.
*/

#include "elimiGaussiana.h"

#define ZERO 1.0E-20
#define true 1
#define false 0

/* Absolute Value Function */
double absval(double val)
{
	if (val >= 0) return val;
	else return -val;
}

//void INPUT(int *, double [][13], int *);

void elimiGaussiana( int N, Matrix& A, Vector& X )
{	
   double C,XM,SUM;
   int M,I,J,ICHG,NN,IP,JJ,K,KK,OK;

      /* STEP 1 */
      NN = N - 1;
      M = N + 1;
      ICHG = 0;
      I = 1;
      while ( (I <= NN)) {
         /* STEP 2 */
         /* use IP instead of p */
         IP = I;
         
		 while ((absval(A[IP-1][I-1]) <= ZERO) && (IP <= N)) 
            IP++;		 

         if (IP == M) OK = false;
         else {
            /* STEP 3 */
            if (IP != I) {
               for (JJ=1; JJ<=M; JJ++) {
                  C = A[I-1][JJ-1];
                  A[I-1][JJ-1] = A[IP-1][JJ-1];
                  A[IP-1][JJ-1] = C;
               }
               ICHG = ICHG + 1;
            }  
            /* STEP 4 */
            JJ = I + 1; 
            for (J=JJ; J<=N; J++) {
               /* STEP 5 */
               /* use XM in place of m(J,I) */
               XM = A[J-1][I-1] / A[I-1][I-1]; 
               /* STEP 6 */
               for (K=JJ; K<=M; K++)
                  A[J-1][K-1] = A[J-1][K-1] - XM * A[I-1][K-1];
               /*  Multiplier XM could be saved in A[J,I].  */
               A[J-1][I-1] = 0.0;
            }  
         }
         I = I + 1;
      }
     
         /* STEP 7 */
         if (absval(A[N-1][N-1]) <= ZERO) OK = false;
         else {
            /* STEP 8 */
            /* start backward substitution */
            X[N-1] = A[N-1][M-1] / A[N-1][N-1];
            /* STEP 9 */
            for (K=1; K<=NN; K++) {
               I = NN - K + 1;
               JJ = I + 1;
               SUM = 0.0;
               for (KK=JJ ; KK<=N; KK++)  
                  SUM = SUM - A[I-1][KK-1] * X[KK-1];
               X[I-1] = (A[I-1][M-1]+SUM) / A[I-1][I-1];
            }
            /* STEP 10 */

         }  
      
}





/*
prueba::

	Vector X(3,0.0);
	Matrix a;
	matrixInit( a, X.size(), X.size()+1, 0 );
	a[0][0]=1; a[0][1]=2;	a[0][2]=3;  a[0][3] = 1;
	a[1][0] = 4; a[1][1] = 5; a[1][2] = 6; a[1][3] = -2;
	a[2][0]=7; a[2][1]=8; a[2][2]=10; a[2][3]=5;

	elimiGaussiana (X.size(), a, X);

	double pepe = X[1];*/

/*
	int N = N2;  // 0 <= I < N AND 0 <= J < M=N+1
      
	int M = N+1;

	for (int i=0; i < N-1; ++i)
	{

		int IP = i;				 //USE IP IN PLACE OF P

		double absolute = abs(A[IP][i]);

		while (absolute  < 1e-20)
			IP+=1 ;
			
		if(IP > N)      
			std::cout <<  "error" ;

		if(IP != i)    //  ne !=
		{
			for (int JJ=0; JJ < M; ++JJ)   //XXX     //aca hago el cambio de fila Ep con Ei
			{
				double aux = A[i][JJ];
				A[i][JJ] = A[IP][JJ];
				A[IP][JJ] = aux;
			}
		}

		for (int j=i+1; j < N; ++j) //XX
		{
			double XM = A[j][i]/A[i][i];					//USE XM IN PLACE OF M(J,I)

			for (int k=i+1; k < M; ++k)  //XX
			{
				A[j][k] = A[j][k]-XM*A[i][k];
			}

			A[j][i] = 0;

		}

	}

	if(abs(A[N][N]) < 1E-20) 	
		std::cout << "error";
			
	else {

//		START BACKWARD SUBSTITUTION
				
		X[N] = A[N][N+1]/A[N][N];  
		
		int L = N-1;
		
		for (int k=0; k < L; ++k) //XX
		{
			int i   = L-k;// +1; //este +1 está mal!
			int JJ  = i+1;

			double SUM = 0.0;

			for (int KK=JJ; KK < N; ++KK)
			{			
				SUM = SUM-A[i][KK]*X[KK];
			}

			X[i] = (A[i][N+1]+SUM)/A[i][i];
		}

	}

}*/ 