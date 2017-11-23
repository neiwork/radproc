//C***********************************************************************
//C                                                                      *
//C                        BISECTION ALGORITHM 2.1                       *
//C                                                                      *
//C***********************************************************************
//C
//C     TO FIND A SOLUTION TO F(X) = 0 GIVEN THE CONTINOUS FUNCTION
//C     F ON THE INTERVAL <A, B>, WHERE F(A) AND F(B) HAVE
//C     OPPOSITE SIGNS :
//C
//C     INPUT : ENDPOINTS A, B; TOLERANCE TOL;
//C              MAXIMUM INTERATIONS N0.
//C
//C     OUTPUT : APPROXIMATE SOLUTION P OR A
//C              MESSAGE THAT THE ALGORITHM FAILS.
//C

#include "bisection.h"

#define ZERO 1.0E-20


double bisection(double A, double B, fun1 fun)
{
	double FA = fun(A);
	double FB = fun(B);

	if ((FA * FB) > 0.0) {
		std::cout << "F(A) and F(B) have same sign";
	}

	double TOL = 1.0e-2;  //Tolerancia
	int N0 = 20; // 'Input maximum number of iterations '

	double P;
//     STEP 1
	int I = 1;

	//     STEP 2

	//016   IF(I.GT.N0) GOTO 020
	while (I <= N0) {
		//STEP 3
		// COMPUTE P(I)
		P = A + (B - A) / 2.0;
		double FP = fun(P);

		//          STEP 4
				//if((ABS(FP) <= 1.0E-20 && (B - A) / 2 < TOL)) {
					// PROCEDURE COMPLETED SUCCESSFULLY
					//deberia salir
				//}
		if (!(abs(FP) <= 1.0E-20 && (B - A) / 2.0 < TOL)) {
			//          STEP 5
			I = I + 1;
			//          STEP 6
			//          COMPUTE A(I) AND B(I)
			if (FA*FP > 0) {
				A = P;
				FA = FP;
			}
			else {
				B = P;
				FB = FP;
			}
		}
		//020   CONTINUE cierro el do while
	}

	//     STEP 7
	//     PROCEDURE COMPLETED UNSUCCESSFULLY
	return P;
}