/*
 *	This file is part of qpDUNES.
 *
 *	qpDUNES -- A DUal NEwton Strategy for convex quadratic programming.
 *	Copyright (C) 2012 by Janick Frasch, Hans Joachim Ferreau, et al.
 *	All rights reserved.
 *
 *	qpDUNES is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	qpDUNES is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with qpDUNES; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


/**
 *	\file examples/nmpcPrototype.c
 *	\author Janick Frasch
 *	\version 1.0beta
 *	\date 2012
 *
 *	Example for NMPC coupling
 *
 *	In the context of NMPC qpDUNES solves a problem
 *
 *	min  sum_{k=0..nI} z_k'*H_k*z_k + g_k'*z_k
 *	s.t. x_{k+1} = C_k * z_k + c_k				for k=0..nI-1
 *	     dLow <= D * z_k <= dUpp				for k=0..nI
 *
 *	where x_k is implicitly defined by z_k = [ x_k  u_k ] as the first nX variables of z_k
 *
 *	It holds
 *	z_k  \in R^nZ  for k=0..nI-1
 *	z_nI \in R*nX
 *
 *	nX < nZ
 *	nU = nZ - nX
 *
 *
 *	USEFUL OPTIONS:
 *	qpOptions.logLevel = QPDUNES_LOG_OFF;	// switches off logging completely (can decrease memory consumption significantly)
 *
 *	COMPILER FLAGS:
 *  -D__SUPPRESS_ALL_WARNINGS__				// suppress warnings
 *  -D__SUPPRESS_ALL_OUTPUT__				// do not print anything at all
 *  -U__MEASURE_TIMINGS__					// exclude code for detailed runtime profiling in qpDUNES
 *  -U__USE_ASSERTS__						// switch some safty checks off
 *
 */


#include <qpDUNES.h>


#define INFTY 1.0e12

//#define USE_D
#undef USE_D


int main( )
{
	return_t statusFlag;

	int iter, ii, kk;

	/** define problem dimensions */
	const unsigned int nI = 1;			/* number of control intervals */
	const int nX = 16;					/* number of states */
	const int nU = 8;					/* number of controls */
	const unsigned int nZ = nX+nU;		/* number of stage variables */
#ifdef USE_D
	unsigned int nD[nI+1];  			/* number of constraints */
	for ( kk=0; kk<nI; ++kk ) {
		nD[kk] = nX+nU;
	}
	nD[nI] = nX;
#else
	unsigned int* nD = 0;	  			/* number of constraints */
#endif /* USE_D */


	/** define problem data */
	const double H[] = 
	{
		0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 96.2461, 96.2361, -2.29185e-018, 96.2361, 0.00803458, -190.945, 0, 0, 0, 0, 0, 0, 0, 0, -1.19191, 0.000633004, -35.832, -23.0887, 0.000550931, -14.0715, -0.000672604, 14.4235,
		0, 0, 96.2361, 96.2461, -1.43167e-018, 96.2361, 0.00803458, -190.945, 0, 0, 0, 0, 0, 0, 0, 0, -1.19191, 0.000633004, -35.832, -23.0887, 0.000550931, -14.0715, -0.000672604, 14.4235,
		0, 0, -2.29185e-018, -1.43167e-018, 23.3474, -0.0247478, 28.8497, 0.0118034, 0, 0, 0, 0, 0, 0, 0, 0, 4.79509, 6.85805, 6.35749e-017, 2.15011e-016, -1.6793, 0.00420693, -2.6386, -0.000901843,
		0, 0, 96.2361, 96.2361, -0.0247478, 96.2461, -0.0225587, -190.945, 0, 0, 0, 0, 0, 0, 0, 0, -1.19699, -0.00663952, -35.832, -23.0887, 0.00233172, -14.0715, 0.00212546, 14.4235,
		0, 0, 0.00803458, 0.00803458, 28.8497, -0.0225587, 35.674, -0.00135022, 0, 0, 0, 0, 0, 0, 0, 0, 5.92759, 8.47792, -0.00299155, -0.00192763, -2.07595, 0.0040258, -3.26183, 8.93338e-005,
		0, 0, -190.945, -190.945, 0.0118034, -190.945, -0.00135022, 378.868, 0, 0, 0, 0, 0, 0, 0, 0, 2.36732, 0.00221266, 71.0953, 45.8109, -0.00194246, 27.9196, -1.86218e-016, -28.618,
		0, 0, 0, 0, 0, 0, 0, 0, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 1.01, -9.15134e-018, 1.257e-018, 0.870344, 0.000295467, 0.793355, 1.76847e-005, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, -9.15134e-018, 1.01, 1, 9.15473e-018, 1, 8.34882e-005, -1.98413, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 1.257e-018, 1, 1.01, 1.62088e-017, 1, 8.34882e-005, -1.98413, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0.870344, 9.15473e-018, 1.62088e-017, 1.01, 2.0383e-017, 0.990272, 0.000138043, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0.000295467, 1, 1, 2.0383e-017, 1.01, 1.25225e-017, -1.98413, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0.793355, 8.34882e-005, 8.34882e-005, 0.990272, 1.25225e-017, 1.01, -4.66628e-019, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 1.76847e-005, -1.98413, -1.98413, 0.000138043, -1.98413, -4.66628e-019, 3.94676, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, -1.19191, -1.19191, 4.79509, -1.19699, 5.92759, 2.36732, 0, 0, 0, 0, 0, 0, 0, 0, 1.01, 1.4091, 0.443788, 0.285959, -0.345049, 0.175143, -0.542138, -0.178823,
		0, 0, 0.000633004, 0.000633004, 6.85805, -0.00663952, 8.47792, 0.00221266, 0, 0, 0, 0, 0, 0, 0, 0, 1.4091, 2.02534, -0.000235689, -0.000151869, -0.493487, 0.00114371, -0.775391, -0.000170148,
		0, 0, -35.832, -35.832, 6.35749e-017, -35.832, -0.00299155, 71.0953, 0, 0, 0, 0, 0, 0, 0, 0, 0.443788, -0.000235689, 14.1972, 9.39594, -0.000234353, 5.61835, 0.000245957, -5.15167,
		0, 0, -23.0887, -23.0887, 2.15011e-016, -23.0887, -0.00192763, 45.8109, 0, 0, 0, 0, 0, 0, 0, 0, 0.285959, -0.000151869, 9.39594, 6.30466, -0.000159793, 3.7342, 0.000157138, -3.25377,
		0, 0, 0.000550931, 0.000550931, -1.6793, 0.00233172, -2.07595, -0.00194246, 0, 0, 0, 0, 0, 0, 0, 0, -0.345049, -0.493487, -0.000234353, -0.000159793, 0.130838, -0.000396372, 0.189866, 0.000139909,
		0, 0, -14.0715, -14.0715, 0.00420693, -14.0715, 0.0040258, 27.9196, 0, 0, 0, 0, 0, 0, 0, 0, 0.175143, 0.00114371, 5.61835, 3.7342, -0.000396372, 2.2374, -0.000379307, -2.01095,
		0, 0, -0.000672604, -0.000672604, -2.6386, 0.00212546, -3.26183, -1.86218e-016, 0, 0, 0, 0, 0, 0, 0, 0, -0.542138, -0.775391, 0.000245957, 0.000157138, 0.189866, -0.000379307, 0.308327, 2.83298e-017,
		0, 0, 14.4235, 14.4235, -0.000901843, 14.4235, 8.93338e-005, -28.618, 0, 0, 0, 0, 0, 0, 0, 0, -0.178823, -0.000170148, -5.15167, -3.25377, 0.000139909, -2.01095, 2.83298e-017, 2.22829,

		0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0.01, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0.01, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0.01, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.01, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.01, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.01, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.01, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.01, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.01,
	};
	
	const double g[] = 
	{
		0, 0, -9.79127, -9.79127, 0.000182849, -9.79127, -0.000591417, 19.4271, 0, 0, 0, 0, 0, 0, 0, 0, 0.121305, -1.06702e-005, 3.64581, 2.34926, -6.92165e-005, 1.43174, 4.77577e-005, -1.46743,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	};

	const double C[] =
	{	
		1, 0, 0, 0, 0, 0, 0, 0, 0.05, 0, 0, 0, 0, 0, 0, 0, 0.00125, 0, 0, 0, 0, 0, 0, 0,
		0, 1, 0, 0, 0, 0, 0, 0, 0, 0.05, 0, 0, 0, 0, 0, 0, 0, 0.00125, 0, 0, 0, 0, 0, 0,
		0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0.05, 0, 0, 0, 0, 0, 0, 0, 0.00125, 0, 0, 0, 0, 0,
		0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0.05, 0, 0, 0, 0, 0, 0, 0, 0.00125, 0, 0, 0, 0,
		0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0.05, 0, 0, 0, 0, 0, 0, 0, 0.00125, 0, 0, 0,
		0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0.05, 0, 0, 0, 0, 0, 0, 0, 0.00125, 0, 0,
		0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0.05, 0, 0, 0, 0, 0, 0, 0, 0.00125, 0,
		0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0.05, 0, 0, 0, 0, 0, 0, 0, 0.00125,
		0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0.05, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0.05, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0.05, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0.05, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0.05, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0.05, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0.05, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0.05,
	};
	
	const double c[] =
	{ 
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	};

#ifdef USE_D
	const double* zLow = 0;
	const double* zUpp = 0;
#else
	const double zLow[] =
	{
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1.078, -1.6762, -1.1973, -2.1893, -0.5644, -1.6249, -1.317, -0.98,
		-4.6078, -1e+012, -0.6128, -1.23, -3.0165, -1.0123, -3.0154, -0.2291, -1.47, -1.1802, -0.9749, -1.1802, -1.2999, -1.2999, -2.0525, -0.245,
	};
	const double zUpp[] =
	{ 
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.078, 1.6762, 1.1973, 2.1893, 0.5644, 1.6249, 1.317, 0.98,
		4.4922, 1e+012, 0.6592, 1.0655, 3.0153, 0.7331, 3.0164, 0.3612, 1.47, 1.1802, 0.9749, 1.1802, 1.2999, 1.2999, 2.0525, 0.245,
	};
#endif /* USE_D */


#ifdef USE_D
	const double D[] =
		{	1.0, 	0.0, 	0.0,
			0.0, 	1.0, 	0.0,
			0.0, 	0.0,	1.0,

			1.0, 	0.0, 	0.0,
			0.0, 	1.0, 	0.0,
			0.0, 	0.0,	1.0,

			1.0, 	0.0, 	0.0,
			0.0, 	1.0, 	0.0,
			0.0, 	0.0,	1.0,

			1.0, 	0.0,
			0.0, 	1.0
		};
	const double dLow[] =
		{	-1.9, -3.0, -30.0,
			-1.9, -3.0, -30.0,
			-1.9, -3.0, -30.0,
			-1.9, -3.0
		};

	const double dUpp[] =
		{	 1.9,  3.0,  30.0,
			 1.9,  3.0,  30.0,
			 1.9,  3.0,  30.0,
			 1.9,  3.0
		};
#else
	const double* D = 0;
	const double* dLow = 0;
	const double* dUpp = 0;
#endif /* USE_D */



	/** define simulation environment */
	double x0[nX] = { 
		4.8078,
		0.1218,
		- 1.5319,
		0.476,
		0.0006,
		0.1396,
		- 0.0005,
		0.7991,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0 
	};	/* initial value */
#ifdef USE_D
	double* z0Low = 0;
	double* z0Upp = 0;
	double d0Low[nX+nU];
	double d0Upp[nX+nU];
#else
	double z0Low[nX+nU];			/* auxiliary variables */
	double z0Upp[nX+nU];
	double* d0Low = 0;
	double* d0Upp = 0;
#endif /* USE_D */

	double zOpt[nI*nZ+nX];			/* primal solution */
	double lambdaOpt[nI*nX];		/* dual solution */
	double muOpt[2*nI*nZ+2*nX];

	const unsigned int nSteps = 2;	/* number of simulation steps */


	/** (1) set up a new qpDUNES problem */
	qpData_t qpData;


	/** (2) set qpDUNES options */
	qpOptions_t qpOptions = qpDUNES_setupDefaultOptions();
// 	qpOptions.maxIter    = 100;
// 	qpOptions.printLevel = 2;
// 	qpOptions.stationarityTolerance = 1.e-6;


	/** (3) allocate data for qpDUNES and set options */
	qpDUNES_setup( &qpData, nI, nX, nU, nD, &(qpOptions) );
	

	/** (4) set sparsity of primal Hessian and local constraint matrix */
	for ( kk=0; kk<nI+1; ++kk ) {
		//qpData.intervals[kk]->H.sparsityType = QPDUNES_DIAGONAL;
		//qpData.intervals[kk]->D.sparsityType = QPDUNES_IDENTITY;
	}
	
	/** (5) initial MPC data setup: components not given here are set to zero (if applicable)
	 *      instead of passing g, D, zLow, zUpp, one can also just pass NULL pointers (0) */
	statusFlag = qpDUNES_init( &qpData, H, g, C, c, zLow,zUpp, D,dLow,dUpp );	/* todo: add constraint vectors, make non-trivial example */
	if (statusFlag != QPDUNES_OK) {
		printf( "Data init failed.\n" );
		return (int)statusFlag;
	}
	

	/** MAIN MPC SIMULATION LOOP */
 	for ( iter=0; iter<nSteps; ++iter ) {
 		/** (1) embed current initial value */
		for ( ii=0; ii<nX; ++ii ) {
#ifdef USE_D
			d0Low[ii] = x0[ii];
			d0Upp[ii] = x0[ii];
#else
			z0Low[ii] = x0[ii];
			z0Upp[ii] = x0[ii];
#endif /* USE_D */
		}
		for ( ii=nX; ii<nX+nU; ++ii ) {
#ifdef USE_D
			d0Low[ii] = dLow[ii];
			d0Upp[ii] = dUpp[ii];
#else
			z0Low[ii] = zLow[ii];
			z0Upp[ii] = zUpp[ii];
#endif /* USE_D */
		}
		statusFlag = qpDUNES_updateIntervalData( &qpData, qpData.intervals[0], 0, 0, 0, 0, z0Low,z0Upp, D,d0Low,d0Upp, 0 );
		if (statusFlag != QPDUNES_OK) {
			printf( "Initial value embedding failed.\n" );
			return (int)statusFlag;
		}

#ifdef USE_D
		qpDUNES_printMatrixData( d0Low, 1, nZ, "d0Low@it%d", iter);
		qpDUNES_printMatrixData( d0Upp, 1, nZ, "d0Upp@it%d", iter);
#else
		qpDUNES_printMatrixData( z0Low, 1, nZ, "z0Low@it%d", iter);
		qpDUNES_printMatrixData( z0Upp, 1, nZ, "z0Upp@it%d", iter);
#endif /* USE_D */


		/** (2) solve QP */
		statusFlag = qpDUNES_solve( &qpData );
		if (statusFlag != QPDUNES_SUCC_OPTIMAL_SOLUTION_FOUND) {
			printf( "QP solution %d failed.\n", iter );
			return (int)statusFlag;
		}
		
		// up to here in fdb step


		/** (3) obtain primal and dual optimal solution */
		qpDUNES_getPrimalSol( &qpData, zOpt );
		qpDUNES_getDualSol( &qpData, lambdaOpt, muOpt );
		qpDUNES_printMatrixData( lambdaOpt, 1, nI*nX, "dual lambda" );
		qpDUNES_printMatrixData( muOpt, 1, 2*nI*nZ + 2*nX, "dual mu" );
		/// ...
		

 		/** (4) prepare QP for next solution */
		qpDUNES_shiftLambda( &qpData );			/* shift multipliers */
		qpDUNES_shiftIntervals( &qpData );		/* shift intervals (particulary important when using qpOASES for underlying local QPs) */
		
		// optional
		
		/// H = ...
		/// g = ...
		/// C = ...
		/// c = ...
		/// zLow = ...
		/// zUpp = ...
		statusFlag = qpDUNES_updateData( &qpData, H, g, C, c, zLow,zUpp, D,dLow,dUpp );		/* data update: components not given here keep their previous value */
		if (statusFlag != QPDUNES_OK) {
			printf( "Data update failed.\n" );
			return (int)statusFlag;
		}
		
		// mandatory


		/** (5) simulate next initial value */
		for (ii=0; ii<nX; ++ii) {
			/// x0 = ...
			x0[ii] = zOpt[1*nZ+ii];
 		}
	}


	/** cleanup of allocated data */
	qpDUNES_cleanup( &qpData );


	return 0;
}


/*
 *	end of file
 */
