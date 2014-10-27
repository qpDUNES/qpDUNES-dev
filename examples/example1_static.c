/*
 *	This file is part of qpDUNES.
 *
 *	qpDUNES -- A DUal NEwton Strategy for convex quadratic programming.
 *	Copyright (C) 2012 by Janick Frasch, Hans Joachim Ferreau et al.
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
 *	\file examples/example1_static.c
 *	\author Janick Frasch
 *	\version 1.0beta
 *	\date 2014
 *
 *	Very simple example for testing the static memory version of qpDUNES.
 */


#define __STATIC_MEMORY__

#include <qpDUNES.h>

#define INFTY 1.0e12


int main( )
{
	int i;
	boolean_t isLTI;
	
	unsigned int nI = 2;
	unsigned int nX = 3;
	unsigned int nU = 2;
	unsigned int* nD = 0;
	
	double Q[3*3] =
		{	1.0, 0.0, 0.0, 
			0.0, 1.0, 0.0,
			0.0, 0.0, 1.0	};
			
	double x0[3] = 
		{ 2.0, 3.0, 4.0 };
		
		
	double R[2*2] =
		{	1.0, 0.0,
			0.0, 1.0	};
	double *S=0;
	
	double* P = Q;
	
	double A[3*3] =
		{	1.0, 0.0, 0.0, 
			0.0, 1.0, 0.0,
			0.0, 0.0, 1.0	};
	double B[3*2] =
		{	1.0, 0.0,
			0.0, 1.0,
			1.0, 1.0	};
	double c[3] = 
		{	5.0,
			5.0,
			5.0	};
	
	double xLow[3] = { -INFTY, -INFTY, -INFTY };
	double xUpp[3] = {  INFTY,  INFTY,  INFTY };
	double uLow[2] = { -INFTY, -INFTY };
	double uUpp[2] = {  INFTY,  INFTY };
	
	qpData_t qpData;


	qpDUNES_setupStatic( &qpData, nD, 0 );	/* passing 0 in the last argument sets the default QP options */

//	for( i=0; i<nI; ++i )
//	{
//		qpDUNES_setupSimpleBoundedInterval(  &qpData, qpData.intervals[i],Q,R,S, A,B,c, xLow,xUpp,uLow,uUpp );
//	}
//	qpDUNES_setupSimpleBoundedInterval(  &qpData, qpData.intervals[nI], P,0,0, 0,0,0, xLow,xUpp,0,0 );
//
//	qpDUNES_setupAllLocalQPs( &qpData, isLTI=QPDUNES_TRUE );	/* determine local QP solvers and set up auxiliary data */
//
//
//	qpDUNES_solve( &qpData );
//
//	for( i=0; i<nI; ++i )
//	{
//		qpDUNES_printMatrixData( qpData.intervals[i]->z.data, 1, nX+nU, "z[%d]:", i );
//	}
//	qpDUNES_printMatrixData( qpData.intervals[nI]->z.data, 1, nX, "z[%d]:", i );
	
	qpData.intervals[1]->yPrev;
	
	printf( "example1 done.\n" );
	
	return 0;
}


/*
 *	end of file
 */
