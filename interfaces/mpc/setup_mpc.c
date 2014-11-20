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
 *	\file interfaces/mpc/mpcDUNES.cpp
 *	\author Jancik Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 *
 *	Interface for comfortable use of qpDUNES in an MPC context
 *
 */


#include <setup_mpc.h>



/* ----------------------------------------------
 * allocate memory
 * 
 # >>>>>>                                           */
return_t mpcDUNES_setup(	mpcProblem_t* const mpcProblem,
							uint_t nI,
							uint_t nX,
							uint_t nU,
							uint_t* nD,
							qpOptions_t* qpOptions
							)
{
	int_t ii;
	int_t nZ = nX+nU;

	/* allocate qpData struct */
	qpDUNES_setup( &(mpcProblem->qpData), nI, nX, nU, nD, qpOptions );
	
	/* allocate xRef, uRef, xOpt, uOpt, lambdaOpt */
	mpcProblem->zRef  = (real_t*)calloc( nI*nZ+nX,sizeof(real_t) );
//	mpcProblem->xRef  = (real_t*)calloc( (nI+1)*nX,sizeof(real_t) );
//	mpcProblem->uRef  = (real_t*)calloc( nI*nU,sizeof(real_t) );
	mpcProblem->xOpt  = (real_t*)calloc( (nI+1)*nX,sizeof(real_t) );
	mpcProblem->uOpt  = (real_t*)calloc( nI*nU,sizeof(real_t) );
	mpcProblem->lambdaOpt  = (real_t*)calloc( nI*nX,sizeof(real_t) );
	
	/* allocate workspace */
	mpcProblem->xnTmp   = (real_t*)calloc( nI*nX,sizeof(real_t) );
	mpcProblem->zn1Tmp   = (real_t*)calloc( nI*nZ+nX,sizeof(real_t) );
	mpcProblem->zn1Tmp2   = (real_t*)calloc( nI*nZ+nX,sizeof(real_t) );
	mpcProblem->xn1Tmp   = (real_t*)calloc( (nI+1)*nX,sizeof(real_t) );
	mpcProblem->xn1Tmp2  = (real_t*)calloc( (nI+1)*nX,sizeof(real_t) );
	mpcProblem->xn1Tmp3  = (real_t*)calloc( (nI+1)*nX,sizeof(real_t) );
	mpcProblem->unTmp   = (real_t*)calloc( nI*nU,sizeof(real_t) );
	mpcProblem->unTmp2  = (real_t*)calloc( nI*nU,sizeof(real_t) );
	
	mpcProblem->z0LowOrig = (real_t*)calloc( nZ,sizeof(real_t) );	/* workspace to save bounds before initial value embedding */
	mpcProblem->z0UppOrig = (real_t*)calloc( nZ,sizeof(real_t) );

	/* initalize solution variables */
	for ( ii=0; ii<(nI+1)*nX; ++ii ) {
		mpcProblem->xOpt[ii] = -mpcProblem->qpData.options.QPDUNES_INFTY;
	}
	for ( ii=0; ii<nI*nU; ++ii ) {
		mpcProblem->uOpt[ii] = -mpcProblem->qpData.options.QPDUNES_INFTY;
	}
	for ( ii=0; ii<nI*nX; ++ii ) {
		mpcProblem->lambdaOpt[ii] = -mpcProblem->qpData.options.QPDUNES_INFTY;
	}
	mpcProblem->optObjVal = -mpcProblem->qpData.options.QPDUNES_INFTY;
	mpcProblem->exitFlag = QPDUNES_UNTERMINATED;

	return QPDUNES_OK;
}
/*<<< END OF mpcDUNES_setup */


/* ----------------------------------------------
 * free memory
 * 
 # >>>>>>           						*/
return_t mpcDUNES_cleanup(	mpcProblem_t* const mpcProblem
							)
{
	qpDUNES_cleanup( &(mpcProblem->qpData) );
	
	qpDUNES_free( &(mpcProblem->zRef) );
	qpDUNES_free( &(mpcProblem->xOpt) );
	qpDUNES_free( &(mpcProblem->uOpt) );
	qpDUNES_free( &(mpcProblem->lambdaOpt) );
	
	qpDUNES_free( &(mpcProblem->xnTmp) );
	qpDUNES_free( &(mpcProblem->zn1Tmp) );
	qpDUNES_free( &(mpcProblem->zn1Tmp2) );
	qpDUNES_free( &(mpcProblem->xn1Tmp) );
	qpDUNES_free( &(mpcProblem->xn1Tmp2) );
	qpDUNES_free( &(mpcProblem->xn1Tmp3) );
	qpDUNES_free( &(mpcProblem->unTmp) );
	qpDUNES_free( &(mpcProblem->unTmp2) );

	qpDUNES_free( &(mpcProblem->z0LowOrig) );
	qpDUNES_free( &(mpcProblem->z0UppOrig) );
	
	return QPDUNES_OK;
}
/*<<< END OF qpDUNES_cleanup */


/* ----------------------------------------------
 * set up an linear time invariant (LTI) MPC problem with simple bounds in xu-style
 * 
 # >>>>>>                                           */
return_t mpcDUNES_initLtiSb_xu(	mpcProblem_t* const mpcProblem,
								const real_t* const Q,
								const real_t* const R,
								const real_t* const S,
								const real_t* const P,
								const real_t* const A,
								const real_t* const B,
								const real_t* const c,
								const real_t* const xLow,
								const real_t* const xUpp,
								const real_t* const uLow,
								const real_t* const uUpp,
								const real_t* const xRef,
								const real_t* const uRef
								)
{
	int_t ii,jj,kk;
	
	return_t statusFlag;

	mpcProblem->isLTI = QPDUNES_TRUE;


	qpData_t* qpData = &(mpcProblem->qpData);
//	int_t nI = qpData->nI;
//	int_t nZ = qpData->nZ;
//	int_t nX = qpData->nX;
//	int_t nU = qpData->nU;
	
	real_t* cMod = mpcProblem->xnTmp;	/* temporary variables for bound setup from workspace */
	real_t* zLowMod = mpcProblem->zn1Tmp;
	real_t* zUppMod = mpcProblem->zn1Tmp2;
	
	/** check existence of data */
	if ( !Q ) {
		qpDUNES_printError( qpData, __FILE__, __LINE__, "Q matrix missing" );
		return QPDUNES_ERR_INVALID_ARGUMENT;
	}
	if ( !R ) {
		qpDUNES_printError( qpData, __FILE__, __LINE__, "R matrix missing" );
		return QPDUNES_ERR_INVALID_ARGUMENT;
	}
	if ( !A ) {
		qpDUNES_printError( qpData, __FILE__, __LINE__, "A matrix missing" );
		return QPDUNES_ERR_INVALID_ARGUMENT;
	}
	if ( !B ) {
		qpDUNES_printError( qpData, __FILE__, __LINE__, "B matrix missing" );
		return QPDUNES_ERR_INVALID_ARGUMENT;
	}


	/** get data */

	/* get reference */
	for ( kk=0; kk<_NI_; ++kk ) {
		if ( xRef != 0 ) {
			for ( ii=0; ii<_NX_; ++ii ) {
				mpcProblem->zRef[kk*_NZ_+ii] = xRef[kk*_NX_+ii];
			}
		}
		else {
			for ( ii=0; ii<_NX_; ++ii ) {
				mpcProblem->zRef[kk*_NZ_+ii] = 0.;
			}
		}
		if ( uRef != 0 ) {
			for ( ii=0; ii<_NU_; ++ii ) {
				mpcProblem->zRef[kk*_NZ_+_NX_+ii] = uRef[kk*_NU_+ii];
			}
		}
		else {
			for ( ii=0; ii<_NU_; ++ii ) {
				mpcProblem->zRef[kk*_NZ_+_NX_+ii] = 0.;
			}
		}
	}
	

	/* translate variables to account for reference zQP := (z-zRef) */
	/* account for variable substitution in dynamics: c += [A B] * [xRef' uRef']' */
	if ( c != 0 ) {
		for ( kk=0; kk<_NI_; ++kk ) {
			for ( ii=0; ii<_NX_; ++ii )  cMod[kk*_NX_+ii] = c[ii];
		}
	}
	else {
		for ( ii=0; ii<_NI_*_NX_; ++ii )  cMod[ii] = 0.;
	}
	for ( kk=0; kk<_NI_; ++kk ) {
		for ( ii=0; ii<_NX_; ++ii ) {
			/* x part */
			if ( xRef != 0 ) {
				for ( jj=0; jj<_NX_; ++jj ) {
					cMod[kk*_NX_+ii] += A[ii*_NX_ + jj] * mpcProblem->zRef[kk*_NZ_+jj] - mpcProblem->zRef[(kk+1)*_NZ_+jj];
				}
			}
			/* u part */
			if ( uRef != 0 ) {
				for ( jj=0; jj<_NU_; ++jj ) {
					cMod[kk*_NX_+ii] += B[ii*_NU_ + jj] * mpcProblem->zRef[kk*_NZ_+_NX_+jj];
				}
			}
		}
	}

	/* account for variable substitution in constraints: [x{l,u} u{l,u}] -= [xRef uRef] */
	/* WARNING: currently no support for general constraint matrices! */
	/* lower bounds */
	/* first nI intervals */
	for ( kk=0; kk<_NI_; ++kk ) {
		/* x part */
		if ( xLow != 0 ) {
			for ( ii=0; ii<_NX_; ++ii )  zLowMod[kk*_NZ_+ii] = xLow[kk*_NX_+ii] - mpcProblem->zRef[kk*_NZ_+ii];
		}
		else {
			for ( ii=0; ii<_NX_; ++ii )  zLowMod[kk*_NZ_+ii] = -qpData->options.QPDUNES_INFTY;
		}
		/* u part */
		if ( uLow != 0 ) {
			for ( ii=0; ii<_NU_; ++ii )  zLowMod[kk*_NZ_+_NX_+ii] = uLow[kk*_NU_+ii] - mpcProblem->zRef[kk*_NZ_+_NX_+ii];
		}
		else {
			for ( ii=0; ii<_NU_; ++ii )  zLowMod[kk*_NZ_+_NX_+ii] = -qpData->options.QPDUNES_INFTY;
		}
	}
	/* last interval */
	if ( xLow != 0 ) {
		for ( ii=0; ii<_NX_; ++ii )  zLowMod[_NI_*_NZ_+ii] = xLow[_NI_*_NX_+ii] - mpcProblem->zRef[_NI_*_NZ_+ii];
	}
	else {
		for ( ii=0; ii<_NX_; ++ii )  zLowMod[_NI_*_NZ_+ii] = -qpData->options.QPDUNES_INFTY;
	}

	/* upper bounds */
	/* first nI intervals */
	for ( kk=0; kk<_NI_; ++kk ) {
		/* x part */
		if ( xUpp != 0 ) {
			for ( ii=0; ii<_NX_; ++ii )  zUppMod[kk*_NZ_+ii] = xUpp[kk*_NX_+ii] - mpcProblem->zRef[kk*_NZ_+ii];
		}
		else {
			for ( ii=0; ii<_NX_; ++ii )  zUppMod[kk*_NZ_+ii] = qpData->options.QPDUNES_INFTY;
		}
		/* u part */
		if ( uUpp != 0 ) {
			for ( ii=0; ii<_NU_; ++ii )  zUppMod[kk*_NZ_+_NX_+ii] = uUpp[kk*_NU_+ii] - mpcProblem->zRef[kk*_NZ_+_NX_+ii];
		}
		else {
			for ( ii=0; ii<_NU_; ++ii )  zUppMod[kk*_NZ_+_NX_+ii] = qpData->options.QPDUNES_INFTY;
		}
	}
	/* last interval */
	if ( xUpp != 0 ) {
		for ( ii=0; ii<_NX_; ++ii )  zUppMod[_NI_*_NZ_+ii] = xUpp[_NI_*_NX_+ii] - mpcProblem->zRef[_NI_*_NZ_+ii];
	}
	else {
		for ( ii=0; ii<_NX_; ++ii )  zUppMod[_NI_*_NZ_+ii] = qpData->options.QPDUNES_INFTY;
	}
	
	
	/* setup regular intervals */
	for( kk=0; kk<_NI_; ++kk )
	{
		/* TODO: do this more efficiently, let qp solver know that intervals are identical... factorizations, etc. */
		statusFlag = qpDUNES_setupRegularInterval( qpData, qpData->intervals[kk],
												0, Q, R, S, 0,
												0, A, B, &(cMod[kk*_NX_]),
												&(zLowMod[kk*_NZ_]), &(zUppMod[kk*_NZ_]), 0,0,0,0,
												0,0,0 );
		if (statusFlag != QPDUNES_OK) {
			qpDUNES_printError(qpData, __FILE__, __LINE__, "Setup of interval %d of %d failed. Bailing out.", kk, _NI_ );
			return statusFlag;
		}
	}
	/* set up final interval */
	if ( P != 0 ) {
		statusFlag = qpDUNES_setupFinalInterval( qpData, qpData->intervals[_NI_],
											  P, 0,
											  &(zLowMod[_NI_*_NZ_]), &(zUppMod[_NI_*_NZ_]),
											  0,0,0 );
	}
	else {
		statusFlag = qpDUNES_setupFinalInterval( qpData, qpData->intervals[_NI_],
											  Q, 0,
											  &(zLowMod[_NI_*_NZ_]), &(zUppMod[_NI_*_NZ_]),
											  0,0,0 );
	}
	if (statusFlag != QPDUNES_OK) {
		qpDUNES_printError(qpData, __FILE__, __LINE__, "Setup of interval %d of %d failed. Bailing out.", _NI_, _NI_ );
		return statusFlag;
	}
	
	/* determine local QP solvers and set up auxiliary data */
	statusFlag = qpDUNES_setupAllLocalQPs( qpData, mpcProblem->isLTI );
	if (statusFlag != QPDUNES_OK) {
		qpDUNES_printError(qpData, __FILE__, __LINE__, "Local QP setup failed. Bailing out.", kk, _NI_ );
		return statusFlag;
	}

	return QPDUNES_OK;
}
/*<<< END OF mpcDUNES_initLtiSb_xu */



/* ----------------------------------------------
 * set up an linear time invariant (LTI) MPC problem with simple bounds in z-style
 *
 # >>>>>>                                           */
return_t mpcDUNES_initLtiSb(	mpcProblem_t* const mpcProblem,
							const real_t* const H_,
							const real_t* const P_,
							const real_t* const g_,
							const real_t* const C_,
							const real_t* const c_,
							const real_t* const zLow_,
							const real_t* const zUpp_,
							const real_t* const zRef_
							)
{
	int_t ii,jj,kk;

	return_t statusFlag;

	mpcProblem->isLTI = QPDUNES_TRUE;


	qpData_t* qpData = &(mpcProblem->qpData);
//	int_t nI = qpData->nI;
//	int_t _NZ_ = qpData->_NZ_;
//	int_t _NX_ = qpData->_NX_;

	real_t* cMod = mpcProblem->xnTmp;	/* temporary variables for bound setup from workspace */
	real_t* zLowMod = mpcProblem->zn1Tmp;
	real_t* zUppMod = mpcProblem->zn1Tmp2;

	/** check existence of data */
	if ( !H_ ) {
		qpDUNES_printError( qpData, __FILE__, __LINE__, "H matrix missing" );
		return QPDUNES_ERR_INVALID_ARGUMENT;
	}
	if ( !C_ ) {
		qpDUNES_printError( qpData, __FILE__, __LINE__, "C matrix missing" );
		return QPDUNES_ERR_INVALID_ARGUMENT;
	}


	/** (2) get data */

	/** (2a) get reference */
	if ( zRef_ != 0 ) {
		for ( ii=0; ii<_NI_*_NZ_+_NX_; ++ii ) {
			mpcProblem->zRef[ii] = zRef_[ii];
		}
	}
	else {
		for ( ii=0; ii<_NI_*_NZ_+_NX_; ++ii ) {
			mpcProblem->zRef[ii] = 0.;
		}
	}

	/* translate variables to account for reference zQP := (z-zRef) */
	/** (2b) account for variable substitution in dynamics: c += [A B] * [xRef' uRef']' */
	if ( c_ != 0 )	{
		for ( kk=0; kk<_NI_; ++kk ) {
			for ( ii=0; ii<_NX_; ++ii ) {
				cMod[kk*_NX_+ii] = c_[ii];
			}
		}
	}
	else {
		for ( ii=0; ii<_NI_*_NX_; ++ii )  cMod[ii] = 0.;
	}
	if ( zRef_ != 0 ) {
		for ( kk=0; kk<_NI_; ++kk ) {
			for ( ii=0; ii<_NX_; ++ii ) {
				for ( jj=0; jj<_NZ_; ++jj ) {
					cMod[kk*_NX_+ii] += C_[ii*_NZ_ + jj] * zRef_[kk*_NZ_+jj] - zRef_[(kk+1)*_NZ_+ii];
				}
			}
		}
	}


	/** (2c) account for variable substitution in stage constraints: [z{l,u}] -= [zRef] */
	/* WARNING: currently no support for general constraint matrices!
	 * to support this in general use z{l,u} -= D*zRef */
	if ( zLow_ != 0 ) {
		if ( zRef_ != 0 ) {
			for ( ii=0; ii<_NI_*_NZ_+_NX_; ++ii )  zLowMod[ii] = zLow_[ii] - zRef_[ii];
		}
		else {
			for ( ii=0; ii<_NI_*_NZ_+_NX_; ++ii )  zLowMod[ii] = zLow_[ii];
		}
	}
	else {
		for ( ii=0; ii<_NI_*_NZ_+_NX_; ++ii )  zLowMod[ii] = -qpData->options.QPDUNES_INFTY;
	}

	if ( zUpp_ != 0 ) {
		if ( zRef_ != 0 ) {
			for ( ii=0; ii<_NI_*_NZ_+_NX_; ++ii )  zUppMod[ii] = zUpp_[ii] - zRef_[ii];
		}
		else {
			for ( ii=0; ii<_NI_*_NZ_+_NX_; ++ii )  zUppMod[ii] = zUpp_[ii];
		}
	}
	else {
		for ( ii=0; ii<_NI_*_NZ_+_NX_; ++ii )  zUppMod[ii] = qpData->options.QPDUNES_INFTY;
	}


	/** (3) setup regular intervals */
	for( kk=0; kk<_NI_; ++kk )
	{
		/* TODO: do this more efficiently, let qp solver know that intervals are identical... factorizations, etc. */
		statusFlag = qpDUNES_setupRegularInterval( qpData, qpData->intervals[kk],
												H_, 0, 0, 0, g_,
												C_, 0, 0, &(cMod[kk*_NX_]),
												&(zLowMod[kk*_NZ_]), &(zUppMod[kk*_NZ_]), 0,0,0,0,
												0, 0,0 );
		if (statusFlag != QPDUNES_OK) {
			qpDUNES_printError(qpData, __FILE__, __LINE__, "Setup of interval %d of %d failed. Bailing out.", kk, _NI_ );
			return statusFlag;
		}
	}
	/* set up final interval */
	statusFlag = qpDUNES_setupFinalInterval( qpData, qpData->intervals[_NI_],
										  P_, g_,
										  &(zLowMod[_NI_*_NZ_]), &(zUppMod[_NI_*_NZ_]),
										  0, 0,0 );
	if (statusFlag != QPDUNES_OK) {
		qpDUNES_printError(qpData, __FILE__, __LINE__, "Setup of interval %d of %d failed. Bailing out.", _NI_, _NI_ );
		return statusFlag;
	}

	/* determine local QP solvers and set up auxiliary data */
	statusFlag = qpDUNES_setupAllLocalQPs( qpData, mpcProblem->isLTI );
	if (statusFlag != QPDUNES_OK) {
		qpDUNES_printError(qpData, __FILE__, __LINE__, "Local QP setup failed. Bailing out.", kk, _NI_ );
		return statusFlag;
	}

	return QPDUNES_OK;
}
/*<<< END OF mpcDUNES_initLtiSb */




/* ----------------------------------------------
 * set up an linear time varying (LTV) MPC problem with simple bounds
 *
 # >>>>>>                                           */
return_t mpcDUNES_initLtvSb(	mpcProblem_t* const mpcProblem,
								const real_t* const H_,
								const real_t* const g_,
								const real_t* const C_,
								const real_t* const c_,
								const real_t* const zLow_,
								const real_t* const zUpp_,
								const real_t* const zRef_
								)
{
	int_t kk, ii, jj;

//	static const boolean_t isLTI = QPDUNES_FALSE;
	mpcProblem->isLTI = QPDUNES_FALSE;

	qpData_t* qpData = &(mpcProblem->qpData);
//	int_t nI = qpData->nI;
//	int_t _NZ_ = qpData->_NZ_;
//	int_t _NX_ = qpData->_NX_;

	real_t* cMod = mpcProblem->xnTmp;	/* temporary variables for bound setup from workspace */
	real_t* zLowMod = mpcProblem->zn1Tmp;
	real_t* zUppMod = mpcProblem->zn1Tmp2;

	/** (1) check existence of data */
	if ( !H_ ) {
		qpDUNES_printError( qpData, __FILE__, __LINE__, "H matrix missing" );
		return QPDUNES_ERR_INVALID_ARGUMENT;
	}
	if ( !C_ ) {
		qpDUNES_printError( qpData, __FILE__, __LINE__, "C matrix missing" );
		return QPDUNES_ERR_INVALID_ARGUMENT;
	}


	/** (2) get data */

	/** (2a) get reference */
	if ( zRef_ != 0 ) {
		for ( ii=0; ii<_NI_*_NZ_+_NX_; ++ii ) {
			mpcProblem->zRef[ii] = zRef_[ii];
		}
	}
	else {
		for ( ii=0; ii<_NI_*_NZ_+_NX_; ++ii ) {
			mpcProblem->zRef[ii] = 0.;
		}
	}

	/* translate variables to account for reference zQP := (z-zRef) */
	/** (2b) account for variable substitution in dynamics: c += [A B] * [xRef' uRef']' */
	if ( c_ != 0 ) {
		for ( ii=0; ii<_NI_*_NX_; ++ii )   cMod[ii] = c_[ii];
	}
	else {
		for ( ii=0; ii<_NI_*_NX_; ++ii )  cMod[ii] = 0.;
	}
	if ( zRef_ != 0 ) {
		for ( kk=0; kk<_NI_; ++kk ) {
			for ( ii=0; ii<_NX_; ++ii ) {
				for ( jj=0; jj<_NZ_; ++jj ) {
					/*                   blocks       rows    cols                  */
					cMod[kk*_NX_+ii] += C_[kk*(_NX_*_NZ_) + ii*_NZ_ + jj] * zRef_[kk*_NZ_+jj] - zRef_[(kk+1)*_NZ_+ii];
				}
			}
		}
	}


	/** (2c) account for variable substitution in stage constraints: [z{l,u}] -= [zRef] */
	/* WARNING: currently no support for general constraint matrices!
	 * to support this in general use z{l,u} -= D*zRef */
	if ( zLow_ != 0 ) {
		if ( zRef_ != 0 ) {
			for ( ii=0; ii<_NI_*_NZ_+_NX_; ++ii )  zLowMod[ii] = zLow_[ii] - zRef_[ii];
		}
		else {
			for ( ii=0; ii<_NI_*_NZ_+_NX_; ++ii )  zLowMod[ii] = zLow_[ii];
//			zLowMod = zLow_;
		}
	}
	else {
		for ( ii=0; ii<_NI_*_NZ_+_NX_; ++ii )  zLowMod[ii] = -qpData->options.QPDUNES_INFTY;
	}

	if ( zUpp_ != 0 ) {
		if ( zRef_ != 0 ) {
			for ( ii=0; ii<_NI_*_NZ_+_NX_; ++ii )  zUppMod[ii] = zUpp_[ii] - zRef_[ii];
		}
		else {
			for ( ii=0; ii<_NI_*_NZ_+_NX_; ++ii )  zUppMod[ii] = zUpp_[ii];
//			zUppMod = zUpp_;
		}
	}
	else {
		for ( ii=0; ii<_NI_*_NZ_+_NX_; ++ii )  zUppMod[ii] = qpData->options.QPDUNES_INFTY;
	}


	/** (3) setup intervals */
	/** (3a) setup regular intervals */
	for( kk=0; kk<_NI_; ++kk )
	{
		if (g_ != 0) {
			qpDUNES_setupRegularInterval( qpData, qpData->intervals[kk],
									   &(H_[kk*_NZ_*_NZ_]), 0, 0, 0, &(g_[kk*_NZ_]),
									   &(C_[kk*_NX_*_NZ_]), 0, 0, &(cMod[kk*_NX_]),
									   &(zLowMod[kk*_NZ_]), &(zUppMod[kk*_NZ_]), 0,0,0,0,
									   0, 0,0 );
		}
		else {
			qpDUNES_setupRegularInterval( qpData, qpData->intervals[kk],
									   &(H_[kk*_NZ_*_NZ_]), 0, 0, 0, 0,
									   &(C_[kk*_NX_*_NZ_]), 0, 0, &(cMod[kk*_NX_]),
									   &(zLowMod[kk*_NZ_]), &(zUppMod[kk*_NZ_]), 0,0,0,0,
									   0, 0,0 );
		}
	}
	/** (3b) set up final interval */
	if (g_ != 0)	qpDUNES_setupFinalInterval( qpData, qpData->intervals[_NI_], &(H_[_NI_*_NZ_*_NZ_]), &(g_[_NI_*_NZ_]), &(zLowMod[_NI_*_NZ_]), &(zUppMod[_NI_*_NZ_]), 0, 0,0 );
	else			qpDUNES_setupFinalInterval( qpData, qpData->intervals[_NI_], &(H_[_NI_*_NZ_*_NZ_]), 0, &(zLowMod[_NI_*_NZ_]), &(zUppMod[_NI_*_NZ_]), 0, 0,0 );



	/** (4) determine local QP solvers and set up auxiliary data */
	qpDUNES_setupAllLocalQPs( qpData, mpcProblem->isLTI );


	return QPDUNES_OK;
}
/*<<< END OF mpcDUNES_initLtvSb */



/* ----------------------------------------------
 * set up an linear time varying (LTV) MPC problem
 * with affine constraints
 *
 * min sum ( (z[k]-zRef[k])' * H[k] * (z[k]-zRef[k])  +  g[k]' * z[k] )
 * s.t.   x[k+1]  = C[k] * z[k] + c[k]
 *       dLow[k] <= D[k] * z[k] <= dUpp[k]
 *
 # >>>>>>                                           */
return_t mpcDUNES_initLtv(	mpcProblem_t* const mpcProblem,
							const real_t* const H_,
							const real_t* const g_,
							const real_t* const C_,
							const real_t* const c_,
							const real_t* const zLow_,
							const real_t* const zUpp_,
							const real_t* const D_,
							const real_t* const dLow_,
							const real_t* const dUpp_,
							const real_t* const zRef_
							)
{
	int_t kk, ii, jj;

	mpcProblem->isLTI = QPDUNES_FALSE;

	qpData_t* qpData = &(mpcProblem->qpData);
//	int_t nI = qpData->nI;
//	int_t _NZ_ = qpData->_NZ_;
//	int_t _NX_ = qpData->_NX_;
	int_t nVttl = _NI_*_NZ_ + _NX_;

	int_t nDoffset;

//	real_t* cMod = mpcProblem->xnTmp;	/* temporary variables for bound setup from workspace */
//	real_t* zLowMod = mpcProblem->zn1Tmp;
//	real_t* zUppMod = mpcProblem->zn1Tmp2;
	real_t* gMod = mpcProblem->zn1Tmp;
	real_t* tmp_HzRef = mpcProblem->zn1Tmp2;
	real_t* qpObjConst = mpcProblem->qpObjConst;

	/** (1) check existence of data */
	if ( !H_ ) {
		qpDUNES_printError( qpData, __FILE__, __LINE__, "H matrix missing" );
		return QPDUNES_ERR_INVALID_ARGUMENT;
	}
	if ( !C_ ) {
		qpDUNES_printError( qpData, __FILE__, __LINE__, "C matrix missing" );
		return QPDUNES_ERR_INVALID_ARGUMENT;
	}


	/** (2) get data */

	/** (2a) get reference */
	if ( zRef_ != 0 ) {
		for ( ii=0; ii<_NI_*_NZ_+_NX_; ++ii ) {
			mpcProblem->zRef[ii] = zRef_[ii];
		}
	}
	else {
		for ( ii=0; ii<_NI_*_NZ_+_NX_; ++ii ) {
			mpcProblem->zRef[ii] = 0.;
		}
	}

	/* translate objective to account for reference
	 * .5 * (z-zRef)' * H * (z-zRef)  =  .5 * z' * H * z  +  (zRef' * H) * z  +  (zRef * H * zRef) */
	/* get H*zRef */
	if ( zRef_ != 0 ) {
		for( kk=0; kk<_NI_; ++kk ) {
			for( ii=0; ii<_NZ_; ++ii ) {
				tmp_HzRef[kk*_NZ_+ii] = 0.;
				for( jj = 0; jj < _NZ_; ++jj ) {
					/*                        blocks     rows    cols        blocks  rows */
					tmp_HzRef[kk*_NZ_+ii] += H_[kk*_NZ_*_NZ_ + ii*_NZ_ + jj] * zRef_[kk*_NZ_ + jj];
				}
			}
		}
		for( ii=0; ii<_NZ_; ++ii ) {
			tmp_HzRef[_NI_*_NZ_+ii] = 0.;
			for( jj=0; jj<_NZ_; ++jj ) {
				/*                        blocks     rows    cols        blocks  rows */
				tmp_HzRef[_NI_*_NZ_+ii] += H_[_NI_*_NZ_*_NZ_ + ii*_NX_ + jj] * zRef_[_NI_*_NZ_ + jj];
			}
		}
	}
	else {
		tmp_HzRef = 0;
	}
	/* set gMod */
	if ( g_ != 0 ) {
		if (tmp_HzRef != 0)	{
			for ( ii=0; ii<nVttl; ++ii )   gMod[ii] = g_[ii] + tmp_HzRef[ii];
		}
		else {
			for ( ii=0; ii<nVttl; ++ii )   gMod[ii] = g_[ii];
		}
	}
	else {
		if (tmp_HzRef != 0)	{
			for ( ii=0; ii<nVttl; ++ii )   gMod[ii] = tmp_HzRef[ii];
		}
		else	{
			gMod = 0;
		}
	}
	/* set constant QP part */
	*qpObjConst = 0.;
	if (zRef_ != 0)	{
		for ( ii=0; ii<nVttl; ++ii )   *qpObjConst += tmp_HzRef[ii] * zRef_[ii];
	}


	/** (3) setup intervals */
	/** (3a) setup regular intervals */
	nDoffset = 0;
	for( kk=0; kk<_NI_; ++kk )
	{
		if (qpData->intervals[kk]->nD > 0)	{
			qpDUNES_setupRegularInterval( qpData, qpData->intervals[kk],
										   offsetArray(H_, kk*_NZ_*_NZ_), 0, 0, 0, offsetArray(gMod, kk*_NZ_),
										   offsetArray(C_, kk*_NX_*_NZ_), 0, 0, offsetArray(c_, kk*_NX_),
										   offsetArray(zLow_, kk*_NZ_), offsetArray(zUpp_, kk*_NZ_), 0, 0, 0, 0,
										   offsetArray(D_, nDoffset*_NZ_), offsetArray(dLow_, nDoffset), offsetArray(dUpp_, nDoffset) );
			nDoffset += qpData->intervals[kk]->nD;
		}
		else {
			qpDUNES_setupRegularInterval( qpData, qpData->intervals[kk],
										   offsetArray(H_, kk*_NZ_*_NZ_), 0, 0, 0, offsetArray(gMod, kk*_NZ_),
										   offsetArray(C_, kk*_NX_*_NZ_), 0, 0, offsetArray(c_, kk*_NX_),
										   offsetArray(zLow_, kk*_NZ_), offsetArray(zUpp_, kk*_NZ_), 0, 0, 0, 0,
										   0, 0,0 );
		}
	}
	/** (3b) set up final interval */
	if (qpData->intervals[_NI_]->nD > 0)	{
		qpDUNES_setupFinalInterval( qpData, qpData->intervals[_NI_],
									 offsetArray(H_, _NI_*_NZ_*_NZ_), offsetArray(gMod, _NI_*_NZ_),
									 offsetArray(zLow_, _NI_*_NZ_), offsetArray(zUpp_, _NI_*_NZ_),
									 offsetArray(D_, nDoffset*_NZ_), offsetArray(dLow_, nDoffset), offsetArray(dUpp_, nDoffset) );
	}
	else {
		qpDUNES_setupFinalInterval( qpData, qpData->intervals[_NI_],
									 offsetArray(H_, _NI_*_NZ_*_NZ_), offsetArray(gMod, _NI_*_NZ_),
									 offsetArray(zLow_, _NI_*_NZ_), offsetArray(zUpp_, _NI_*_NZ_),
									 0, 0,0 );
	}


	/** (4) determine local QP solvers and set up auxiliary data */
	qpDUNES_setupAllLocalQPs( qpData, mpcProblem->isLTI );


	return QPDUNES_OK;
}
/*<<< END OF mpcDUNES_initLtv */




/* ----------------------------------------------
 * solve an MPC problem for a given initial value
 * 
 # >>>>>>                                           */
return_t mpcDUNES_solve(	mpcProblem_t* const mpcProblem,
							const real_t* const x0
							)
{
	/* TODO: make shift in data optional */

	int_t kk, ii;
	int_t nQpoasesIter;
	interval_t* interval;
	qpData_t* qpData = &(mpcProblem->qpData);

	/* (0) save bounds on first interval before initial value embedding for recovery afterwards */
	qpDUNES_copyArray( mpcProblem->z0LowOrig, qpData->intervals[0]->zLow.data, _NZ_ );
	qpDUNES_copyArray( mpcProblem->z0UppOrig, qpData->intervals[0]->zUpp.data, _NZ_ );

	/* (1) change bounds to embed initial value x0 */
	if ( x0 != 0 ) {
		for ( ii=0; ii<_NX_; ++ii ) {
//			qpData->intervals[0]->zLow.data[ii] = x0[ii] - mpcProblem->zRef[0+ii];
//			qpData->intervals[0]->zUpp.data[ii] = x0[ii] - mpcProblem->zRef[0+ii];
			qpData->intervals[0]->zLow.data[ii] = x0[ii];
			qpData->intervals[0]->zUpp.data[ii] = x0[ii];
		}
	}

	/* (2) resolve first stage QP */
	/* update stage QPs */
	switch ( qpData->intervals[0]->qpSolverSpecification ) {
		case QPDUNES_STAGE_QP_SOLVER_CLIPPING:
			qpDUNES_setupClippingSolver( qpData, qpData->intervals[0], QPDUNES_FALSE );
			break;

		case QPDUNES_STAGE_QP_SOLVER_QPOASES:
//			qpDUNES_updateQpoases( qpData, interval, QPDUNES_FALSE, g_changed, QPDUNES_TRUE, QPDUNES_TRUE, QPDUNES_FALSE, QPDUNES_FALSE, QPDUNES_FALSE );
			qpDUNES_updateQpoases( qpData, qpData->intervals[0], QPDUNES_FALSE, QPDUNES_FALSE, QPDUNES_TRUE, QPDUNES_TRUE, QPDUNES_FALSE, QPDUNES_FALSE, QPDUNES_FALSE, &nQpoasesIter );
			break;

		default:
			qpDUNES_printError( qpData, __FILE__, __LINE__, "Unknown stage QP solver selected." );
			return QPDUNES_ERR_INVALID_ARGUMENT;

	}



	/* (3) solve QP */
	mpcProblem->exitFlag = qpDUNES_solve( qpData );
	

	/* (4) recover MPC solution */
	/*  - primal: recover zMPC := zQP + zRef */
	for ( kk=0; kk<_NI_; ++kk ) { /* regular intervals */
		for ( ii=0; ii<_NX_; ++ii ) {
			mpcProblem->xOpt[kk*_NX_+ii] = qpData->intervals[kk]->z.data[ii];	// + mpcProblem->zRef[kk*_NZ_+ii];
		}
		for ( ii=0; ii<_NU_; ++ii ) {
			mpcProblem->uOpt[kk*_NU_+ii] = qpData->intervals[kk]->z.data[_NX_+ii]; // + mpcProblem->zRef[kk*_NZ_+_NX_+ii];
		}
	}
	for ( ii=0; ii<_NX_; ++ii ) {		/* last interval */
		mpcProblem->xOpt[_NI_*_NX_+ii] = qpData->intervals[_NI_]->z.data[ii]; // + mpcProblem->zRef[_NI_*_NZ_+ii];
	}
	
	/*  - dual */
	for ( ii=0; ii<_NI_*_NX_; ++ii ) {
		mpcProblem->lambdaOpt[ii] = qpData->lambda.data[ii];
	}
	
	/*  - objective value */
	mpcProblem->optObjVal = qpDUNES_computeObjectiveValue( qpData );
	

	/* (5) prepare next QP solution */
	/*  - reset original variable bounds (important in LTI case with shift, might be redundant if data is updated) */
	qpDUNES_updateIntervalData( qpData,
								mpcProblem->qpData.intervals[0],
								0, 0, 0, 0,
								mpcProblem->z0LowOrig,mpcProblem->z0UppOrig,
								0, 0,0, 0 );
	/*  - shift variables */
	qpDUNES_shiftLambda( qpData );			/* shift multipliers */
	/* shift intervals (particulary important when using qpOASES for underlying local QPs) */
	if (mpcProblem->isLTI == QPDUNES_TRUE)	{
		qpDUNES_shiftIntervalsLTI( qpData );		/* do not force hessian refactorization */
	}
	else {
		qpDUNES_shiftIntervals( qpData );
	}


	/* re-setup stage solvers to account for lambda shift.
	 * TODO: this is a bit inefficient, but makes sure the local clipping QPs are initialized correctly (qStep) for the next iteration
	 * if lambda and intervals are always shifted alongside, if would be sufficient to only update first and second-but-last interval... */
	for( kk=0; kk<_NI_+1; ++kk ) {
		interval = qpData->intervals[kk];
		switch ( interval->qpSolverSpecification ) {
			case QPDUNES_STAGE_QP_SOLVER_CLIPPING:
				qpDUNES_setupClippingSolver( qpData, interval, QPDUNES_FALSE );
				break;

			case QPDUNES_STAGE_QP_SOLVER_QPOASES:
				qpDUNES_updateQpoases( qpData, interval,
										QPDUNES_FALSE,
										QPDUNES_TRUE,
										QPDUNES_FALSE,
										QPDUNES_FALSE,
										QPDUNES_FALSE,
										QPDUNES_FALSE,
										QPDUNES_FALSE,
										&nQpoasesIter );
				break;

			default:
				qpDUNES_printError( qpData, __FILE__, __LINE__, "MPC shift failed. Unknown stage QP solver selected." );
				return QPDUNES_ERR_INVALID_ARGUMENT;

		}
	}

	return mpcProblem->exitFlag;
}
/*<<< END OF mpcDUNES_solve */




/*
 *	end of file
 */
