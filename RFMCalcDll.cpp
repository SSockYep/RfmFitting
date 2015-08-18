// RFMCalcDll.cpp : Defines the exported functions for the DLL application.
//


#include "RFMCalcDll.h"
#include "Mat.h"

//#include "GaussRand.h"

#include <vector>
#include <iostream>
using namespace std;
// This is an example of an exported variable


//Calc A or C val
double calcHeadVal( double* pHeadCoef, double dP, double dL, double dH )
{
	double dRes = pHeadCoef[0] + pHeadCoef[1]*dP + pHeadCoef[2]*dL + pHeadCoef[3]*dH + pHeadCoef[4]*dP*dL + pHeadCoef[5]*dP*dH + pHeadCoef[6]*dL*dH + pHeadCoef[7]*dP*dP + \
		pHeadCoef[8]*dL*dL + pHeadCoef[9]*dH*dH + pHeadCoef[10]*dP*dL*dH + pHeadCoef[11]*dP*dP*dP + pHeadCoef[12]*dP*dL*dL + pHeadCoef[13]*dP*dH*dH + \
		pHeadCoef[14]*dP*dP*dL + pHeadCoef[15]*dL*dL*dL + pHeadCoef[16]*dL*dH*dH + pHeadCoef[17]*dP*dP*dH + pHeadCoef[18]*dL*dL*dH + pHeadCoef[19]*dH*dH*dH;

	return dRes;
}
//Calc B or D val
double calcTailVal( double* pTailCoef, double dP, double dL, double dH )
{
	double dRes = 1 + pTailCoef[0]*dP + pTailCoef[1]*dL + pTailCoef[2]*dH + pTailCoef[3]*dP*dL + pTailCoef[4]*dP*dH + pTailCoef[5]*dL*dH + pTailCoef[6]*dP*dP + \
		pTailCoef[7]*dL*dL + pTailCoef[8]*dH*dH + pTailCoef[9]*dP*dL*dH + pTailCoef[10]*dP*dP*dP + pTailCoef[11]*dP*dL*dL + pTailCoef[12]*dP*dH*dH + \
		pTailCoef[13]*dP*dP*dL + pTailCoef[14]*dL*dL*dL + pTailCoef[15]*dL*dH*dH + pTailCoef[16]*dP*dP*dH + pTailCoef[17]*dL*dL*dH + pTailCoef[18]*dH*dH*dH;
	return dRes;
}
//Coefs in normal equation,=39

int getFullRPCACoefs( double* pCoefs, double dP, double dL, double dH, double dPara)
{
	pCoefs[0] = 1;
	pCoefs[1] = dP;
	pCoefs[2] = dL;
	pCoefs[3] = dH;
	pCoefs[4] = dP*dL;
	pCoefs[5] = dP*dH;
	pCoefs[6] = dL*dH;
	pCoefs[7] = dP*dP;
	pCoefs[8] = dL*dL;
	pCoefs[9] = dH*dH;
	pCoefs[10] = dP*dL*dH;
	pCoefs[11] = dP*dP*dP;
	pCoefs[12] = dP*dL*dL;
	pCoefs[13] = dP*dH*dH;
	pCoefs[14] = dP*dP*dL;
	pCoefs[15] = dL*dL*dL;
	pCoefs[16] = dL*dH*dH;
	pCoefs[17] = dP*dP*dH;
	pCoefs[18] = dL*dL*dH;
	pCoefs[19] = dH*dH*dH;

	for ( int i = 0; i < 19; ++i )
	{
		pCoefs[i+20] = -dPara * pCoefs[i+1];
	}
	return 1;
}
RFMCALC_EXTERN_C int rfmTest( double *pImgPts, int nImgPts, double *pGeoPts, int nGeoPts, \
						   RPCparas *pRpcs, RFMCalcPara para, RFMType rfmType )
{
	//CCassiniPDSData clsCasPDS;
	//
	//clsCasPDS.ReadBIDRPDS( "H:\\Cassini\\T16\\1.IMG" );
	//CMat mat;
	////MatrixXd matA1(4,4);
	////matA1 << 94.61, -22.11, -11.45, -6.95,
	////	     -22.11, 70.51, -6.95, -8.42,
	////		 -11.45, -6.95, 96.09, -20.21,
	////		 -6.96, -8.42, -20.21, 66.63;

	////VectorXd vecL1( 4 );
	////vecL1 << -43.52, 178.81, -120.11, -30.07;
	////MatrixXd matW1(4,4);
	////matW1.setIdentity();

	////mat.SetMat( matA1, vecL1, matW1 );
	////VectorXd vecRes1;
	////mat.RobustSolve( vecRes1 );

	////int k1 = 1;

	////////////////////////////////////////////////////////////////////////////

	////MatrixXd matA2(4,4);
	////matA2 << 2.122, -0.0269, 0.0048, 0.0019,
	////	     -0.0269, 1.9949, -0.0432, -0.0168,
	////		 0.0048, -0.0432, 0.8892, 1.0178,
	////		 0.0019, -0.0168, 1.0178, 1.1656;

	////VectorXd vecL2( 4 );
	////vecL2 << 2.090, 3.768, 6.656, 7.684;
	////MatrixXd matW2(4,4);
	////matW2.setIdentity();

	////mat.SetMat( matA2, vecL2, matW2 );
	////VectorXd vecRes2;
	////mat.RobustSolve( vecRes2 );

	////int k2 = 1;
	////////////////////////////////////////////////////////////////////////////

	//MatrixXd matCur;//(10,5);

	//matCur.resize( 10, 5 );


	//matCur <<  2.0 ,  -5.0 ,  1.0 ,  1.0 ,  -9.5
	//	,  -2.0 ,  4.0 ,  1.0 ,  -1.05 ,  8.5
	//	,  -2.0 ,  1.0 ,  1.0 ,  -1.0 ,  2.4
	//	,  -1.0 ,  2.5 ,  4.0 ,  -0.5 ,  7.0
	//	,  -1.0 ,  3.2 ,  4.0 ,  -.5 ,  8.4
	//	,  1.0 ,  1.0 ,  -3.0 ,  .4 ,  .49
	//	,  3.0 ,  7.0 ,  -3.0 ,  1.5 ,  12.7
	//	,  5.0 ,  -1.0 ,  -2.0 ,  2.5 ,  -3.0
	//	,  4.0 ,  2.0 ,  -2.0 ,  2.01 ,  3.0
	//	,  4.0 ,  3.0 ,  -2.0 ,  2.0 ,  5.0;

	//CGaussRand clsRand( 0, 1 );
	//VectorXd vecX(matCur.cols());
	//for ( int i = 0; i < matCur.cols(); ++i )
	//{
	//	vecX(i) = 1 ;//+ clsRand.GetRand();
	//}
	//VectorXd vecL = matCur*vecX;
	//vecL(0) = vecL(0) + (-1.08906);
	//vecL(1) = vecL(1) + (0.0325575);
	//vecL(2) = vecL(2) + (0.552527);
	//vecL(3) = vecL(3) + ( 1.10061);
	//vecL(4) = vecL(4) + (1.54421);

	//vecL(5) = vecL(5) + ( 0.0859311);
	//vecL(6) = vecL(6) + ( -1.49159);
	//vecL(7) = vecL(7) + (-0.742302);
	//vecL(8) = vecL(8) + (-1.06158);
	//vecL(9) = vecL(9) + (2.35046 );

	//MatrixXd matVal = matCur.transpose() * matCur;

	//double dRes = mat.GetConditionNum( matVal );


	////double dCoef = mat.estimateKParaUseL_Curve( matCur, vecL );
	//MatrixXd matW(matCur.rows(),matCur.rows());
	//matW.setIdentity();

	//mat.SetMat( matCur, vecL, matW );
	//VectorXd vecXRes, vecResids;
	//mat.RobustSolve( vecXRes, vecResids, para.nL_Curve_Points, para.dICCVThres );
	return 0;
}

RFMCALC_EXTERN_C int rfmCalc( double *pImgPts, int nImgPts, double *pGeoPts, int nGeoPts, \
							 RPCparas *pRpcs, RFMCalcPara &para, double* pResErr, RFMType rfmType )
{
	double first_sqrt_mean;
	cout << "rfmCalc" << endl;
	if(nImgPts != nGeoPts)
	{
		printf("error occurred\n");
		return -1;
	}
	MatrixXd M_mat(nImgPts, 39);
	MatrixXd M_mat_cpy(nImgPts, 39);
	MatrixXd N_mat(nImgPts, 39);
	VectorXd row_M(nImgPts), col_N(nImgPts);

	M_mat.setZero();
	N_mat.setZero();
	row_M.setZero();
	col_N.setZero();

	cout << "init Row and Col Values" << endl;
	//init Row and Col Values
	for(int i=0; i<nImgPts; i++)
	{
		double cX, cY, cZ;
		double iX, iY;
		cX = (pGeoPts[3*i] - pRpcs->dObjXOff) / pRpcs->dObjXScale;
		cY = (pGeoPts[3*i+1] - pRpcs->dObjYOff) / pRpcs->dObjYScale ;
		cZ = (pGeoPts[3*i+2] - pRpcs->dObjZOff) / pRpcs->dObjZScale;
		iX = (pImgPts[2*i] - pRpcs->dXOff) / pRpcs->dXScale;
		iY = (pImgPts[2*i] - pRpcs->dYOff) / pRpcs->dYScale;

		row_M(i) = iX;
		col_N(i) = iY;

		setFormerMat(M_mat, i, cX, cY, cZ);
		setFormerMat(N_mat, i, cX, cY, cZ);

		for(int j=20; j<39; j++)
		{
			M_mat(i,j) = (-1.0*iX) * M_mat(i, j-19);
			N_mat(i,j) = (-1.0*iY) * N_mat(i, j-19);
		}
	}

	VectorXd r_coeff(39);
	VectorXd c_coeff(39);
	//gsl_vector* r_coeff = gsl_vector_alloc(39);
	//gsl_vector* c_coeff = gsl_vector_alloc(39);
	MatrixXd workMat(M_mat.cols(), M_mat.rows());
	MatrixXd weight_M (M_mat.rows() , M_mat.rows());
	MatrixXd weight_N (M_mat.rows() , M_mat.rows());
	MatrixXd not_squared_weight_M (M_mat.rows() , M_mat.rows());
	MatrixXd M_weighted(M_mat.cols(), M_mat.cols());
	MatrixXd N_weighted(N_mat.cols(), N_mat.cols());
	VectorXd row_M_weighted(M_mat.cols());
	VectorXd col_N_weighted(N_mat.cols());
	VectorXd error_vector(M_mat.rows());
	weight_M.setZero();
	weight_N.setZero();
	
	cout << "first calculate" << endl;

	//calculate weighted M matrix (N x N), row_M vector

	for(int i = 0; i < weight_M.rows(); i++)
	{
		weight_M(i,i) = 1.0;
		weight_N(i,i) = 1.0;
		//gsl_matrix_set(weight_M,i,i,1.0);
		//gsl_matrix_set(weight_N,i,i,1.0);
	}
	workMat.setZero();
	error_vector.setZero();

	workMat = M_mat.transpose() * weight_M;
	row_M_weighted = workMat * row_M;
	M_weighted = workMat * M_mat;

	//calculate weighted N matrix (N x N)
	workMat.setZero();
	workMat = N_mat.transpose() * weight_N;
	col_N_weighted = workMat * col_N;
	N_weighted = workMat * N_mat;

	r_coeff = solveMat(M_weighted, row_M_weighted);
	c_coeff = solveMat(N_weighted, col_N_weighted);

	error_vector = weight_M * row_M;
	M_mat_cpy = weight_M * M_mat;
	error_vector = (M_mat_cpy * r_coeff) - error_vector;

	double sum = 0;
	for(int i = 0; i < nImgPts; i++)
		sum += pow(error_vector(i), 2.0);
	first_sqrt_mean = sqrt(sum/nImgPts);
	/*
	 * a: samp_num_coeff
	 * b: samp_den_coeff
	 * c: line_num_coeff
	 * d: line_den_coeff
	 * */
	for(int i=0; i<20; i++)
		pRpcs->c[i] = r_coeff(i);
	for(int i=0; i<20; i++)
	{
		if(i==0) pRpcs->d[i] = 1.0;
		else pRpcs->d[i] = r_coeff(i+19);
	}
	for(int i=0; i<20; i++)
			pRpcs->a[i] = c_coeff(i);
	for(int i=0; i<20; i++)
	{
		if(i==0) pRpcs->b[i] = 1.0;
		else pRpcs->b[i] = c_coeff(i+19);
	}

	workMat.setZero();
	weight_M.setZero();
	weight_N.setZero();
	M_weighted.setZero();
	N_weighted.setZero();
	row_M_weighted.setZero();
	col_N_weighted.setZero();

	r_coeff.setZero();
	c_coeff.setZero();

	//iterative Routine
	int iter_n = 1;
	cout << "iteration start" << endl;

	int count = 0;
	do
	{
		weight_M.setZero();
		weight_N.setZero();

		M_weighted.setZero();
		N_weighted.setZero();
		row_M_weighted.setZero();
		col_N_weighted.setZero();

		setupWeightMatrix(weight_M, pRpcs->d, nGeoPts, pGeoPts);
		setupWeightMatrix(weight_N, pRpcs->b, nGeoPts, pGeoPts);
		not_squared_weight_M = weight_M;

		weight_M = weight_M * weight_M;
		weight_N = weight_N * weight_N;
		// (NxM) * (MxM) => (NxM)

		workMat.setZero();
		error_vector.setZero();

		workMat = M_mat.transpose() * weight_M;
		row_M_weighted = workMat * row_M;
		M_weighted = workMat * M_mat;

		workMat.setZero();
		workMat = N_mat.transpose() * weight_N;
		col_N_weighted = workMat * col_N;
		N_weighted = workMat * N_mat;

		r_coeff = solveMat(M_weighted, row_M_weighted);
		c_coeff = solveMat(N_weighted, col_N_weighted);

		error_vector = weight_M * row_M;
		M_mat_cpy = weight_M * M_mat;
		error_vector = (M_mat_cpy * r_coeff) - error_vector;

		double sum = 0;
		double sqrt_mean;
		for(int i = 0; i < nImgPts; i++)
			sum += pow(error_vector(i), 2.0);
		sqrt_mean = sqrt(sum/nImgPts);

		for(int i=0; i<20; i++)
			pRpcs->c[i] = r_coeff(i);
		for(int i=0; i<20; i++)
		{
			if(i==0) pRpcs->d[i] = 1.0;
			else pRpcs->d[i] = r_coeff(i+19);
		}
		for(int i=0; i<20; i++)
				pRpcs->a[i] = c_coeff(i);
		for(int i=0; i<20; i++)
		{
			if(i==0) pRpcs->b[i] = 1.0;
			else pRpcs->b[i] = c_coeff(i+19);
		}

		printf("%d: %+.5E\n", iter_n, sqrt_mean);
		if(sqrt_mean < para.dICCVThres)
			count++;
		else
			count = 0;
		if(count >= 5)
			break;
		/* CHECK END CONDITION */
		//printf("%dth iteration done\n", iter_n);
	} while(++iter_n <= para.nMaxItera);

	return 1;
}
RFMCALC_EXTERN_C int rfmCalcSimple( double *pImgPts, int nImgPts, double *pGeoPts, int nGeoPts, \
							 RPCparas *pRpcs, RFMCalcPara &para, double* pResErr, RFMType rfmType )
{
	if ( nImgPts != nGeoPts )
	{
		return false;
	}

	switch ( rfmType )
	{
	case Three_Order_Coefs:
		if ( nImgPts < 39 ) return -1;
		break;
	}

	vector<double> vecImgX, vecImgY, vecObjX, vecObjY, vecObjZ;
	for( int i = 0; i < nImgPts; ++i )
	{
		vecImgX.push_back( ( pImgPts[2*i] - pRpcs->dXOff )/pRpcs->dXScale );
		vecImgY.push_back( ( pImgPts[2*i+1] - pRpcs->dYOff )/pRpcs->dYScale );

		vecObjX.push_back( ( pGeoPts[3*i] - pRpcs->dObjXOff )/pRpcs->dObjXScale );
		vecObjY.push_back( ( pGeoPts[3*i+1] - pRpcs->dObjYOff )/pRpcs->dObjYScale );
		vecObjZ.push_back( ( pGeoPts[3*i+2] - pRpcs->dObjZOff )/pRpcs->dObjZScale );

		//vecImgX.push_back(  pImgPts[2*i]  );
		//vecImgY.push_back(  pImgPts[2*i+1]  );

		//vecObjX.push_back(  pGeoPts[3*i]  );
		//vecObjY.push_back(  pGeoPts[3*i+1]  );
		//vecObjZ.push_back(  pGeoPts[3*i+2]  );
	}
	int nCoefNum = 78;
	MatrixXd matA( nImgPts*2, nCoefNum );

	VectorXd vecL( nImgPts*2 );

	if ( rfmType == Three_Order_Coefs )
	{
		matA.setZero();

		vecL.setZero();

		for ( int i = 0; i < nImgPts; ++i )
		{
			double dCoefsHead[39], dCoefsTail[39];
			getFullRPCACoefs( dCoefsHead, vecObjX[i], vecObjY[i], vecObjZ[i], vecImgX[i] );
			getFullRPCACoefs( dCoefsTail, vecObjX[i], vecObjY[i], vecObjZ[i], vecImgY[i] );
			for ( int j = 0; j < 39; ++j )
			{
				matA(i,j) = dCoefsHead[j];
				matA(i+nImgPts,j+39) = dCoefsTail[j];
			}
			vecL(i) = vecImgX[i];
			vecL(i+nImgPts) = vecImgY[i];
		}

		VectorXd vecPreX(78);
		vecPreX.setZero();

		int nItera = 0;
		VectorXd vecRes, vecResids, vecResidVals;

		CMat matCalc;
		if ( 1 != matCalc.RobustSolveSimple( matA, vecL, vecRes, vecResids, para.nMaxItera, para.nL_Curve_Points, para.dICCVThres ) )
		{
			return -1;
		}

		for ( int i = 0; i < 20; ++i )
		{
			pRpcs->a[i] = vecRes(i);
			pRpcs->c[i] = vecRes(i+39);
			if ( i == 0 )
			{
				pRpcs->b[i] = 1;
				pRpcs->d[i] = 1;
			}
			else
			{
				pRpcs->b[i] = vecRes(i+19);
				pRpcs->d[i] = vecRes(i+58);

			}
		}
	}
	if ( pResErr != NULL )
	{
		for ( int i = 0; i < nImgPts; ++i )
		{
			double dCurTmp = 0;
			double dResX, dResY;
			if ( -1 != rfmGetxy( pGeoPts[3*i], pGeoPts[3*i+1], pGeoPts[3*i+2], pRpcs, &dResX, &dResY ) )
			{
				dResX = pImgPts[2*i] - dResX;
				dResY = pImgPts[2*i+1] - dResY;
				pResErr[i] = sqrt( dResX*dResX + dResY*dResY );
			}
			else
			{
				pResErr[i] = 0;
			}
		}
	}

	return 1;
}
RFMCALC_EXTERN_C int rfmGetNormalizedCoefs( double *pImgPts, int nImgPts, double *pGeoPts, int nGeoPts, RPCparas *pRpcs )
{
	if ( nImgPts != nGeoPts )
	{
		return -1;
	}
	vector<double> vecImgX, vecImgY, vecObjX, vecObjY, vecObjZ;
	for ( int i = 0; i < nImgPts; ++i )
	{
		pRpcs->dXOff += pImgPts[2*i];
		pRpcs->dYOff += pImgPts[2*i+1];
		pRpcs->dObjXOff += pGeoPts[3*i];
		pRpcs->dObjYOff += pGeoPts[3*i+1];
		pRpcs->dObjZOff += pGeoPts[3*i+2];

		vecImgX.push_back( pImgPts[2*i] );
		vecImgY.push_back( pImgPts[2*i+1] );

		vecObjX.push_back( pGeoPts[3*i] );
		vecObjY.push_back( pGeoPts[3*i+1] );
		vecObjZ.push_back( pGeoPts[3*i+2] );
	}
	pRpcs->dXOff /= nImgPts;
	pRpcs->dYOff /= nImgPts;
	pRpcs->dObjXOff /= nGeoPts;
	pRpcs->dObjYOff /= nGeoPts;
	pRpcs->dObjZOff /= nGeoPts;

	double dImgXMax = *max_element( vecImgX.begin(), vecImgX.end() );
	double dImgXMin = *min_element( vecImgX.begin(), vecImgX.end() );

	double dImgYMax = *max_element( vecImgY.begin(), vecImgY.end() );
	double dImgYMin = *min_element( vecImgY.begin(), vecImgY.end() );

	double dObjXMax = *max_element( vecObjX.begin(), vecObjX.end() );
	double dObjXMin = *min_element( vecObjX.begin(), vecObjX.end() );
	double dObjYMax = *max_element( vecObjY.begin(), vecObjY.end() );
	double dObjYMin = *min_element( vecObjY.begin(), vecObjY.end() );
	double dObjZMax = *max_element( vecObjZ.begin(), vecObjZ.end() );
	double dObjZMin = *min_element( vecObjZ.begin(), vecObjZ.end() );


	pRpcs->dXScale = max(fabs(dImgXMax-pRpcs->dXOff), fabs(dImgXMin-pRpcs->dXOff));
	pRpcs->dYScale = max(fabs(dImgYMax-pRpcs->dYOff), fabs(dImgYMin-pRpcs->dYOff));
	pRpcs->dObjXScale = max(fabs(dObjXMax-pRpcs->dObjXOff), fabs(dObjXMin-pRpcs->dObjXOff));
	pRpcs->dObjYScale = max(fabs(dObjYMax-pRpcs->dObjYOff), fabs(dObjYMin-pRpcs->dObjYOff));
	pRpcs->dObjZScale = max(fabs(dObjZMax-pRpcs->dObjZOff), fabs(dObjZMin-pRpcs->dObjZOff));

	return 1;
}

RFMCALC_EXTERN_C int rfmGetxy( double dX, double dY, double dZ, RPCparas *pRpcs, double *px, double *py )
{
	double P = ( dX - pRpcs->dObjXOff ) / pRpcs->dObjXScale;
	double L = ( dY - pRpcs->dObjYOff ) / pRpcs->dObjYScale;
	double H = ( dZ - pRpcs->dObjZOff ) / pRpcs->dObjZScale;


	double dNewImgX = calcHeadVal( pRpcs->a, P, L, H ) / calcTailVal( (pRpcs->b+1), P, L, H );
	double dNewImgY = calcHeadVal( pRpcs->c, P, L, H ) / calcTailVal( (pRpcs->d+1), P, L, H );

	*px = dNewImgX*pRpcs->dXScale + pRpcs->dXOff;
	*py = dNewImgY*pRpcs->dYScale + pRpcs->dYOff;

	return 1;
}

RFMCALC_EXTERN_C int rfmReadGeoPtsAndNormalizedCoefs( const char* strGeoFile, double  *&pImgPts, int &nImgPts, double *&pGeoPts, int &nGeoPts, \
													 RPCparas *pRpcs )
{
	FILE *pF = NULL;
	if( NULL == ( pF = fopen( strGeoFile, "r") ) )
	{
		return -1;
	}
	double dCurVal = 0;

	fscanf( pF, "SAMPLE_OFF:	%lf\n", &pRpcs->dXOff );
	fscanf( pF, "LINE_OFF:	%lf\n", &pRpcs->dYOff );
	fscanf( pF, "LON_OFF:  %lf\n", &pRpcs->dObjXOff );
	fscanf( pF, "LAT_OFF:    %lf\n", &pRpcs->dObjYOff );
	fscanf( pF, "HEIGHT_OFF:    %lf\n", &pRpcs->dObjZOff );

	fscanf( pF, "SAMPLE_SCALE:	%lf\n", &pRpcs->dXScale );
	fscanf( pF, "LINE_SCALE:	%lf\n", &pRpcs->dYScale );
	fscanf( pF, "LON_SCALE:   %lf\n", &pRpcs->dObjXScale );
	fscanf( pF, "LAT_SCALE:  %lf\n", &pRpcs->dObjYScale );
	fscanf( pF, "HEIGHT_SCALE:  %lf\n", &pRpcs->dObjZScale );

	double dTmp1, dTmp2, dTmp3, dTmp4, dTmp5;
	dTmp1 = dTmp2 = dTmp3 = dTmp4 = dTmp5 = 0;

	fscanf( pF, "%lf %lf %lf %lf %lf\n", &dTmp1, &dTmp2, &dTmp3, &dTmp4, &dTmp5 );
	fscanf( pF, "%lf %lf %lf %lf %lf\n", &dTmp1, &dTmp2, &dTmp3, &dTmp4, &dTmp5 );
	fscanf( pF, "%lf %lf %lf %lf %lf\n", &dTmp1, &dTmp2, &dTmp3, &dTmp4, &dTmp5 );
	fscanf( pF, "%lf %lf %lf %lf %lf\n", &dTmp1, &dTmp2, &dTmp3, &dTmp4, &dTmp5 );

	int nTmp1, nTmp2;
	vector<double> vecRes;
	do
	{
		fscanf( pF, "%lf %lf %lf %lf %lf", &dTmp1, &dTmp2, &dTmp3, &dTmp4, &dTmp5  );
		vecRes.push_back( dTmp1 );
		vecRes.push_back( dTmp2 );
		vecRes.push_back( dTmp3 );
		vecRes.push_back( dTmp4 );
		vecRes.push_back( dTmp5 );
	} while ( fgetc(pF) != EOF );
	fclose(pF);

	nImgPts = vecRes.size() / 5;
	nGeoPts = vecRes.size() / 5;

	pImgPts = new double[nImgPts*2 ];
	pGeoPts = new double[nImgPts*3 ];

	for ( int i = 0; i < nImgPts; ++i )
	{
		pImgPts[2*i] = vecRes[5*i];
		pImgPts[2*i+1] = vecRes[5*i+1];

		pGeoPts[3*i] = vecRes[5*i+2];
		pGeoPts[3*i+1] = vecRes[5*i+3];
		pGeoPts[3*i+2] = vecRes[5*i+4];
	}
	return 1;

}

RFMCALC_EXTERN_C int rfmWriteRPCs( const char* strRPCsFile, RPCparas *pRpcs )
{
	FILE * fp;

	if((fp = fopen(strRPCsFile, "w")) == NULL)
	{
		printf("Failed");
		return 0;
	}

	char str[50];
	int i;


	// Store IKONOS-II RPC
	printf("store a Corrected IKONOS RPC file.\n");
	sprintf(str, "LINE_OFF: %+#010.2f pixels\n", pRpcs->dYOff);
	fputs(str, fp);

	sprintf(str, "SAMP_OFF: %+#010.2f pixels\n", pRpcs->dXOff);
	fputs(str, fp);

	sprintf(str, "LAT_OFF: %+#012.8f degrees\n", pRpcs->dObjYOff);
	fputs(str, fp);

	sprintf(str, "LONG_OFF: %+#013.8f degrees\n", pRpcs->dObjXOff);
	fputs(str, fp);

	sprintf(str, "HEIGHT_OFF: %+#09.3f meters\n", pRpcs->dObjZOff);
	fputs(str, fp);

	sprintf(str, "LINE_SCALE: %+#010.2f pixels\n", pRpcs->dYScale);
	fputs(str, fp);


	sprintf(str, "SAMP_SCALE: %+#010.2f pixels\n", pRpcs->dXScale);
	fputs(str, fp);

	sprintf(str, "LAT_SCALE: %+#012.8f degrees\n", pRpcs->dObjYScale);
	fputs(str, fp);

	sprintf(str, "LONG_SCALE: %+#013.8f degrees\n", pRpcs->dObjXScale);
	fputs(str, fp);

	sprintf(str, "HEIGHT_SCALE: %+#09.3f meters\n", pRpcs->dObjZScale);
	fputs(str, fp);

	for(i = 0; i < 20; i++)
	{
		sprintf(str, "LINE_NUM_COEFF_%d: %+.15E\n", i+1, pRpcs->c[i]);
		fputs(str, fp);
	}
	for(i = 0; i < 20; i++)
	{
		sprintf(str, "LINE_DEN_COEFF_%d: %+.15E\n", i+1, pRpcs->d[i]);
		fputs(str, fp);
	}
	for(i = 0; i < 20; i++)
	{
		sprintf(str, "SAMP_NUM_COEFF_%d: %+.15E\n", i+1, pRpcs->a[i]);
		fputs(str, fp);
	}
	for(i = 0; i < 20; i++)
	{
		sprintf(str, "SAMP_DEN_COEFF_%d: %+.15E\n", i+1, pRpcs->b[i]);
		fputs(str, fp);
	}




	fclose(fp);
	return 1;
}
double extractnumber(const char * str)
{
	int count = 0;
	int k = 0;
	char value[30];

	memset(value, 0, 30);

	while(str[count] != ' ')
	{

		count = count + 1;
	}
	while(str[count+1] != ' ' && str[count+1] != '\n')
	{
		value[k] = str[count+1];
		k++;
		count++;
	}

	return atof(value);
}
RFMCALC_EXTERN_C int rfmReadRPCs(const char* strRPCsFile, RPCparas *pRpcs )
{
	FILE * fp = NULL;

	if((fp = fopen(strRPCsFile, "r")) == NULL)
	{
		return 0;
	}

	double num = 0;
	int i;
	char str[50];


	//read a IKONOS RPC file.
	printf("read a IKONOS RPC file.\n");
	fgets(str, 50, fp);
	pRpcs->dYOff = extractnumber(str);
	fgets(str, 50, fp);
	pRpcs->dXOff = extractnumber(str);
	fgets(str, 50, fp);
	pRpcs->dObjYOff = extractnumber(str);

	fgets(str, 50, fp);
	pRpcs->dObjXOff = extractnumber(str);
	fgets(str, 50, fp);
	pRpcs->dObjZOff = extractnumber(str);
	fgets(str, 50, fp);
	pRpcs->dYScale = extractnumber(str);
	fgets(str, 50, fp);
	pRpcs->dXScale = extractnumber(str);
	fgets(str, 50, fp);
	pRpcs->dObjYScale = extractnumber(str);
	fgets(str, 50, fp);
	pRpcs->dObjXScale = extractnumber(str);
	fgets(str, 50, fp);
	pRpcs->dObjZScale = extractnumber(str);

	for(i = 0; i < 20; i++)
	{
		fgets(str, 50, fp);
		pRpcs->c[i] = extractnumber(str);
	}
	for(i = 0; i < 20; i++)
	{
		fgets(str, 50, fp);
		pRpcs->d[i] = extractnumber(str);
	}
	for(i = 0; i < 20; i++)
	{
		fgets(str, 50, fp);
		pRpcs->a[i] = extractnumber(str);
	}
	for(i = 0; i < 20; i++)
	{
		fgets(str, 50, fp);
		pRpcs->b[i] = extractnumber(str);
	}



	fclose(fp);
	return 1;
}

double getXEquationVal( double dX, double dY, double dZ, RPCparas *pPara, double dP )
{
	double dXNorm = ( dX - pPara->dObjXOff ) / pPara->dObjXScale;
	double dYNorm = ( dY - pPara->dObjYOff ) / pPara->dObjYScale;
	double dZNorm = ( dZ - pPara->dObjZOff ) / pPara->dObjZScale;
	double dTmp1 = calcHeadVal( pPara->a, dXNorm, dYNorm, dZNorm );
	double dTmp2 = calcTailVal( (pPara->b+1), dXNorm, dYNorm, dZNorm );
	return -(dTmp1-dP*dTmp2);
}
double getYEquationVal( double dX, double dY, double dZ, RPCparas *pPara, double dP )
{
	double dXNorm = ( dX - pPara->dObjXOff ) / pPara->dObjXScale;
	double dYNorm = ( dY - pPara->dObjYOff ) / pPara->dObjYScale;
	double dZNorm = ( dZ - pPara->dObjZOff ) / pPara->dObjZScale;
	double dTmp1 = calcHeadVal( pPara->c, dXNorm, dYNorm, dZNorm );
	double dTmp2 = calcTailVal( (pPara->d+1), dXNorm, dYNorm, dZNorm );
	return -(dTmp1-dP*dTmp2);
}

double gettodXCoef(  double dXNorm, double dYNorm, double dZNorm, double *para, double dP  )
{
	double dRes1 = para[1] + para[4]*dYNorm + para[5]*dZNorm + para[7]*2.0*dXNorm + para[10]*dYNorm*dZNorm + para[11]*3.0*dXNorm*dXNorm + \
		           para[12]*dYNorm*dYNorm + para[13]*dZNorm*dZNorm + para[14]*2.0*dXNorm*dYNorm + para[17]*2.0*dXNorm*dZNorm;
	return (dRes1*dP);
}
double gettodYCoef(  double dXNorm, double dYNorm, double dZNorm, double *para, double dP  )
{
	double dRes1 = para[2] + para[4]*dXNorm + para[6]*dZNorm + para[8]*2.0*dYNorm + para[10]*dXNorm*dZNorm + para[12]*2.0*dXNorm*dYNorm + \
		para[14]*dXNorm*dXNorm + para[15]*3.0*dYNorm*dYNorm + para[16]*2.0*dZNorm*dZNorm + para[18]*2.0*dYNorm*dZNorm;
	return (dRes1*dP);
}
double gettodZCoef(  double dXNorm, double dYNorm, double dZNorm, double *para, double dP  )
{
	double dRes1 = para[3] + para[5]*dXNorm + para[6]*dYNorm + para[9]*2.0*dZNorm + para[10]*dXNorm*dYNorm + para[13]*2.0*dXNorm*dZNorm + \
		para[16]*2.0*dYNorm*dZNorm + para[17]*dXNorm*dXNorm + para[18]*dYNorm*dYNorm + para[17]*3.0*dZNorm*dZNorm;
	return (dRes1*dP);
}
bool getxToCoef( double dX, double dY, double dZ, RPCparas para, double dP, double *pCoef )
{
	double dXNorm = ( dX - para.dObjXOff ) / para.dObjXScale;
	double dYNorm = ( dY - para.dObjYOff ) / para.dObjYScale;
	double dZNorm = ( dZ - para.dObjZOff ) / para.dObjZScale;
	double dRes1 = gettodXCoef( dXNorm, dYNorm, dZNorm, para.a, 1 );
	double dRes2 = gettodXCoef( dXNorm, dYNorm, dZNorm, para.b, dP );
	pCoef[0] = (dRes1-dRes2) / para.dObjXScale;

	double dRes3 = gettodYCoef( dXNorm, dYNorm, dZNorm, para.a, 1 );
	double dRes4 = gettodYCoef( dXNorm, dYNorm, dZNorm, para.b, dP );
	pCoef[1] = (dRes3-dRes4) / para.dObjYScale;

	double dRes5 = gettodZCoef( dXNorm, dYNorm, dZNorm, para.a, 1 );
	double dRes6 = gettodZCoef( dXNorm, dYNorm, dZNorm, para.b, dP );
	pCoef[2] = (dRes5-dRes6) / para.dObjZScale;
	return true;
}
bool getyToCoef( double dX, double dY, double dZ, RPCparas para, double dP, double *pCoef )
{
	double dXNorm = ( dX - para.dObjXOff ) / para.dObjXScale;
	double dYNorm = ( dY - para.dObjYOff ) / para.dObjYScale;
	double dZNorm = ( dZ - para.dObjZOff ) / para.dObjZScale;
	double dRes1 = gettodXCoef( dXNorm, dYNorm, dZNorm, para.c, 1 );
	double dRes2 = gettodXCoef( dXNorm, dYNorm, dZNorm, para.d, dP );
	pCoef[0] = (dRes1-dRes2) / para.dObjXScale;

	double dRes3 = gettodYCoef( dXNorm, dYNorm, dZNorm, para.c, 1 );
	double dRes4 = gettodYCoef( dXNorm, dYNorm, dZNorm, para.d, dP );
	pCoef[1] = (dRes3-dRes4) / para.dObjYScale;

	double dRes5 = gettodZCoef( dXNorm, dYNorm, dZNorm, para.c, 1 );
	double dRes6 = gettodZCoef( dXNorm, dYNorm, dZNorm, para.d, dP );
	pCoef[2] = (dRes5-dRes6) / para.dObjZScale;
	return true;
}

RFMCALC_EXTERN_C int rfmIntersection( RPCparas pParasL, RPCparas pParasR, double dLx, double dLy, double dRx, double dRy, double &dX, double &dY, double &dZ)
{
	double dNormLx = ( dLx - pParasL.dXOff ) / pParasL.dXScale;
	double dNormLy = ( dLy - pParasL.dYOff ) / pParasL.dYScale;

	double dNormRx = ( dRx - pParasR.dXOff ) / pParasR.dXScale;
	double dNormRy = ( dRy - pParasR.dYOff ) / pParasR.dYScale;

	MatrixXd matA(4,3);
	VectorXd vecL(4);

	matA(0,0) = (pParasL.a[1]-dNormLx*pParasL.b[1])/pParasL.dObjXScale;
	matA(0,1) = (pParasL.a[2]-dNormLx*pParasL.b[2])/pParasL.dObjYScale;
	matA(0,2) = (pParasL.a[3]-dNormLx*pParasL.b[3])/pParasL.dObjZScale;
	matA(1,0) = (pParasR.a[1]-dNormRx*pParasR.b[1])/pParasR.dObjXScale;
	matA(1,1) = (pParasR.a[2]-dNormRx*pParasR.b[2])/pParasR.dObjYScale;
	matA(1,2) = (pParasR.a[3]-dNormRx*pParasR.b[3])/pParasR.dObjZScale;

	matA(2,0) = (pParasL.c[1]-dNormLy*pParasL.d[1])/pParasL.dObjXScale;
	matA(2,1) = (pParasL.c[2]-dNormLy*pParasL.d[2])/pParasL.dObjYScale;
	matA(2,2) = (pParasL.c[3]-dNormLy*pParasL.d[3])/pParasL.dObjZScale;
	matA(3,0) = (pParasR.c[1]-dNormRy*pParasR.d[1])/pParasR.dObjXScale;
	matA(3,1) = (pParasR.c[2]-dNormRy*pParasR.d[2])/pParasR.dObjYScale;
	matA(3,2) = (pParasR.c[3]-dNormRy*pParasR.d[3])/pParasR.dObjZScale;

	vecL(0) = dNormLx - pParasL.a[0] + (pParasL.a[1]-dNormLx*pParasL.b[1])*pParasL.dObjXOff/pParasL.dObjXScale + (pParasL.a[2]-dNormLx*pParasL.b[2])*pParasL.dObjYOff/pParasL.dObjYScale + (pParasL.a[3]-dNormLx*pParasL.b[3])*pParasL.dObjZOff/pParasL.dObjZScale;
	vecL(1) = dNormRx - pParasR.a[0] + (pParasR.a[1]-dNormRx*pParasR.b[1])*pParasR.dObjXOff/pParasR.dObjXScale + (pParasR.a[2]-dNormRx*pParasR.b[2])*pParasR.dObjYOff/pParasR.dObjYScale + (pParasR.a[3]-dNormRx*pParasR.b[3])*pParasR.dObjZOff/pParasR.dObjZScale;
	vecL(2) = dNormLy - pParasL.c[0] + (pParasL.c[1]-dNormLy*pParasL.d[1])*pParasL.dObjXOff/pParasL.dObjXScale + (pParasL.c[2]-dNormLy*pParasL.d[2])*pParasL.dObjYOff/pParasL.dObjYScale + (pParasL.c[3]-dNormLy*pParasL.d[3])*pParasL.dObjZOff/pParasL.dObjZScale;
	vecL(3) = dNormRy - pParasR.c[0] + (pParasR.c[1]-dNormRy*pParasR.d[1])*pParasR.dObjXOff/pParasR.dObjXScale + (pParasR.c[2]-dNormRy*pParasR.d[2])*pParasR.dObjYOff/pParasR.dObjYScale + (pParasR.c[3]-dNormRy*pParasR.d[3])*pParasR.dObjZOff/pParasR.dObjZScale;;

	VectorXd vecInitRes = matA.jacobiSvd(ComputeThinU | ComputeThinV).solve(vecL);

	double dTmpX = vecInitRes(0);
	double dTmpY = vecInitRes(1);
	double dTmpZ = vecInitRes(2);
	//////////////////////////////////////////////////////////////////////////
	double dThres = 1.0e-5;
	int nItera = 0;
	do
	{
		double dResX = dTmpX;
		double dResY = dTmpY;
		double dResZ = dTmpZ;
		double dLXCoef[3], dRXCoef[3], dLYCoef[3], dRYCoef[3];

		getxToCoef( dResX, dResY, dResZ, pParasL, dNormLx, dLXCoef );
		getxToCoef( dResX, dResY, dResZ, pParasR, dNormRx, dRXCoef );
		getyToCoef( dResX, dResY, dResZ, pParasL, dNormLy, dLYCoef );
		getyToCoef( dResX, dResY, dResZ, pParasR, dNormRy, dRYCoef );

		matA(0,0) = dLXCoef[0];
		matA(0,1) = dLXCoef[1];
		matA(0,2) = dLXCoef[2];

		matA(1,0) = dRXCoef[0];
		matA(1,1) = dRXCoef[1];
		matA(1,2) = dRXCoef[2];

		matA(2,0) = dLYCoef[0];
		matA(2,1) = dLYCoef[1];
		matA(2,2) = dLYCoef[2];

		matA(3,0) = dRYCoef[0];
		matA(3,1) = dRYCoef[1];
		matA(3,2) = dRYCoef[2];

		vecL(0) = getXEquationVal( dResX, dResY, dResZ, &pParasL, dNormLx );
		vecL(1) = getXEquationVal( dResX, dResY, dResZ, &pParasR, dNormRx );
		vecL(2) = getYEquationVal( dResX, dResY, dResZ, &pParasL, dNormLy );
		vecL(3) = getYEquationVal( dResX, dResY, dResZ, &pParasR, dNormRy );

		VectorXd vecCurRes = matA.jacobiSvd(ComputeThinU | ComputeThinV).solve(vecL);
		double dResiErrX = vecCurRes(0);
		double dResiErrY = vecCurRes(1);
		double dResiErrZ = vecCurRes(2);
		dTmpX = dResX + vecCurRes(0);
		dTmpY = dResY + vecCurRes(1);
		dTmpZ = dResZ + vecCurRes(2);
		if ( (fabs( dResiErrX )<dThres)&&(fabs(dResiErrY)<dThres)&&(fabs(dResiErrZ)< dThres) )
		{
			break;
		}
	} while ( nItera++ < 20 );
	dX = dTmpX;
	dY = dTmpY;
	dZ = dTmpZ;
	return 1;
}
RFMCALC_EXTERN_C int rfmFreeData( double * &pData )
{
	if ( pData != NULL )
	{
		delete[] pData;
	}
	else
	{
		return -1;
	}
	return 1;
}

void setFormerMat(MatrixXd& mat, int idx, double x, double y, double z){
		mat(idx, 0) = 1.0;

		mat(idx, 1) = z;
		mat(idx, 2) = y;
		mat(idx, 3) = x;

		mat(idx, 4) = z*y;
		mat(idx, 5) = z*x;
		mat(idx, 6) = y*x;

		mat(idx, 7) = z*z;
		mat(idx, 8) = y*y;
		mat(idx, 9) = x*x;


		mat(idx, 10) = z*y*x;
		mat(idx, 11) = z*z*y;
		mat(idx, 12) = z*z*x;
		mat(idx, 13) = y*y*z;
		mat(idx, 14) = y*y*x;
		mat(idx, 15) = z*x*x;
		mat(idx, 16) = y*y*x;
		mat(idx, 17) = z*z*z;
		mat(idx, 18) = y*y*y;
		mat(idx, 19) = x*x*x;

}

VectorXd solveMat(MatrixXd& SRC_A, VectorXd& B){
	JacobiSVD<MatrixXd> svd(SRC_A, ComputeThinU | ComputeThinV);
	VectorXd X = svd.solve(B);
	return X;

/*
	gsl_matrix* U = gsl_matrix_alloc(SRC_A->size1, SRC_A->size2);
	gsl_matrix_set_all(U, 0.0);
	gsl_matrix_memcpy(U, SRC_A);

	gsl_matrix* V = gsl_matrix_alloc(U->size2, U->size2);
	gsl_vector* S = gsl_vector_alloc(U->size2);

	gsl_matrix_set_all(V, 0.0);
	gsl_vector_set_all(S, 0.0);
	gsl_vector_set_all( X , 0.0 );

	// gsl_linalg_SV_decomp(U,V,S,X);
	gsl_vector_set_all( X , 0.0 );

	gsl_linalg_SV_decomp_jacobi(U, V, S);
	gsl_linalg_SV_solve(U, V, S, B, X);

	gsl_matrix_free(U);
	gsl_matrix_free(V);
	gsl_vector_free(S);
*/
}

void setupWeightMatrix(MatrixXd& wMat, double* denominator, int nGeoPts,  double *pGeoPts)
{
		for(int i=0; i<nGeoPts; i++){
		double x, y, z;
		x = pGeoPts[3*i];
		y = pGeoPts[3*i+1];
		z = pGeoPts[3*i+2];
		double v = denominator[0]* 1.0;
		v += denominator[1]* z;
		v += denominator[2]* y;
		v += denominator[3]* x;

		v += denominator[4]* z*y;
		v += denominator[5]* z*x;
		v += denominator[6]* y*x;
		v += denominator[7]* z*z;
		v += denominator[8]* y*y;
		v += denominator[9]* x*x;

		v += denominator[10] * z*y*x;
		v += denominator[11] * z*z*y;
		v += denominator[12] * z*z*x;
		v += denominator[13] * y*y*z;
		v += denominator[14] * y*y*x;
		v += denominator[15] * z*x*x;
		v += denominator[16] * y*y*x;
		v += denominator[17] * z*z*z;
		v += denominator[18] * y*y*y;
		v += denominator[19] * x*x*x;
		wMat(i, i) = 1.0 / v;
		//printf("weight: %+.5E\n", v);
	}

}

void squareMat(MatrixXd& mat)
{
	for(int i=0; i<mat.rows(); i++){
		double v = mat(i,i);
		v = v*v;
		mat(i,i)=v;
	}
}
