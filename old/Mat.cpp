#include "Mat.h"
#include <stdio.h>
#include <complex>
#include "Eigen/SVD"
CMat::CMat(void)
{
}

CMat::~CMat(void)
{
}
bool CMat::SetMat( MatrixXd matA, VectorXd matL, MatrixXd matW  )
{
	m_MatA = matA;
	m_VecL = matL;
	m_MatW = matW;
	return true;
}
int CMat::RobustSolve( VectorXd &vecRes, VectorXd &vecResids, int nIteraMax, int nCurvePts, double dICCVThres  )
{
	if ( m_MatA.rows() != m_VecL.size() )
	{
		return -1;
	}
	double dCoef = estimateKParaUseL_Curve( m_MatA, m_VecL, nCurvePts );
	
	MatrixXd matInter(m_MatA.cols(),m_MatA.cols() );
	matInter.setIdentity();

	MatrixXd matTmp = m_MatA.transpose() * m_MatW * m_MatA +  matInter* dCoef;
	VectorXd vecTmpX = matTmp.inverse() * m_MatA.transpose() * m_MatW * m_VecL;
	
	VectorXd vecXRes;
	if ( 1 == solveByICCV( vecTmpX, vecXRes, 20, dICCVThres ) )
	{
		vecRes = vecXRes;

		vecResids = m_MatW * m_MatA * vecRes - m_MatW * m_VecL;
		return 1;
	}
	else
	{
		return -1;
	}

}
int CMat::RobustSolveSimple( MatrixXd matA, MatrixXd matL, VectorXd &vecRes, VectorXd &vecResids, int nIteraMax, int nCurvePts, double dICCVThres  )
{
	if ( matA.rows() != matL.size() )
	{
		return -1;
	}
	double dCoef = estimateKParaUseL_Curve( matA, matL, nCurvePts );

	MatrixXd matInter(matA.cols(),matA.cols() );
	matInter.setIdentity();

	MatrixXd matTmp = matA.transpose() * matA +  matInter* dCoef;
	VectorXd vecTmpX = matTmp.inverse() * matA.transpose() * matL;

	VectorXd vecXRes;
	if ( 1 == solveByICCV( matA, matL, vecTmpX, vecXRes, nIteraMax, dICCVThres ) )
	{
		vecRes = vecXRes;

		vecResids = matA * vecRes - matL;
		return 1;
	}
	else
	{
		return -1;
	}

}
double CMat::GetConditionNum( MatrixXd mat )
{
	if ( mat.rows() != mat.cols() )
	{
		return -1;
	}
	
	double dRes = GetEuclidean2NormMat( mat ) * GetEuclidean2NormMat( mat.inverse() );
	
	return dRes;
}
double CMat::GetEuclidean2NormMat( MatrixXd mat )
{
	if ( ( mat.rows() == 0 ) || ( mat.cols() == 0 ) )
	{
		return -1;
	}
	if ( mat.cols() == 1 )
	{
		return GetEuclidean2NormVec( mat.col(0) );
	}

	MatrixXd tmpRes = mat.transpose() * mat;

	VectorXcd eivals = tmpRes.eigenvalues();
	
	double dMaxVal = -10e37;

	for ( int i = 0; i < eivals.cols(); ++i )
	{
		std::complex<double> curVal = eivals[i];
		if ( fabs( curVal.imag() ) < 1.0e-7 )
		{
			if ( dMaxVal < curVal.real() )
			{
				dMaxVal = curVal.real();
			}
		}
	}
	return sqrt( dMaxVal );
}
double CMat::GetEuclidean2NormVec( VectorXd vec )
{
	return vec.norm();
	double dRes = 0;
	for ( int i = 0; i < vec.size(); ++i )
	{
		dRes += vec(i)*vec(i);
	}
	//dRes /= vec.size();
	return sqrt( dRes );
}

double CMat::estimateKParaUseL_Curve( MatrixXd matA, VectorXd vecL, int nPoints )
{
	JacobiSVD<MatrixXd> svd( matA, ComputeThinU | ComputeThinV);

	MatrixXd matU = svd.matrixU();
	VectorXd vecS = svd.singularValues();//col vector
	MatrixXd matV = svd.matrixV();

	VectorXd vecB = matU.transpose() * vecL;

	double dB2 = pow( GetEuclidean2NormVec( vecL), 2 ) - pow( GetEuclidean2NormVec( vecB ), 2 );

	VectorXd vecParas( nPoints );
	
	
	vecParas( vecParas.size()-1 ) = std::max( vecS(vecS.size()-1), 16*vecS(0)*eps );
	double dRatio = pow ( ( vecS(0) / vecParas( vecParas.size() - 1 ) ), (1.0 / ( nPoints-1) ) );

	for( int i = nPoints-2; i >= 0; --i )
	{
		vecParas(i) = dRatio * vecParas(i+1);
	}
	
	VectorXd vecXNorm(nPoints);
	VectorXd vecAXLNorm(nPoints);

	for ( int i = 0; i < nPoints; ++i )
	{
		//VectorXd vecF( vecS.size() );
		double dCurPara = vecParas(i);
		//VectorXd vecX( vecS.size() );
		//VectorXd vecAXL( vecS.size() );
		double dF, dX, dAXL;
		dF = dX = dAXL = 0;
		double dInDev = 0;

		for ( int j = 0; j < vecS.size(); ++j )
		{
			double dCurSS = vecS(j)*vecS(j);
			dF = dCurSS / ( dCurSS + dCurPara*dCurPara );
			
			dX += pow( (dF * vecB(j) / vecS(j)), 2 );

			dAXL += pow( (( 1.0 - dF ) * vecB(j)), 2 );

			dInDev += ( 1-dF) * dF * dF * vecB(j) * vecB(j) / ( vecS(j) * vecS(j) );
		}
		double dXNormS = sqrt(dX) ;
		double dAXLNormS = sqrt(dAXL);

		vecXNorm(i) = dXNormS;
		vecAXLNorm(i) = sqrt( dAXLNormS * dAXLNormS /*+ dB2*/ );

	}
	VectorXd vecK( nPoints );
	for ( int i = 0; i < nPoints; ++i )
	{
		VectorXd vecF( vecS.size() ),vecCF( vecS.size() ), vecF1( vecS.size() ), vecF2( vecS.size() );
		VectorXd vecPhi( vecS.size() ), vecPsi( vecS.size() ), vecdPhi( vecS.size() ), vecdPsi( vecS.size() );
		for ( int j = 0; j < vecS.size(); ++j )
		{
			vecF(j) = vecS(j)*vecS(j) / (vecS(j)*vecS(j) + vecParas(i)*vecParas(i) );
			vecCF(j) = 1 - vecF(j);
			vecF1(j) = -2.0 * vecF(j) * vecCF(j) / vecParas(i);
			vecF2(j) = -vecF1(j) * ( 3 - 4*vecF(j) ) / vecParas(i) ;

			vecPhi(j) = vecF(j)*vecF1(j)*( vecB(j) * vecB(j) ) / ( vecS(j)*vecS(j) );
			vecPsi(j) = vecCF(j)*vecF1(j)*vecB(j)*vecB(j);
			vecdPhi(j) = (vecF1(j)*vecF1(j)+vecF(j)*vecF2(j)) * ( vecB(j) * vecB(j) ) / ( vecS(j)*vecS(j) );
			vecdPsi(j) = ( vecCF(j)*vecF2(j) - vecF1(j)*vecF1(j) ) * vecB(j)*vecB(j);
		}
		double dPhi = vecPhi.sum();
		double dPsi = vecPsi.sum();
		double ddPhi = vecdPhi.sum();
		double ddPsi = vecdPsi.sum();

		double dDeta = dPhi / vecXNorm(i);
		double dRho = -dPsi / vecAXLNorm(i);
		double ddDeta = ddPhi / vecXNorm(i) - dDeta * ( dDeta / vecXNorm(i) );
		double ddRho = -ddPsi / vecAXLNorm(i) - dRho * ( dRho / vecAXLNorm(i) );

		double dLogeta = dDeta / vecXNorm(i);
		double dLogrho = dRho / vecAXLNorm(i);
		double ddLogeta = ddDeta / vecXNorm(i) - dLogeta * dLogeta;
		double ddLogrho = ddRho / vecAXLNorm(i) - dLogrho * dLogrho;

		vecK(i) = -(dLogrho*ddLogeta - ddLogrho * dLogeta ) / pow( (dLogrho*dLogrho+dLogeta*dLogeta), 1.5);
		 
	}
	
	double dMinVal = 1.0e22;
	int nFlag = -1;
	for ( int i = 0; i < nPoints; ++i )
	{
		if ( vecK(i) < dMinVal )
		{
			nFlag = i;
			dMinVal = vecK(i);
		}
	}
	
	if( nFlag  != -1 )
	{
		if ( dMinVal < 0 )
		{
			return vecParas(nFlag);
		}
		else
		{
			return vecParas( vecParas.size() - 1 );
		}
		
	}
	else
	{
		return -1.0e22;
	}
}
int CMat::solveByICCV( MatrixXd matA, VectorXd vecL, VectorXd vecInitX, VectorXd &vecX, int nMaxItera , double dThres  )
{
	MatrixXd matI( vecInitX.size(), vecInitX.size() );
	matI.setIdentity();
	MatrixXd matTmp = matA.transpose() * matA + matI;
	MatrixXd matQ = matTmp.inverse();

	int nIter = 0;
	VectorXd vecTmpX = vecInitX;
	while( nIter++ < nMaxItera )
	{
		VectorXd vecTmpRes = matQ * ( matA.transpose() * vecL + vecTmpX );
		bool bIsOk = true;
		for ( int i = 0; i < vecTmpRes.size(); ++i )
		{
			double dResid = fabs( vecTmpRes(i)-vecTmpX(i) );
			if ( dResid > dThres )
			{
				bIsOk = false;
			}
		}
		vecTmpX = vecTmpRes;
		if ( bIsOk == true )
		{
			break;
		}
	}
	vecX = vecTmpX;

	return 1;
}
int CMat::solveByICCV( VectorXd vecInitX, VectorXd &vecX, int nMaxItera, double dThres )
{
	MatrixXd matI( vecInitX.size(), vecInitX.size() );
	matI.setIdentity();
	MatrixXd matTmp = m_MatA.transpose() * m_MatW * m_MatA + matI;
	MatrixXd matQ = matTmp.inverse();
	
	int nIter = 0;
	VectorXd vecTmpX = vecInitX;
	while( nIter++ < nMaxItera )
	{
		VectorXd vecTmpRes = matQ * ( m_MatA.transpose() * m_MatW * m_VecL + vecTmpX );
		bool bIsOk = true;
		for ( int i = 0; i < vecTmpRes.size(); ++i )
		{
			double dResid = fabs( vecTmpRes(i)-vecTmpX(i) );
			if ( dResid > dThres )
			{
				bIsOk = false;
			}
		}
		vecTmpX = vecTmpRes;
		if ( bIsOk == true )
		{
			break;
		}
	}
	vecX = vecTmpX;

	return 1;
}
VectorXd CMat::dotMulti( VectorXd v1, VectorXd v2 )
{
	if ( v1.size() != v2.size() )
	{
		return v1;
	}
	VectorXd vRes(v1.size());
	for ( int i = 0; i < v1.size(); ++i )
	{
		vRes(i) = v1(i) * v2(i);
	}
	return vRes;
}
VectorXd CMat::dotAdd( VectorXd v1, VectorXd v2 )
{
	if ( v1.size() != v2.size() )
	{
		return v1;
	}
	VectorXd vRes(v1.size());
	for ( int i = 0; i < v1.size(); ++i )
	{
		vRes(i) = v1(i) + v2(i);
	}
	return vRes;
}
