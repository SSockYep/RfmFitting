#ifndef _MAT_H
#define _MAT_H


#include "Eigen/Dense"
#include "Eigen/Core"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics.h>
using namespace Eigen;

#ifndef eps
#define eps  2.2204e-016
#endif
class CMat
{
public:
	CMat(void);
	~CMat(void);
public:
	bool SetMat( MatrixXd matA, VectorXd vecL, MatrixXd matW );
	//L_Curve + ICCV
	int RobustSolve( VectorXd &vecRes, VectorXd &vecResids,  int nIteraMax, int nCurvePts, double dICCVThres );

	int RobustSolveSimple( MatrixXd matA, MatrixXd matL, VectorXd &vecRes, VectorXd &vecResids,  int nIteraMax, int nCurvePts, double dICCVThres );

	//Get Matrix Condition Count Num
	double GetConditionNum( MatrixXd mat );

	double GetEuclidean2NormVec( VectorXd vec );

	//Get MatEuclidean2Norm of Matrix
	double GetEuclidean2NormMat( MatrixXd mat );


private:
	double estimateKParaUseL_Curve( MatrixXd matA, VectorXd vecL, int nPoints = 200 );

	int solveByICCV( VectorXd vecInitX, VectorXd &vecX, int nMaxItera = 20, double dThres = 1.0e-7 );

	int solveByICCV( MatrixXd matA, VectorXd vecL, VectorXd vecInitX, VectorXd &vecX, int nMaxItera = 20, double dThres = 1.0e-7 );

	VectorXd dotMulti( VectorXd v1, VectorXd v2 );
	VectorXd dotAdd( VectorXd v1, VectorXd v2 );

	MatrixXd m_MatA;
	MatrixXd m_MatW;
	VectorXd m_VecL;
};








#endif
