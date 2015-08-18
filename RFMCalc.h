#ifndef _RFM_CALC_H_
#define _RFM_CALC_H_

#include <vector>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics.h>
#include <Eigen/Dense>
#include <Eigen/Core>

using namespace Eigen;
using namespace std;

typedef struct tagRPCCoefs
{
	double a[20], b[20], c[20], d[20];
	double dXScale, dYScale, dXOff, dYOff;
	double dObjXScale, dObjYScale, dObjZScale, dObjXOff, dObjYOff, dObjZOff;
	tagRPCCoefs()
	{
		for ( int i = 0; i < 20; ++i )
		{
			a[i] = b[i] = c[i] = d[i] = 0;
		}
		dXScale = dYScale = dXOff = dYOff = 0;
		dObjXScale = dObjYScale = dObjZScale = dObjXOff = dObjYOff = dObjZOff = 0;
	}
}RPCparas;

typedef enum tagRFMType
{
	Three_Order_Coefs = 0,
	Two_Order_Coefs = 1,
}RFMType;

typedef struct tagRFMCalcPara
{
	int nL_Curve_Points;
	double dICCVThres;
	tagRFMCalcPara()
	{
		nL_Curve_Points = 200;
		dICCVThres = 1.0e-4;
		nMaxItera = 20;
		dThres - 1.0e-6;
	}
}RFMCalcPara;

void setFormerMat(MatrixXd& mat, int idx, double x, double y, double z);
void setFormerMatQuad(MatrixXd& mat, int idx, double x, double y, double z);
VectorXd solveMat(MatrixXd& SRC_A, VectorXd& B);
void setupWeightMatrix(MatrixXd& wMat, double* denominator, vector<double> geoPts);
void setupWeightMatrixQuad(MatrixXd& wMat, double* denominator, vector<double> geoPts);
//void squareMat(MatrixXd& mat);

int rfmCalc(vector<double> imgPts, vector<double> geoPts, \
					RPCparas& rpcs, RFMCalcPara para, RFMType rfmType = Three_Order_Coefs);

int rfmGetNormalizedCoefs(vector<double> imgPts, vector<double> geoPts, RPCparas& rpcs );
int rfmReadGeoPtsAndNormalizedCoefs(const char* strGeoFile, vector<double>& imgPts, vector<double>& geoPts, \
					RPCparas& rpcs, RFMType rfmType = Three_Order_Coefs);
int rfmWriteRPCs(const char* strRPCsFile, PRCparas& rpcs, RFMType rfmType = Three_Order_Coefs);
int rfmReadRPCs(const char* strRPCsFile, RPCparas& rpcs, RFMType rfmType = Three_Order_Coefs);

#endif
