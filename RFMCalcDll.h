// The following ifdef block is the standard way of creating macros which make exporting
// from a DLL simpler. All files within this DLL are compiled with the RFMCALCDLL_EXPORTS
// symbol defined on the command line. this symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see
// RFMCALCDLL_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.

#ifndef _RFMCLACHEAD
#define _RFMCLACHEAD

#if defined _WIN32 || defined _WIN64
#include "targetver.h"
#include <windows.h>
#endif

#if defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64
#define RFMCALC_CDECL __cdecl
#define RFMCALC_STDCALL __stdcall
#else
#define RFMCALC_CDECL
#define RFMCALC_STDCALL
#include <stddef.h>
#include <stdio.h>
#endif

#ifndef RFMCALC_EXTERN_C
#ifdef __cplusplus
#define RFMCALC_EXTERN_C extern "C"
#define RFMCALC_DEFAULT(val) = val
#else
#define RFMCALC_EXTERN_C
#define RFMCALC_DEFAULT(val)
#endif
#endif

#if ( defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64 ) && defined RFMCALCDLL_EXPORTS
#define RFMCALC_EXPORTS __declspec(dllexport)
#else
#define RFMCALC_EXPORTS
#endif

#ifndef RFMCALC_API
#define RFMCALC_API(rettype) RFMCALC_EXTERN_C RFMCALC_EXPORTS rettype RFMCALC_CDECL
#endif

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics.h>
#include "Mat.h"

//////////////////////////////////////////////////////////////////////////
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
	Three_Order_Coefs = 1,
	Two_Order_Coefs = 0

}RFMType;

typedef struct tagRFMCalcPara
{
	int nL_Curve_Points;
	double dICCVThres;

	int nMaxItera;
	double dThres;
	tagRFMCalcPara()
	{
		nL_Curve_Points = 200;
		dICCVThres = 1.0e-4;
		nMaxItera = 20;
		dThres = 1.0e-6;
	}
}RFMCalcPara;
void setFormerMat(MatrixXd& mat, int idx, double x, double y, double z);
VectorXd solveMat(MatrixXd& SRC_A, VectorXd& B);
void setupWeightMatrix(MatrixXd& wMat, double* denominator, int nGeoPts,  double *pGeoPts);
void squareMat(MatrixXd& mat);


RFMCALC_API( int ) rfmTest( double *pImgPts, int nImgPts, double *pGeoPts, int nGeoPts, \
						    RPCparas *pRpcs, RFMCalcPara para, RFMType rfmType = Three_Order_Coefs );


RFMCALC_API( int ) rfmCalc( double *pImgPts, int nImgPts, double *pGeoPts, int nGeoPts, \
						   RPCparas *pRpcs, RFMCalcPara &para, double* pResErr = NULL, RFMType rfmType = Three_Order_Coefs );

RFMCALC_API( int ) rfmCalcSimple( double *pImgPts, int nImgPts, double *pGeoPts, int nGeoPts, \
						   RPCparas *pRpcs, RFMCalcPara &para, double* pResErr = NULL, RFMType rfmType = Three_Order_Coefs );

RFMCALC_API( int ) rfmGetxy( double dX, double dY, double dZ, RPCparas *pRpcs, double *px, double *py );

RFMCALC_API( int ) rfmGetNormalizedCoefs( double *pImgPts, int nImgPts, double *pGeoPts, int nGeoPts, \
						   RPCparas *pRpcs );


RFMCALC_API( int ) rfmReadGeoPtsAndNormalizedCoefs( const char* strGeoFile, double  *&pImgPts, int &nImgPts, double *&pGeoPts, int &nGeoPts, \
										 RPCparas *pRpcs );

RFMCALC_API( int ) rfmWriteRPCs( const char* strRPCsFile, RPCparas *pRpcs );

RFMCALC_API( int ) rfmFreeData( double * &pData );

RFMCALC_API( int ) rfmReadRPCs(const char* strRPCsFile, RPCparas *pRpcs );

RFMCALC_API( int ) rfmIntersection( RPCparas pParasL, RPCparas pParasR, double dLx, double dLy, double dRx, double dRy, double &dX, double &dY, double &dZ);

#endif
