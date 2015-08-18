#include <string.h>
#include "RFMCalcDll.h"

int main(int argc, char** argv) //gcpname writerpc maxiteration sur
{
	int nImgPts = 0;
	int nGeoPts = 0;
	double * pImgPts, *pGeoPts, *pErr;
	char* gcp_name, *w_file_name;
	int max_iter;
	double s;
	pImgPts = pGeoPts = 0;
	RPCparas rpcs;
	RFMCalcPara rfmCalcPara;
	sscanf(argv[3], "%d", &max_iter);
	sscanf(argv[4], "%lf", &s);
	rfmCalcPara.dICCVThres = s;
	rfmCalcPara.nMaxItera = max_iter;
	gcp_name = (char*)malloc(sizeof(char)*(strlen(argv[1])+1));
	w_file_name = (char*)malloc(sizeof(char)*(strlen(argv[2])+1));
	strcpy(gcp_name, argv[1]);
	strcpy(w_file_name, argv[2]);
printf("reading file....\n");
	rfmReadGeoPtsAndNormalizedCoefs( gcp_name, pImgPts, nImgPts, pGeoPts, nGeoPts, &rpcs );

	pErr = new double[nImgPts];
//	rfmCalcSimple( pImgPts, nImgPts, pGeoPts, nGeoPts, &rpcs, rfmCalcPara, pErr );
printf("calculating...\n");
	rfmCalc( pImgPts, nImgPts, pGeoPts, nGeoPts, &rpcs, rfmCalcPara, pErr );
printf("writing file...\n");
        rfmWriteRPCs( w_file_name, &rpcs );
	rfmFreeData( pImgPts );
	rfmFreeData( pGeoPts );
	rfmFreeData( pErr );

	free(gcp_name);
	free(w_file_name);
	
}
