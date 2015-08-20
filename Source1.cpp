#include "RFMCalc.h"
#include <cstring>

int main(int argc, char** argv)
{
	vector<double> imgPts, geoPts;

	char* gcp_name, * w_file_name;
	int max_iter;
	double s;
	RPCparas rpcs;
	RFMCalcPara rfmCalcPara;

	sscanf(argv[3], "%d", &max_iter);
	sscanf(argv[4], "%lf", &s);
	rfmCalcPara.nMaxItera = max_iter;
	gcp_name = (char*)malloc(sizeof(char)*(strlen(argv[1])+1));
	w_file_name = (char*)malloc(sizeof(char)*(strlen(argv[2])+1));

	strcpy(gcp_name, argv[1]);
	strcpy(w_file_name, argv[2]);

	rfmReadGeoPtsAndNormalizedCoefs(gcp_name, imgPts, geoPts, rpcs);
	cout << "after dObjZScale: " << rpcs.dObjZScale << endl;

	rfmCalc(imgPts, geoPts, rpcs, rfmCalcPara);

	rfmWriteRPCs(w_file_name, rpcs);

	free(gcp_name);
	free(w_file_name);
	
	return 0;
}
