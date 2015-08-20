#include "RFMCalc.h"
#include <cstring>

int main(int argc, char** argv)
{
	vector<double> imgPts, geoPts;
	char* gcp_name, * w_file_name;
	int max_iter;
	double thre;
	RPCparas rpcs;
	RFMCalcPara rfmCalcPara;
	RFMType rfmType = Three_Order_Coefs;
	bool standard_out = false;

	for (int i = 1; i < argc; i++)
	{
		if(!strcmp(argv[i], "-infile"))
		{
			gcp_name = (char*)malloc(sizeof(char)*(strlen(argv[i+1])+1));
			strcpy(gcp_name, argv[i+1]);
			cout << "infile: " << gcp_name << endl;
		}
		else if(!strcmp(argv[i], "-outfile"))
		{
			w_file_name = (char*)malloc(sizeof(char)*(strlen(argv[i+1])+1));
			strcpy(w_file_name, argv[i+1]);
			cout << "outfile: " << w_file_name << endl;
		}
		else if(!strcmp(argv[i], "-order"))
		{
			int tmp;
			sscanf(argv[i+1], "%d", &tmp);
			if(tmp == 2)
				rfmType = Two_Order_Coefs;
			else if (tmp == 3)
				rfmType = Three_Order_Coefs;
			else
			{
				cout << "wrong order number" << endl;
				return -1;
			}
		}
		else if(!strcmp(argv[i], "-max_iter"))
		{
			sscanf(argv[i+1], "%d", &max_iter);
			rfmCalcPara.nMaxItera = max_iter;
		}
		
		else if(!strcmp(argv[i], "-threshold"))
		{
			sscanf(argv[i+1], "%lf", &thre);
			rfmCalcPara.dICCVThres = thre;
		}
		else if(!strcmp(argv[i], "-out_type"))
		{
			if(!strcmp(argv[i+1], "standard"))
				standard_out = true;
		}
	}

	rfmReadGeoPtsAndNormalizedCoefs(gcp_name, imgPts, geoPts, rpcs);

	rfmCalc(imgPts, geoPts, rpcs, rfmCalcPara, rfmType);

	if(standard_out)
		rfmWriteRPCs(w_file_name, rpcs, Three_Order_Coefs);
	else
		rfmWriteRPCs(w_file_name, rpcs, rfmType);

	free(gcp_name);
	free(w_file_name);
	
	return 0;
}
