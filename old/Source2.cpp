#include "RFMCalc.h"
#include <cstring>

int main(int argc, char** argv)
{
	vector<double> imgPts, geoPts;
	char* gcp_name, * w_file_name;
	int max_iter = 10;
	double thre = 0.1;
	RPCparas rpcs;
	RFMCalcPara rfmCalcPara;
	RFMType rfmType = Three_Order_Coefs;

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
	}

	rfmReadGeoPtsAndNormalizedCoefs(gcp_name, imgPts, geoPts, rpcs);

	rfmCalc(imgPts, geoPts, rpcs, rfmCalcPara, Two_Order_Coefs);

	rfmWriteRPCs(w_file_name, rpcs, Two_Order_Coefs);

	free(gcp_name);
	free(w_file_name);
	
	return 0;
}
