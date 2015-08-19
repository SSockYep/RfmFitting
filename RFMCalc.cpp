#include "RFMCalc.h"

int rfmCalc(vector<double> imgPts, vector<double> geoPts, \
					RPCparas& rpcs, RFMCalcPara para, RFMType rfmType)
{
	if(rfmType == Three_Order_Coefs)
	{
		double first_sqrt_mean;
		if(imgPts.size()/2 != geoPts.size()/3)
		{
			cout << "error occured" << endl;
			return -1;
		}

		MatrixXd M_mat(imgPts.size()/2, 39);
		MatrixXd M_mat_cpy(imgPts.size()/2, 39);
		MatrixXd N_mat(imgPts.size()/2, 39);
		VectorXd row_M(imgPts.size()/2), col_N(imgPts.size()/2);
		
		M_mat.setZero();
		M_mat_cpy.setZero();
		N_mat.setZero();
		row_M.setZero();
		col_N.setZero();

		//init first values
		for(int i=0; i<imgPts.size()/2; i++)
		{
			double cX, cY, cZ;
			double iX, iY;
			cX = (geoPts[3*i] - rpcs.dObjXOff) / rpcs.dObjXScale;
			cY = (geoPts[3*i+1] - rpcs.dObjYOff) / rpcs.dObjYScale ;
			cZ = (geoPts[3*i+2] - rpcs.dObjZOff) / rpcs.dObjZScale;
			iX = (imgPts[2*i] - rpcs.dXOff) / rpcs.dXScale;
			iY = (imgPts[2*i] - rpcs.dYOff) / rpcs.dYScale;

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

		//calculate weighted M matrix (N x N), row_M vector

		for(int i = 0; i < weight_M.rows(); i++)
		{
			weight_M(i,i) = 1.0;
			weight_N(i,i) = 1.0;
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
		for(int i = 0; i < imgPts.size()/2; i++)
			sum += pow(error_vector(i), 2.0);
		first_sqrt_mean = sqrt(sum/imgPts.size()/2);
		/*
		 * a: samp_num_coeff
		 * b: samp_den_coeff
		 * c: line_num_coeff
		 * d: line_den_coeff
		 * */
		for(int i=0; i<20; i++)
			rpcs.c[i] = r_coeff(i);
		for(int i=0; i<20; i++)
		{
			if(i==0) rpcs.d[i] = 1.0;
			else rpcs.d[i] = r_coeff(i+19);
		}
		for(int i=0; i<20; i++)
				rpcs.a[i] = c_coeff(i);
		for(int i=0; i<20; i++)
		{
			if(i==0) rpcs.b[i] = 1.0;
			else rpcs.b[i] = c_coeff(i+19);
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

		int count = 0;
		do
		{
			weight_M.setZero();
			weight_N.setZero();

			M_weighted.setZero();
			N_weighted.setZero();
			row_M_weighted.setZero();
			col_N_weighted.setZero();

			setupWeightMatrix(weight_M, rpcs.d, geoPts);
			setupWeightMatrix(weight_N, rpcs.b, geoPts);
			not_squared_weight_M = weight_M;

			weight_M = weight_M * weight_M;
			weight_N = weight_N * weight_N;
			//squareMat(weight_M);
			//squareMat(weight_N);
			// (NxM) * (MxM) => (NxM)

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
			double sqrt_mean;
			for(int i = 0; i < imgPts.size()/2; i++)
				sum += pow(error_vector(i), 2.0);
			sqrt_mean = sqrt(sum/imgPts.size()/2);

			for(int i=0; i<20; i++)
				rpcs.c[i] = r_coeff(i);
			for(int i=0; i<20; i++)
			{
				if(i==0) rpcs.d[i] = 1.0;
				else rpcs.d[i] = r_coeff(i+19);
			}
			for(int i=0; i<20; i++)
					rpcs.a[i] = c_coeff(i);
			for(int i=0; i<20; i++)
			{
				if(i==0) rpcs.b[i] = 1.0;
				else rpcs.b[i] = c_coeff(i+19);
			}

			printf("%d: %+.5E\n", iter_n, sqrt_mean);
			if(sqrt_mean < para.dICCVThres)
				count++;
			else
				count = 0;
			if(count >= 5)
				break;
			/* CHECK END CONDITION */
		} while(++iter_n <= para.nMaxItera);
	}
	
	else if(rfmType == Two_Order_Coefs)
	{
		double first_sqrt_mean;
		if(imgPts.size()/2 != geoPts.size()/3)
		{
			cout << "error occured" << endl;
			return -1;
		}

		MatrixXd M_mat(imgPts.size()/2, 19);
		MatrixXd M_mat_cpy(imgPts.size()/2, 19);
		MatrixXd N_mat(imgPts.size()/2, 19);
		VectorXd row_M(imgPts.size()/2), col_N(imgPts.size()/2);
		
		M_mat.setZero();
		M_mat_cpy.setZero();
		N_mat.setZero();
		row_M.setZero();
		col_N.setZero();

		//init first values
		for(int i=0; i<imgPts.size()/2; i++)
		{
			double cX, cY, cZ;
			double iX, iY;
			cX = (geoPts[3*i] - rpcs.dObjXOff) / rpcs.dObjXScale;
			cY = (geoPts[3*i+1] - rpcs.dObjYOff) / rpcs.dObjYScale ;
			cZ = (geoPts[3*i+2] - rpcs.dObjZOff) / rpcs.dObjZScale;
			iX = (imgPts[2*i] - rpcs.dXOff) / rpcs.dXScale;
			iY = (imgPts[2*i] - rpcs.dYOff) / rpcs.dYScale;

			row_M(i) = iX;
			col_N(i) = iY;

			setFormerMatQuad(M_mat, i, cX, cY, cZ);
			setFormerMatQuad(N_mat, i, cX, cY, cZ);

			for(int j=10; j<19; j++)
			{
				M_mat(i,j) = (-1.0*iX) * M_mat(i, j-9);
				N_mat(i,j) = (-1.0*iY) * N_mat(i, j-9);
			}
		}

		VectorXd r_coeff(19);
		VectorXd c_coeff(19);

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

		//calculate weighted M matrix (N x N), row_M vector

		for(int i = 0; i < weight_M.rows(); i++)
		{
			weight_M(i,i) = 1.0;
			weight_N(i,i) = 1.0;
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
		for(int i = 0; i < imgPts.size()/2; i++)
			sum += pow(error_vector(i), 2.0);
		first_sqrt_mean = sqrt(sum/imgPts.size()/2);
		/*
		 * a: samp_num_coeff
		 * b: samp_den_coeff
		 * c: line_num_coeff
		 * d: line_den_coeff
		 * */
		for(int i=0; i<10; i++)
			rpcs.c[i] = r_coeff(i);
		for(int i=0; i<10; i++)
		{
			if(i==0) rpcs.d[i] = 1.0;
			else rpcs.d[i] = r_coeff(i+9);
		}
		for(int i=0; i<10; i++)
				rpcs.a[i] = c_coeff(i);
		for(int i=0; i<10; i++)
		{
			if(i==0) rpcs.b[i] = 1.0;
			else rpcs.b[i] = c_coeff(i+9);
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

		int count = 0;
		do
		{
			weight_M.setZero();
			weight_N.setZero();

			M_weighted.setZero();
			N_weighted.setZero();
			row_M_weighted.setZero();
			col_N_weighted.setZero();

			setupWeightMatrixQuad(weight_M, rpcs.d, geoPts);
			setupWeightMatrixQuad(weight_N, rpcs.b, geoPts);
			not_squared_weight_M = weight_M;

			weight_M = weight_M * weight_M;
			weight_N = weight_N * weight_N;
			//squareMat(weight_M);
			//squareMat(weight_N);
			// (NxM) * (MxM) => (NxM)

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
			double sqrt_mean;
			for(int i = 0; i < imgPts.size()/2; i++)
				sum += pow(error_vector(i), 2.0);
			sqrt_mean = sqrt(sum/imgPts.size()/2);

			for(int i=0; i<10; i++)
				rpcs.c[i] = r_coeff(i);
			for(int i=0; i<10; i++)
			{
				if(i==0) rpcs.d[i] = 1.0;
				else rpcs.d[i] = r_coeff(i+9);
			}
			for(int i=0; i<10; i++)
					rpcs.a[i] = c_coeff(i);
			for(int i=0; i<10; i++)
			{
				if(i==0) rpcs.b[i] = 1.0;
				else rpcs.b[i] = c_coeff(i+9);
			}

			printf("%d: %+.5E\n", iter_n, sqrt_mean);
			if(sqrt_mean < para.dICCVThres)
				count++;
			else
				count = 0;
			if(count >= 5)
				break;
			/* CHECK END CONDITION */
		} while(++iter_n <= para.nMaxItera);
	}
}

int rfmGetNormalizedCoefs(vector<double> imgPts, vector<double> geoPts, RPCparas& rpcs)
{
	if ( imgPts.size()/2 != geoPts.size()/3 )
	{
		return -1;
	}
	vector<double> vecImgX, vecImgY, vecObjX, vecObjY, vecObjZ;
	for ( int i = 0; i < imgPts.size()/2; ++i )
	{
		rpcs.dXOff += imgPts[2*i];
		rpcs.dYOff += imgPts[2*i+1];
		rpcs.dObjXOff += geoPts[3*i];
		rpcs.dObjYOff += geoPts[3*i+1];
		rpcs.dObjZOff += geoPts[3*i+2];

		vecImgX.push_back( imgPts[2*i] );
		vecImgY.push_back( imgPts[2*i+1] );

		vecObjX.push_back( geoPts[3*i] );
		vecObjY.push_back( geoPts[3*i+1] );
		vecObjZ.push_back( geoPts[3*i+2] );
	}
	rpcs.dXOff /= imgPts.size()/2;
	rpcs.dYOff /= imgPts.size()/2;
	rpcs.dObjXOff /= geoPts.size()/3;
	rpcs.dObjYOff /= geoPts.size()/3;
	rpcs.dObjZOff /= geoPts.size()/3;

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


	rpcs.dXScale = max(fabs(dImgXMax-rpcs.dXOff), fabs(dImgXMin-rpcs.dXOff));
	rpcs.dYScale = max(fabs(dImgYMax-rpcs.dYOff), fabs(dImgYMin-rpcs.dYOff));
	rpcs.dObjXScale = max(fabs(dObjXMax-rpcs.dObjXOff), fabs(dObjXMin-rpcs.dObjXOff));
	rpcs.dObjYScale = max(fabs(dObjYMax-rpcs.dObjYOff), fabs(dObjYMin-rpcs.dObjYOff));
	rpcs.dObjZScale = max(fabs(dObjZMax-rpcs.dObjZOff), fabs(dObjZMin-rpcs.dObjZOff));

	return 1;
}

VectorXd solveMat(MatrixXd& SRC_A, VectorXd& B)
{
	JacobiSVD<MatrixXd> svd(SRC_A, ComputeThinU | ComputeThinV);
	VectorXd X = svd.solve(B);
	return X;
}

void setFormerMat(MatrixXd& mat, int idx, double x, double y, double z)
{
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

void setFormerMatQuad(MatrixXd& mat, int idx, double x, double y, double z)
{
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

}

void setupWeightMatrix(MatrixXd& wMat, double* denominator, vector<double> geoPts)
{
		for(int i=0; i<geoPts.size()/3; i++){
		double x, y, z;
		x = geoPts[3*i];
		y = geoPts[3*i+1];
		z = geoPts[3*i+2];
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

void setupWeightMatrixQuad(MatrixXd& wMat, double* denominator, vector<double> geoPts)
{
		for(int i=0; i<geoPts.size()/3; i++){
		double x, y, z;
		x = geoPts[3*i];
		y = geoPts[3*i+1];
		z = geoPts[3*i+2];
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

		wMat(i, i) = 1.0 / v;
		//printf("weight: %+.5E\n", v);
	}

}

int rfmReadGeoPtsAndNormalizedCoefs(const char* strGeoFile, vector<double>& imgPts, vector<double>& geoPts, \
					RPCparas& rpcs)
{
	FILE *pF = NULL;
	if(NULL == (pF = fopen(strGeoFile, "r")) )
		return -1;

	double dCurVal = 0;

	fscanf( pF, "SAMPLE_OFF:	%lf\n", &rpcs.dXOff );
	fscanf( pF, "LINE_OFF:	%lf\n", &rpcs.dYOff );
	fscanf( pF, "LON_OFF:  %lf\n", &rpcs.dObjXOff );
	fscanf( pF, "LAT_OFF:    %lf\n", &rpcs.dObjYOff );
	fscanf( pF, "HEIGHT_OFF:    %lf\n", &rpcs.dObjZOff );

	fscanf( pF, "SAMPLE_SCALE:	%lf\n", &rpcs.dXScale );
	fscanf( pF, "LINE_SCALE:	%lf\n", &rpcs.dYScale );
	fscanf( pF, "LON_SCALE:   %lf\n", &rpcs.dObjXScale );
	fscanf( pF, "LAT_SCALE:  %lf\n", &rpcs.dObjYScale );
	fscanf( pF, "HEIGHT_SCALE:  %lf\n", &rpcs.dObjZScale );

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

	imgPts.resize(vecRes.size() / 5 * 2);
	geoPts.resize(vecRes.size() / 5 * 3);

	for ( int i = 0; i < imgPts.size()/2; ++i )
	{
		imgPts[2*i] = vecRes[5*i];
		imgPts[2*i+1] = vecRes[5*i+1];

		geoPts[3*i] = vecRes[5*i+2];
		geoPts[3*i+1] = vecRes[5*i+3];
		geoPts[3*i+2] = vecRes[5*i+4];
	}
	return 1;
}

int rfmWriteRPCs(const char* strRPCsFile, RPCparas& rpcs, RFMType rfmType)
{
	FILE* fp;

	if((fp = fopen(strRPCsFile, "w")) == NULL)
	{
		printf("failed");
		return -1;
	}

	char str[50];
	int i;

	printf("store a Corrected IKONOS RPC file.\n");
	sprintf(str, "LINE_OFF: %+#010.2f pixels\n", rpcs.dYOff);
	fputs(str, fp);

	sprintf(str, "SAMP_OFF: %+#010.2f pixels\n", rpcs.dXOff);
	fputs(str, fp);

	sprintf(str, "LAT_OFF: %+#012.8f degrees\n", rpcs.dObjYOff);
	fputs(str, fp);

	sprintf(str, "LONG_OFF: %+#013.8f degrees\n", rpcs.dObjXOff);
	fputs(str, fp);

	sprintf(str, "HEIGHT_OFF: %+#09.3f meters\n", rpcs.dObjZOff);
	fputs(str, fp);

	sprintf(str, "LINE_SCALE: %+#010.2f pixels\n", rpcs.dYScale);
	fputs(str, fp);


	sprintf(str, "SAMP_SCALE: %+#010.2f pixels\n", rpcs.dXScale);
	fputs(str, fp);

	sprintf(str, "LAT_SCALE: %+#012.8f degrees\n", rpcs.dObjYScale);
	fputs(str, fp);

	sprintf(str, "LONG_SCALE: %+#013.8f degrees\n", rpcs.dObjXScale);
	fputs(str, fp);

	sprintf(str, "HEIGHT_SCALE: %+#09.3f meters\n", rpcs.dObjZScale);
	fputs(str, fp);

	if (rfmType == Three_Order_Coefs)
	{
		for(i = 0; i < 20; i++)
		{
			sprintf(str, "LINE_NUM_COEFF_%d: %+.15E\n", i+1, rpcs.c[i]);
			fputs(str, fp);
		}
		for(i = 0; i < 20; i++)
		{
			sprintf(str, "LINE_DEN_COEFF_%d: %+.15E\n", i+1, rpcs.d[i]);
			fputs(str, fp);
		}
		for(i = 0; i < 20; i++)
		{
			sprintf(str, "SAMP_NUM_COEFF_%d: %+.15E\n", i+1, rpcs.a[i]);
			fputs(str, fp);
		}
		for(i = 0; i < 20; i++)
		{
			sprintf(str, "SAMP_DEN_COEFF_%d: %+.15E\n", i+1, rpcs.b[i]);
			fputs(str, fp);
		}
	}
	else
	{
		for(i = 0; i < 10; i++)
		{
			sprintf(str, "LINE_NUM_COEFF_%d: %+.15E\n", i+1, rpcs.c[i]);
			fputs(str, fp);
		}
		for(i = 0; i < 10; i++)
		{
			sprintf(str, "LINE_DEN_COEFF_%d: %+.15E\n", i+1, rpcs.d[i]);
			fputs(str, fp);
		}
		for(i = 0; i < 10; i++)
		{
			sprintf(str, "SAMP_NUM_COEFF_%d: %+.15E\n", i+1, rpcs.a[i]);
			fputs(str, fp);
		}
		for(i = 0; i < 10; i++)
		{
			sprintf(str, "SAMP_DEN_COEFF_%d: %+.15E\n", i+1, rpcs.b[i]);
			fputs(str, fp);
		}
	}

	return 1;
}


double extractnumber(const char* str)
{
	int count = 0;
	int k = 0;
	char value[30];

	memset(value, 0, 30);

	while(str[count] != ' ')
		count++;

	while(str[count+1] != ' ' && str[count+1] != '\n')
	{
		value[k] = str[count+1];
		k++;
		count++;
	}

	return atof(value);
}

int rfmReadRpcs(const char* strRPCsFile, RPCparas& rpcs, RFMType rfmType)
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
	rpcs.dYOff = extractnumber(str);
	fgets(str, 50, fp);
	rpcs.dXOff = extractnumber(str);
	fgets(str, 50, fp);
	rpcs.dObjYOff = extractnumber(str);

	fgets(str, 50, fp);
	rpcs.dObjXOff = extractnumber(str);
	fgets(str, 50, fp);
	rpcs.dObjZOff = extractnumber(str);
	fgets(str, 50, fp);
	rpcs.dYScale = extractnumber(str);
	fgets(str, 50, fp);
	rpcs.dXScale = extractnumber(str);
	fgets(str, 50, fp);
	rpcs.dObjYScale = extractnumber(str);
	fgets(str, 50, fp);
	rpcs.dObjXScale = extractnumber(str);
	fgets(str, 50, fp);
	rpcs.dObjZScale = extractnumber(str);

	for(i = 0; i < 20; i++)
	{
		fgets(str, 50, fp);
		rpcs.c[i] = extractnumber(str);
	}
	for(i = 0; i < 20; i++)
	{
		fgets(str, 50, fp);
		rpcs.d[i] = extractnumber(str);
	}
	for(i = 0; i < 20; i++)
	{
		fgets(str, 50, fp);
		rpcs.a[i] = extractnumber(str);
	}
	for(i = 0; i < 20; i++)
	{
		fgets(str, 50, fp);
		rpcs.b[i] = extractnumber(str);
	}


	fclose(fp);
	return 1;	
}


