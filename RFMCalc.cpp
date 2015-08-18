#include "RFMCalc.h"

int rfmCalc(vector<double> imgPts, vector<double> geoPts, \
					RPCparas& rpcs, RFMCalcPara para, RFMType rfmType)
{
	if(rfmType == Three_Order_Coefs)
	{
		if(imgPts.size() != geoPts.size)
		{
			cout << "error occured" << endl;
			return -1;
		}

		MatrixXd M_mat(imgPts.size(), 19);
		MatrixXd M_mat_cpy(imgPts.size(), 19);
		MatrixXd N_mat(imgPts.size(), 19);
		VectorXd row_M(imgPts.size()), col_M(imgPts.size());
		
		M_mat.setZero();
		M_mat_cpy.setZero();
		N_mat.setZero();
		row_M.setZero();
		col_N.setZero();

		//init first values
		for(int i=0; i<nImgPts; i++)
		{
			double cX, cY, cZ;
			double iX, iY;
			cX = (geoPts[3*i] - rpcs.dObjXOff) / rpcs.dObjXScale;
			cY = (geoPts[3*i+1] - rpcs.dObjYOff) / rpcs.dObjYScale ;
			cZ = (geoPts[3*i+2] - rpcs.dObjZOff) / rpcs.dObjZScale;
			iX = (pImgPts[2*i] - rpcs.dXOff) / rpcs.dXScale;
			iY = (pImgPts[2*i] - rpcs.dYOff) / rpcs.dYScale;

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
			for(int i = 0; i < nImgPts; i++)
				sum += pow(error_vector(i), 2.0);
			sqrt_mean = sqrt(sum/nImgPts);

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
		if(imgPts.size() != geoPts.size)
		{
			cout << "error occured" << endl;
			return -1;
		}

		MatrixXd M_mat(imgPts.size(), 19);
		MatrixXd M_mat_cpy(imgPts.size(), 19);
		MatrixXd N_mat(imgPts.size(), 19);
		VectorXd row_M(imgPts.size()), col_M(imgPts.size());
		
		M_mat.setZero();
		M_mat_cpy.setZero();
		N_mat.setZero();
		row_M.setZero();
		col_N.setZero();

		//init first values
		for(int i=0; i<nImgPts; i++)
		{
			double cX, cY, cZ;
			double iX, iY;
			cX = (geoPts[3*i] - rpcs.dObjXOff) / rpcs.dObjXScale;
			cY = (geoPts[3*i+1] - rpcs.dObjYOff) / rpcs.dObjYScale ;
			cZ = (geoPts[3*i+2] - rpcs.dObjZOff) / rpcs.dObjZScale;
			iX = (pImgPts[2*i] - rpcs.dXOff) / rpcs.dXScale;
			iY = (pImgPts[2*i] - rpcs.dYOff) / rpcs.dYScale;

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
		for(int i = 0; i < nImgPts; i++)
			sum += pow(error_vector(i), 2.0);
		first_sqrt_mean = sqrt(sum/nImgPts);
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
			for(int i = 0; i < nImgPts; i++)
				sum += pow(error_vector(i), 2.0);
			sqrt_mean = sqrt(sum/nImgPts);

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
	if ( imgPts.size() != geoPts.size() )
	{
		return -1;
	}
	vector<double> vecImgX, vecImgY, vecObjX, vecObjY, vecObjZ;
	for ( int i = 0; i < nImgPts; ++i )
	{
		rpcs.dXOff += pImgPts[2*i];
		rpcs.dYOff += pImgPts[2*i+1];
		rpcs.dObjXOff += pGeoPts[3*i];
		rpcs.dObjYOff += pGeoPts[3*i+1];
		rpcs.dObjZOff += pGeoPts[3*i+2];

		vecImgX.push_back( pImgPts[2*i] );
		vecImgY.push_back( pImgPts[2*i+1] );

		vecObjX.push_back( pGeoPts[3*i] );
		vecObjY.push_back( pGeoPts[3*i+1] );
		vecObjZ.push_back( pGeoPts[3*i+2] );
	}
	rpcs.dXOff /= nImgPts;
	rpcs.dYOff /= nImgPts;
	rpcs.dObjXOff /= nGeoPts;
	rpcs.dObjYOff /= nGeoPts;
	rpcs.dObjZOff /= nGeoPts;

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
		mat(idx, 6) = y*x;f

		mat(idx, 7) = z*z;
		mat(idx, 8) = y*y;
		mat(idx, 9) = x*x;

}

void setupWeightMatrix(MatrixXd& wMat, double* denominator, vector<double> geoPts)
{
		for(int i=0; i<nGeoPts; i++){
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
		for(int i=0; i<nGeoPts; i++){
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

