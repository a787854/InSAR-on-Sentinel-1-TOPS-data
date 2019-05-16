#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#define __CUDA_INTERNAL_COMPILATION__
#include <math_functions.h>
#include <math_constants.h>
#include <device_functions.h>
#include <cublas_v2.h>
#include "cusolverDn.h"
//#include "helper_cusolver.h"
#include "cuda_occupancy.h"
#include "cuConstants.cuh"
#define __CUDA_INTERNAL_COMPILATION__
#include <math.h>
#include <iostream>
#include<time.h>
#include <stdio.h>


using namespace std;


extern "C"	void GeoCoreg_warp_cuda(
	int m_burstIdx,//burst index
	int s_burstIdx,
	int xmin,
	int xmax,
	int ymin,
	int ymax,
	double dem_DeltaLat,
	double dem_DeltaLon,
	int DemLines,
	int DemPixels,
	double *AzCpm,
	double *RgCpm,
	double lat0,
	double lon0,
	int demnodata,
	double m_azimuthTimeInterval,
	double s_azimuthTimeInterval,
	double wavelength,
	double m_linesPerBurst,
	double s_linesPerBurst,
	double m_burstFirstLineTime,
	double s_burstFirstLineTime,
	int m_numOfSamples,
	int s_numOfSamples,
	double *m_orbitAzTime,
	double *s_orbitAzTime,
	int    orb_m_numberofpoints,
	int    orb_m_numbercoefs,
	int    orb_s_numberofpoints,
	int    orb_s_numbercoefs,
	double m_slrTimeToFirstPixel,
	double s_slrTimeToFirstPixel,
	double m_rangeSpacing,
	double s_rangeSpacing,
	bool m_nearRangeOnLeft,
	bool s_nearRangeOnLeft,
	int* inputDEM,
	double *m_x_coef,
	double *m_y_coef,
	double *m_z_coef,
	double *s_x_coef,
	double *s_y_coef,
	double *s_z_coef,
	double m_aziInitial,
	double s_aziInitial
	);



__device__ inline double sqr(const double x)
{
	return x*x;
}



__device__ inline double Distance(double* Vec3D_A,
	double* Vec3D_B)
{
	return sqrt(sqr(Vec3D_A[0] - Vec3D_B[0]) +
		sqr((Vec3D_A[1] - Vec3D_B[1])) +
		sqr(Vec3D_A[2] - Vec3D_B[2]));
}

__device__ inline void llh2xyz(double phi, double lambda, float height, double Res[3])
{
	double sincosLam[2];
	double sincosPhi[2];
	sincos(lambda, sincosLam, sincosLam + 1);
	sincos(phi, sincosPhi, sincosPhi + 1);

	double Nph = b_ell_a / sqrt(1.0 - b_ell_e2*sqr(sincosPhi[0])) + height;
	Res[0] = Nph * sincosPhi[1] * sincosLam[1];
	Res[1] = Nph * sincosPhi[1] * sincosLam[0];
	Res[2] = (Nph - b_ell_e2*(Nph - height))   *sincosPhi[0];

}

__device__ inline double SolveNewtonEq(double postar[3], double possat[3], double velsat[3], double accsat[3])
{
	possat[0] = postar[0] - possat[0];
	possat[1] = postar[1] - possat[1];
	possat[2] = postar[2] - possat[2];


	return( -(velsat[0] * possat[0] + velsat[1] * possat[1] + velsat[2] * possat[2]) /
		((accsat[0] * possat[0] + accsat[1] * possat[1] + accsat[2] * possat[2])
		- velsat[0] * velsat[0] - velsat[1] * velsat[1] - velsat[2] * velsat[2]));

}



inline __device__ __host__  double deg2rad(const double x)
{
	return x * CUDART_PI / 180.0;
}



inline __device__ __host__ void getPos(double t, double* possat, double* coef_x,
	double *coef_y, double* coef_z, int n)
{

	possat[0] = 0.0;
	possat[1] = 0.0;
	possat[2] = 0.0;
	for (int i = n - 1; i >= 0; --i)
	{
		possat[0] *= t;
		possat[0] += coef_x[i];

		possat[1] *= t;
		possat[1] += coef_y[i];

		possat[2] *= t;
		possat[2] += coef_z[i];

	}


}

inline __device__ __host__ void getVelo(double t, double* velsat, double* coef_x,
	double *coef_y, double* coef_z, int n)
{
	velsat[0] = coef_x[1];
	velsat[1] = coef_y[1];
	velsat[2] = coef_z[1];


	double powt;


	for (int i = 2; i < n; ++i)
	{
		powt = double(i)*pow(t, i - 1);


		velsat[0] += coef_x[i] * powt;
		velsat[1] += coef_y[i] * powt;
		velsat[2] += coef_z[i] * powt;

	}

	velsat[0] /= 10.0;
	velsat[1] /= 10.0;
	velsat[2] /= 10.0;



}

inline __device__ __host__ void getAcc(double t, double*accsat, double* coef_x,
	double *coef_y, double* coef_z, int n)
{

	
	accsat[0] = 0.0;
	accsat[1] = 0.0;
	accsat[2] = 0.0;


	double powt1;

	for (int i = 2; i < n; ++i)
	{

		powt1 = double((i - 1)*i)*pow(t, i - 2);


		accsat[0] += coef_x[i] * powt1;
		accsat[1] += coef_y[i] * powt1;
		accsat[2] += coef_z[i] * powt1;


	}

	accsat[0] /= 100.0;
	accsat[1] /= 100.0;
	accsat[2] /= 100.0;

}





__global__ void ell2lpGeo
(   /*input*/
int *d_DEM,
int    orb_numberofpoints,
int    orb_numberofcoefs,
double   midtime,
double *  orb_coef_x,
double *  orb_coef_y,
double *  orb_coef_z,
double *Linebuffer,
double *Pixelbuffer,
double azimuthTimeInterval,
double firstLineTime,
double slrTimeToFirstPixel,
double     aziInitial,
int demNoData,
bool nearRangeOnLeft,
int LineOffset,
double rangeSpacing
)
{
	unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;

	if (tid < b_nrows*b_ncols)
	{
		unsigned int row = tid / b_ncols;
		unsigned int col = tid%b_ncols;

		int ncoefs = orb_numberofcoefs;


		



		float height = d_DEM[tid];

		

		if (height != demNoData)
		{


			double targetpos[3];
			llh2xyz(deg2rad(b_lat0 - b_DEMdeltalat*row), deg2rad(b_lon0 + b_DEMdeltalon*col),
			height, targetpos);


			double sol = 0.0;
			double azimuth_tmp = aziInitial;

			double possat[3];
			double velsat[3];
			double accsat[3];



			double midTime = midtime;


			double t_tmp;
			for (int i = 0; i <= 10; i++)
			{

				t_tmp = (azimuth_tmp - midTime) / 10.0;
				getPos(t_tmp, possat, orb_coef_x, orb_coef_y, orb_coef_z, ncoefs);

				getVelo(t_tmp, velsat, orb_coef_x, orb_coef_y, orb_coef_z, ncoefs);

				getAcc(t_tmp, accsat, orb_coef_x, orb_coef_y, orb_coef_z, ncoefs);

				sol=SolveNewtonEq(targetpos, possat, velsat, accsat);

				
		


				azimuth_tmp += sol;

				if (abs(sol) < 1e-010)
					break;


			}

			t_tmp = (azimuth_tmp - midTime) / 10.0;

			getPos(t_tmp, possat, orb_coef_x, orb_coef_y, orb_coef_z, ncoefs);


	
			Linebuffer[tid] = LineOffset + (azimuth_tmp - firstLineTime) / azimuthTimeInterval;

			double range_tmp = sqrt(sqr(targetpos[0] - possat[0]) + 
				sqr(targetpos[1] - possat[1]) + sqr(targetpos[2] - possat[2]));

			Pixelbuffer[tid] = (range_tmp - slrTimeToFirstPixel*SOL_cuda) / rangeSpacing;




		}
		else
		{
			Linebuffer[tid] = -99999999;
			Pixelbuffer[tid] = -99999999;


		}

	}

}

__device__  inline double getDopplerFreq(double *earthPos, double* sensorPos, double* sensorVel, double wavelength)
{
	double xDiff = earthPos[0] - sensorPos[0];
	double yDiff = earthPos[1] - sensorPos[1];
	double zDiff = earthPos[2] - sensorPos[2];
	double distance = sqrt(xDiff * xDiff + yDiff * yDiff + zDiff * zDiff);

	return 2.0 * (sensorVel[0] * xDiff + sensorVel[1] * yDiff + sensorVel[2] * zDiff) / (distance * wavelength);

}

__global__ void ell2lpGeo_BinarySearch
(   /*input*/
int *d_DEM,
int    orb_numberofpoints,
int    orb_numberofcoefs,
double   midtime,
double *  orb_coef_x,
double *  orb_coef_y,
double *  orb_coef_z,
double *Linebuffer,
double *Pixelbuffer,
double azimuthTimeInterval,
double firstLineTime,
double slrTimeToFirstPixel,
double     aziInitial,
int demNoData,
bool nearRangeOnLeft,
int LineOffset,
double rangeSpacing,
double waveLength,
double* azimuthTime,
int numberofAzimuthTime
//int *NI
)
{

	unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;

	if (tid < b_nrows*b_ncols)
	{
		unsigned int row = tid / b_ncols;
		unsigned int col = tid%b_ncols;

		int ncoefs = orb_numberofcoefs;


		double phi = deg2rad(b_lat0 - b_DEMdeltalat*row);
		double lambda = deg2rad(b_lon0 + b_DEMdeltalon*col);



		double height = d_DEM[tid];



		if (height != demNoData)
		{

			
			double earthPoint[3];
			llh2xyz(deg2rad(b_lat0 - b_DEMdeltalat*row), deg2rad(b_lon0 + b_DEMdeltalon*col),
				height, earthPoint);

			double sol = 0.0;
			//double azimuth_tmp = aziInitial;//rowd_azi[col];

			double possat[3];
			double velsat[3];
			//double satacc[3];


			double delta_x;
			double delta_y;
			double delta_z;

			double firstVecTime = 0.0;
			double secondVecTime = 0.0;
			double firstVecFreq = 0.0;
			double secondVecFreq = 0.0;

			double midAziTime = midtime;

			double t_tmp = (azimuthTime[0] - midAziTime) / 10.0;

			getPos(t_tmp, possat, orb_coef_x, orb_coef_y, orb_coef_z, ncoefs);
			getVelo(t_tmp, velsat, orb_coef_x, orb_coef_y, orb_coef_z, ncoefs);

			double currentFreq = getDopplerFreq(earthPoint, possat, velsat, waveLength);


			firstVecFreq = currentFreq;
			firstVecTime = azimuthTime[0];

			for (int i = 1; i < numberofAzimuthTime; i++)
			{
				t_tmp = (azimuthTime[i] - midAziTime) / 10.0;
				getPos(t_tmp, possat, orb_coef_x, orb_coef_y, orb_coef_z, ncoefs);
				getVelo(t_tmp, velsat, orb_coef_x, orb_coef_y, orb_coef_z, ncoefs);

				currentFreq = getDopplerFreq(earthPoint, possat, velsat, waveLength);

				if (firstVecFreq*currentFreq > 0)
				{
					firstVecTime = azimuthTime[i];
					firstVecFreq = currentFreq;

				}
				else
				{
					secondVecTime = azimuthTime[i];
					secondVecFreq = currentFreq;
					break;

				}

			}

			


			if (firstVecFreq * secondVecFreq < 0.0)
			{
				double lowerBoundTime = firstVecTime;
				double upperBoundTime = secondVecTime;
				double lowerBoundFreq = firstVecFreq;
				double upperBoundFreq = secondVecFreq;
				double midTime = 0.0;
				double midFreq = 0.0;
				int nv = 8;
				int i0, iN;
				double weight;
				double time2;



				double diffTime = abs(upperBoundTime - lowerBoundTime);

				double absLineTimeInterval = abs(azimuthTimeInterval / 10);
				int totalIterations = (int)(diffTime / absLineTimeInterval);
				int numIterations = 0;

				while (diffTime > absLineTimeInterval&&numIterations <= totalIterations)
				{

					possat[0] = possat[1] = possat[2] = 0.0;
					velsat[0] = velsat[1] = velsat[2] = 0.0;
					midTime = (upperBoundTime + lowerBoundTime) / 2.0;

					t_tmp = (midTime - midAziTime) / 10.0;

					getPos(t_tmp, possat, orb_coef_x, orb_coef_y, orb_coef_z, ncoefs);
					getVelo(t_tmp, velsat, orb_coef_x, orb_coef_y, orb_coef_z, ncoefs);


					midFreq = getDopplerFreq(earthPoint, possat, velsat,
						waveLength);



					if (midFreq * lowerBoundFreq > 0.0)
					{
						lowerBoundTime = midTime;
						lowerBoundFreq = midFreq;
					}
					else if (midFreq * upperBoundFreq > 0.0)
					{
						upperBoundTime = midTime;
						upperBoundFreq = midFreq;
					}
					else if (midFreq == 0.0)
					{
						break;
					}

					diffTime = abs(upperBoundTime - lowerBoundTime);
					numIterations++;


				}

				midTime = lowerBoundTime - lowerBoundFreq * (upperBoundTime - lowerBoundTime) / (upperBoundFreq - lowerBoundFreq);

				t_tmp = (midTime - midAziTime) / 10.0;
				getPos(t_tmp, possat, orb_coef_x, orb_coef_y, orb_coef_z, ncoefs);


				delta_x = earthPoint[0] - possat[0];
				delta_y = earthPoint[1] - possat[1];
				delta_z = earthPoint[2] - possat[2];

				Linebuffer[tid] = LineOffset + (midTime - firstLineTime) / azimuthTimeInterval;
				//NI[tid] = numIterations;

				


				double range_tmp = sqrt(delta_x*delta_x + delta_y*delta_y + delta_z*delta_z);

				Pixelbuffer[tid] = (range_tmp - slrTimeToFirstPixel*SOL_cuda) / rangeSpacing;


			}



		}
		else
		{
			Linebuffer[tid] = -99999999;
			Pixelbuffer[tid] = -99999999;


		}

	}

}





__global__ void WeighthSet
(   /*input*/
int *d_DEM,
double *Linebuffer,
double *Pixelbuffer,
int demNoData,
double* weight,
int xmin,
int xmax,
int ymin,
int ymax
)
{

	int tid = blockIdx.x*blockDim.x + threadIdx.x;
	if (tid<b_ncols*b_nrows)
	{
		int height = d_DEM[tid];
		if (height != demNoData && Linebuffer[tid] >= ymin && Linebuffer[tid] <= ymax &&
			Pixelbuffer[tid] >= xmin && Pixelbuffer[tid] <= xmax)
		{
			weight[tid] = 1.0;
		}

	}

}

__global__ void AMatrixSet
(   /*input*/
double *Linebuffer,
double *Pixelbuffer,
double *A,
int xmin,
int xmax,
int ymin,
int ymax
)
{
	int tid = blockIdx.x*blockDim.x + threadIdx.x;

	if (tid < b_nrows* b_ncols)
	{
		double indexL = normalizeWarp(Linebuffer[tid], ymin, ymax);
		double indexP = normalizeWarp(Pixelbuffer[tid], xmin, xmax);


		int offset = tid * 6;

		A[offset] = 1;
		A[offset + 1] = indexL;
		A[offset + 2] = indexP;
		A[offset + 3] = indexL*indexL;
		A[offset + 4] = indexL*indexP;
		A[offset + 5] = indexP*indexP;
	}

}




void GeoCoreg_warp_cuda(
	int m_burstIdx,//burst index
	int s_burstIdx,
	int xmin,
	int xmax,
	int ymin,
	int ymax,
	double dem_DeltaLat,
	double dem_DeltaLon,
	int DemLines,
	int DemPixels,
	double *AzCpm,
	double *RgCpm,
	double lat0,
	double lon0,
	int demnodata,
	double m_azimuthTimeInterval,
	double s_azimuthTimeInterval,
	double wavelength,
	double m_linesPerBurst,
	double s_linesPerBurst,
	double m_burstFirstLineTime,
	double s_burstFirstLineTime,
	int m_numOfSamples,
	int s_numOfSamples,
	double *m_orbitAzTime,
	double *s_orbitAzTime,
	int    orb_m_numberofpoints,
	int    orb_m_numbercoefs,
	int    orb_s_numberofpoints,
	int    orb_s_numbercoefs,
	double m_slrTimeToFirstPixel,
	double s_slrTimeToFirstPixel,
	double m_rangeSpacing,
	double s_rangeSpacing,
	bool m_nearRangeOnLeft,
	bool s_nearRangeOnLeft,
	int* inputDEM,
	double *m_x_coef,
	double *m_y_coef,
	double *m_z_coef,
	double *s_x_coef,
	double *s_y_coef,
	double *s_z_coef,
	double m_aziInitial,
	double s_aziInitial
	)
{

	cublasHandle_t Cublashandle;
	cublasStatus_t status = cublasCreate(&Cublashandle);
	cusolverDnHandle_t CuSolverhandle = NULL;
	cusolverDnCreate(&CuSolverhandle);

	double alpha = 1.0f;
	double beta = 0.0f;
	double b = 6356752.3142451794975639665996337;
	double a = 6378137.0;
	double e2 = 1.0 - (b*b / a / a);

	double s_min4picdivlam = (-4.0*CUDART_PI*b_SOL_cuda) / wavelength;
	int m_lineOffset = m_burstIdx*m_linesPerBurst;
	int s_lineOffset = s_burstIdx*s_linesPerBurst;
	
	
	int numPixels = DemPixels;
	
	int numLines = DemLines;
	int offset = numLines*numPixels;
	double m_midTime = m_orbitAzTime[orb_m_numberofpoints / 2];
	double s_midTime = s_orbitAzTime[orb_s_numberofpoints / 2];



	cudaHostRegister(inputDEM, numLines*numPixels*sizeof(short), cudaHostRegisterDefault);
	cudaHostRegister(AzCpm, 6 * sizeof(double), cudaHostRegisterDefault);
	cudaHostRegister(RgCpm, 6 * sizeof(double), cudaHostRegisterDefault);




	int *d_inputDEM;
	size_t d_pitch;
	cudaMalloc((void**)&d_inputDEM, numPixels*numLines*sizeof(int));

	
	double *d_masterAZ, *d_masterRG, *d_slaveAZ, *d_slaveRG;
	cudaMalloc((void**)&d_masterAZ, numLines*numPixels*sizeof(double));
	cudaMalloc((void**)&d_masterRG, numLines*numPixels*sizeof(double));
	cudaMalloc((void**)&d_slaveAZ, numLines*numPixels*sizeof(double));
	cudaMalloc((void**)&d_slaveRG, numLines*numPixels*sizeof(double));

	bool BinaryOn = false;
	double * d_m_azimuthTime, *d_s_azimuthTime;
	if (BinaryOn)
	{
		
		cudaMalloc((void**)&d_m_azimuthTime, orb_m_numberofpoints*sizeof(double));
		cudaMalloc((void**)&d_s_azimuthTime, orb_s_numberofpoints*sizeof(double));

		cudaMemcpy(d_m_azimuthTime, m_orbitAzTime, orb_m_numberofpoints*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(d_s_azimuthTime, s_orbitAzTime, orb_s_numberofpoints*sizeof(double), cudaMemcpyHostToDevice);
	}

	double  *d_m_x_coef, *d_m_y_coef, *d_m_z_coef;
	double  *d_s_x_coef, *d_s_y_coef, *d_s_z_coef;

	cudaMalloc((void**)&d_s_x_coef,  orb_s_numbercoefs*sizeof(double));
	cudaMalloc((void**)&d_s_y_coef,  orb_s_numbercoefs*sizeof(double));
	cudaMalloc((void**)&d_s_z_coef,  orb_s_numbercoefs*sizeof(double));



	cudaMalloc((void**)&d_m_x_coef,  orb_m_numbercoefs*sizeof(double));
	cudaMalloc((void**)&d_m_y_coef,  orb_m_numbercoefs*sizeof(double));
	cudaMalloc((void**)&d_m_z_coef,  orb_m_numbercoefs*sizeof(double));



	double *d_weight = NULL;
	double *d_A = NULL;
	double *d_WeightA = NULL;
	double *d_N = NULL;
	double *d_AzCpm = NULL;
	double *d_RgCpm = NULL;
	cudaMalloc((void**)&d_weight, numLines*numPixels*sizeof(double));
	cudaMalloc((void**)&d_A, numLines*numPixels * 6 * sizeof(double));
	cudaMalloc((void**)&d_WeightA, numLines*numPixels * 6 * sizeof(double));
	cudaMalloc((void**)&d_N, 36 * sizeof(double));
	cudaMalloc((void**)&d_AzCpm, 12 * sizeof(double));
	d_RgCpm = d_AzCpm + 6;


	//for (int i = 0; i < 10; i++)

	cudaEvent_t g_start, g_stop;
	float time_cost;
	cudaEventCreate(&g_start);
	cudaEventCreate(&g_stop);
	cudaEventRecord(g_start, 0);

	cudaMemset(d_AzCpm, 0, 6 * sizeof(double));
	cudaMemset(d_RgCpm, 0, 6 * sizeof(double));
	cudaMemset(d_weight, 0, numLines*numPixels*sizeof(double));


	cudaMemcpyToSymbol(b_ell_e2, &e2, sizeof(double), 0, cudaMemcpyHostToDevice);

	unsigned int threads = 256;
	unsigned int blocksNum = (numLines*numPixels + 255) / 256;




	cudaMemcpyToSymbol(b_nrows, &numLines, sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(b_ncols, &numPixels, sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(b_lat0, &lat0, sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(b_lon0, &lon0, sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(b_CorrPitch, &d_pitch, sizeof(size_t), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(b_DEMdeltalat, &dem_DeltaLat, sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(b_DEMdeltalon, &dem_DeltaLon, sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(b_s_min4picdivlam, &s_min4picdivlam, sizeof(double), 0, cudaMemcpyHostToDevice);



	

	cudaMemcpy(d_m_x_coef, m_x_coef, orb_m_numbercoefs*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_m_y_coef, m_y_coef, orb_m_numbercoefs*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_m_z_coef, m_z_coef, orb_m_numbercoefs*sizeof(double), cudaMemcpyHostToDevice);


	cudaMemcpy(d_s_x_coef, s_x_coef, orb_s_numbercoefs*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_s_y_coef, s_y_coef, orb_s_numbercoefs*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_s_z_coef, s_z_coef, orb_s_numbercoefs*sizeof(double), cudaMemcpyHostToDevice);



	cudaMemcpy(d_inputDEM, inputDEM, numLines*numPixels*sizeof(short), cudaMemcpyHostToDevice);



	if (BinaryOn)
	{
		ell2lpGeo_BinarySearch << <blocksNum, threads >> >(d_inputDEM, orb_m_numberofpoints, orb_m_numbercoefs, m_midTime,
			d_m_x_coef, d_m_y_coef, d_m_z_coef, d_masterAZ, d_masterRG, m_azimuthTimeInterval, m_burstFirstLineTime, m_slrTimeToFirstPixel,
			m_aziInitial, demnodata, m_nearRangeOnLeft, m_lineOffset, m_rangeSpacing, wavelength, d_m_azimuthTime, orb_m_numberofpoints);

		ell2lpGeo_BinarySearch << <blocksNum, threads >> >(d_inputDEM, orb_s_numberofpoints, orb_s_numbercoefs, s_midTime,
			d_s_x_coef, d_s_y_coef, d_s_z_coef, d_slaveAZ, d_slaveRG, s_azimuthTimeInterval, s_burstFirstLineTime, s_slrTimeToFirstPixel,
			s_aziInitial, demnodata, s_nearRangeOnLeft, s_lineOffset, s_rangeSpacing, wavelength, d_s_azimuthTime, orb_s_numberofpoints);
	}

	ell2lpGeo << <blocksNum, threads >> >(d_inputDEM, orb_m_numberofpoints, orb_m_numbercoefs, m_midTime,
		d_m_x_coef, d_m_y_coef, d_m_z_coef, d_masterAZ, d_masterRG, m_azimuthTimeInterval, m_burstFirstLineTime, m_slrTimeToFirstPixel,
		m_aziInitial, demnodata, m_nearRangeOnLeft, m_lineOffset, m_rangeSpacing);


	ell2lpGeo << <blocksNum, threads >> >(d_inputDEM, orb_s_numberofpoints, orb_s_numbercoefs, s_midTime,
		d_s_x_coef, d_s_y_coef, d_s_z_coef, d_slaveAZ, d_slaveRG, s_azimuthTimeInterval, s_burstFirstLineTime, s_slrTimeToFirstPixel,
		s_aziInitial, demnodata, s_nearRangeOnLeft, s_lineOffset, s_rangeSpacing);

	WeighthSet << <blocksNum, threads >> > (d_inputDEM, d_masterAZ, d_masterRG, demnodata, d_weight, xmin, xmax,
		ymin, ymax);

	AMatrixSet << <blocksNum, threads >> > (d_masterAZ, d_masterRG, d_A, xmin, xmax, ymin, ymax);



	int SizeN = numPixels*numLines;


	cublasDdgmm(Cublashandle, CUBLAS_SIDE_RIGHT, 6, SizeN, d_A, 6, d_weight, 1, d_WeightA, 6);

	

	 cublasDgemm(Cublashandle, CUBLAS_OP_N, CUBLAS_OP_T, 6, 6, SizeN, &alpha, d_WeightA, 6, d_A, 6, &beta, d_N, 6);
	

	status = cublasDgemm(Cublashandle, CUBLAS_OP_N, CUBLAS_OP_T, 6, 1, SizeN, &alpha, d_WeightA, 6, d_slaveAZ, 1, &beta, d_AzCpm, 6);
	status = cublasDgemm(Cublashandle, CUBLAS_OP_N, CUBLAS_OP_T, 6, 1, SizeN, &alpha, d_WeightA, 6, d_slaveRG, 1, &beta, d_RgCpm, 6);


	int bufferSize = 0;
	int *info = NULL;
	double *buffer = NULL;
	double *A = NULL;
	int h_info = 0;
	cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;
	cusolverStatus_t output;
	output = cusolverDnDpotrf_bufferSize(CuSolverhandle, uplo, 6, d_N, 6, &bufferSize);
	cudaMalloc(&info, sizeof(int));
	cudaMalloc(&buffer, sizeof(double)*bufferSize);
	h_info = cusolverDnDpotrf(CuSolverhandle, uplo, 6, d_N, 6, buffer, bufferSize, info);
	//checkCudaErrors(cudaMemcpy(&h_info, info, sizeof(int), cudaMemcpyDeviceToHost));
	if ( h_info !=0 ){
		fprintf(stderr, "Error: Cholesky factorization failed\n");
	}
	else
	{
		cusolverDnDpotrs(CuSolverhandle, uplo, 6, 1, d_N, 6, d_AzCpm, 6, info);
		cusolverDnDpotrs(CuSolverhandle, uplo, 6, 1, d_N, 6, d_RgCpm, 6, info);

		cudaMemcpy(AzCpm, d_AzCpm, 6 * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(RgCpm, d_RgCpm, 6 * sizeof(double), cudaMemcpyDeviceToHost);
	}
	
	cudaEventRecord(g_stop, 0);
	cudaDeviceSynchronize();
	cudaEventElapsedTime(&time_cost, g_start, g_stop);
	cout << " Geometric Coregistration Elapsed Time :" << time_cost << endl;


	//ofstream GTout("E:\\GTtimer.txt", ios::app);
	//GTout << time_cost << endl;
	//GTout.close();
	cublasDestroy(Cublashandle);
	cusolverDnDestroy(CuSolverhandle);
	cudaFree(buffer);
	cudaFree(info);
	




	cudaHostUnregister(AzCpm);
	cudaHostUnregister(RgCpm);
	cudaHostUnregister(inputDEM);
	cudaFree(d_inputDEM);
	cudaFree(d_masterAZ);
	cudaFree(d_masterRG);
	cudaFree(d_slaveAZ);
	cudaFree(d_slaveRG);
	cudaFree(d_m_x_coef);
	cudaFree(d_m_y_coef);
	cudaFree(d_m_z_coef);
	cudaFree(d_s_x_coef);
	cudaFree(d_s_y_coef);
	cudaFree(d_s_z_coef);
	cudaFree(d_A);
	cudaFree(d_WeightA);
	cudaFree(d_N);
	cudaFree(d_AzCpm);
	if (BinaryOn)
	{
		cudaFree(d_m_azimuthTime);
		cudaFree(d_s_azimuthTime);
	}

	cudaDeviceReset();




}




