
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#define __CUDA_INTERNAL_COMPILATION__
#include <math_functions.h>
#include <math_constants.h>
#include <device_functions.h>
#include <cuComplex.h>
//#include <helper_cuda.h>
#include <cuda_occupancy.h>
#include "cuConstants.cuh"
#include "cuda_profiler_api.h"
#undef  __CUDA_INTERNAL_COMPILATION__
#include <math.h>
#include<time.h>
#include<complex>
#include <iostream>





#include <stdio.h>

using namespace std;

extern "C" void DerampDemodResample(
	complex<short>*SlaveArray,
	double *CpmAz,
	double *CpmRg,
	double AzimuthShift,
	complex<float>* output,
	float *KernelAz,
	float *KernelRg,
	int sBurstIdx,
	int slave_pixels,
	int slave_lines,
	int MasterBox[4],
	int SlaveBox[4],
	int S_linesPerBurst,
	int S_SamplesPerBurst,
	double azimuthTimeInterval,
	double* dopplerRate,
	double* referenceTime,
	double* dopplerCentroid,
	int Npoints
	);

extern "C" cuComplex* DerampDemodResample_ESD(
	complex<short>*SlaveArray,
	double *CpmAz,
	double *CpmRg,
	double AzimuthShift,
	complex<float>* output,
	float *KernelAz,
	float *KernelRg,
	int sBurstIdx,
	int slave_pixels,
	int slave_lines,
	int MasterBox[4],
	int SlaveBox[4],
	int linesPerBurst,
	int SamplesPerBurst,
	double azimuthTimeInterval,
	double* dopplerRate,
	double* referenceTime,
	double* dopplerCentroid,
	int Npoints
	);

extern "C" cuComplex* ResampleFirstBurst(
	complex<float>*SlaveArray,
	int ww,
	int hh
	);

__device__  cuComplex PartialCompute(cuComplex* Input, cuComplex *kernel)
{
	cuComplex tempAdd = make_cuComplex(0.0f, 0.0f);
	cuComplex tempMul;

	for (int i = 0; i < npoints; i++)
	{
		tempMul = cuCmulf(Input[i], kernel[i]);
		tempAdd = cuCaddf(tempAdd, tempMul);
	}

	return tempAdd;

}

__device__ inline double d_sqr(double x)
{
	return x*x;
}
__device__  inline cuComplex CmulfFloat(cuComplex Left, float Right)
{
	cuComplex Res;


	Res.x = Left.x*Right;
	Res.y = Left.y*Right;

	return Res;

}

__device__ inline void GetIndexes(double x, int Indexed[2])
{
	Indexed[0] = trunc(x);
	Indexed[1] = ((x - Indexed[0]) * c_Interval + 0.5);
}




__global__ void DerampDemod_Shared(float *PhaseArray, cuComplex* SlaveArray,
	double* dopplerRate, double *referenceTime, double* dopplerCentroid,
	int x0, int y0, int numLines, int numPixels, int firstLineInBurst,
	double azimuthTimeInterval, size_t d_pitch1,
	size_t d_pitch2, size_t d_pitch3, short2* SlaveArrayS2, size_t d_pitchS2)
{
	__shared__ double tile_Doppler[16];
	__shared__ double tile_ReferenceTime[16];
	__shared__ double tile_DopplerCentroid[16];


	const int row = blockIdx.y*blockDim.y + threadIdx.y;
	const int col = blockIdx.x*blockDim.x + threadIdx.x;

	int y = row + y0;
	int x = col + x0;


	if (threadIdx.y == 0)
	{

		tile_Doppler[threadIdx.x] = dopplerRate[x];
		tile_ReferenceTime[threadIdx.x] = referenceTime[x];
		tile_DopplerCentroid[threadIdx.x] = dopplerCentroid[x];

	}



	__syncthreads();

	if (row < numLines&&col < numPixels)
	{

		double ta = (y - firstLineInBurst)*azimuthTimeInterval;



		double kt = tile_Doppler[threadIdx.x];
		double deramp = -CUDART_PI*kt*d_sqr(ta - tile_ReferenceTime[threadIdx.x]);
		double demod = -2.0 * CUDART_PI*tile_DopplerCentroid[threadIdx.x] * ta;
		double phase = deramp + demod;




		float* rowPhaseArray = (float *)((char*)PhaseArray + row*d_pitch1);
		rowPhaseArray[col] = phase;


		cuComplex* rowSlaveArrayOutput = (cuComplex *)((char*)SlaveArray + row*d_pitch2);
		short2* rowSlaveArrayInput = (short2 *)((char*)SlaveArrayS2 + row*d_pitchS2);

		double SinCos[2];
		sincos(phase, SinCos, SinCos + 1);


		rowSlaveArrayOutput[col] = cuCmulf(make_cuComplex(rowSlaveArrayInput[col].x, rowSlaveArrayInput[col].y),
			make_cuComplex(SinCos[1], SinCos[0]));

		/*if (row == 50 && col == 2000)
		{
			printf("Deramped Slave:(%lf,%lf)\n", rowSlaveArrayOutput[col].x, rowSlaveArrayOutput[col].y);
		}*/
	}

}


__global__ void
resample_texture_kernel_6p(cuComplex *output, double* SlaveAz,
double* SlaveRg, size_t CorrPitch_1)
{
	__shared__ double SslaveRg[16][16];
	__shared__ double SslaveAz[16][16];

	const int row = blockIdx.y*blockDim.y + threadIdx.y;
	const int col = blockIdx.x*blockDim.x + threadIdx.x;


	//const int index = row*(overlap_pixehi - overlap_pixelo + 1) + col;
	cuComplex *rowoutput = (cuComplex *)((char*)output + row*CorrPitch_1);


	// adjust pCorr to point to row
	//output = (float *)((char*)output + row*CorrPitch);

	int Npoints2m1 = 6 / 2 - 1;



	if (row < c_mLines && col < c_mPixels)
	{


		double *rowSlaveRg = (double *)((char*)SlaveRg + row*CorrPitch_1);
		double *rowSlaveAz = (double *)((char*)SlaveAz + row*CorrPitch_1);


		SslaveRg[threadIdx.y][threadIdx.x] = rowSlaveRg[col];
		SslaveAz[threadIdx.y][threadIdx.x] = rowSlaveAz[col] + d_AzimuthShift;

		__syncthreads();


		const int  fl_interpL = int(SslaveAz[threadIdx.y][threadIdx.x]);
		const int fl_interpP = int(SslaveRg[threadIdx.y][threadIdx.x]);




		if (fl_interpL > c_sYmax || fl_interpL<c_sY0
			|| fl_interpP>c_sXmax || fl_interpP < c_sX0)
		{
			rowoutput[col] = make_cuComplex(0.0f, 0.0f);

		}

		else if (fl_interpL>c_sYmax - Npoints2m1 - 2 || fl_interpL<c_sY0 + Npoints2m1
			|| fl_interpP>c_sXmax - Npoints2m1 - 2 || fl_interpP < c_sX0 + Npoints2m1)
		{
			int indexL = fl_interpL - c_sY0;
			int indexP = fl_interpP - c_sX0;



			cuComplex SlaveArray = tex2D(tex_slave, indexP, indexL);
			float sampleI = cuCrealf(SlaveArray);
			float sampleQ = cuCimagf(SlaveArray);
			float samplePhase = tex2D(tex_PhaseArray, indexP, indexL);

			double ReSampleI = sampleI*cos(samplePhase) + sampleQ*sin(samplePhase);
			double ReSampleQ = -sampleI*sin(samplePhase) + sampleQ*cos(samplePhase);

			rowoutput[col] = make_cuComplex(ReSampleI, ReSampleQ);





		}
		else
		{
			const int indexL = fl_interpL - c_sY0;
			const int indexP = fl_interpP - c_sX0;

			const float interpLdec = SslaveAz[threadIdx.y][threadIdx.x] - fl_interpL; //shared[offset1];//fl_interpL;            // e.g. .35432
			const float interpPdec = SslaveRg[threadIdx.y][threadIdx.x] - fl_interpP;//shared[offset1 + 1];//fl_interpP; // e.g. .5232
			//float* rowPhaseArray = (float*)((char*)PhaseArray + indexL*Spitch);


			const int kernelnoL = int(interpLdec * 2047 + 0.5); // lookup table index
			const int kernelnoP = int(interpPdec * 2047 + 0.5); // lookup table index


			cuComplex kernelL[6], kernelP[6];
#pragma unroll 2
			for (int i = 0; i < 6; ++i)
			{
				kernelL[i].x = tex2D(tex_kernelAz, i, kernelnoL);
				kernelP[i].x = tex2D(tex_kernelRg, i, kernelnoP);
				//Axis[i] = tex2D(tex_Axis, i, kernelnoL);
			}

			//azimuth_shift(kernelL, Axis, interpP, d_rsr2x, d_f_DC_a0, d_f_DC_a1, d_f_DC_a2, d_slave_prf, npoints);





			cuComplex B[36];
			cuComplex Phase[36];
			float temp;
			//B[0] = tex2D(tex_slave, fl_interpP-slave_pixelo-3, fl_interpL-slave_linelo-3);
#pragma unroll 2
			for (int j = 0; j < 6; ++j)
			{
				for (int i = 0; i < 6; ++i)
				{
					B[j * 6 + i] = tex2D(tex_slave, indexP - Npoints2m1 + i, indexL - Npoints2m1 + j);
					temp = tex2D(tex_PhaseArray, indexP - Npoints2m1 + i, indexL - Npoints2m1 + j);
					Phase[j * 6 + i] = make_cuComplex(temp, 0.0f);

				}
			}


			cuComplex sum = make_cuComplex(0.0f, 0.0f);
			cuComplex sum1 = make_cuComplex(0.0f, 0.0f);
			cuComplex Result[6];
			cuComplex Result1[6];

#pragma unroll 2
			for (int jj = 0; jj < 6; ++jj)
			{
				for (int k = 0; k < 6; ++k)
				{
					sum = cuCaddf(sum, cuCmulf(kernelP[k], B[jj * 6 + k]));
					sum1 = cuCaddf(sum1, cuCmulf(kernelP[k], Phase[jj * 6 + k]));

				}
				Result[jj] = sum;
				Result1[jj] = sum1;
				sum = make_cuComplex(0.0f, 0.0f);// complex requires this
				sum1 = make_cuComplex(0.0f, 0.0f);
			}

			/*行卷积获得最终结果*/
			sum = cuCaddf(cuCaddf(cuCaddf(cuCaddf(cuCaddf(cuCmulf(Result[0], kernelL[0]), cuCmulf(Result[1], kernelL[1])), cuCmulf(Result[2], kernelL[2])), cuCmulf(Result[3], kernelL[3])),
				cuCmulf(Result[4], kernelL[4])), cuCmulf(Result[5], kernelL[5]));
			sum1 = cuCaddf(cuCaddf(cuCaddf(cuCaddf(cuCaddf(cuCmulf(Result1[0], kernelL[0]), cuCmulf(Result1[1], kernelL[1])), cuCmulf(Result1[2], kernelL[2])), cuCmulf(Result1[3], kernelL[3])),
				cuCmulf(Result1[4], kernelL[4])), cuCmulf(Result1[5], kernelL[5]));


			float sampleI = cuCrealf(sum);
			float sampleQ = cuCimagf(sum);
			float samplePhase = cuCrealf(sum1);

			float ReSampleI = sampleI*cos(samplePhase) + sampleQ*sin(samplePhase);
			float ReSampleQ = -sampleI*sin(samplePhase) + sampleQ*cos(samplePhase);

			rowoutput[col] = make_cuComplex(ReSampleI, ReSampleQ);


		}
	}

}



__global__ void
Interpolation_12p
(cuComplex *output, int LineOffset, int Lines, size_t CorrPitch_1)
{

	const unsigned int row = blockIdx.y*blockDim.y + threadIdx.y;
	const unsigned int col = blockIdx.x*blockDim.x + threadIdx.x;

	cuComplex * rowoutput = (cuComplex *)((char*)output + row*CorrPitch_1);
	
	


	if (row < Lines && col < c_mPixels)
	{
		
		double Temp[2];
		Temp[0] = normalizeWarp(col, c_MasterBox[0], c_MasterBox[1]);
		Temp[1] = normalizeWarp(row +LineOffset + c_mY0, c_MasterBox[2], c_MasterBox[3]);

		

	
		
		int IndexesL[2];//[interger index, decimal index]
		int IndexesP[2];
		

		GetIndexes(my_polyval(Temp[1], Temp[0], c_CpmAz) + d_AzimuthShift, IndexesL);
		GetIndexes(my_polyval(Temp[1], Temp[0], c_CpmRg), IndexesP);

		
	


		if (IndexesL[0] > c_sYmax || IndexesL[0]<c_sY0
			|| IndexesP[0]>c_sXmax || IndexesP[0] < c_sX0)
		{
			rowoutput[col] = make_cuComplex(0.0f, 0.0f);

		}

		else if (IndexesL[0]>c_sYmax - c_Npoints2m1 - 2 || IndexesL[0]<c_sY0 + c_Npoints2m1
			|| IndexesP[0]>c_sXmax - c_Npoints2m1 - 2 || IndexesP[0] < c_sX0 + c_Npoints2m1)
		{
			IndexesP[0] -= c_sX0;
			IndexesL[0] -= c_sY0;
			cuComplex SlaveArray = tex2D(tex_slave, IndexesP[0], IndexesL[0]);
			float samplePhase = tex2D(tex_PhaseArray, IndexesP[0] , IndexesL[0] );

			sincos(samplePhase, Temp, Temp + 1);
			rowoutput[col] = cuCmulf(SlaveArray, make_cuComplex(Temp[1], -Temp[0]));



		}
		else
		{
			
			
			//Read look-up tables of interpolation convolution kernels 
			float kernelL[12];
			float kernelP[12];

			
			kernelL[0] = tex2D(tex_kernelAz, 0, IndexesL[1]);
			kernelL[1] = tex2D(tex_kernelAz, 1, IndexesL[1]);
			kernelL[2] = tex2D(tex_kernelAz, 2, IndexesL[1]);
			kernelL[3] = tex2D(tex_kernelAz, 3, IndexesL[1]);
			kernelL[4] = tex2D(tex_kernelAz, 4, IndexesL[1]);
			kernelL[5] = tex2D(tex_kernelAz, 5, IndexesL[1]);
			kernelL[6] = tex2D(tex_kernelAz, 6, IndexesL[1]);
			kernelL[7] = tex2D(tex_kernelAz, 7, IndexesL[1]);
			kernelL[8] = tex2D(tex_kernelAz, 8, IndexesL[1]);
			kernelL[9] = tex2D(tex_kernelAz, 9, IndexesL[1]);
			kernelL[10] = tex2D(tex_kernelAz, 10, IndexesL[1]);
			kernelL[11] = tex2D(tex_kernelAz, 11, IndexesL[1]);

			kernelP[0] = tex2D(tex_kernelRg, 0, IndexesP[1]);
			kernelP[1] = tex2D(tex_kernelRg, 1, IndexesP[1]);
			kernelP[2] = tex2D(tex_kernelRg, 2, IndexesP[1]);
			kernelP[3] = tex2D(tex_kernelRg, 3, IndexesP[1]);
			kernelP[4] = tex2D(tex_kernelRg, 4, IndexesP[1]);
			kernelP[5] = tex2D(tex_kernelRg, 5, IndexesP[1]);
			kernelP[6] = tex2D(tex_kernelRg, 6, IndexesP[1]);
			kernelP[7] = tex2D(tex_kernelRg, 7, IndexesP[1]);
			kernelP[8] = tex2D(tex_kernelRg, 8, IndexesP[1]);
			kernelP[9] = tex2D(tex_kernelRg, 9, IndexesP[1]);
			kernelP[10] = tex2D(tex_kernelRg, 10, IndexesP[1]);
			kernelP[11] = tex2D(tex_kernelRg, 11, IndexesP[1]);


		

			cuComplex tempComplex;


			cuComplex sum = make_cuComplex(0.0f, 0.0f);
		
			double sum1 = 0.0;
			double tempFloat;

			IndexesP[0] -=(c_sX0+ c_Npoints2m1);
			IndexesL[0] -=(c_sY0+ c_Npoints2m1);

			
			

			// Partially unroll the loop to reduce register pressure
			// Interpolate the slave image at first
			for (int j = 0; j < 12; ++j)
			{
				IndexesL[1] = IndexesL[0] + j;
				tempComplex =
					cuCaddf(CmulfFloat(tex2D(tex_slave, IndexesP[0] + 11,IndexesL[1]), kernelP[11]),
					cuCaddf(CmulfFloat(tex2D(tex_slave, IndexesP[0] + 10,IndexesL[1]), kernelP[10]),
					cuCaddf(CmulfFloat(tex2D(tex_slave, IndexesP[0] + 9,IndexesL[1]), kernelP[9]),
					cuCaddf(CmulfFloat(tex2D(tex_slave, IndexesP[0] + 8,IndexesL[1]), kernelP[8]),
					cuCaddf(CmulfFloat(tex2D(tex_slave, IndexesP[0] + 7,IndexesL[1]), kernelP[7]),
					cuCaddf(CmulfFloat(tex2D(tex_slave, IndexesP[0] + 6,IndexesL[1]), kernelP[6]),
					cuCaddf(CmulfFloat(tex2D(tex_slave, IndexesP[0] + 5,IndexesL[1]), kernelP[5]),
					cuCaddf(CmulfFloat(tex2D(tex_slave, IndexesP[0] + 4,IndexesL[1]), kernelP[4]),
					cuCaddf(CmulfFloat(tex2D(tex_slave, IndexesP[0] + 3,IndexesL[1]), kernelP[3]),
					cuCaddf(CmulfFloat(tex2D(tex_slave, IndexesP[0] + 2,IndexesL[1]), kernelP[2]),
					cuCaddf(CmulfFloat(tex2D(tex_slave, IndexesP[0],IndexesL[1]), kernelP[0]),
					CmulfFloat(tex2D(tex_slave, IndexesP[0] + 1,IndexesL[1]), kernelP[1])
					)))))))))));

				sum = cuCaddf(sum, CmulfFloat(tempComplex, kernelL[j]));
			}

		

			//Interpolate the deramping phase
			for (int j = 0; j < 12; ++j)
			{
				IndexesL[1] = IndexesL[0] + j;
				tempFloat = (double)tex2D(tex_PhaseArray, IndexesP[0],IndexesL[1])* kernelP[0]
					+ (double)tex2D(tex_PhaseArray, IndexesP[0] + 1,IndexesL[1])*kernelP[1]
					+ (double)tex2D(tex_PhaseArray, IndexesP[0] + 2,IndexesL[1])*kernelP[2]
					+ (double)tex2D(tex_PhaseArray, IndexesP[0] + 3,IndexesL[1])*kernelP[3]
					+ (double)tex2D(tex_PhaseArray, IndexesP[0] + 4,IndexesL[1])*kernelP[4]
					+ (double)tex2D(tex_PhaseArray, IndexesP[0] + 5,IndexesL[1])*kernelP[5]
					+ (double)tex2D(tex_PhaseArray, IndexesP[0] + 6,IndexesL[1])*kernelP[6]
					+ (double)tex2D(tex_PhaseArray, IndexesP[0] + 7,IndexesL[1])*kernelP[7]
					+ (double)tex2D(tex_PhaseArray, IndexesP[0] + 8,IndexesL[1])*kernelP[8]
					+ (double)tex2D(tex_PhaseArray, IndexesP[0] + 9,IndexesL[1])*kernelP[9]
					+ (double)tex2D(tex_PhaseArray, IndexesP[0] + 10,IndexesL[1])*kernelP[10]
					+ (double)tex2D(tex_PhaseArray, IndexesP[0] + 11,IndexesL[1])*kernelP[11];

				
				
				sum1 = sum1 + tempFloat * kernelL[j];
				//sum1 = __fma_rn(tempFloat, kernelL[j], sum1);


			}

		
			//Compute Sin(sum1) and Cos(sum1) and store them into the Temp Array
			sincos(sum1, Temp, Temp+1);
			
			//Rereamp and output aligned signals
			rowoutput[col] = cuCmulf(sum, make_cuComplex(Temp[1], -Temp[0]));

		
	
		

		}
	}

}

__global__ void
Interpolation_12p_ESDshift
(cuComplex *output, int LineOffset, int Lines, size_t CorrPitch_1)
{

	const unsigned int row = blockIdx.y*blockDim.y + threadIdx.y;
	const unsigned int col = blockIdx.x*blockDim.x + threadIdx.x;

	cuComplex * rowoutput = (cuComplex *)((char*)output + row*CorrPitch_1);




	if (row < Lines && col < c_mPixels)
	{

		double Temp[2];
		Temp[0] = normalizeWarp(col, c_MasterBox[0], c_MasterBox[1]);
		Temp[1] = normalizeWarp(row + LineOffset + c_mY0, c_MasterBox[2], c_MasterBox[3]);





		int IndexesL[2];//[interger index, decimal index]
		int IndexesP[2];


		GetIndexes(my_polyval(Temp[1], Temp[0], c_CpmAz) + d_AzimuthShift, IndexesL);
		GetIndexes(my_polyval(Temp[1], Temp[0], c_CpmRg), IndexesP);





		if (IndexesL[0] > c_sYmax || IndexesL[0]<c_sY0
			|| IndexesP[0]>c_sXmax || IndexesP[0] < c_sX0)
		{
			rowoutput[col] = make_cuComplex(0.0f, 0.0f);

		}

		else if (IndexesL[0]>c_sYmax - c_Npoints2m1 - 2 || IndexesL[0]<c_sY0 + c_Npoints2m1
			|| IndexesP[0]>c_sXmax - c_Npoints2m1 - 2 || IndexesP[0] < c_sX0 + c_Npoints2m1)
		{
			IndexesP[0] -= c_sX0;
			IndexesL[0] -= c_sY0;
			cuComplex SlaveArray = tex2D(tex_slave, IndexesP[0], IndexesL[0]);
			float samplePhase = tex2D(tex_PhaseArray, IndexesP[0], IndexesL[0]);

			sincos(samplePhase, Temp, Temp + 1);
			rowoutput[col] = cuCmulf(SlaveArray, make_cuComplex(Temp[1], -Temp[0]));



		}
		else
		{


			//Read look-up tables of interpolation convolution kernels 
			float kernelL[12];
			float kernelP[12];


			kernelL[0] = tex2D(tex_kernelAz, 0, IndexesL[1]);
			kernelL[1] = tex2D(tex_kernelAz, 1, IndexesL[1]);
			kernelL[2] = tex2D(tex_kernelAz, 2, IndexesL[1]);
			kernelL[3] = tex2D(tex_kernelAz, 3, IndexesL[1]);
			kernelL[4] = tex2D(tex_kernelAz, 4, IndexesL[1]);
			kernelL[5] = tex2D(tex_kernelAz, 5, IndexesL[1]);
			kernelL[6] = tex2D(tex_kernelAz, 6, IndexesL[1]);
			kernelL[7] = tex2D(tex_kernelAz, 7, IndexesL[1]);
			kernelL[8] = tex2D(tex_kernelAz, 8, IndexesL[1]);
			kernelL[9] = tex2D(tex_kernelAz, 9, IndexesL[1]);
			kernelL[10] = tex2D(tex_kernelAz, 10, IndexesL[1]);
			kernelL[11] = tex2D(tex_kernelAz, 11, IndexesL[1]);

			kernelP[0] = tex2D(tex_kernelRg, 0, IndexesP[1]);
			kernelP[1] = tex2D(tex_kernelRg, 1, IndexesP[1]);
			kernelP[2] = tex2D(tex_kernelRg, 2, IndexesP[1]);
			kernelP[3] = tex2D(tex_kernelRg, 3, IndexesP[1]);
			kernelP[4] = tex2D(tex_kernelRg, 4, IndexesP[1]);
			kernelP[5] = tex2D(tex_kernelRg, 5, IndexesP[1]);
			kernelP[6] = tex2D(tex_kernelRg, 6, IndexesP[1]);
			kernelP[7] = tex2D(tex_kernelRg, 7, IndexesP[1]);
			kernelP[8] = tex2D(tex_kernelRg, 8, IndexesP[1]);
			kernelP[9] = tex2D(tex_kernelRg, 9, IndexesP[1]);
			kernelP[10] = tex2D(tex_kernelRg, 10, IndexesP[1]);
			kernelP[11] = tex2D(tex_kernelRg, 11, IndexesP[1]);




			cuComplex tempComplex;


			cuComplex sum = make_cuComplex(0.0f, 0.0f);

			double sum1 = 0.0;
			double tempFloat;

			IndexesP[0] -= (c_sX0 + c_Npoints2m1);
			IndexesL[0] -= (c_sY0 + c_Npoints2m1);




			// Partially unroll the loop to reduce register pressure
			// Interpolate the slave image at first
			for (int j = 0; j < 12; ++j)
			{
				IndexesL[1] = IndexesL[0] + j;
				tempComplex =
					cuCaddf(CmulfFloat(tex2D(tex_slave, IndexesP[0] + 11, IndexesL[1]), kernelP[11]),
					cuCaddf(CmulfFloat(tex2D(tex_slave, IndexesP[0] + 10, IndexesL[1]), kernelP[10]),
					cuCaddf(CmulfFloat(tex2D(tex_slave, IndexesP[0] + 9, IndexesL[1]), kernelP[9]),
					cuCaddf(CmulfFloat(tex2D(tex_slave, IndexesP[0] + 8, IndexesL[1]), kernelP[8]),
					cuCaddf(CmulfFloat(tex2D(tex_slave, IndexesP[0] + 7, IndexesL[1]), kernelP[7]),
					cuCaddf(CmulfFloat(tex2D(tex_slave, IndexesP[0] + 6, IndexesL[1]), kernelP[6]),
					cuCaddf(CmulfFloat(tex2D(tex_slave, IndexesP[0] + 5, IndexesL[1]), kernelP[5]),
					cuCaddf(CmulfFloat(tex2D(tex_slave, IndexesP[0] + 4, IndexesL[1]), kernelP[4]),
					cuCaddf(CmulfFloat(tex2D(tex_slave, IndexesP[0] + 3, IndexesL[1]), kernelP[3]),
					cuCaddf(CmulfFloat(tex2D(tex_slave, IndexesP[0] + 2, IndexesL[1]), kernelP[2]),
					cuCaddf(CmulfFloat(tex2D(tex_slave, IndexesP[0], IndexesL[1]), kernelP[0]),
					CmulfFloat(tex2D(tex_slave, IndexesP[0] + 1, IndexesL[1]), kernelP[1])
					)))))))))));

				sum = cuCaddf(sum, CmulfFloat(tempComplex, kernelL[j]));
			}



			//Interpolate the deramping phase
			for (int j = 0; j < 12; ++j)
			{
				IndexesL[1] = IndexesL[0] + j;
				tempFloat = (double)tex2D(tex_PhaseArray, IndexesP[0], IndexesL[1])* kernelP[0]
					+ (double)tex2D(tex_PhaseArray, IndexesP[0] + 1, IndexesL[1])*kernelP[1]
					+ (double)tex2D(tex_PhaseArray, IndexesP[0] + 2, IndexesL[1])*kernelP[2]
					+ (double)tex2D(tex_PhaseArray, IndexesP[0] + 3, IndexesL[1])*kernelP[3]
					+ (double)tex2D(tex_PhaseArray, IndexesP[0] + 4, IndexesL[1])*kernelP[4]
					+ (double)tex2D(tex_PhaseArray, IndexesP[0] + 5, IndexesL[1])*kernelP[5]
					+ (double)tex2D(tex_PhaseArray, IndexesP[0] + 6, IndexesL[1])*kernelP[6]
					+ (double)tex2D(tex_PhaseArray, IndexesP[0] + 7, IndexesL[1])*kernelP[7]
					+ (double)tex2D(tex_PhaseArray, IndexesP[0] + 8, IndexesL[1])*kernelP[8]
					+ (double)tex2D(tex_PhaseArray, IndexesP[0] + 9, IndexesL[1])*kernelP[9]
					+ (double)tex2D(tex_PhaseArray, IndexesP[0] + 10, IndexesL[1])*kernelP[10]
					+ (double)tex2D(tex_PhaseArray, IndexesP[0] + 11, IndexesL[1])*kernelP[11];



				sum1 = sum1 + tempFloat * kernelL[j];
				//sum1 = __fma_rn(tempFloat, kernelL[j], sum1);


			}


			//Compute Sin(sum1) and Cos(sum1) and store them into the Temp Array
			sincos(sum1, Temp, Temp + 1);

			//Rereamp and output aligned signals
			rowoutput[col] = cuCmulf(sum, make_cuComplex(Temp[1], -Temp[0]));





		}
	}

}


__global__ void
resample_texture_kernel_12p_warpFunction_raw
(cuComplex *output,
cuComplex* slaveArray,
float *phaseArray, 
float *kernelAz, 
float *kernelRg, 
size_t d_pitch1,
size_t d_pitch2,
int Lines,
size_t CorrPitch_1
)
{


	const int row = blockIdx.y*blockDim.y + threadIdx.y;
	const int col = blockIdx.x*blockDim.x + threadIdx.x;

	cuComplex * rowoutput = (cuComplex *)((char*)output + row*CorrPitch_1);

	// adjust pCorr to point to row
	//output = (float *)((char*)output + row*CorrPitch);

	int Npoints2m1 = 12 / 2 - 1;
	//int ncols = SlavePitch / 8;



	if (row < Lines && col < c_mPixels)
	{
	

		double normalizeP = normalizeWarp((double)col, c_MasterBox[0], c_MasterBox[1]);
		double normalizeL = normalizeWarp((double)row, c_MasterBox[2], c_MasterBox[3]);

		double SlavePos_x = my_polyval(normalizeL, normalizeP, c_CpmRg);
		double SlavePos_y = my_polyval(normalizeL, normalizeP, c_CpmAz) + d_AzimuthShift;
		

		int  fl_interpL = int(SlavePos_y);
		int fl_interpP = int(SlavePos_x);


		int indexL = fl_interpL - c_sY0;
		int indexP = fl_interpP - c_sX0;



		if (fl_interpL > c_sYmax || fl_interpL<c_sY0
			|| fl_interpP>c_sXmax || fl_interpP < c_sX0)
		{

			rowoutput[col] = make_cuComplex(0.0f, 0.0f);

		}

		else if (fl_interpL>c_sYmax - Npoints2m1 - 2 || fl_interpL<c_sY0 + Npoints2m1
			|| fl_interpP>c_sXmax - Npoints2m1 - 2 || fl_interpP < c_sX0 + Npoints2m1)
		{

		
			cuComplex* rowSlave = (cuComplex *)((char*)slaveArray + indexL*d_pitch2);
			float* rowPhase = (float *)((char*)phaseArray + indexL*d_pitch1);
			cuComplex Slave = rowSlave[indexP];//tex2D(tex_slave, indexP, indexL);//SlaveArray[indexL*ncols + indexP];
			float samplePhase = rowPhase[indexP];  //tex2D(tex_PhaseArray, indexP, indexL);

			rowoutput[col] = cuCmulf(Slave, make_cuComplex(cos(samplePhase), -sin(samplePhase)));



		}
		else
		{

			const int kernelnoL = int((SlavePos_y - fl_interpL) * 2047 + 0.5); // lookup table index
			const int kernelnoP = int((SlavePos_x - fl_interpP) * 2047 + 0.5); // lookup table index

		/*	float kernelL[12], kernelP[12];


			for (int i = 0; i < 12; i++)
			{
			
				kernelL[i] = kernelAz[kernelnoL * 12 + i];
				kernelP[i] = kernelRg[kernelnoP * 12 + i];

			}*/





			/*stupid but neccessary for experiment*/

			//cuComplex B[144];
			//float Phase[144];
			//cuComplex * rowB;
			//float* rowPhase;
			//


			//for (int j = 0; j < 12; ++j)
			//{
			//	rowB = (cuComplex *)((char*)slaveArray + (indexL - Npoints2m1 + j)*d_pitch2);
			//	rowPhase = (float *)((char*)phaseArray + (indexL - Npoints2m1 + j)*d_pitch1);
			//	for (int i = 0; i < 12; ++i)
			//	{
			//		
			//		B[j * 12 + i] = rowB[indexP - Npoints2m1 + i]; //tex2D(tex_slave, indexP - Npoints2m1 + i, indexL - Npoints2m1 + j);
			//		Phase[j * 12 + i] = rowPhase[indexP - Npoints2m1 + i];//tex2D(tex_PhaseArray, indexP - Npoints2m1 + i, indexL - Npoints2m1 + j);
			//		

			//	}
			//}


			//cuComplex sum = make_cuComplex(0.0f, 0.0f);
			//float sum1 = 0.0f;
			//cuComplex Result[12];
			//float Result1[12];


			//for (int jj = 0; jj < 12; ++jj)
			//{
			//	for (int k = 0; k < 12; ++k)
			//	{
			//		sum = cuCaddf(sum, CmulfFloat(B[jj * 12 + k],kernelP[k]));
			//		sum1 += Phase[jj * 12 + k] * kernelP[k];

			//	}
			//	Result[jj] = sum;
			//	Result1[jj] = sum1;
			//	sum = make_cuComplex(0.0f, 0.0f);// complex requires this
			//	sum1 = 0.0f;
			//}

			//for (int iter = 0; iter < 12; iter++)
			//{
			//	sum = cuCaddf(sum, CmulfFloat(Result[iter], kernelL[iter]));
			//	sum1 += Result1[iter]* kernelL[iter];
			//}

			/*stupid but neccessary for experiment*/

			/*clever but not enough */


			cuComplex sum = make_cuComplex(0.0f, 0.0f);
			cuComplex tempComplex = make_cuComplex(0.0f, 0.0f);
			double tempFloat=0.0;

			float sum1 = 0.0f;
			cuComplex * rowB;
			float* rowPhase;

			indexL -= Npoints2m1;
			indexP -= Npoints2m1;

//			for (int j = 0; j < 12; ++j)
//			{
//				rowB = (cuComplex *)((char*)slaveArray + (indexL + j)*d_pitch2);
//
//
//
//
//				tempComplex = make_cuComplex(0.0f, 0.0f);
//				for (int i = 0; i < 12; i++)
//				{
//
//					tempComplex = cuCaddf(tempComplex, CmulfFloat(rowB[indexP+i], kernelRg[kernelnoP * 12 + i]));
//					/*tempComplex = cuCaddf(tempComplex, CmulfFloat(rowB[indexP + 1], kernelP[1]));
//					tempComplex = cuCaddf(tempComplex, CmulfFloat(rowB[indexP + 2], kernelP[2]));
//					tempComplex = cuCaddf(tempComplex, CmulfFloat(rowB[indexP + 3], kernelP[3]));
//					tempComplex = cuCaddf(tempComplex, CmulfFloat(rowB[indexP + 4], kernelP[4]));
//					tempComplex = cuCaddf(tempComplex, CmulfFloat(rowB[indexP + 5], kernelP[5]));
//					tempComplex = cuCaddf(tempComplex, CmulfFloat(rowB[indexP + 6], kernelP[6]));
//					tempComplex = cuCaddf(tempComplex, CmulfFloat(rowB[indexP + 7], kernelP[7]));
//					tempComplex = cuCaddf(tempComplex, CmulfFloat(rowB[indexP + 8], kernelP[8]));
//					tempComplex = cuCaddf(tempComplex, CmulfFloat(rowB[indexP + 9], kernelP[9]));
//					tempComplex = cuCaddf(tempComplex, CmulfFloat(rowB[indexP + 10], kernelP[10]));
//					tempComplex = cuCaddf(tempComplex, CmulfFloat(rowB[indexP + 11], kernelP[11]));
//*/
//
//				}
//
//				//tempComplex = CmulfFloat(tempComplex, kernelL[j]);
//				tempComplex = CmulfFloat(tempComplex, kernelAz[kernelnoL * 12 + j]);
//
//				sum = cuCaddf(sum, tempComplex);
//
//
//			}
//
//			for (int j = 0; j < 12; ++j)
//			{
//				rowPhase = (float *)((char*)phaseArray + (indexL + j)*d_pitch1);
//				tempFloat = 0.0;
//
//				for (int i = 0; i < 12; i++)
//				{
//
//					tempFloat += (double)rowPhase[indexP+i] * kernelRg[kernelnoP * 12 + i];
//				}
//				/*tempFloat += (double)rowPhase[indexP + 1] * kernelP[1];
//				tempFloat += (double)rowPhase[indexP + 2] * kernelP[2];
//				tempFloat += (double)rowPhase[indexP + 3] * kernelP[3];
//				tempFloat += (double)rowPhase[indexP + 4] * kernelP[4];
//				tempFloat += (double)rowPhase[indexP + 5] * kernelP[5];
//				tempFloat += (double)rowPhase[indexP + 6] * kernelP[6];
//				tempFloat += (double)rowPhase[indexP + 7] * kernelP[7];
//				tempFloat += (double)rowPhase[indexP + 8] * kernelP[8];
//				tempFloat += (double)rowPhase[indexP + 9] * kernelP[9];
//				tempFloat += (double)rowPhase[indexP + 10] * kernelP[10];
//				tempFloat += (double)rowPhase[indexP + 11] * kernelP[11];*/
//
//
//
//				tempFloat *= kernelAz[kernelnoL * 12 + j];
//
//				//Result[j].y = 0.0f;
//				sum1 += tempFloat;
//
//			}
			
			for (int i = 0; i < 12; i++)
			{
				rowB = (cuComplex *)((char*)slaveArray + (indexL + i)*d_pitch2);
				rowPhase = (float *)((char*)phaseArray + (indexL + i)*d_pitch1);

				//tempComplex = make_cuComplex(0.0f, 0.0f);
				//tempFloat = 0.0f;

				//float tmpAz = kernelAz[kernelnoL * 12 + i];
				for (int j = 0; j < 12; j++)
				{
					
					float tmpRg = kernelRg[kernelnoP * 12 + j];
				
					tempComplex = cuCaddf(tempComplex, 
						CmulfFloat(CmulfFloat(rowB[indexP + j], tmpRg), kernelAz[kernelnoL * 12 + i]));
					tempFloat += (double)rowPhase[indexP + j] * tmpRg * kernelAz[kernelnoL * 12 + i];

					//tempComplex = cuCaddf(tempComplex,
						//CmulfFloat(rowB[indexP + j], kernelRg[kernelnoP * 12 + j]));
					//tempFloat += (double)rowPhase[indexP + j] * kernelRg[kernelnoP * 12 + j];


				}
				//sum = cuCaddf(sum, CmulfFloat(tempComplex, kernelAz[kernelnoL * 12 + i]));
				//sum1 += tempFloat*kernelAz[kernelnoL * 12 + i];
			}
			sum = tempComplex;
			sum1 = tempFloat;
		
		
			rowoutput[col] = cuCmulf(sum, make_cuComplex(cos(sum1), -sin(sum1)));
	
		}
	}

}


__global__ void
resample_kernel_12p_textureTest
(cuComplex *output, int LineOffset, int Lines, size_t CorrPitch_1)
{

	const int row = blockIdx.y*blockDim.y + threadIdx.y;
	const int col = blockIdx.x*blockDim.x + threadIdx.x;

	cuComplex * rowoutput = (cuComplex *)((char*)output + row*CorrPitch_1);

	// adjust pCorr to point to row


	int Npoints2m1 = 12 / 2 - 1;



	if (row < Lines && col < c_mPixels)
	{

		double normalizeP = normalizeWarp((double)col, c_MasterBox[0], c_MasterBox[1]);
		double normalizeL = normalizeWarp((double)(row + LineOffset + c_mY0), c_MasterBox[2], c_MasterBox[3]);



		double SlavePos_x = my_polyval(normalizeL, normalizeP, c_CpmRg);
		double SlavePos_y = my_polyval(normalizeL, normalizeP, c_CpmAz) + d_AzimuthShift;




		int  fl_interpL = int(SlavePos_y);
		int fl_interpP = int(SlavePos_x);


		int indexL = fl_interpL - c_sY0;
		int indexP = fl_interpP - c_sX0;



		if (fl_interpL > c_sYmax || fl_interpL<c_sY0
			|| fl_interpP>c_sXmax || fl_interpP < c_sX0)
		{
			rowoutput[col] = make_cuComplex(0.0f, 0.0f);

		}

		else if (fl_interpL>c_sYmax - Npoints2m1 - 2 || fl_interpL<c_sY0 + Npoints2m1
			|| fl_interpP>c_sXmax - Npoints2m1 - 2 || fl_interpP < c_sX0 + Npoints2m1)
		{

			//float* rowPhaseArray = (float*)((char*)PhaseArray + indexL*Spitch);


			cuComplex SlaveArray = tex2D(tex_slave, indexP, indexL);
			float samplePhase = tex2D(tex_PhaseArray, indexP, indexL);


			/*double ReSampleI = sampleI*cos(samplePhase) + sampleQ*sin(samplePhase);
			double ReSampleQ = -sampleI*sin(samplePhase) + sampleQ*cos(samplePhase);*/

			rowoutput[col] = cuCmulf(SlaveArray, make_cuComplex(cos(samplePhase), -sin(samplePhase)));



		}
		else
		{
			//const float interpLdec = SlavePos_y - fl_interpL; //shared[offset1];//fl_interpL;            // e.g. .35432
			//const float interpPdec = SlavePos_x - fl_interpP;//shared[offset1 + 1];//fl_interpP; // e.g. .5232


			const int kernelnoL = int((SlavePos_y - fl_interpL) * 2047 + 0.5); // lookup table index
			const int kernelnoP = int((SlavePos_x - fl_interpP) * 2047 + 0.5); // lookup table index
			//float kernelL_real[6], /*kernelL_img[6],*/ Axis[6];
			//float kernelP_real[6];
			//float Axis[4];
			/*float kernelL[12], kernelP[12];*/

			//for (int i = 0; i < 12; ++i)
			//{
			/*kernelL[0] = tex2D(tex_kernelAz, 0, kernelnoL);
			kernelL[1] = tex2D(tex_kernelAz, 1, kernelnoL);
			kernelL[2] = tex2D(tex_kernelAz, 2, kernelnoL);
			kernelL[3] = tex2D(tex_kernelAz, 3, kernelnoL);
			kernelL[4] = tex2D(tex_kernelAz, 4, kernelnoL);
			kernelL[5] = tex2D(tex_kernelAz, 5, kernelnoL);
			kernelL[6] = tex2D(tex_kernelAz, 6, kernelnoL);
			kernelL[7] = tex2D(tex_kernelAz, 7, kernelnoL);
			kernelL[8] = tex2D(tex_kernelAz, 8, kernelnoL);
			kernelL[9] = tex2D(tex_kernelAz, 9, kernelnoL);
			kernelL[10] = tex2D(tex_kernelAz, 10, kernelnoL);
			kernelL[11] = tex2D(tex_kernelAz, 11, kernelnoL);

			kernelP[0] = tex2D(tex_kernelRg, 0, kernelnoP);
			kernelP[1] = tex2D(tex_kernelRg, 1, kernelnoP);
			kernelP[2] = tex2D(tex_kernelRg, 2, kernelnoP);
			kernelP[3] = tex2D(tex_kernelRg, 3, kernelnoP);
			kernelP[4] = tex2D(tex_kernelRg, 4, kernelnoP);
			kernelP[5] = tex2D(tex_kernelRg, 5, kernelnoP);
			kernelP[6] = tex2D(tex_kernelRg, 6, kernelnoP);
			kernelP[7] = tex2D(tex_kernelRg, 7, kernelnoP);
			kernelP[8] = tex2D(tex_kernelRg, 8, kernelnoP);
			kernelP[9] = tex2D(tex_kernelRg, 9, kernelnoP);
			kernelP[10] = tex2D(tex_kernelRg, 10, kernelnoP);
			kernelP[11] = tex2D(tex_kernelRg, 11, kernelnoP);*/
			//Axis[i] = tex2D(tex_Axis, i, kernelnoL);
			//}



			//cuComplex Result[12];
			cuComplex tempComplex;


			cuComplex sum = make_cuComplex(0.0f, 0.0f);
			//cuComplex sum1=make_cuComplex(0.0f,0.0f);
			double sum1 = 0.0;
			double tempFloat;
			//cuComplex temp[12];
			//cuComplex temp[144];
			indexP -= Npoints2m1;
			indexL -= Npoints2m1;


			for (int i = 0; i < 12; i++)
			{
				for (int j = 0; j < 12; j++)
				{
					sum = cuCaddf(sum,
						CmulfFloat(CmulfFloat(tex2D(tex_slave, indexP + j, indexL + i), tex2D(tex_kernelAz, i, kernelnoL)), tex2D(tex_kernelRg, j, kernelnoP)));
					
				}
			}

			for (int i = 0; i < 12; i++)
			{
				for (int j = 0; j < 12; j++)
				{

					sum1 += (double)tex2D(tex_PhaseArray, indexP + j, indexL + i)* tex2D(tex_kernelAz, i, kernelnoL) * tex2D(tex_kernelRg, j, kernelnoP);
				}
			}
			

			rowoutput[col] = cuCmulf(sum, make_cuComplex(cos(sum1), -sin(sum1)));

		}
	}

}

    void DerampDemodResample(
	complex<short>*SlaveArray,
	double *CpmAz,
	double *CpmRg,
	double AzimuthShift,
	complex<float>* output,
	float *KernelAz,
	float *KernelRg,
	int sBurstIdx,
	int slave_pixels,
	int slave_lines,
	int MasterBox[4],
	int SlaveBox[4],
	int linesPerBurst,
	int SamplesPerBurst,
	double azimuthTimeInterval,
	double* dopplerRate,
	double* referenceTime,
	double* dopplerCentroid,
	int Npoints
	)
{

	int sLines = SlaveBox[3] - SlaveBox[2] + 1;
	int sPixels = SlaveBox[1] - SlaveBox[0] + 1;
	int mLines = MasterBox[3] - MasterBox[2]+1;
	int mPixels = MasterBox[1] - MasterBox[0] + 1;

	

	int Npoints2m1 = Npoints / 2 - 1;


	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 32, 0, 0, cudaChannelFormatKindFloat);
	cudaChannelFormatDesc channelDesc_1 = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);

	int sfirstLineInBurst = sBurstIdx*linesPerBurst;

	cudaHostRegister(SlaveArray, sLines*sPixels*sizeof(short2), cudaHostRegisterDefault);
	cudaHostRegister(dopplerCentroid, SamplesPerBurst*sizeof(double), cudaHostRegisterDefault);
	cudaHostRegister(referenceTime, SamplesPerBurst*sizeof(double), cudaHostRegisterDefault);
	cudaHostRegister(dopplerRate, SamplesPerBurst*sizeof(double), cudaHostRegisterDefault);
	cudaHostRegister(output, mLines*mPixels*sizeof(cuComplex), cudaHostRegisterDefault);
	cudaHostRegister(KernelAz, 2048 * 12 * sizeof(float), cudaHostRegisterDefault);
	cudaHostRegister(KernelRg, 2048 * 12 * sizeof(float), cudaHostRegisterDefault);


	// (cudaHostRegister(PhaseArray, sLines*sPixels*sizeof(float), cudaHostRegisterDefault));
	size_t d_pitch1, d_pitch2, d_pitch3, d_pitchS2;


	// It is worth to use another array to save complex<short>
	float* d_PhaseArray;
	cudaMallocPitch((void**)&d_PhaseArray, &d_pitch1, sPixels*sizeof(float), sLines);

	short2* d_SlaveArrayS2;
	cudaMallocPitch((void**)&d_SlaveArrayS2, &d_pitchS2, sPixels*sizeof(short2), sLines);
	cuComplex* d_SlaveArray;
	cudaMallocPitch((void**)&d_SlaveArray, &d_pitch2, sPixels*sizeof(cuComplex), sLines);


	double* d_dopplerRate, *d_referenceTime, *d_dopplerCentroid;
	cudaMallocPitch((void**)&d_dopplerRate, &d_pitch3, SamplesPerBurst*sizeof(double), 1);
	cudaMallocPitch((void**)&d_referenceTime, &d_pitch3, SamplesPerBurst*sizeof(double), 1);
	cudaMallocPitch((void**)&d_dopplerCentroid, &d_pitch3, SamplesPerBurst*sizeof(double), 1);


	size_t CorrPitch;

	cuComplex * d_resample;
	cudaMallocPitch((void **)&d_resample, &CorrPitch, mPixels*sizeof(cuComplex), mLines);


	cudaArray *KernelAzArray = NULL;
	cudaArray *KernelRgArray = NULL;
	cudaMallocArray(&KernelAzArray, &channelDesc_1, Npoints, 2048);
	cudaMallocArray(&KernelRgArray, &channelDesc_1, Npoints, 2048);


	dim3 threads(16, 16);
	dim3 blocks = dim3((sPixels + 15) / 16, (sLines + 15) / 16);
	cudaStream_t stream[4];
	for (int i = 0; i < 4; i++)
	{
		cudaStreamCreate(&stream[i]);
	}
	cudaEvent_t g_start, g_stop;
	//cudaFuncSetCacheConfig(resample_texture_kernel_12p_overlap_warpFunction_test, cudaFuncCachePreferL1);

	//Memcpy to Constant  Variables
	cudaMemcpyToSymbol(npoints, &Npoints, sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_mLines, &mLines, sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_mPixels, &mPixels, sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_sX0, &SlaveBox[0], sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_sXmax, &SlaveBox[1], sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_sY0, &SlaveBox[2], sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_sYmax, &SlaveBox[3], sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_mY0, &MasterBox[2], sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_AzimuthShift, &AzimuthShift, sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_MasterBox, MasterBox, 4 * sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_Npoints2m1, &Npoints2m1, sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_CpmAz, CpmAz, 6 * sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_CpmRg, CpmRg, 6 * sizeof(double), 0, cudaMemcpyHostToDevice);




	// (cudaMemcpy2DAsync(d_SlaveArray, d_pitch2, SlaveArray, sPixels*sizeof(cuComplex), sPixels*sizeof(cuComplex), sLines, cudaMemcpyHostToDevice));
	// (cudaMemcpy2DAsync(d_dopplerRate, d_pitch3, dopplerRate, SamplesPerBurst*sizeof(double), SamplesPerBurst*sizeof(double), 1, cudaMemcpyHostToDevice));
	// (cudaMemcpy2DAsync(d_referenceTime, d_pitch3, referenceTime, SamplesPerBurst*sizeof(double), SamplesPerBurst*sizeof(double), 1, cudaMemcpyHostToDevice));
	// (cudaMemcpy2DAsync(d_dopplerCentroid, d_pitch3, dopplerCentroid, SamplesPerBurst*sizeof(double), SamplesPerBurst*sizeof(double), 1, cudaMemcpyHostToDevice));

	// (cudaMemcpy2D(d_SlaveArray, d_pitch2, SlaveArray, sPixels*sizeof(cuComplex), sPixels*sizeof(cuComplex), sLines, cudaMemcpyHostToDevice));
	cudaMemcpy2D(d_SlaveArrayS2, d_pitchS2, SlaveArray, sPixels*sizeof(short2), sPixels*sizeof(short2), sLines, cudaMemcpyHostToDevice);
	cudaMemcpy2D(d_dopplerRate, d_pitch3, dopplerRate, SamplesPerBurst*sizeof(double), SamplesPerBurst*sizeof(double), 1, cudaMemcpyHostToDevice);
	cudaMemcpy2D(d_referenceTime, d_pitch3, referenceTime, SamplesPerBurst*sizeof(double), SamplesPerBurst*sizeof(double), 1, cudaMemcpyHostToDevice);
	cudaMemcpy2D(d_dopplerCentroid, d_pitch3, dopplerCentroid, SamplesPerBurst*sizeof(double), SamplesPerBurst*sizeof(double), 1, cudaMemcpyHostToDevice);


	//DerampDemod_Shared << <blocks, threads >> >(d_PhaseArray, d_SlaveArray, d_dopplerRate, d_referenceTime, d_dopplerCentroid, sX0, sY0, sLines, sPixels, sfirstLineInBurst,
	//azimuthTimeInterval, d_pitch1, d_pitch2, d_pitch3);

	DerampDemod_Shared << <blocks, threads >> >(d_PhaseArray, d_SlaveArray, d_dopplerRate, d_referenceTime, d_dopplerCentroid, SlaveBox[0], SlaveBox[2], sLines, sPixels, sfirstLineInBurst,
		azimuthTimeInterval, d_pitch1, d_pitch2, d_pitch3, d_SlaveArrayS2, d_pitchS2);

	//cudaEventRecord(g_stop, 0);
	//cudaEventSynchronize(g_stop);
	//cudaEventElapsedTime(&time_cost1, g_start, g_stop);
	//cout << "DeRamping duration:" << time_cost1 << endl;

	
	



	float *d_KernelAz;
	float *d_KernelRg;
	cudaMalloc((void**)&d_KernelAz, Npoints * 2048 * sizeof(float));
	cudaMalloc((void**)&d_KernelRg, Npoints * 2048 * sizeof(float));
	cudaMemcpy(d_KernelAz, KernelAz, Npoints * 2048 * sizeof(float),cudaMemcpyDeviceToHost);
	cudaMemcpy(d_KernelRg, KernelRg, Npoints * 2048 * sizeof(float), cudaMemcpyDeviceToHost);

	cudaMemcpyToArray(KernelAzArray, 0, 0, KernelAz, Npoints * 2048 * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpyToArray(KernelRgArray, 0, 0, KernelRg, Npoints * 2048 * sizeof(float), cudaMemcpyHostToDevice);
	cudaBindTextureToArray(tex_kernelAz, KernelAzArray, channelDesc_1);
	cudaBindTextureToArray(tex_kernelRg, KernelRgArray, channelDesc_1);
	cudaBindTexture2D(0, tex_PhaseArray, d_PhaseArray, channelDesc_1, sPixels, sLines, d_pitch1);
	cudaBindTexture2D(0, tex_slave, d_SlaveArray, channelDesc, sPixels, sLines, d_pitch2);


	size_t SPitch, MPitch;


	threads = dim3(16, 16);
	blocks = dim3((mPixels + threads.x - 1) / threads.x,
		(mLines + threads.y - 1) / threads.y);


	bool OverlapOrNot = false;

	//for Subsets
	int partMlines = mLines / 4;
	int RemainMlines = mLines % 4;
	dim3 Partblocks = dim3((mPixels + threads.x - 1) / threads.x,
		(partMlines + threads.y - 1) / threads.y);
	dim3 Lastblocks = dim3((mPixels + threads.x - 1) / threads.x,
		(partMlines + RemainMlines + threads.y - 1) / threads.y);

	int PartOffsetD = partMlines*CorrPitch / 8;
	int PartOffsetH = partMlines*mPixels;




	//cudaFuncSetCacheConfig(Interpolation_12p, cudaFuncCachePreferL1); Boosting configuration for kepler architecture 
	int NumTest = 10;

	float time_cost1, time_cost2;
	cudaEventCreate(&g_start);
	cudaEventCreate(&g_stop);
	cudaEventRecord(g_start, 0);


	if (Npoints == 6)
	{

	}
	
	

	if (Npoints == 12)
	{
		//for raw test

		if (!OverlapOrNot)
		{
			/*resample_texture_kernel_12p_warpFunction_raw << <blocks, threads >> >(d_resample, d_SlaveArray, d_PhaseArray,
				d_KernelAz, d_KernelRg, d_pitch1, d_pitch2,mLines, CorrPitch);*/
			//resample_kernel_12p_textureTest << <blocks, threads >> >(d_resample, 0, mLines, CorrPitch);

			 //for whole test
			Interpolation_12p
				<< <blocks, threads >> >(d_resample, 0,
				mLines, CorrPitch);
			
		
		}
		//SubSet1

		//Interpolation_12p
		//<< <Partblocks, threads, 0, stream[0] >> >(d_resample, 0,
		//	partMlines, CorrPitch);

		//cudaMemcpy2DAsync(output, (mPixels)*sizeof(cuComplex), d_resample, CorrPitch, (mPixels)*sizeof(cuComplex), partMlines, cudaMemcpyDeviceToHost, stream[0]);



		////SubSet2
		//Interpolation_12p
		//	<< <Partblocks, threads, 0, stream[1] >>>
		//	(d_resample + PartOffsetD, partMlines,
		//	partMlines, CorrPitch);
		//cudaMemcpy2DAsync(output + PartOffsetH, (mPixels)*sizeof(cuComplex), 
		//	d_resample + PartOffsetD, CorrPitch, (mPixels)*sizeof(cuComplex),
		//	partMlines, cudaMemcpyDeviceToHost, stream[1]);

		////SubSet3
		//
		//Interpolation_12p
		//<< <Partblocks, threads, 0, stream[2] >>>
		//(d_resample + 2 * PartOffsetD, 2 * partMlines,
		//	partMlines, CorrPitch);
		//cudaMemcpy2DAsync(output + 2 * PartOffsetH, (mPixels)*sizeof(cuComplex),
		//	d_resample + 2 * PartOffsetD, CorrPitch, (mPixels)*sizeof(cuComplex),
		//	partMlines, cudaMemcpyDeviceToHost, stream[2]);


		////SubSet4
		//Interpolation_12p
		//	<< <Lastblocks, threads, 0, stream[3] >> >
		//	(d_resample + 3 * PartOffsetD, 3 * partMlines,
		//	partMlines + RemainMlines, CorrPitch);
		//cudaMemcpy2DAsync(output + 3 * PartOffsetH, (mPixels)*sizeof(cuComplex), 
		//	d_resample + 3 * PartOffsetD, CorrPitch, (mPixels)*sizeof(cuComplex), 
		//	(partMlines + RemainMlines), cudaMemcpyDeviceToHost, stream[3]);

		

		

	}
	

	cudaEventRecord(g_stop, 0);
	cudaEventSynchronize(g_stop);
	cudaEventElapsedTime(&time_cost2, g_start, g_stop);
	cout << "kernel duration:" << time_cost2  <<"ms"<< endl;
	cudaEventDestroy(g_start);
	cudaEventDestroy(g_stop);
	
	cudaMemcpy2D(output, mPixels*sizeof(cuComplex),
		d_resample, CorrPitch, mPixels*sizeof(cuComplex),
		mLines, cudaMemcpyDeviceToHost);


	for (int i = 0; i < 4; i++)
	{
		cudaStreamDestroy(stream[i]);
	}



	cudaHostUnregister(SlaveArray);
	cudaHostUnregister(output);
	cudaHostUnregister(dopplerCentroid);
	cudaHostUnregister(dopplerRate);
	cudaHostUnregister(referenceTime);
	cudaHostUnregister(KernelAz);
	cudaHostUnregister(KernelRg);

	cudaUnbindTexture(tex_kernelAz);
	cudaUnbindTexture(tex_kernelRg);
	cudaUnbindTexture(tex_slave);
	cudaUnbindTexture(tex_PhaseArray);

	cudaFree(d_PhaseArray);
	cudaFree(d_SlaveArray);
	cudaFree(d_dopplerRate);
	cudaFree(d_referenceTime);
	cudaFree(d_dopplerCentroid);
	cudaFree(d_resample);
	cudaFreeArray(KernelAzArray);
	cudaFreeArray(KernelRgArray);
	cudaFree(d_SlaveArrayS2);

	

	
	cudaDeviceReset();



	
}



cuComplex* DerampDemodResample_ESD(
	complex<short>*SlaveArray,
	double *CpmAz,
	double *CpmRg,
	double AzimuthShift,
	complex<float>* output,
	float *KernelAz,
	float *KernelRg,
	int sBurstIdx,
	int slave_pixels,
	int slave_lines,
	int MasterBox[4],
	int SlaveBox[4],
	int linesPerBurst,
	int SamplesPerBurst,
	double azimuthTimeInterval,
	double* dopplerRate,
	double* referenceTime,
	double* dopplerCentroid,
	int Npoints
	)
{

	cudaDeviceReset();

	int sLines = SlaveBox[3] - SlaveBox[2] + 1;
	int sPixels = SlaveBox[1] - SlaveBox[0] + 1;
	int mLines = MasterBox[3] - MasterBox[2] + 1;
	int mPixels = MasterBox[1] - MasterBox[0] + 1;




	int Npoints2m1 = Npoints / 2 - 1;


	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 32, 0, 0, cudaChannelFormatKindFloat);
	cudaChannelFormatDesc channelDesc_1 = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);

	int sfirstLineInBurst = sBurstIdx*linesPerBurst;

	cudaHostRegister(SlaveArray, sLines*sPixels*sizeof(short2), cudaHostRegisterDefault);
	cudaHostRegister(dopplerCentroid, SamplesPerBurst*sizeof(double), cudaHostRegisterDefault);
	cudaHostRegister(referenceTime, SamplesPerBurst*sizeof(double), cudaHostRegisterDefault);
	cudaHostRegister(dopplerRate, SamplesPerBurst*sizeof(double), cudaHostRegisterDefault);
	cudaHostRegister(KernelAz, 2048 * 12 * sizeof(float), cudaHostRegisterDefault);
	cudaHostRegister(KernelRg, 2048 * 12 * sizeof(float), cudaHostRegisterDefault);
	cudaHostRegister(output, mLines*mPixels*sizeof(cuComplex), cudaHostRegisterDefault);


	// (cudaHostRegister(PhaseArray, sLines*sPixels*sizeof(float), cudaHostRegisterDefault));
	size_t d_pitch1, d_pitch2, d_pitch3, d_pitchS2;

	size_t TotalBytes = 0;

	// It is worth to use another array to save complex<short>
	float* d_PhaseArray;
	 (cudaMallocPitch((void**)&d_PhaseArray, &d_pitch1, sPixels*sizeof(float), sLines));
	
	

	short2* d_SlaveArrayS2;
	 (cudaMallocPitch((void**)&d_SlaveArrayS2, &d_pitchS2, sPixels*sizeof(short2), sLines));



	cuComplex* d_SlaveArray;
	 (cudaMallocPitch((void**)&d_SlaveArray, &d_pitch2, sPixels*sizeof(cuComplex), sLines));
	

	double* d_dopplerRate, *d_referenceTime, *d_dopplerCentroid;
	 (cudaMallocPitch((void**)&d_dopplerRate, &d_pitch3, SamplesPerBurst*sizeof(double), 1));
	 (cudaMallocPitch((void**)&d_referenceTime, &d_pitch3, SamplesPerBurst*sizeof(double), 1));
	 (cudaMallocPitch((void**)&d_dopplerCentroid, &d_pitch3, SamplesPerBurst*sizeof(double), 1));


	size_t CorrPitch;

	cuComplex * d_resample;
	 (cudaMallocPitch((void **)&d_resample, &CorrPitch, mPixels*sizeof(cuComplex), mLines));
	

	cudaArray *KernelAzArray = NULL;
	cudaArray *KernelRgArray = NULL;
	 (cudaMallocArray(&KernelAzArray, &channelDesc_1, Npoints, 2048));
	 (cudaMallocArray(&KernelRgArray, &channelDesc_1, Npoints, 2048));


	dim3 threads(16, 16);
	dim3 blocks = dim3((sPixels + 15) / 16, (sLines + 15) / 16);
	cudaStream_t stream[4];
	for (int i = 0; i < 4; i++)
	{
		cudaStreamCreate(&stream[i]);
	}
	cudaEvent_t g_start, g_stop;
	//cudaFuncSetCacheConfig(resample_texture_kernel_12p_overlap_warpFunction_test, cudaFuncCachePreferL1);

	//Memcpy to Constant  Variables
	cudaMemcpyToSymbol(npoints, &Npoints, sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_mLines, &mLines, sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_mPixels, &mPixels, sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_sX0, &SlaveBox[0], sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_sXmax, &SlaveBox[1], sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_sY0, &SlaveBox[2], sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_sYmax, &SlaveBox[3], sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_mY0, &MasterBox[2], sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_AzimuthShift, &AzimuthShift, sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_MasterBox, MasterBox, 4 * sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_Npoints2m1, &Npoints2m1, sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_CpmAz, CpmAz, 6 * sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_CpmRg, CpmRg, 6 * sizeof(double), 0, cudaMemcpyHostToDevice);




	

	
	 (cudaMemcpy2D(d_SlaveArrayS2, d_pitchS2, SlaveArray, sPixels*sizeof(short2), sPixels*sizeof(short2), sLines, cudaMemcpyHostToDevice));
	 (cudaMemcpy2D(d_dopplerRate, d_pitch3, dopplerRate, SamplesPerBurst*sizeof(double), SamplesPerBurst*sizeof(double), 1, cudaMemcpyHostToDevice));
	 (cudaMemcpy2D(d_referenceTime, d_pitch3, referenceTime, SamplesPerBurst*sizeof(double), SamplesPerBurst*sizeof(double), 1, cudaMemcpyHostToDevice));
	 (cudaMemcpy2D(d_dopplerCentroid, d_pitch3, dopplerCentroid, SamplesPerBurst*sizeof(double), SamplesPerBurst*sizeof(double), 1, cudaMemcpyHostToDevice));


	//DerampDemod_Shared << <blocks, threads >> >(d_PhaseArray, d_SlaveArray, d_dopplerRate, d_referenceTime, d_dopplerCentroid, sX0, sY0, sLines, sPixels, sfirstLineInBurst,
	//azimuthTimeInterval, d_pitch1, d_pitch2, d_pitch3);

	DerampDemod_Shared << <blocks, threads >> >(d_PhaseArray, d_SlaveArray, d_dopplerRate, d_referenceTime, d_dopplerCentroid, SlaveBox[0], SlaveBox[2], sLines, sPixels, sfirstLineInBurst,
		azimuthTimeInterval, d_pitch1, d_pitch2, d_pitch3, d_SlaveArrayS2, d_pitchS2);

	 (cudaDeviceSynchronize());

	//cudaEventRecord(g_stop, 0);
	//cudaEventSynchronize(g_stop);
	//cudaEventElapsedTime(&time_cost1, g_start, g_stop);
	//cout << "DeRamping duration:" << time_cost1 << endl;







	cudaMemcpyToArray(KernelAzArray, 0, 0, KernelAz, Npoints * 2048 * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpyToArray(KernelRgArray, 0, 0, KernelRg, Npoints * 2048 * sizeof(float), cudaMemcpyHostToDevice);
	cudaBindTextureToArray(tex_kernelAz, KernelAzArray, channelDesc_1);
	cudaBindTextureToArray(tex_kernelRg, KernelRgArray, channelDesc_1);
	cudaBindTexture2D(0, tex_PhaseArray, d_PhaseArray, channelDesc_1, sPixels, sLines, d_pitch1);
	cudaBindTexture2D(0, tex_slave, d_SlaveArray, channelDesc, sPixels, sLines, d_pitch2);

	//float *d_KernelAz, *d_KernelRg;
	//cudaMalloc((void**)&d_KernelAz, 2048 * Npoints*sizeof(float));
	//cudaMalloc((void**)&d_KernelRg, 2048 * Npoints*sizeof(float));
	// (cudaMemcpy(d_KernelAz, KernelAz, 2048 * Npoints*sizeof(float), cudaMemcpyHostToDevice));
	// (cudaMemcpy(d_KernelRg, KernelRg, 2048 * Npoints*sizeof(float), cudaMemcpyHostToDevice));


	size_t SPitch, MPitch;


	threads = dim3(16, 16);
	blocks = dim3((mPixels + threads.x - 1) / threads.x,
		(mLines + threads.y - 1) / threads.y);



	//for Subsets
	int partMlines = mLines / 4;
	int RemainMlines = mLines % 4;
	dim3 Partblocks = dim3((mPixels + threads.x - 1) / threads.x,
		(partMlines + threads.y - 1) / threads.y);
	dim3 Lastblocks = dim3((mPixels + threads.x - 1) / threads.x,
		(partMlines + RemainMlines + threads.y - 1) / threads.y);

	int PartOffsetD = partMlines*CorrPitch / 8;
	int PartOffsetH = partMlines*mPixels;




	//cudaFuncSetCacheConfig(Interpolation_12p, cudaFuncCachePreferL1); Set for more L1 cache. You can comment out to test. This is for Kepler architecture at 


	float time_cost1, time_cost2;
	cudaEventCreate(&g_start);
	cudaEventCreate(&g_stop);
	cudaEventRecord(g_start, 0);

	 (cudaDeviceSynchronize());
	if (Npoints == 6)
	{

	}


	if (Npoints == 12)
	{
		bool  OverlapOrNot = false; //Using concurrent computing or not
		if (!OverlapOrNot)
		{
			/*resample_texture_kernel_12p_warpFunction_raw << <blocks, threads >> >(d_resample, d_SlaveArray, d_PhaseArray,
			d_KernelAz, d_KernelRg, d_pitch1, d_pitch2,mLines, CorrPitch);*/
			//resample_kernel_12p_textureTest << <blocks, threads >> >(d_resample, 0, mLines, CorrPitch);

			//for whole test
			Interpolation_12p
				<< <blocks, threads >> >(d_resample, 0,
				mLines, CorrPitch);

		

		}

		//SubSet1

		if (OverlapOrNot)
		{
			Interpolation_12p
				<< <Partblocks, threads, 0, stream[0] >> >(d_resample, 0,
				partMlines, CorrPitch);

			cudaMemcpy2DAsync(output, (mPixels)*sizeof(cuComplex), d_resample, CorrPitch, (mPixels)*sizeof(cuComplex), partMlines, cudaMemcpyDeviceToHost, stream[0]);



			//SubSet2
			Interpolation_12p
				<< <Partblocks, threads, 0, stream[1] >> >
				(d_resample + PartOffsetD, partMlines,
				partMlines, CorrPitch);
			cudaMemcpy2DAsync(output + PartOffsetH, (mPixels)*sizeof(cuComplex),
				d_resample + PartOffsetD, CorrPitch, (mPixels)*sizeof(cuComplex),
				partMlines, cudaMemcpyDeviceToHost, stream[1]);

			//SubSet3

			Interpolation_12p
				<< <Partblocks, threads, 0, stream[2] >> >
				(d_resample + 2 * PartOffsetD, 2 * partMlines,
				partMlines, CorrPitch);
			cudaMemcpy2DAsync(output + 2 * PartOffsetH, (mPixels)*sizeof(cuComplex),
				d_resample + 2 * PartOffsetD, CorrPitch, (mPixels)*sizeof(cuComplex),
				partMlines, cudaMemcpyDeviceToHost, stream[2]);


			//SubSet4
			Interpolation_12p
				<< <Lastblocks, threads, 0, stream[3] >> >
				(d_resample + 3 * PartOffsetD, 3 * partMlines,
				partMlines + RemainMlines, CorrPitch);
			cudaMemcpy2DAsync(output + 3 * PartOffsetH, (mPixels)*sizeof(cuComplex),
				d_resample + 3 * PartOffsetD, CorrPitch, (mPixels)*sizeof(cuComplex),
				(partMlines + RemainMlines), cudaMemcpyDeviceToHost, stream[3]);

		}



	}
	

	 (cudaDeviceSynchronize());
	cudaEventRecord(g_stop, 0);
	cudaEventSynchronize(g_stop);
	cudaEventElapsedTime(&time_cost2, g_start, g_stop);
	cout << "kernel duration:" << time_cost2<<"ms" << endl;
	cudaEventDestroy(g_start);
	cudaEventDestroy(g_stop);

	cudaMemcpy2D(output, mPixels*sizeof(cuComplex),
		d_resample, CorrPitch, mPixels*sizeof(cuComplex),
		mLines, cudaMemcpyDeviceToHost);


	for (int i = 0; i < 4; i++)
	{
		cudaStreamDestroy(stream[i]);
	}

	 (cudaDeviceSynchronize());

	cudaHostUnregister(SlaveArray);
	cudaHostUnregister(output);
	cudaHostUnregister(dopplerCentroid);
	cudaHostUnregister(dopplerRate);
	cudaHostUnregister(referenceTime);
	cudaHostUnregister(KernelAz);
	cudaHostUnregister(KernelRg);

	cudaUnbindTexture(tex_kernelAz);
	cudaUnbindTexture(tex_kernelRg);
	cudaUnbindTexture(tex_slave);
	cudaUnbindTexture(tex_PhaseArray);

	cudaFree(d_PhaseArray);
	cudaFree(d_SlaveArray);
	cudaFree(d_dopplerRate);
	cudaFree(d_referenceTime);
	cudaFree(d_dopplerCentroid);
	// (cudaFree(d_resample));
	cudaFreeArray(KernelAzArray);
	cudaFreeArray(KernelRgArray);
	cudaFree(d_SlaveArrayS2);



	cuComplex *d_output = d_resample;
	return d_output;

}

cuComplex* ResampleFirstBurst(
	complex<float>*SlaveArray,
	int ww,
	int hh
	)
{
	
	//Page-Locking host Memory
	cudaHostRegister(SlaveArray, ww*hh*sizeof(complex<float>), cudaHostRegisterDefault);
	
	size_t RePitch;
	cuComplex * d_Reslave;
	cudaMallocPitch((void **)&d_Reslave, &RePitch, ww*sizeof(cuComplex), hh);
	cudaMemcpy2D(d_Reslave, RePitch, SlaveArray, 
		ww*sizeof(cuComplex), ww*sizeof(cuComplex), hh, cudaMemcpyHostToDevice);

	cudaHostUnregister(SlaveArray);
	return d_Reslave;
}