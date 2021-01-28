#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#define __CUDA_INTERNAL_COMPILATION__
#include <math_functions.h>
#include <math_constants.h>
#include <device_functions.h>
#include <cuComplex.h>
#include "cuConstants.cuh"
#undef  __CUDA_INTERNAL_COMPILATION__
#include <math.h>
#include <iostream>
#include<omp.h>
#include <complex>
#include <math.h>
#include<time.h>
#include<complex>
#include <iostream>
//#include <helper_cuda.h>






#include <stdio.h>

using namespace std;
extern "C" void EstESDShifts(
	complex<short>*MasterFor, //Input
	complex<short>*MasterBack, //Input
	complex<float>*SlaveFor,// Input
	complex<float>*SlaveBack, //Input
	int* OverlapSizeArray, //Input
	double * ShiftArray,//Output
	int numOverlap,
	int ArrayLines,
	int ArrayPixels,
	int CohWinAz,
	int CohWinRg,
	int OverlapMaxLines,
	float CohThresHold,
	double spectralSeparation,
	double azimuthTimeInterval
	);



int  MultipleOf8_ESD(int inputLines)
{
	int remains = inputLines % 8;
	int multiple = inputLines / 8;
	if (remains == 0)
		return inputLines;
	else{

		return (multiple + 1) * 8;


	}
}

__global__ void DoubleDifference4OneD(short2 *MasterFor, short2 *MasterBack,
	cuComplex * SlaveFor, cuComplex * SlaveBack, unsigned int Lines, unsigned int Pixels)
{
	const int row = blockIdx.y*blockDim.y + threadIdx.y;
	const int col = blockIdx.x*blockDim.x + threadIdx.x;


	if (row<Lines&&col<Pixels)
	{
		unsigned int i = row*Pixels + col;

		short2 Temp1 = MasterFor[i];
		short2 Temp2 = MasterBack[i];


	
		

		SlaveBack[i] = cuCmulf(cuCmulf(make_cuComplex(Temp1.x, Temp1.y), cuConjf(SlaveFor[i])),
			cuConjf(cuCmulf(make_cuComplex(Temp2.x, Temp2.y), cuConjf(SlaveBack[i]))));



	

	
	}


}


__global__ void CohMatrixSet(short2 *MasterFor, cuComplex * SlaveFor,
	cuComplex *MasterCoh, cuComplex *SlaveCoh,
	int inputLines, int inputPixels, int dh, int dw, size_t PitchCoh)
{
	const int row = blockIdx.y*blockDim.y + threadIdx.y;
	const int col = blockIdx.x*blockDim.x + threadIdx.x;


	if (row < inputLines&&col < inputPixels)
	{
		unsigned int id = row*inputPixels + col;


		cuComplex *rowMasterCoh = (cuComplex *)((char*)MasterCoh + (row + dh)*PitchCoh);
		cuComplex *rowSlaveCoh = (cuComplex *)((char*)SlaveCoh + (row + dh)*PitchCoh);

		int regdw = dw;


		cuComplex tempMaster = make_cuComplex(MasterFor[id].x, MasterFor[id].y);
		cuComplex tempSlave = SlaveFor[id];
		rowMasterCoh[col + regdw] = cuCmulf(tempMaster, cuConjf(tempSlave));


		float slavenorm = tempSlave.x*tempSlave.x + tempSlave.y*tempSlave.y;
		float masternorm = tempMaster.x*tempMaster.x + tempMaster.y*tempMaster.y;

		rowSlaveCoh[col + regdw] = make_cuComplex(slavenorm, masternorm);

	

	}


}

__global__ void MaskOut(cuComplex *InputArray, float*Coh, float threshold, int inputLines, int inputPixels,
	size_t FloatPitch)
{
	const int row = blockIdx.y*blockDim.y + threadIdx.y;
	const int col = blockIdx.x*blockDim.x + threadIdx.x;

	if (row < inputLines&&col < inputPixels)
	{
		unsigned int id = row*inputPixels + col;



		float *rowCoh = (float *)((char*)Coh + row*FloatPitch);


		if (rowCoh[col] < threshold || rowCoh[col] >= 1)
		{
			InputArray[id] = make_cuComplex(0.0f, 0.0f);

		}
		else
		{
			double theta = atan2(InputArray[id].y, InputArray[id].x);
			double EstSC[2];
			sincos(theta, EstSC, EstSC + 1);
			InputArray[id] = make_cuComplex(EstSC[1], EstSC[0]);

		}
	}


}

/*Process two pixels in vertical for a thread*/
__global__ void CoherenceBasic42
(float *Result, cuComplex *MasterArray, cuComplex *SlaveArray,
int inputLines, int inputPixels, int FloatPixels, int ComplexPixels,
int winY, int winX)
{


	const int row = blockIdx.y*blockDim.y + threadIdx.y;
	const int col = blockIdx.x*blockDim.x + threadIdx.x;


	if (row * 2 < (inputLines) && col < inputPixels)
	{


		int RwinY = winY;
		int RwinX = winX;
		int RComplexPixels = ComplexPixels;



		cuComplex temp0 = make_cuComplex(0.0f, 0.0f);
		cuComplex sum1 = temp0;
		cuComplex power1 = temp0;
		cuComplex sum1Upper = temp0;
		cuComplex power1Upper = temp0;
		cuComplex sum1Lower = temp0;
		cuComplex power1Lower = temp0;

		for (int iU = col; iU < col + RwinX; iU++)
		{
			sum1Upper = cuCaddf(sum1Upper, MasterArray[row * 2 * RComplexPixels + iU]);
			power1Upper = cuCaddf(power1Upper, SlaveArray[row * 2 * RComplexPixels + iU]);

		}

		for (int j = row * 2 + 1; j < row * 2 + RwinY; j++)
		{
	

			for (int i = col; i <col + RwinX; i++)
			{
				sum1 = cuCaddf(sum1, MasterArray[j*RComplexPixels + i]);
				power1 = cuCaddf(power1, SlaveArray[j*RComplexPixels + i]);



			}

		}

		for (int iL = col; iL < col + RwinX; iL++)
		{
			sum1Lower = cuCaddf(sum1Lower, MasterArray[(row * 2 + RwinY) * RComplexPixels + iL]);
			power1Lower = cuCaddf(power1Lower, SlaveArray[(row * 2 + RwinY) * RComplexPixels + iL]);
		}

		sum1Upper = cuCaddf(sum1Upper, sum1);
		power1Upper = cuCaddf(power1Upper, power1);
		float pU = power1Upper.x*power1Upper.y;
		Result[(row * 2)*FloatPixels + col] = (pU > 0.0f) ?
			cuCabsf(make_cuComplex(sum1Upper.x / sqrtf(pU), sum1Upper.y / sqrtf(pU))
			) : 0.0f;


		sum1Lower = cuCaddf(sum1Lower, sum1);
		power1Lower = cuCaddf(power1Lower, power1);
		float pL = power1Lower.x*power1Lower.y;
		Result[(row * 2 + 1)*FloatPixels + col] = (pL > 0.0f) ?
			cuCabsf(make_cuComplex(sum1Lower.x / sqrtf(pL), sum1Lower.y / sqrtf(pL))
			) : 0.0f;
	}



}



__global__ void Divide2
(float *Result1, float *Result2, size_t FloatPitch, unsigned int Lines, unsigned int Pixels)
{
	const int row = blockIdx.y*blockDim.y + threadIdx.y;
	const int col = blockIdx.x*blockDim.x + threadIdx.x;

	if (row < Lines&&col < Pixels)
	{

		float* rowResult1 = (float *)((char*)Result1 + row*FloatPitch);
		float* rowResult2 = (float *)((char*)Result2 + row*FloatPitch);

		rowResult1[col] = (rowResult1[col] + rowResult2[col]) / 2.0;

	}

}

__global__ void
reduce0(cuComplex *g_idata, cuComplex *g_odata,
unsigned int n)
{
	__shared__ cuComplex sdata[256];


	// load shared mem
	unsigned int tid = threadIdx.x;
	unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;



	sdata[tid] = (i < n) ? g_idata[i] : make_cuComplex(0.0f, 0.0f);


	__syncthreads();

	// do reduction in shared mem
	for (unsigned int s = 1; s < blockDim.x; s *= 2)
	{
		// modulo arithmetic is slow!
		if ((tid % (2 * s)) == 0)
		{
			sdata[tid] = cuCaddf(sdata[tid], sdata[tid + s]);

		}

		__syncthreads();
	}


	// write result for this block to global mem
	if (tid == 0) g_odata[blockIdx.x] = sdata[0], sdata[0];
}

__global__ void CohFirstLineSet_ESD(cuComplex *d_Sum1, cuComplex *d_Power1,
	cuComplex *MasterArray2, cuComplex* SlaveArray2, int winY, int winX,
	int inputPixels, size_t InPitch)
{


	const int col = blockIdx.x*blockDim.x + threadIdx.x;


	if (col < inputPixels)
	{
		int InPixels = InPitch / 8;



		cuComplex sum1 = make_cuComplex(0.0f, 0.0f);
		cuComplex power1 = make_cuComplex(0.0f, 0.0f);

		for (int j = 0; j < winY; j++)
		{


			for (int i = col; i < col + winX; i++)
			{
				sum1 = cuCaddf(sum1,
					MasterArray2[j*InPixels + i]);
				power1 = cuCaddf(power1,
					SlaveArray2[j*InPixels + i]);

			}
		}

		d_Sum1[col] = sum1;
		d_Power1[col] = power1;


	}



}

__global__ void CohDiffer_vertical_ESD(cuComplex *Sum1Array, cuComplex* Power1Array,
	cuComplex *MasterArray2, cuComplex* SlaveArray2, int winY, int winX,
	int inputLines, int inputPixels, size_t SrcPitch, size_t DstPitch)
{


	const int row = blockIdx.y*blockDim.y + threadIdx.y;
	const int col = blockIdx.x*blockDim.x + threadIdx.x;

	if (row < inputLines&&col < inputPixels)
	{

		cuComplex sum1 = make_cuComplex(0.0f, 0.0f);
		cuComplex power1 = make_cuComplex(0.0f, 0.0f);
		int RwinX = winX;
		int RwinY = winY;
		int SrcPixels = SrcPitch / 8;
		int DstPixels = DstPitch / 8;


		for (int i = 0; i < RwinX; i++)
		{
			sum1 = cuCaddf(sum1,
				cuCsubf(MasterArray2[(row + RwinY)*(SrcPixels)+col + i],
				MasterArray2[row *(SrcPixels)+col + i]));

		}
		Sum1Array[(row + 1)*(DstPixels)+col] = sum1;
		for (int i = 0; i < RwinX; i++)
		{
			power1 = cuCaddf(power1, cuCsubf(SlaveArray2[(row + RwinY)*(SrcPixels)+col + i],
				SlaveArray2[row *(SrcPixels)+col + i]));

		}
		Power1Array[(row + 1)*(DstPixels)+col] = power1;





	}

}

//2pixels
__global__ void CohDiffer42_vertical_ESD(cuComplex *Sum1Array, cuComplex* Power1Array,
	cuComplex *MasterArray2, cuComplex* SlaveArray2, int winY, int winX,
	int inputLines, int inputPixels, size_t SrcPitch, size_t DstPitch)
{


	const int row = blockIdx.y*blockDim.y + threadIdx.y;
	const int col = blockIdx.x*blockDim.x + threadIdx.x;

	if (row < inputLines&&col * 2 < inputPixels)
	{

		cuComplex complex0 = make_cuComplex(0.0f, 0.0f);
		cuComplex sum1 = complex0;
		cuComplex power1 = complex0;
		cuComplex sumLeft = complex0;
		cuComplex powerLeft = complex0;
		cuComplex sumRight = complex0;
		cuComplex powerRight = complex0;
		int RwinX = winX;
		int RwinY = winY;
		int SrcPixels = SrcPitch / 8;
		int DstPixels = DstPitch / 8;

		int IdxWinY = row + RwinY;

		sumLeft = cuCsubf(MasterArray2[IdxWinY*SrcPixels + col * 2],
			MasterArray2[row *SrcPixels + col * 2]);
		for (int i = col * 2 + 1; i <col * 2 + RwinX; i++)
		{
			sum1 = cuCaddf(sum1,
				cuCsubf(MasterArray2[IdxWinY*SrcPixels + i],
				MasterArray2[row *SrcPixels + i]));

		}
		sumRight = cuCsubf(MasterArray2[IdxWinY*SrcPixels + col * 2 + RwinX],
			MasterArray2[row *SrcPixels + col * 2 + RwinX]);

		sumLeft = cuCaddf(sum1, sumLeft);
		sumRight = cuCaddf(sum1, sumRight);

		Sum1Array[(row + 1)*(DstPixels)+col * 2] = sumLeft;
		Sum1Array[(row + 1)*(DstPixels)+col * 2 + 1] = sumRight;

		powerLeft = cuCsubf(SlaveArray2[IdxWinY*SrcPixels + col * 2],
			SlaveArray2[row *SrcPixels + col * 2]);
		for (int i = col * 2 + 1; i <col * 2 + RwinX; i++)
		{
			power1 = cuCaddf(power1, cuCsubf(SlaveArray2[IdxWinY*SrcPixels + i],
				SlaveArray2[row *SrcPixels + i]));

		}
		powerRight = cuCsubf(SlaveArray2[IdxWinY*SrcPixels + col * 2 + RwinX],
			SlaveArray2[row *SrcPixels + col * 2 + RwinX]);
		powerLeft = cuCaddf(power1, powerLeft);
		powerRight = cuCaddf(power1, powerRight);

		Power1Array[(row + 1)*(DstPixels)+col * 2] = powerLeft;
		Power1Array[(row + 1)*(DstPixels)+col * 2 + 1] = powerRight;



	}

}



__global__ void shfl_vertical_ESD(cuComplex *Src, int width, int height, int Nstep)
{
	__shared__ float sums_x[32][9];
	__shared__ float sums_y[32][9];
	int tidx = blockIdx.x * blockDim.x + threadIdx.x;


	//int warp_id = threadIdx.x / warpSize ;
	unsigned int lane_id = tidx % 8;
	cuComplex stepSum = make_cuComplex(0.0f, 0.0f);

	sums_x[threadIdx.x][threadIdx.y] = 0;
	sums_y[threadIdx.x][threadIdx.y] = 0;
	__syncthreads();

	for (int step = 0; step < Nstep; step++)
	{
		cuComplex sum = make_cuComplex(0.0f, 0.0f);
		cuComplex *p = Src + (threadIdx.y + step * 8)*width + tidx;

		sum = *p;
		sums_x[threadIdx.x][threadIdx.y] = sum.x;
		sums_y[threadIdx.x][threadIdx.y] = sum.y;
		__syncthreads();

		// place into SMEM
		// shfl scan reduce the SMEM, reformating so the column
		// sums are computed in a warp
		// then read out properly
		float partial_sum_x = 0;
		float partial_sum_y = 0;
		int j = threadIdx.x % 8;
		int k = threadIdx.x / 8 + threadIdx.y * 4;

		partial_sum_x = sums_x[k][j];
		partial_sum_y = sums_y[k][j];

		for (int i = 1; i <= 8; i *= 2)
		{
			float n_x = __shfl_up_sync(0xffffffff, partial_sum_x, i, 32);
			float n_y = __shfl_up_sync(0xffffffff, partial_sum_y, i, 32);// removed for debug
			//float n_x = 0;
			//float n_y = 0;
			if (lane_id >= i)
			{
				partial_sum_x += n_x;
				partial_sum_y += n_y;
			}
		}

		sums_x[k][j] = partial_sum_x;
		sums_y[k][j] = partial_sum_y;
		__syncthreads();

		if (threadIdx.y > 0)
		{
			sum.x += sums_x[threadIdx.x][threadIdx.y - 1];
			sum.y += sums_y[threadIdx.x][threadIdx.y - 1];
		}

		sum.x += stepSum.x;
		sum.y += stepSum.y;
		stepSum.x += sums_x[threadIdx.x][blockDim.y - 1];
		stepSum.y += sums_y[threadIdx.x][blockDim.y - 1];
		__syncthreads();
		*p = sum;
	}

};


__global__ void CpCoh2_ESD(float *Res, cuComplex *Sum1Array, cuComplex *Power1Array, int OutputLines, int OutputPixels,
	size_t InPitch, size_t Outpitch)
{
	int row = blockIdx.y*blockDim.y + threadIdx.y;
	int col = blockIdx.x*blockDim.x + threadIdx.x;


	if (row < OutputLines&&col < OutputPixels)
	{

		cuComplex* rowSum = (cuComplex *)((char*)Sum1Array + row*InPitch);
		cuComplex* rowPower = (cuComplex *)((char*)Power1Array + row*InPitch);
		cuComplex tmpSum = rowSum[col];
		cuComplex tmpPower = rowPower[col];

		float pp = tmpPower.x*tmpPower.y;

		cuComplex TmpRes = (pp > 0.0) ? (cuCdivf(tmpSum, make_cuComplex(sqrt(pp), 0.0f))) : make_cuComplex(0.0f, 0.0f);;

		Res[row*(Outpitch / 4) + col] = cuCabsf(TmpRes);


	}
}






void EstESDShifts_Overlap(
	complex<short>*MasterFor, //Input
	complex<short>*MasterBack, //Input
	complex<float>*SlaveFor,// Input
	complex<float>*SlaveBack, //Input
	int* OverlapSizeArray, //Input
	double * ShiftArray,//Output
	int numOverlap,
	int ArrayLines,
	int ArrayPixels,
	int CohWinAz,
	int CohWinRg,
	int OverlapMaxLines,
	float CohThresHold,
	double spectralSeparation,
	double azimuthTimeInterval
	)
{

	short2 *d_MasterFor, *d_MasterBack;
	cuComplex *d_SlaveFor, *d_SlaveBack;


	unsigned int ArraySize = ArrayLines*ArrayPixels;

	

	cudaHostRegister(MasterFor, ArraySize*sizeof(short2), cudaHostRegisterDefault);
	cudaHostRegister(MasterBack, ArraySize*sizeof(short2), cudaHostRegisterDefault);
	cudaHostRegister(SlaveFor, ArraySize*sizeof(cuComplex), cudaHostRegisterDefault);
	cudaHostRegister(SlaveBack, ArraySize*sizeof(cuComplex), cudaHostRegisterDefault);



	size_t d_pitch4short2, //pitch for short2
		d_pitch4float2; //pitch for float2

	OverlapMaxLines = MultipleOf8_ESD(OverlapMaxLines);
	


	//
	cuComplex *SumOverlap;
	cudaMallocHost((void**)&SumOverlap, numOverlap*sizeof(cuComplex));
	size_t d_pitchCohfloat, d_pitchCohfloat2;
	size_t MasterShortPitch, SlaveComplexPitch;
	cuComplex *d_MasterCoh, *d_SlaveCoh;

	float *d_coh1, *d_coh2;



	int cohw = ArrayPixels + CohWinRg - 1;
	int cohh = OverlapMaxLines + CohWinAz - 1;



	int dw = (CohWinRg - 1) / 2;
	int dh = (CohWinAz - 1) / 2;

	
	
	//for estimating coherence
	(cudaMallocPitch((void**)&d_MasterCoh, &d_pitchCohfloat2, cohw*sizeof(cuComplex), cohh*numOverlap));
	(cudaMallocPitch((void**)&d_SlaveCoh, &d_pitchCohfloat2, cohw*sizeof(cuComplex), cohh*numOverlap));

	(cudaMallocPitch((void**)&d_coh1, &d_pitchCohfloat, ArrayPixels*sizeof(float), OverlapMaxLines*numOverlap));
	(cudaMallocPitch((void**)&d_coh2, &d_pitchCohfloat, ArrayPixels*sizeof(float), OverlapMaxLines*numOverlap));




	int FloatPixels = d_pitchCohfloat / 4;
	int ComplexPixels = d_pitchCohfloat2 / 8;

	//Set Reduce Parameters
	cuComplex *d_blockSum, *d_blockSum1, *d_blockSum2;
	int threadNum = 256;
	int blockNum = (OverlapMaxLines*ArrayPixels + threadNum - 1) / threadNum;

	int blockNum1 = (blockNum + 255) / threadNum;

	(cudaMalloc((void**)&d_blockSum, numOverlap*blockNum*sizeof(cuComplex)));
	(cudaMalloc((void**)&d_blockSum1, numOverlap*blockNum1*sizeof(cuComplex)));
	(cudaMalloc((void**)&d_blockSum2, numOverlap*sizeof(cuComplex)));



	cudaMalloc((void**)&d_MasterFor, ArraySize*sizeof(short2));
	cudaMalloc((void**)&d_MasterBack, ArraySize*sizeof(short2));

	size_t d_pitchCohTmp;
	cuComplex * d_sum1, *d_power1;
	cudaMallocPitch((void**)&d_sum1, &d_pitchCohTmp, ArrayPixels*sizeof(cuComplex), OverlapMaxLines*numOverlap);
	cudaMallocPitch((void**)&d_power1, &d_pitchCohTmp, ArrayPixels*sizeof(cuComplex), OverlapMaxLines*numOverlap);

	//d_MasterBack = d_MasterFor + ArraySize;

	cudaMalloc((void**)&d_SlaveFor, ArraySize*sizeof(cuComplex));
	cudaMalloc((void**)&d_SlaveBack, ArraySize*sizeof(cuComplex));
	//d_SlaveBack = d_SlaveFor + ArraySize;








	dim3 threads(16, 16);
	dim3 blocks((ArrayPixels + 15) / 16, (ArrayLines + 15) / 16);


	int *PrefixOvlines = new int[numOverlap];


	PrefixOvlines[0] = 0;
	for (int i = 1; i < numOverlap; i++)
	{
		PrefixOvlines[i] = 0;
		for (int j = 0; j < i; j++)
		{
			PrefixOvlines[i] += OverlapSizeArray[j];
		}

	}






	dim3 blockSz(32, 8);
	int PitchPixels = d_pitchCohTmp / sizeof(cuComplex);
	//For horizontal
	//int Nstep = PitchPixels / 8;

	if (PitchPixels% blockSz.x != 0)
	{
		cout << "Errors:PitchSize is not multiple of 32!" << endl;
		return;

	}
	int numX = PitchPixels;
	dim3 GridPrefix(numX / blockSz.x, 1);
	int PartSteps = OverlapMaxLines / 8;



	bool OverlaporNot = false;

	int ArrayPixels42 = ArrayPixels % 2 ? ArrayPixels + 1 : ArrayPixels;
	dim3 threadsV(128, 2);
	cudaEvent_t g_start, g_stop;
	float time_cost;
	cudaEventCreate(&g_start);
	cudaEventCreate(&g_stop);
	cudaEventRecord(g_start, 0);

	if (OverlaporNot)
	{



		cudaStream_t stream[12];
		for (int i = 0; i < numOverlap; i++)
		{
			cudaStreamCreate(&stream[i]);
		}

		for (int i = 0; i < numOverlap; i++)
		{
			int overlapLines = OverlapSizeArray[i];
			int overlapN = overlapLines*ArrayPixels;
			int LinesOffset = PrefixOvlines[i];
			int DeviceOffset = LinesOffset*ArrayPixels;
			int CohOffset2 = i*cohh*d_pitchCohfloat2 / sizeof(cuComplex);
			int CohOffset1 = OverlapMaxLines*i*d_pitchCohfloat / sizeof(float);
			int SumOffset = OverlapMaxLines*i*d_pitchCohTmp / sizeof(cuComplex);


			dim3 PartBlocksDiff((ArrayPixels + threadsV.x - 1) / threadsV.x,
				(overlapLines - 1 + threadsV.y - 1) / threadsV.y);
			dim3 PartBlocksDiff42((ArrayPixels42 / 2 + threadsV.x - 1) / threadsV.x,
				(OverlapMaxLines - 1 + threadsV.y - 1) / threadsV.y);

			int tmpBlockNum = (overlapLines*ArrayPixels + threadNum - 1) / threadNum;
			int s = (tmpBlockNum + 255) / 256;
			int s1 = (s + 255) / 256;

			dim3 dimGrid = dim3(tmpBlockNum, 1, 1);
			dim3 dimBlocks = dim3(threadNum, 1, 1);

			blocks = dim3((ArrayPixels + 15) / 16, (overlapLines + 15) / 16);
			dim3 threads42 = dim3(64, 4);
			dim3 blocks42 = dim3((ArrayPixels + threads42.x - 1) / threads42.x, (OverlapMaxLines / 2 + threads42.y - 1) / threads42.y);


			cudaMemcpyAsync(d_MasterFor + DeviceOffset, MasterFor + DeviceOffset,
				overlapN*sizeof(short2), cudaMemcpyHostToDevice, stream[i]);
			cudaMemcpyAsync(d_MasterBack + DeviceOffset, MasterBack + DeviceOffset,
				overlapN*sizeof(short2), cudaMemcpyHostToDevice, stream[i]);
			cudaMemcpyAsync(d_SlaveFor + DeviceOffset, SlaveFor + DeviceOffset,
				overlapN*sizeof(cuComplex), cudaMemcpyHostToDevice, stream[i]);
			cudaMemcpyAsync(d_SlaveBack + DeviceOffset, SlaveBack + DeviceOffset,
				overlapN*sizeof(cuComplex), cudaMemcpyHostToDevice, stream[i]);



			//cudaMemsetAsync(d_MasterCoh + CohOffset2, 0, d_pitchCohfloat2*cohh, stream[i]);
			//cudaMemsetAsync(d_SlaveCoh + CohOffset2, 0, d_pitchCohfloat2*cohh, stream[i]);

			//Estimating the coherence of forward overlap band 
			CohMatrixSet << <blocks, threads, 0, stream[i] >> >(d_MasterFor + DeviceOffset, d_SlaveFor + DeviceOffset,
				d_MasterCoh + CohOffset2, d_SlaveCoh + CohOffset2,
				overlapLines, ArrayPixels, dh, dw, d_pitchCohfloat2);

			CohFirstLineSet_ESD << <(ArrayPixels + 255) / 256, 256, 0, stream[i] >> >(d_sum1 + SumOffset,
				d_power1 + SumOffset, d_MasterCoh + CohOffset2, d_SlaveCoh + CohOffset2, CohWinAz, CohWinRg,
				ArrayPixels, d_pitchCohfloat2);


			CohDiffer_vertical_ESD << <PartBlocksDiff, threadsV, 0, stream[i] >> >(d_sum1 + SumOffset,
				d_power1 + SumOffset, d_MasterCoh + CohOffset2, d_SlaveCoh + CohOffset2, CohWinAz, CohWinRg,
				overlapLines - 1, ArrayPixels,
				d_pitchCohfloat2, d_pitchCohTmp);




			shfl_vertical_ESD << <GridPrefix, blockSz, 0, stream[i] >> >(d_sum1 + SumOffset, numX, OverlapMaxLines, PartSteps);
			shfl_vertical_ESD << <GridPrefix, blockSz, 0, stream[i] >> >(d_power1 + SumOffset, numX, OverlapMaxLines, PartSteps);

			CpCoh2_ESD << <blocks, threads, 0, stream[i] >> >(d_coh1 + CohOffset1,
				d_sum1 + SumOffset, d_power1 + SumOffset, overlapLines, ArrayPixels, d_pitchCohTmp, d_pitchCohfloat);


			//Estimating the coherence of backward overlap band 

			CohMatrixSet << <blocks, threads, 0, stream[i] >> >(d_MasterBack + DeviceOffset, d_SlaveBack + DeviceOffset,
				d_MasterCoh + CohOffset2, d_SlaveCoh + CohOffset2,
				overlapLines, ArrayPixels, dh, dw, d_pitchCohfloat2);


			CohFirstLineSet_ESD << <(ArrayPixels + 255) / 256, 256, 0, stream[i] >> >(d_sum1 + SumOffset,
				d_power1 + SumOffset, d_MasterCoh + CohOffset2, d_SlaveCoh + CohOffset2, CohWinAz, CohWinRg,
				ArrayPixels, d_pitchCohfloat2);

			CohDiffer_vertical_ESD << <PartBlocksDiff, threadsV, 0, stream[i] >> >(d_sum1 + SumOffset,
				d_power1 + SumOffset, d_MasterCoh + CohOffset2, d_SlaveCoh + CohOffset2, CohWinAz, CohWinRg,
				overlapLines - 1, ArrayPixels, d_pitchCohfloat2, d_pitchCohTmp);



			shfl_vertical_ESD << <GridPrefix, blockSz, 0, stream[i] >> >(d_sum1 + SumOffset, numX, OverlapMaxLines, PartSteps);
			shfl_vertical_ESD << <GridPrefix, blockSz, 0, stream[i] >> >(d_power1 + SumOffset, numX, OverlapMaxLines, PartSteps);

			CpCoh2_ESD << <blocks, threads, 0, stream[i] >> >(d_coh2 + CohOffset1,
				d_sum1 + SumOffset, d_power1 + SumOffset, overlapLines, ArrayPixels, d_pitchCohTmp, d_pitchCohfloat);



			//Average the coherence of forward and backward bands
			Divide2 << <blocks, threads, 0, stream[i] >> > (d_coh1 + CohOffset1, d_coh2 + CohOffset1, d_pitchCohfloat, overlapLines, ArrayPixels);
			//cudaStreamSynchronize(stream[i]);
			//Estimte the double differential interferogram
			DoubleDifference4OneD << <blocks, threads, 0, stream[i] >> >
				(d_MasterFor + DeviceOffset, d_MasterBack + DeviceOffset,
				d_SlaveFor + DeviceOffset, d_SlaveBack + DeviceOffset,
				overlapLines, ArrayPixels);

			//Masking out incoherent pixels
			MaskOut << <blocks, threads, 0, stream[i] >> >(d_SlaveBack + DeviceOffset, d_coh1 + CohOffset1, CohThresHold, overlapLines, ArrayPixels,
				d_pitchCohfloat);

			//Using parallel recuding operation to accumulate ESD phases

			reduce0 << <tmpBlockNum, 256, 0, stream[i] >> >(d_SlaveBack + DeviceOffset, d_blockSum + i*blockNum, overlapLines*ArrayPixels);


			reduce0 << <s, 256, 0, stream[i] >> >(d_blockSum + i*blockNum, d_blockSum1 + i*blockNum1, tmpBlockNum);


			reduce0 << <s1, 256, 0, stream[i] >> >(d_blockSum1 + i*blockNum1, d_blockSum2 + i, (tmpBlockNum + 255) / 256);


			//cudaStreamSynchronize(stream[i]);

		}

		for (int i = 0; i < numOverlap; i++)
		{
			cudaStreamDestroy(stream[i]);
		}
	}

	

	cudaMemcpy(SumOverlap, d_blockSum2, numOverlap * sizeof(cuComplex), cudaMemcpyDeviceToHost);

	cudaEventRecord(g_stop, 0);
	cudaEventSynchronize(g_stop);
	cudaEventElapsedTime(&time_cost, g_start, g_stop);
	cout << "ESDCorretion duration:" << time_cost<<"ms" << endl;
	cudaEventDestroy(g_start);
	cudaEventDestroy(g_stop);
	//ofstream ShiftOut("E:\\ESDTimer.txt", ios::app);
	//ShiftOut << time_cost << endl;
	//ShiftOut.close();

	for (int i = 0; i < numOverlap; i++)
	{
	
		ShiftArray[i] = atan2((double)SumOverlap[i].y, (double)SumOverlap[i].x) / (2 * CUDART_PI*spectralSeparation*azimuthTimeInterval);
	}

	




	cudaHostUnregister(MasterFor);
	cudaHostUnregister(MasterBack);
	cudaHostUnregister(SlaveFor);
	cudaHostUnregister(SlaveBack);

	cudaFree(d_MasterFor);
	cudaFree(d_MasterBack);
	cudaFree(d_SlaveFor);
	cudaFree(d_SlaveBack);



	cudaFree(d_sum1);
	cudaFree(d_power1);
	cudaFree(d_MasterCoh);
	cudaFree(d_SlaveCoh);
	cudaFree(d_coh1);
	cudaFree(d_coh2);
	cudaFree(d_blockSum);



	//delete[] SumOverlap;
	cudaFreeHost(SumOverlap);
	
	cudaDeviceSynchronize();
	int cuda_last_error_flag = cudaGetLastError();

	if (cuda_last_error_flag != 0)
	{
		printf("Returns an error code (%d) from GPU execution!\n Please debug it or report it to developer!\n", cuda_last_error_flag);
		exit(0);
	}
	
	cudaDeviceReset();
	delete[] PrefixOvlines;
	

}

void EstESDShifts_NonOverlap(
	complex<short>*MasterFor, //Input
	complex<short>*MasterBack, //Input
	complex<float>*SlaveFor,// Input
	complex<float>*SlaveBack, //Input
	int* OverlapSizeArray, //Input
	double * ShiftArray,//Output
	int numOverlap,
	int ArrayLines,
	int ArrayPixels,
	int CohWinAz,
	int CohWinRg,
	int OverlapMaxLines,
	float CohThresHold,
	double spectralSeparation,
	double azimuthTimeInterval
	)
{
	short2 *d_MasterFor, *d_MasterBack;
	cuComplex *d_SlaveFor, *d_SlaveBack;


	unsigned int ArraySize = ArrayLines*ArrayPixels;



	cudaHostRegister(MasterFor, ArraySize*sizeof(short2), cudaHostRegisterDefault);
	cudaHostRegister(MasterBack, ArraySize*sizeof(short2), cudaHostRegisterDefault);
	cudaHostRegister(SlaveFor, ArraySize*sizeof(cuComplex), cudaHostRegisterDefault);
	cudaHostRegister(SlaveBack, ArraySize*sizeof(cuComplex), cudaHostRegisterDefault);



	size_t d_pitch4short2, //pitch for short2
		d_pitch4float2; //pitch for float2

	OverlapMaxLines = MultipleOf8_ESD(OverlapMaxLines);



	//
	cuComplex *SumOverlap;
	cudaMallocHost((void**)&SumOverlap, numOverlap*sizeof(cuComplex));
	size_t d_pitchCohfloat, d_pitchCohfloat2;
	size_t MasterShortPitch, SlaveComplexPitch;
	cuComplex *d_MasterCoh, *d_SlaveCoh;

	float *d_coh1, *d_coh2;



	int cohw = ArrayPixels + CohWinRg - 1;
	int cohh = OverlapMaxLines + CohWinAz - 1;



	int dw = (CohWinRg - 1) / 2;
	int dh = (CohWinAz - 1) / 2;



	//for estimating coherence
	cudaMallocPitch((void**)&d_MasterCoh, &d_pitchCohfloat2, cohw*sizeof(cuComplex), cohh);
	(cudaMallocPitch((void**)&d_SlaveCoh, &d_pitchCohfloat2, cohw*sizeof(cuComplex), cohh));

	(cudaMallocPitch((void**)&d_coh1, &d_pitchCohfloat, ArrayPixels*sizeof(float), OverlapMaxLines));
	(cudaMallocPitch((void**)&d_coh2, &d_pitchCohfloat, ArrayPixels*sizeof(float), OverlapMaxLines));




	int FloatPixels = d_pitchCohfloat / 4;
	int ComplexPixels = d_pitchCohfloat2 / 8;

	//Set Reduce Parameters
	cuComplex *d_blockSum, *d_blockSum1, *d_blockSum2;
	int threadNum = 256;
	int blockNum = (OverlapMaxLines*ArrayPixels + threadNum - 1) / threadNum;

	int blockNum1 = (blockNum + 255) / threadNum;

	(cudaMalloc((void**)&d_blockSum, blockNum*sizeof(cuComplex)));
	(cudaMalloc((void**)&d_blockSum1, blockNum1*sizeof(cuComplex)));
	(cudaMalloc((void**)&d_blockSum2, sizeof(cuComplex)));


	int Maxize = OverlapMaxLines*ArrayPixels;
	cudaMalloc((void**)&d_MasterFor, Maxize*sizeof(short2));
	cudaMalloc((void**)&d_MasterBack, Maxize*sizeof(short2));

	size_t d_pitchCohTmp;
	cuComplex * d_sum1, *d_power1;
	cudaMallocPitch((void**)&d_sum1, &d_pitchCohTmp, ArrayPixels*sizeof(cuComplex), OverlapMaxLines);
	cudaMallocPitch((void**)&d_power1, &d_pitchCohTmp, ArrayPixels*sizeof(cuComplex), OverlapMaxLines);

	//d_MasterBack = d_MasterFor + ArraySize;

	cudaMalloc((void**)&d_SlaveFor, Maxize*sizeof(cuComplex));
	cudaMalloc((void**)&d_SlaveBack, Maxize*sizeof(cuComplex));
	//d_SlaveBack = d_SlaveFor + ArraySize;








	dim3 threads(16, 16);
	dim3 blocks((ArrayPixels + 15) / 16, (ArrayLines + 15) / 16);


	int *PrefixOvlines = new int[numOverlap];


	PrefixOvlines[0] = 0;
	for (int i = 1; i < numOverlap; i++)
	{
		PrefixOvlines[i] = 0;
		for (int j = 0; j < i; j++)
		{
			PrefixOvlines[i] += OverlapSizeArray[j];
		}

	}






	dim3 blockSz(32, 8);
	int PitchPixels = d_pitchCohTmp / sizeof(cuComplex);
	//For horizontal
	//int Nstep = PitchPixels / 8;

	if (PitchPixels% blockSz.x != 0)
	{
		cout << "Errors:PitchSize is not multiple of 32!" << endl;
		return;

	}
	int numX = PitchPixels;
	dim3 GridPrefix(numX / blockSz.x, 1);
	int PartSteps = OverlapMaxLines / 8;



	bool OverlaporNot = false;

	int ArrayPixels42 = ArrayPixels % 2 ? ArrayPixels + 1 : ArrayPixels;
	dim3 threadsV(128, 2);
	cudaEvent_t g_start, g_stop;
	float time_cost;
	cudaEventCreate(&g_start);
	cudaEventCreate(&g_stop);
	cudaEventRecord(g_start, 0);
	if (!OverlaporNot)
	{
		for (int i = 0; i < numOverlap; i++)
		{
			int overlapLines = OverlapSizeArray[i];
			int overlapN = overlapLines*ArrayPixels;
			int LinesOffset = PrefixOvlines[i];
			int DeviceOffset = LinesOffset*ArrayPixels;
			int CohOffset2 = i*cohh*d_pitchCohfloat2 / sizeof(cuComplex);
			int CohOffset1 = OverlapMaxLines*i*d_pitchCohfloat / sizeof(float);
			int SumOffset = OverlapMaxLines*i*d_pitchCohTmp / sizeof(cuComplex);


			dim3 PartBlocksDiff((ArrayPixels + threadsV.x - 1) / threadsV.x,
				(overlapLines - 1 + threadsV.y - 1) / threadsV.y);
			dim3 PartBlocksDiff42((ArrayPixels42 / 2 + threadsV.x - 1) / threadsV.x,
				(OverlapMaxLines - 1 + threadsV.y - 1) / threadsV.y);

			int tmpBlockNum = (overlapLines*ArrayPixels + threadNum - 1) / threadNum;
			int s = (tmpBlockNum + 255) / 256;
			int s1 = (s + 255) / 256;

			dim3 dimGrid = dim3(tmpBlockNum, 1, 1);
			dim3 dimBlocks = dim3(threadNum, 1, 1);

			blocks = dim3((ArrayPixels + 15) / 16, (overlapLines + 15) / 16);
			dim3 threads42 = dim3(64, 4);
			dim3 blocks42 = dim3((ArrayPixels + threads42.x - 1) / threads42.x, (OverlapMaxLines / 2 + threads42.y - 1) / threads42.y);

			cudaMemset(d_SlaveBack, 0, Maxize*sizeof(cuComplex));

			cudaMemcpy(d_MasterFor , MasterFor + DeviceOffset,
				overlapN*sizeof(short2), cudaMemcpyHostToDevice);
			cudaMemcpy(d_MasterBack , MasterBack + DeviceOffset,
				overlapN*sizeof(short2), cudaMemcpyHostToDevice);
			cudaMemcpy(d_SlaveFor , SlaveFor + DeviceOffset,
				overlapN*sizeof(cuComplex), cudaMemcpyHostToDevice);
			cudaMemcpy(d_SlaveBack , SlaveBack + DeviceOffset,
				overlapN*sizeof(cuComplex), cudaMemcpyHostToDevice);



			cudaMemset(d_MasterCoh , 0, d_pitchCohfloat2*cohh);
			cudaMemset(d_SlaveCoh , 0, d_pitchCohfloat2*cohh);

			//Estimating the coherence of forward overlap band 
			CohMatrixSet << <blocks, threads >> >(d_MasterFor , d_SlaveFor ,
				d_MasterCoh , d_SlaveCoh, overlapLines, ArrayPixels, dh, dw, d_pitchCohfloat2);

			CohFirstLineSet_ESD << <(ArrayPixels + 255) / 256, 256 >> >(d_sum1 ,
				d_power1 , d_MasterCoh , d_SlaveCoh , CohWinAz, CohWinRg,
				ArrayPixels, d_pitchCohfloat2);


			CohDiffer_vertical_ESD << <PartBlocksDiff, threadsV >> >(d_sum1 ,
				d_power1 , d_MasterCoh , d_SlaveCoh , CohWinAz, CohWinRg,
				overlapLines - 1, ArrayPixels,
				d_pitchCohfloat2, d_pitchCohTmp);




			shfl_vertical_ESD << <GridPrefix, blockSz >> >(d_sum1, numX, OverlapMaxLines, PartSteps);
			shfl_vertical_ESD << <GridPrefix, blockSz >> >(d_power1, numX, OverlapMaxLines, PartSteps);

			CpCoh2_ESD << <blocks, threads >> >(d_coh1 ,
				d_sum1 , d_power1 , overlapLines, ArrayPixels, d_pitchCohTmp, d_pitchCohfloat);


			//Estimating the coherence of backward overlap band 
			cudaMemset(d_MasterCoh, 0, d_pitchCohfloat2*cohh);
			cudaMemset(d_SlaveCoh, 0, d_pitchCohfloat2*cohh);

			CohMatrixSet << <blocks, threads >> >(d_MasterBack , d_SlaveBack ,
				d_MasterCoh , d_SlaveCoh ,overlapLines, ArrayPixels, dh, dw, d_pitchCohfloat2);


			CohFirstLineSet_ESD << <(ArrayPixels + 255) / 256, 256 >> >(d_sum1 ,
				d_power1 , d_MasterCoh , d_SlaveCoh , CohWinAz, CohWinRg,
				ArrayPixels, d_pitchCohfloat2);

			CohDiffer_vertical_ESD << <PartBlocksDiff, threadsV >> >(d_sum1 ,
				d_power1 , d_MasterCoh , d_SlaveCoh , CohWinAz, CohWinRg,
				overlapLines - 1, ArrayPixels, d_pitchCohfloat2, d_pitchCohTmp);



			shfl_vertical_ESD << <GridPrefix, blockSz >> >(d_sum1 , numX, OverlapMaxLines, PartSteps);
			shfl_vertical_ESD << <GridPrefix, blockSz >> >(d_power1 , numX, OverlapMaxLines, PartSteps);

			CpCoh2_ESD << <blocks, threads >> >(d_coh2 ,
				d_sum1 , d_power1 , overlapLines, ArrayPixels, d_pitchCohTmp, d_pitchCohfloat);



			//Average the coherence of forward and backward bands
			Divide2 << <blocks, threads >> > (d_coh1 , d_coh2 , d_pitchCohfloat, overlapLines, ArrayPixels);
			//cudaStreamSynchronize(stream[i]);
			//Estimte the double differential interferogram
			
			DoubleDifference4OneD << <blocks, threads >> >
				(d_MasterFor , d_MasterBack ,
				d_SlaveFor , d_SlaveBack ,
				overlapLines, ArrayPixels);
			
			//Masking out incoherent pixels
			MaskOut << <blocks, threads >> >(d_SlaveBack , d_coh1 , CohThresHold, overlapLines, ArrayPixels,
				d_pitchCohfloat);

			//Using parallel recuding operation to accumulate ESD phases
			cudaMemset(d_blockSum,0, blockNum*sizeof(cuComplex));
			cudaMemset(d_blockSum,0, blockNum1*sizeof(cuComplex));
			cudaMemset(d_blockSum,0, sizeof(cuComplex));

			reduce0 << <tmpBlockNum, 256 >> >(d_SlaveBack , d_blockSum , overlapLines*ArrayPixels);


			reduce0 << <s, 256 >> >(d_blockSum , d_blockSum1 , tmpBlockNum);


			reduce0 << <s1, 256 >> >(d_blockSum1 , d_blockSum2 , (tmpBlockNum + 255) / 256);

			cudaMemcpy(SumOverlap+i, d_blockSum2, 1 * sizeof(cuComplex), cudaMemcpyDeviceToHost);
			//cudaStreamSynchronize(stream[i]);

		}

	}


	

	cudaEventRecord(g_stop, 0);
	cudaEventSynchronize(g_stop);
	cudaEventElapsedTime(&time_cost, g_start, g_stop);
	cout << "ESDCorretion duration:" << time_cost << "ms" << endl;
	cudaEventDestroy(g_start);
	cudaEventDestroy(g_stop);
	//ofstream ShiftOut("E:\\ESDTimer.txt", ios::app);
	//ShiftOut << time_cost << endl;
	//ShiftOut.close();

	for (int i = 0; i < numOverlap; i++)
	{

		ShiftArray[i] = atan2((double)SumOverlap[i].y, (double)SumOverlap[i].x) / (2 * CUDART_PI*spectralSeparation*azimuthTimeInterval);
	}






	cudaHostUnregister(MasterFor);
	cudaHostUnregister(MasterBack);
	cudaHostUnregister(SlaveFor);
	cudaHostUnregister(SlaveBack);

	cudaFree(d_MasterFor);
	cudaFree(d_MasterBack);
	cudaFree(d_SlaveFor);
	cudaFree(d_SlaveBack);



	cudaFree(d_sum1);
	cudaFree(d_power1);
	cudaFree(d_MasterCoh);
	cudaFree(d_SlaveCoh);
	cudaFree(d_coh1);
	cudaFree(d_coh2);
	cudaFree(d_blockSum);
	
	cudaDeviceSynchronize();
	int cuda_last_error_flag = cudaGetLastError();

	if (cuda_last_error_flag != 0)
	{
		printf("Returns an error code (%d) from GPU execution!\n Please debug it or report it to developer!\n", cuda_last_error_flag);
		exit(0);
	}
}


void EstESDShifts(
	complex<short>*MasterFor, //Input
	complex<short>*MasterBack, //Input
	complex<float>*SlaveFor,// Input
	complex<float>*SlaveBack, //Input
	int* OverlapSizeArray, //Input
	double * ShiftArray,//Output
	int numOverlap,
	int ArrayLines,
	int ArrayPixels,
	int CohWinAz,
	int CohWinRg,
	int OverlapMaxLines,
	float CohThresHold,
	double spectralSeparation,
	double azimuthTimeInterval
	)
{
	bool OverlaporNot = false;

	if (!OverlaporNot)
	{
		 EstESDShifts_NonOverlap(MasterFor, MasterBack, SlaveFor,SlaveBack, 
			 OverlapSizeArray, ShiftArray,numOverlap, ArrayLines,ArrayPixels,
			 CohWinAz,CohWinRg,OverlapMaxLines,CohThresHold,spectralSeparation,
			 azimuthTimeInterval);
	}

	if (OverlaporNot)
	{
		EstESDShifts_Overlap(MasterFor, MasterBack, SlaveFor, SlaveBack,
			OverlapSizeArray, ShiftArray, numOverlap, ArrayLines, ArrayPixels,
			CohWinAz, CohWinRg, OverlapMaxLines, CohThresHold, spectralSeparation,
			azimuthTimeInterval);
	}









}
