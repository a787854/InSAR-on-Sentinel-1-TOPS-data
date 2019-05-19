#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#define __CUDA_INTERNAL_COMPILATION__
#include <math_functions.h>
#include <math_constants.h>
#include <device_functions.h>
#include <cuComplex.h>
//#include <helper_cuda.h>
#include "cuConstants.cuh"
#undef  __CUDA_INTERNAL_COMPILATION__
#include <math.h>
#include <iostream>
#include<omp.h>
#include <complex>


using namespace std;


#define THREADBLOCK_SIZE 256




int  MultipleOf8(int inputLines)
{
	int remains = inputLines % 8;
	int multiple = inputLines / 8;
	if (remains == 0)
		return inputLines;
	else{

		return (multiple + 1) * 8;


	}
}


int  MultipleOf32(int inputLines)
{
	int remains = inputLines % 32;
	int multiple = inputLines / 32;
	if (remains == 0)
		return inputLines;
	else{

		return (multiple + 1) * 32;


	}
}


static int iDivUp(int dividend, int divisor)
{
	return ((dividend % divisor) == 0) ? (dividend / divisor) : (dividend / divisor + 1);
}




__global__ void CohMatSet(cuComplex *DstMaster, cuComplex *DstSlave, cuComplex *SrcMaster,
	cuComplex *SrcSlave, int SrcLines, int SrcPixels,
	size_t DstPitch, size_t SrcPitch, int dx, int dy)
{

	const int row = blockIdx.y*blockDim.y + threadIdx.y;
	const int col = blockIdx.x*blockDim.x + threadIdx.x;



	if (row < SrcLines&&col < SrcPixels)
	{


		cuComplex* rowDstMaster = (cuComplex *)((char*)DstMaster + (row + dy)*DstPitch);
		cuComplex* rowDstSlave = (cuComplex *)((char*)DstSlave + (row + dy)*DstPitch);
		cuComplex* rowSrcMaster = (cuComplex *)((char*)SrcMaster + row*SrcPitch);
		cuComplex* rowSrcSlave = (cuComplex *)((char*)SrcSlave + row*SrcPitch);

		cuComplex MasterTemp = rowSrcMaster[col];
		cuComplex SlaveTemp = rowSrcSlave[col];
		rowDstMaster[col + dx] = cuCmulf(MasterTemp, cuConjf(SlaveTemp));



		rowDstSlave[col + dx] =
			make_cuComplex(SlaveTemp.x*SlaveTemp.x + SlaveTemp.y*SlaveTemp.y,
			MasterTemp.x*MasterTemp.x + MasterTemp.y*MasterTemp.y);




	}



}

__global__ void CohMatSet
(cuComplex *DstMaster, cuComplex *DstSlave, short2 *SrcMaster,
cuComplex *SrcSlave, int SrcLines, int SrcPixels,
size_t DstPitch, size_t SrcSlavePitch,size_t SrcMasterPitch, int dx, int dy)
{

	const int row = blockIdx.y*blockDim.y + threadIdx.y;
	const int col = blockIdx.x*blockDim.x + threadIdx.x;



	if (row < SrcLines&&col < SrcPixels)
	{
	

		cuComplex* rowDstMaster = (cuComplex *)((char*)DstMaster + (row + dy)*DstPitch);
		cuComplex* rowDstSlave = (cuComplex *)((char*)DstSlave + (row + dy)*DstPitch);
		short2* rowSrcMaster = (short2 *)((char*)SrcMaster + row*SrcMasterPitch);
		cuComplex* rowSrcSlave = (cuComplex *)((char*)SrcSlave + row*SrcSlavePitch);

		cuComplex MasterTemp = make_cuComplex(rowSrcMaster[col].x, rowSrcMaster[col].y);
		cuComplex SlaveTemp = rowSrcSlave[col];



		rowDstMaster[col + dx] = cuCmulf(MasterTemp, cuConjf(SlaveTemp));;

    
		rowDstSlave[col + dx] = make_cuComplex(SlaveTemp.x*SlaveTemp.x + SlaveTemp.y*SlaveTemp.y,
			MasterTemp.x*MasterTemp.x + MasterTemp.y*MasterTemp.y);

		

	
	
		
	}

}












__global__ void CuCplCoherence_vertical(cuComplex *Sum1Array, cuComplex* Power1Array,
	int inputLines, int inputPixels, size_t DstPitch, int LinesOffset)
{


	const int row = blockIdx.y*blockDim.y + threadIdx.y;
	const int col = blockIdx.x*blockDim.x + threadIdx.x;




	if (row < inputLines&&col < inputPixels)
	{
		cuComplex sum1 = make_cuComplex(0.0f, 0.0f);
		cuComplex power1 = make_cuComplex(0.0f, 0.0f);
		int Roffset = LinesOffset;


		for (int i = 0; i < c_winX; i++)
		{
			sum1 = cuCaddf(sum1, cuCsubf(tex2D(tex_input, col + i, row + Roffset + c_winY),
				tex2D(tex_input, col + i, row + Roffset)));
		}
		Sum1Array[(row + 1)*(DstPitch / 8) + col] = sum1;
		for (int i = 0; i < c_winX; i++)
		{
			power1 = cuCaddf(power1, cuCsubf(tex2D(tex_norms, col + i, row + Roffset + c_winY),
				tex2D(tex_norms, col + i, row + Roffset)));
		}
		Power1Array[(row + 1)*(DstPitch / 8) + col] = power1;

	

	}

}





__global__ void CuCplCoherence_basic(float *Result, cuComplex* JudgeArray, int inputLines, int inputPixels, size_t DstPitch, size_t judgePitch)
{


	const int row = blockIdx.y*blockDim.y + threadIdx.y;
	const int col = blockIdx.x*blockDim.x + threadIdx.x;




	if (row < inputLines&&col < inputPixels)
	{
		int winY = c_winY;
		int winX = c_winX;

		cuComplex sum1 = make_cuComplex(0.0f, 0.0f);
		cuComplex power1 = make_cuComplex(0.0f, 0.0f);
		for (int j = 0; j < winY; j++)
		{
			for (int i = col; i < col + winX; i++)
			{

				sum1 = cuCaddf(sum1, tex2D(tex_input, i, row + j));
				power1 = cuCaddf(power1, tex2D(tex_norms, i, row + j));


			}
		}

		float p = power1.x*power1.y;

		float* rowResult = (float *)((char*)Result + row*DstPitch);



		rowResult[col] = (p > 0.0f) ? cuCabsf((cuCdivf(sum1, make_cuComplex(sqrt(p), 0.0f)))) : 0.0f;

	}



}


__global__ void CohFirstLineSet(cuComplex *d_Sum1, cuComplex *d_Power1, int inputPixels, size_t InPitch, int Offset)
{


	const int col = blockIdx.x*blockDim.x + threadIdx.x;


	if (col < inputPixels)
	{

		int winY = c_winY;
		int winX = c_winX;

		cuComplex sum1 = make_cuComplex(0.0f, 0.0f);
		cuComplex power1 = make_cuComplex(0.0f, 0.0f);

		for (int j = Offset; j <Offset+ winY; j++)
		{


			for (int i = col; i < col + winX; i++)
			{
				sum1 = cuCaddf(sum1, tex2D(tex_input, i, j));
				power1 = cuCaddf(power1, tex2D(tex_norms, i, j));

			}
		}


		d_Sum1[col] = sum1;
		d_Power1[col] = power1;

		

		


	}



}






__global__ void shfl_vertical(cuComplex *Src, int width, int height, int Nstep)
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
			float n_x = __shfl_up(partial_sum_x, i, 32);
			float n_y = __shfl_up(partial_sum_y, i, 32);// removed for debug
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





__global__ void CpCoh2(float *Res, cuComplex *Sum1Array, cuComplex *Power1Array, cuComplex*JudgeArray, size_t Pitch,
	size_t Pitch2, size_t Pitch3, int inputLines)
{
	int row = blockIdx.y*blockDim.y + threadIdx.y;
	int col = blockIdx.x*blockDim.x + threadIdx.x;


	if (row < inputLines&&col < c_Dw)
	{

		size_t ComplexPitch = Pitch;
		size_t FloatPitch = Pitch2;
		cuComplex* rowSum = (cuComplex *)((char*)Sum1Array + row*Pitch3);
		cuComplex* rowPower = (cuComplex *)((char*)Power1Array + row*Pitch3);
		float* rowRes = (float *)((char*)Res + row*FloatPitch);
		cuComplex* rowJudgeArray = (cuComplex *)((char*)JudgeArray + row*ComplexPitch);

		cuComplex tmpSum = rowSum[col];
		cuComplex tmpPower = rowPower[col];


		//if (row == 1000 && col == 5000)
		//{
		////printf("Judge Array:(%lf,%lf)\n", rowJudgeArray[col].x, rowJudgeArray[col].y);
		//	printf("TmpSum:(%lf,%lf)\n", tmpSum.x, tmpSum.y);
		//	printf("TmpPower:(%lf,%lf)\n", tmpPower.x, tmpPower.y);
		//}

		float pp = tmpPower.x*tmpPower.y;

		



		if (cuCabsf(rowJudgeArray[col])<0.001)
		{
			rowRes[col] = 0.0f;
		}
		else
		{
		

			cuComplex TmpRes = (pp > 0.0) ? (cuCdivf(tmpSum, make_cuComplex(sqrt(pp), 0.0f))) : make_cuComplex(0.0f, 0.0f);


			rowRes[col] = cuCabsf(TmpRes);
		}


	}


}

__global__ void Cohbasic(float *Result, int inputLines, int inputPixels, size_t DstPitch, cuComplex*JudgeArray,size_t complexpitch)
{


	const int row = blockIdx.y*blockDim.y + threadIdx.y;
	const int col = blockIdx.x*blockDim.x + threadIdx.x;


	if (row < inputLines&&col < inputPixels)
	{
		int winY = c_winY;
		int winX = c_winX;

		cuComplex sum1 = make_cuComplex(0.0f, 0.0f);
		cuComplex power1 = make_cuComplex(0.0f, 0.0f);
		for (int j = 0; j < winY; j++)
		{
			for (int i = col; i < col + winX; i++)
			{



				sum1 = cuCaddf(sum1, tex2D(tex_input, i, row + j));
				power1 = cuCaddf(power1, tex2D(tex_norms, i, row + j));
				

			}
		}

		float p = power1.x*power1.y;

		float* rowResult = (float *)((char*)Result + row*DstPitch);
		cuComplex* rowJudgeArray = (cuComplex *)((char*)JudgeArray + row*complexpitch);
		if (cuCabsf(rowJudgeArray[col])<0.001)
		{
			rowResult[col] = 0.0f;
		}
		else
		{

			rowResult[col] = (p > 0.0f) ? cuCabsf((cuCdivf(sum1, make_cuComplex(sqrt(p), 0.0f)))) : 0.0f;
		}
	}
}





extern "C" void CoherenceEst_ESD
(int Dx0,
int Dy0,
int Dw,
int Dh,
int dx,
int dy,
int cohWinRg,
int cohWinAz,
complex<short>* masterArray,
cuComplex* d_slaveArray,
float* cohdata
)
{
	int dataBytes = Dw*Dh*sizeof(complex<float>);

	cudaHostRegister(cohdata, dataBytes / 2.0, cudaHostRegisterDefault);
	cudaHostRegister(masterArray, dataBytes / 2.0, cudaHostRegisterDefault);


	size_t Pitch3;

	cuComplex*d_masterArray = NULL;
	cuComplex* d_masterArray2 = NULL;
	cuComplex *d_slaveArray2 = NULL;
	cuComplex *d_Sum1, *d_Power1;
	short2* d_TempMasterArray;

	float *d_CohMat_Amp = NULL;
	size_t FloatPitch;
	size_t Pitch;
	size_t FloatPitch3;


	
	cudaMallocPitch((void**)&d_TempMasterArray, &FloatPitch, Dw*sizeof(short2), Dh);// For copy in complex<short> data

	

	//for vertical prefix sum arrays
	int Lines = MultipleOf32(Dh);
	int cohw = Dw + cohWinRg - 1;
	int cohh = Lines + cohWinAz - 1;
	cudaMallocPitch((void**)&d_CohMat_Amp, &FloatPitch, Dw*sizeof(float), Lines);
	
	

	

	cudaMallocPitch((void**)&d_Sum1, &Pitch, Dw*sizeof(cuComplex), Lines);
	cudaMallocPitch((void**)&d_Power1, &Pitch, Dw*sizeof(cuComplex), Lines);
	

	int Nstep = Lines / 8;

	dim3 blockSz(32, 8);
	int PitchPixels = Pitch / sizeof(cuComplex);
	//For horizontal
	//int Nstep = PitchPixels / 8;
	if (PitchPixels% blockSz.x != 0)
	{
		cout << "Warning:PitchSize is not multiple of 32!" << endl;
		cout << "Do not worry! we can handle it!" << endl;
		cudaFree(d_Sum1);
		cudaFree(d_Power1);
		Pitch = MultipleOf32(PitchPixels)*sizeof(cuComplex);
		cudaMalloc((void**)&d_Sum1, Pitch*Lines);
		cudaMalloc((void**)&d_Power1, Pitch*Lines);
		PitchPixels = Pitch / sizeof(cuComplex);

	}
	int numX = PitchPixels;

	dim3 Grid(numX / blockSz.x, 1);


	size_t Pitch2;
	cudaMallocPitch((void**)&d_masterArray2, &Pitch2, cohw*sizeof(cuComplex), cohh);
	cudaMallocPitch((void**)&d_slaveArray2, &Pitch2, cohw*sizeof(cuComplex), cohh);
;


	cudaMemcpyToSymbol(c_Dy0, &Dy0, sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_Dw, &Dw, sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_Dh, &Dh, sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_winY, &cohWinAz, sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_winX, &cohWinRg, sizeof(int), 0, cudaMemcpyHostToDevice);

	cudaChannelFormatDesc channelDesc_Complex = cudaCreateChannelDesc(32, 32, 0, 0, cudaChannelFormatKindFloat);


	dim3 threads(16, 16);
	dim3 Outblocks = dim3((Dw + threads.x-1) /threads.x , (Dh + threads.y-1) / threads.y);
	dim3 Cohblocks = dim3((cohw + threads.x-1) / threads.x, (cohh + threads.y-1) / threads.y);




	cudaStream_t stream[4];
	for (int i = 0; i<4; i++)
		cudaStreamCreate(&stream[i]);




	cudaEvent_t g_start, g_stop;
	float time_cost1, time_cost2;
	cudaEventCreate(&g_start);
	cudaEventCreate(&g_stop);
	cudaEventRecord(g_start, 0);

	cudaMemcpy2DAsync(d_TempMasterArray, FloatPitch, masterArray, Dw*sizeof(short2),
		Dw*sizeof(short2), Dh, cudaMemcpyHostToDevice, stream[0]);

	cudaDeviceSynchronize();

	


	//Windows does not need to memset the array with zero values
	// (cudaMemsetAsync(d_masterArray2, 0, Pitch2*cohh, stream1));
	// (cudaMemsetAsync(d_slaveArray2, 0, Pitch2*cohh, stream1));



	CohMatSet << <Outblocks, threads, 0, stream[0] >> >
		(d_masterArray2, d_slaveArray2, d_TempMasterArray, d_slaveArray, Dh,
		Dw, Pitch2, Pitch,FloatPitch, dx, dy);


	//BindTexture
	cudaBindTexture2D(0, tex_input, d_masterArray2, channelDesc_Complex, cohw, cohh, Pitch2);
	cudaBindTexture2D(0, tex_norms, d_slaveArray2, channelDesc_Complex, cohw, cohh, Pitch2);



	dim3 threadsV(16, 16);
	dim3 blocks = dim3((Dw + threadsV.x - 1) / threadsV.x, (Dh - 1 + threadsV.y - 1) / threadsV.y);

	

	/*Hide the time for coping back the coherence map*/
	int PartLines = Lines / 4;
	int PartOffset = PartLines*numX;
	int HostOffset = PartLines*Dw;

	//
	bool BasicOrNot = false;
	bool OverlapOrNot = false;

	
	if (OverlapOrNot)
	{
		/*Overlap transfer*/

		int PartNstep = PartLines / 8;
		dim3 PartBlocksDiff((Dw + threads.x - 1) / threads.x, (PartLines - 1 + threads.y - 1) / threads.y);
		dim3 PartBlocks((Dw + threads.x - 1) / threads.x, (PartLines + threads.y - 1) / threads.y);
		
		int outputLines[4];
		outputLines[0] = PartLines;
		outputLines[1] = PartLines;
		outputLines[2] = PartLines;
		outputLines[3] = Dh - 3 * PartLines;



		for (int N = 0; N < 4; N++)
		{
			CohFirstLineSet << <(Dw + 255) / 256, 256, 0, stream[N] >> >(d_Sum1 + N* PartOffset,
				d_Power1 + N*PartOffset, Dw, Pitch2, N*PartLines);

			CuCplCoherence_vertical << <PartBlocksDiff, threads, 0, stream[N] >> >
				(d_Sum1 + N* PartOffset, d_Power1 + N*PartOffset, PartLines - 1, Dw, Pitch, N*PartLines);


			shfl_vertical << <Grid, blockSz, 0, stream[N] >> >(d_Sum1 + N*PartOffset, PitchPixels, Lines, PartNstep);
			shfl_vertical << <Grid, blockSz, 0, stream[N] >> >(d_Power1 + N*PartOffset, PitchPixels, Lines, PartNstep);


			CpCoh2 << <PartBlocks, threads, 0, stream[N] >> >
				(d_CohMat_Amp + N*PartOffset, d_Sum1 + N*PartOffset, d_Power1 + N*PartOffset,
				d_slaveArray + N*PartOffset,
				Pitch, FloatPitch, Pitch, PartLines);

			cudaStreamSynchronize(stream[N]);

			cudaMemcpy2DAsync(cohdata + N*HostOffset, Dw*sizeof(float),
				d_CohMat_Amp + N*PartOffset, FloatPitch, Dw*sizeof(float), outputLines[N], cudaMemcpyDeviceToHost, stream[N]);


		}
	}

	
	if (!OverlapOrNot)
	{
		
		int threads1D = 256;
		dim3 Blocks_vertical = dim3((Dw + threads.x - 1) / threads.x, (Dh - 1 + threads.y - 1) / threads.y);
		CohFirstLineSet << <(Dw + threads1D - 1) / threads1D, threads1D,0,stream[2] >> >(d_Sum1, d_Power1, Dw, Pitch2, 0);
		CuCplCoherence_vertical << <Blocks_vertical, threads, 0, stream[2] >> >(d_Sum1, d_Power1, Dh - 1, Dw, Pitch, 0);
		shfl_vertical << <Grid, blockSz, 0, stream[2] >> >(d_Sum1, PitchPixels, Lines, Nstep);
		shfl_vertical << <Grid, blockSz, 0, stream[2] >> >(d_Power1, PitchPixels, Lines, Nstep);

		CpCoh2 << <Outblocks, threads,0, stream[2] >> >
			(d_CohMat_Amp , d_Sum1 , d_Power1 ,d_slaveArray ,
			Pitch, FloatPitch, Pitch, Dh);
		cudaMemcpy2D(cohdata , Dw*sizeof(float),
			d_CohMat_Amp, FloatPitch, Dw*sizeof(float), Dh, cudaMemcpyDeviceToHost);
	}

	
	if (BasicOrNot)
	{
		Cohbasic << <Outblocks, threads, 0, stream[3] >> >(d_CohMat_Amp, Dh, Dw, FloatPitch,d_slaveArray,Pitch);

		cudaMemcpy2D(cohdata, Dw*sizeof(float),
			d_CohMat_Amp, FloatPitch, Dw*sizeof(float), Dh, cudaMemcpyDeviceToHost);
	}

	

	for (int i = 0; i<4; i++)
		cudaStreamDestroy(stream[i]);
	cudaEventRecord(g_stop, 0);

	

	cudaEventSynchronize(g_stop);
	cudaEventElapsedTime(&time_cost2, g_start, g_stop);
	cout << "coherence kernel duration:" << time_cost2<<"ms"<< endl;
	cudaEventDestroy(g_start);
	cudaEventDestroy(g_stop);



	cudaHostUnregister(cohdata);
	cudaHostUnregister(masterArray);







	cudaFree(d_slaveArray);


	cudaUnbindTexture(tex_input);
	cudaUnbindTexture(tex_norms);


	cudaFree(d_CohMat_Amp);

	cudaFree(d_masterArray2);
	cudaFree(d_slaveArray2);
	cudaFree(d_Power1);
	cudaFree(d_Sum1);
	cudaFree(d_TempMasterArray);

	cudaDeviceReset();



}



extern "C" void CoherenceEst_Basic
(int Dx0,
int Dy0,
int Dw,
int Dh,
int dx,
int dy,
int cohWinRg,
int cohWinAz,
complex<short>* masterArray,
cuComplex* d_slaveArray,
float* cohdata
)
{
	int dataBytes = Dw*Dh*sizeof(complex<float>);

	cudaHostRegister(cohdata, dataBytes / 2.0, cudaHostRegisterDefault);
	cudaHostRegister(masterArray, dataBytes / 2.0, cudaHostRegisterDefault);


	size_t Pitch3;

	cuComplex*d_masterArray = NULL;
	cuComplex* d_masterArray2 = NULL;
	cuComplex *d_slaveArray2 = NULL;
	short2* d_TempMasterArray;

	float *d_CohMat_Amp = NULL;
	size_t FloatPitch;
	size_t Pitch;
	size_t FloatPitch3;


	size_t TotalBytes = 0;
	cudaMallocPitch((void**)&d_TempMasterArray, &FloatPitch, Dw*sizeof(short2), Dh);// For copy in complex<short> data
	


	//for vertical prefix sum arrays
	int Lines = MultipleOf32(Dh);
	int cohw = Dw + cohWinRg - 1;
	int cohh = Lines + cohWinAz - 1;
	

	cuComplex* PitchArray;
	cudaMallocPitch((void**)&PitchArray, &Pitch, Dw*sizeof(cuComplex), 1);
	cudaFree(PitchArray);

	size_t Pitch2;
	cudaMallocPitch((void**)&d_masterArray2, &Pitch2, cohw*sizeof(cuComplex), cohh);
	cudaMallocPitch((void**)&d_slaveArray2, &Pitch2, cohw*sizeof(cuComplex), cohh);



	cudaMemcpyToSymbol(c_Dy0, &Dy0, sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_Dw, &Dw, sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_Dh, &Dh, sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_winY, &cohWinAz, sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_winX, &cohWinRg, sizeof(int), 0, cudaMemcpyHostToDevice);

	cudaChannelFormatDesc channelDesc_Complex = cudaCreateChannelDesc(32, 32, 0, 0, cudaChannelFormatKindFloat);


	dim3 threads(16, 16);
	dim3 Outblocks = dim3((Dw + threads.x - 1) / threads.x, (Dh + threads.y - 1) / threads.y);
	dim3 Cohblocks = dim3((cohw + threads.x - 1) / threads.x, (cohh + threads.y - 1) / threads.y);









	
	 (cudaMemcpy2D(d_TempMasterArray, FloatPitch, masterArray, Dw*sizeof(short2),
		Dw*sizeof(short2), Dh, cudaMemcpyHostToDevice));





	//Windows does not need to memset the array with zero values
	// (cudaMemsetAsync(d_masterArray2, 0, Pitch2*cohh, stream1));
	// (cudaMemsetAsync(d_slaveArray2, 0, Pitch2*cohh, stream1));



	CohMatSet << <Outblocks, threads >> >
		(d_masterArray2, d_slaveArray2, d_TempMasterArray, d_slaveArray, Dh,
		Dw, Pitch2, Pitch, FloatPitch, dx, dy);

	
	
	

	cudaMallocPitch((void**)&d_CohMat_Amp, &FloatPitch, Dw*sizeof(float), Lines);


	//BindTexture
	cudaBindTexture2D(0, tex_input, d_masterArray2, channelDesc_Complex, cohw, cohh, Pitch2);
	cudaBindTexture2D(0, tex_norms, d_slaveArray2, channelDesc_Complex, cohw, cohh, Pitch2);



	dim3 threadsV(16, 16);
	dim3 blocks = dim3((Dw + threadsV.x - 1) / threadsV.x, (Dh - 1 + threadsV.y - 1) / threadsV.y);



	cudaEvent_t g_start, g_stop;
	float time_cost1;
	cudaEventCreate(&g_start);
	cudaEventCreate(&g_stop);
	cudaEventRecord(g_start, 0);



	//Basic Coherence estimation


	Cohbasic << <Outblocks, threads >> >(d_CohMat_Amp, Dh, Dw, FloatPitch, d_slaveArray, Pitch);
	
	 (cudaMemcpy2D(cohdata, Dw*sizeof(float),
			d_CohMat_Amp, FloatPitch, Dw*sizeof(float), Dh, cudaMemcpyDeviceToHost));





	
	cudaEventRecord(g_stop, 0);
	cudaEventSynchronize(g_stop);
	 (cudaEventElapsedTime(&time_cost1, g_start, g_stop));
	cout << "Basic coherence kernel duration:" << time_cost1 << "ms" << endl;
	cudaEventDestroy(g_start);
	cudaEventDestroy(g_stop);
	 (cudaDeviceSynchronize());

	cudaHostUnregister(cohdata);
	cudaHostUnregister(masterArray);









	cudaUnbindTexture(tex_input);
	cudaUnbindTexture(tex_norms);


	 (cudaFree(d_CohMat_Amp));

	 (cudaFree(d_masterArray2));
	 (cudaFree(d_slaveArray2));
	
	 (cudaFree(d_slaveArray));
	 (cudaFree(d_TempMasterArray));

	cudaDeviceReset();



}


extern "C" void CoherenceEst
(int Dx0,
int Dy0,
int Dw,
int Dh,
int dx,
int dy,
int cohWinRg,
int cohWinAz,
complex<short>* masterArray,
cuComplex* d_slaveArray,
float* cohdata
)
{
	bool GPUMemSmall = true; 

	//if the total memory is less than 2Gb, go to coherence basic estimation version.
	if (GPUMemSmall)
	{
		CoherenceEst_Basic(Dx0, Dy0, Dw, Dh, dx, dy, cohWinRg, cohWinAz, masterArray,
			d_slaveArray, cohdata);
	}
	else
	{
		CoherenceEst_ESD(Dx0, Dy0, Dw, Dh, dx, dy, cohWinRg, cohWinAz, masterArray,
			d_slaveArray, cohdata);

	}
	



}






