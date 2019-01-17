#ifndef CUCONSTANTS_CUH
#define CUCONSTANTS_CUH
#include "cuComplex.h"

#define b_ell_a  6378137.0
#define b_ell_b  6356752.3142451794975639665996337
#define b_SOL_cuda  299792458.0

//GT
__constant__  int   b_nrows;
__constant__  int   b_ncols;
__constant__ double SOL_cuda = 299792458.0;
__constant__ double b_ell_e2;
__constant__ size_t b_CorrPitch;
__constant__ double b_lat0;
__constant__ double b_lon0;
__constant__ double b_DEMdeltalat;
__constant__ double b_DEMdeltalon;
__constant__ double b_s_min4picdivlam;
__constant__ double d_AzimuthShift;

//DR
__constant__  int npoints;
__constant__ int c_sX0;
__constant__ int c_sY0;
__constant__ int c_sXmax;
__constant__ int c_sYmax;
__constant__ int c_mLines;
__constant__ int c_mPixels;
__constant__ int c_mY0;
__constant__ int c_MasterBox[4];
__constant__ double c_CpmAz[6];
__constant__ double c_CpmRg[6];
__constant__ int c_Npoints2m1;
__constant__ int c_Interval=2047;

texture<cuComplex, 2, cudaReadModeElementType> tex_slave;
texture<float, 2, cudaReadModeElementType> tex_kernelRg;
texture<float, 2, cudaReadModeElementType> tex_kernelAz;
texture<float, 2, cudaReadModeElementType>tex_PhaseArray;

//DC
__constant__ double c_coef[48];
__constant__ int c_Dy0;
__constant__ int c_Dw;
__constant__ int c_Dh;
__constant__ int c_winY;
__constant__ int c_winX;
texture<cuComplex, 2, cudaReadModeElementType> tex_input;
texture<cuComplex, 2, cudaReadModeElementType> tex_norms;

inline __device__ double normalizeWarp(double data, double min, double max)
{
	return (data - (0.5*(min + max))) / (0.25*(max - min));
}

__device__ __host__ inline double my_polyval(double x, double y, double *coeff)
{
	double sum = coeff[0];

	sum += (coeff[1] * x
		+ coeff[2] * y
		+ coeff[3] * x*x
		+ coeff[4] * x*y
		+ coeff[5] * y*y);

	return sum;
}

#endif