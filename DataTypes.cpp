#include"DataTypes.h"



/*****************************************************************/
/*Class:RefDem*/
/*****************************************************************/
void RefDem::Init(const char* DEMPath)
{
	GDALAllRegister();
	pData_dem = GDALOpen(DEMPath, GA_ReadOnly);
	if (pData_dem == NULL)
	{
		cout << "can not open DEM: " << DEMPath << " !\n";
		system("pause");
		exit(0);
	}
	//获取DEM地理坐标信息
	double adfGeoTransform[6];
	if (GDALGetGeoTransform(pData_dem, adfGeoTransform) == CE_None)
	{
		lon_min = adfGeoTransform[0];//pixel_left_up
		deltaLon = fabs(adfGeoTransform[1]);//pixel_interval
		lat_max = adfGeoTransform[3];//line_left_up
		deltaLat = fabs(adfGeoTransform[5]);//line_interval
	}
	Lines = GDALGetRasterYSize(pData_dem);
	Pixels = GDALGetRasterXSize(pData_dem);

	lat_min = lat_max - (Lines - 1)*deltaLat;
	lon_max = lon_min + (Pixels - 1)*deltaLon;
}

RefDem:: ~RefDem()
{
	GDALClose(pData_dem);
}
void RefDem::getData(int x0, int y0, int ww, int hh, short* demBuffer)
{
	GDALRasterBandH hBand_dem;
	hBand_dem = GDALGetRasterBand(pData_dem, 1);
	GDALRasterIO(hBand_dem, GF_Read, x0, y0, ww, hh, demBuffer, ww, hh, GDT_Int16, 0, 0);
}

void RefDem::getData(double lat_min, double lat_max, double lon_min, double lon_max, double extralat,
	double extralon, short*&demBuffer, int& Lines, int& Pixels)
{
	lat_min -= extralat;
	lat_max += extralat;
	lon_min -= extralon;
	lon_max += extralon;

	double UpperLeft[2];//[LatIndex, LonIndex]
	double LowerRight[2];

	//
	getIndex(lat_max, lon_min, UpperLeft);
	getIndex(lat_min, lon_max, LowerRight);
	UpperLeft[0] = floor(UpperLeft[0]);
	UpperLeft[1] = floor(UpperLeft[1]);
	LowerRight[0] = ceil(LowerRight[0]);
	LowerRight[1] = ceil(LowerRight[1]);

	Lines= LowerRight[0] - UpperLeft[0] + 1;
	Pixels = LowerRight[1]-UpperLeft[1] +1;

	demBuffer = new short[Lines*Pixels];
	GDALRasterBandH hBand_dem;
	hBand_dem = GDALGetRasterBand(pData_dem, 1);
	GDALRasterIO(hBand_dem, GF_Read, UpperLeft[1], UpperLeft[0], Pixels, Lines, demBuffer, Pixels, Lines, GDT_Int16, 0, 0);
	
}

void RefDem::getIndex(double lat, double lon, double Res[2])
{
	Res[0] = (lat_max - lat) / deltaLat;
	Res[1] = (lon - lon_min) / deltaLon;
	
}
/*****************************************************************/

/*****************************************************************
*               Class:TransFormCoef                              *
*****************************************************************/

void TransFormCoef::getBurstCoeff(int BurstId, double CoeffAz[6], double CoeffRg[6])
{
	if (CpmAz == NULL|| CpmRg==NULL)
	{
		cout << " coefficients are not initialized!" << endl;
		return;
	}

	int burstPos = BurstId - burst0;

	for (int i = 0; i < 6; i++)
	{

		CoeffAz[i] = CpmAz[burstPos * 6 + i];
		CoeffRg[i] = CpmRg[burstPos * 6 + i];
		
	}


}

void TransFormCoef::setBurstCoeff(int BurstId, double CoeffAz[6], double CoeffRg[6])
{
	if (CpmAz == NULL || CpmRg == NULL)
	{
		cout << " coefficients are not initialized!" << endl;
		return;
	}

	int burstPos = BurstId - burst0;

	for (int i = 0; i < 6; i++)
	{

		CpmAz[burstPos * 6 + i] = CoeffAz[i];
		CpmRg[burstPos * 6 + i] = CoeffRg[i];

	}


}

void TransFormCoef::Init(int BurstBeign, int BurstEnd)
{
	if (CpmAz != NULL || CpmRg != NULL)
	{
		cout << "The coefficients arrays have been Initialized!" << endl;
		return;
	}
	
	int NumB = BurstEnd - BurstBeign + 1;
	burst0 = BurstBeign;
	Nbust = NumB;

	CpmAz = new double[6 * NumB];
	CpmRg = new double[6 * NumB];

}

double* TransFormCoef::getAzCoeff(int BurstId)
{
	if (CpmAz == NULL || CpmRg == NULL)
	{
		cout << " coefficients are not initialized!" << endl;
		return NULL;
	}
	int burstPos = BurstId - burst0;
	return CpmAz+burstPos * 6;

}

double* TransFormCoef::getRgCoeff(int BurstId)
{
	if (CpmAz == NULL || CpmRg == NULL)
	{
		cout << " coefficients are not initialized!" << endl;
		return NULL;
	}
	int burstPos = BurstId - burst0;
	return CpmRg + burstPos * 6;

}


void TransFormCoef::clear()
{
	if (CpmAz != NULL)
	{
		delete[] CpmAz;
		CpmAz = NULL;
	}

	if (CpmRg != NULL)
	{
		delete[] CpmRg;
		CpmRg = NULL;
	}
}

/*****************************************************************
*               Class:TiffRead                                   *
*****************************************************************/

void TiffRead::Init(const char* TiffIn)
{
	GDALAllRegister();
	pData_In = GDALOpen(TiffIn, GA_ReadOnly);
	if (pData_In == NULL)
	{
		cout << "can not open " << TiffIn << " !\n";
		system("pause");
		exit(0);
	}

}

void TiffRead::ReadFloat(int x0, int y0, int Lines, int Pixels, float *Buffer)
{

	GDALRasterBandH hBand = GDALGetRasterBand(pData_In, 1);
	GDALRasterIO(hBand, GF_Read, x0, y0, Pixels, Lines,
		Buffer, Pixels, Lines, GDT_Float32, 0, 0);

}

void TiffRead::ReadDouble(int x0, int y0, int Lines, int Pixels, double *Buffer)
{

	GDALRasterBandH hBand = GDALGetRasterBand(pData_In, 1);
	GDALRasterIO(hBand, GF_Read, x0, y0, Pixels, Lines,
		Buffer, Pixels, Lines, GDT_Float64, 0, 0);

}

void TiffRead::ReadCpxShort(int x0,int y0,int Lines, int Pixels, complex<short> *Buffer)
{

	GDALRasterBandH hBand = GDALGetRasterBand(pData_In, 1);
	GDALRasterIO(hBand, GF_Read, x0, y0, Pixels, Lines,
		Buffer, Pixels, Lines, GDT_CInt16, 0, 0);

}
void TiffRead::ReadCpxShort(complex<short>*dataFor, complex<short>*dataBack,
	const int Numoverlap, int* overlapSizeArray, int linesPerBurst,
	int x0, int y0, int width)
{
	GDALRasterBandH hBand = GDALGetRasterBand(pData_In, 1);
	int arrayOffset = 0;

	for (int i = 0; i < Numoverlap; i++)
	{
		int overlapsize = overlapSizeArray[i];
		int OffsetPerOverlapFor = y0 + (i + 1)*linesPerBurst - overlapsize;
		int OffsetPerOverlapBack = y0 + (i + 1)*linesPerBurst;


	

		GDALRasterIO(hBand, GF_Read, x0, OffsetPerOverlapFor,
			width, overlapsize, &dataFor[arrayOffset*width], width, overlapsize, GDT_CInt16, 0, 0);
		

		GDALRasterIO(hBand, GF_Read, x0, OffsetPerOverlapBack, width, overlapsize
			, &dataBack[arrayOffset*width], width, overlapsize, GDT_CInt16, 0, 0);
		
		arrayOffset += overlapsize;

	}
}
void TiffRead::ReadCpxFloat(int x0, int y0, int Lines, int Pixels, complex<float> *Buffer)
{

	GDALRasterBandH hBand = GDALGetRasterBand(pData_In, 1);
	GDALRasterIO(hBand, GF_Read, x0, y0, Pixels, Lines,
		Buffer, Pixels, Lines, GDT_CFloat32, 0, 0);

}
void TiffRead::ReadCpxFloat(complex<float>*dataFor, complex<float>*dataBack,
	const int Numoverlap, int* overlapSizeArray, int linesPerBurst,
	int x0, int y0, int width)
{
	GDALRasterBandH hBand = GDALGetRasterBand(pData_In, 1);
	int arrayOffset = 0;

	for (int i = 0; i < Numoverlap; i++)
	{
		int overlapsize = overlapSizeArray[i];
		int OffsetPerOverlapFor = y0 + (i + 1)*linesPerBurst - overlapsize;
		int OffsetPerOverlapBack = y0 + (i + 1)*linesPerBurst;

		GDALRasterIO(hBand, GF_Read, x0, OffsetPerOverlapFor,
			width, overlapsize, &dataFor[arrayOffset*width], width, overlapsize, GDT_CFloat32, 0, 0);

		GDALRasterIO(hBand, GF_Read, x0, OffsetPerOverlapBack, width, overlapsize
			, &dataBack[arrayOffset*width], width, overlapsize, GDT_CFloat32, 0, 0);

		arrayOffset += overlapsize;

	}
}

void TiffRead::Close()
{
	GDALClose(pData_In);
}


/*****************************************************************
*               Class:TiffWrite                                   *
*****************************************************************/
void TiffWrite::Init(const char* TiffOut, int types, int Pixels, int Lines)
{
	dataTypes = types;


	if (Pixels <= 0 || Lines <= 0)
	{
		cout << "gdal create image error ! width =" << Pixels << "  height=" << Lines << " \n";
		system("pause");
		exit(0);
	}


	GDALAllRegister();
	GDALDriverH hDriver = GDALGetDriverByName("GTiff");
	if (hDriver == NULL)
		cout << "gdal driver error!\n";
	char **papszParmList = NULL;

	if (dataTypes == GDT_CInt16)
	{
		pData_Out = GDALCreate(hDriver, TiffOut, Pixels, Lines, 1, GDT_CInt16, papszParmList);
		double afxgeo[6];
		afxgeo[0] = 0.0;
		afxgeo[1] = 1.0;
		afxgeo[2] = 0.0;
		afxgeo[3] = 0.0;
		afxgeo[4] = 0.0;
		afxgeo[5] = 1.0;
		GDALSetGeoTransform(pData_Out, afxgeo);

	}
	else if (dataTypes == GDT_CFloat32)
	{
		pData_Out = GDALCreate(hDriver, TiffOut, Pixels, Lines, 1, GDT_CFloat32, papszParmList);
		double afxgeo[6];
		afxgeo[0] = 0.0;
		afxgeo[1] = 1.0;
		afxgeo[2] = 0.0;
		afxgeo[3] = 0.0;
		afxgeo[4] = 0.0;
		afxgeo[5] = 1.0;
		GDALSetGeoTransform(pData_Out, afxgeo);

	}
	else if (dataTypes == GDT_Float32)
	{
		pData_Out = GDALCreate(hDriver, TiffOut, Pixels, Lines, 1, GDT_Float32, papszParmList);
		//double afxgeo[6];
		//afxgeo[0] = 0.0;
		//afxgeo[1] = 1.0;
		//afxgeo[2] = 0.0;
		//afxgeo[3] = 0.0;
		//afxgeo[4] = 0.0;
		//afxgeo[5] = 1.0;
		//GDALSetGeoTransform(pData_Out, afxgeo);
	}
	else if (dataTypes == GDT_Float64)
	{
		pData_Out = GDALCreate(hDriver, TiffOut, Pixels, Lines, 1, GDT_Float64, papszParmList);
		double afxgeo[6];
		afxgeo[0] = 0.0;
		afxgeo[1] = 1.0;
		afxgeo[2] = 0.0;
		afxgeo[3] = 0.0;
		afxgeo[4] = 0.0;
		afxgeo[5] = 1.0;
		GDALSetGeoTransform(pData_Out, afxgeo);
	}
	else
	{
		cout << "unsupported gdal data type \n";
		system("pause");
		exit(0);

	}

	if (pData_Out == NULL)
	{
		cout << "can not create " << TiffOut << " !\n";
		system("pause");
		exit(0);
	}
}
void TiffWrite::WriteFloat(int x0, int y0, int Lines, int Pixels, float *Buffer)
{
	if (pData_Out == NULL)
	{
		cout << "please initilize createUpdataImage first !\n";
		system("pause");
		exit(0);
	}
	GDALRasterBandH hBand = GDALGetRasterBand(pData_Out, 1);


	GDALRasterIO(hBand, GF_Write, x0, y0, Pixels, Lines, Buffer, Pixels, Lines, GDT_Float32, 0, 0);
	
}
void TiffWrite::WriteDouble(int x0, int y0, int Lines, int Pixels, double *Buffer)
{
	if (pData_Out == NULL)
	{
		cout << "please initilize createUpdataImage first !\n";
		system("pause");
		exit(0);
	}
	GDALRasterBandH hBand = GDALGetRasterBand(pData_Out, 1);


	GDALRasterIO(hBand, GF_Write, x0, y0, Pixels, Lines, Buffer, Pixels, Lines, GDT_Float64, 0, 0);
}
void TiffWrite::WriteCpxFloat(int x0, int y0, int Lines, int Pixels, complex<float> *Buffer)
{
	if (pData_Out == NULL)
	{
		cout << "please initilize createUpdataImage first !\n";
		system("pause");
		exit(0);
	}
	GDALRasterBandH hBand = GDALGetRasterBand(pData_Out, 1);


	GDALRasterIO(hBand, GF_Write, x0, y0, Pixels, Lines, Buffer, Pixels, Lines, GDT_CFloat32, 0, 0);
}

void TiffWrite::Close()
{
	GDALClose(pData_Out);
}