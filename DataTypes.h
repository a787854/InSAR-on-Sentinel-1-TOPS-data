
#ifndef DATATYPES_H
#define DATATYPES_H

#include <string>
#include <iostream>
#include <complex>
#include <fstream>
#include "ogr_api.h"
#include "ogr_core.h"
#include "ogr_srs_api.h"
#include "gdal.h"

using namespace std;

struct Polynomial{
	double aztime;
	double t0;
	double c0;
	double c1;
	double c2;
};

struct ConfigSet
{
	string process_dir;
	string masterpath;
	string slavepath;
	int firstsubswath;
	int lastsubswath;
	string polarisation;
	string preciseOrbitMaster;
	string preciseOrbitSlave;
	int multi_az;
	int multi_rg;
	int burst0;
	int burstN;
	bool SpecificDemorNo;
	string SpecificDemPath;
};

struct LatLonMaxMin
{
	double lat_max;
	double lat_min;
	double lon_max;
	double lon_min;

};

class SubSwathInfo
{
public:
	//subswath info
	string subSwathName;
	string ImgPath;
	string XmlPath;
	string polarisation;
	string firstLineUTC;
	string lastLineUTC;
	int numOfLines;
	int numOfSamples;
	double firstLineTime;
	double lastLineTime;
	double firstValidLineTime;
	double lastValidLineTime;
	double slrTimeToFirstPixel;//one way [s]
	double slrTimeToLastPixel;
	double slrTimeToFirstValidPixel;
	double slrTimeToLastValidPixel;
	double azimuthTimeInterval;
	double rangePixelSpacing;
	double azimuthPixelSpacing;
	double radarFrequency;
	double wavelength; // in m
	double rangeSamplingRate;
	double azimuthSteeringRate; //[degrees / s]
	double ascendingNodeTime;
	int firstValidPixel;
	int lastValidPixel;
	double range_bandwidth;
	double azimuth_bandwidth;
	double prf;

	double headingAngle;



	// bursts info
	int numOfBursts;
	int linesPerBurst;
	int samplesPerBurst;
	double burstFirstLineTime[16];//16 is enough to store
	double burstLastLineTime[16];
	double burstFirstValidLineTime[16];
	double burstLastValidLineTime[16];
	int firstValidLine[16];
	int lastValidLine[16];

	int *firstValidSample;
	int *lastValidSample;
	double* rangeDependDopplerRate;
	double* dopplerRate;
	double* referenceTime;
	double* dopplerCentroid;


	// GeoLocationGridPoint
	int numOfGeoLines;
	int numOfGeoPointsPerLine;
	double *azimuthTime;
	double *slantRangeTime; //one way [s]
	double *latitude;
	double *longitude;
	double *incidenceAngle;

	bool nearRangeOnLeft;
	

	LatLonMaxMin latlon;
	LatLonMaxMin *latlonBurst;

	//Coarse Orbit Information attached in S1 leader file
	int numOfOrbit;
	double* orbitAzTime = NULL;
	double* x_pos = NULL;
	double* y_pos  = NULL;
	double* z_pos = NULL;
	double* x_vel = NULL;
	double* y_vel = NULL;
	double* z_vel = NULL;


	int clear();



};
class SentinelTOPS
{
public:
	SubSwathInfo *SubSwath;
	string MissionID;
	double azimuthSpacing;
	double rangeSpacing;
	int NumSubSwath;

	double lat_max;
	double lat_min;
	double lon_max;
	double lon_min;

	double getPRF();//Get azimuth sampling rate
	double getABW();//Get azimuth bandwidth
	double getRSR2X();// Get two times of range sampling rate
	void clear();

};

class S1PreciseOrbit
{
public:
	int NumPoints;
	int NumCoeff;
	double TimeInerval;
	double* orbitAzTime;
	double* x_pos;
	double* y_pos;
	double* z_pos;
	double* x_vel;
	double* y_vel;
	double* z_vel;

	double *coef_x;
	double *coef_y;
	double *coef_z;

	void clear();
};

class ellipsoid_WGS84
{
public:
	double e2;
	double e2b;
	double a;
	double b;

	ellipsoid_WGS84();
};

class  ResampleTable
{
public:
	int Npoints;
	float* KernelAz=NULL;
	float* KernelRg=NULL;


	void clear();

};

class RefDem
{
	public:

	double deltaLat;
	double deltaLon;
	double lon_min;
	double lat_max; 
	double lat_min; 
	double lon_max;
	int Lines;
	int Pixels;

	GDALDatasetH pData_dem;
	void Init(const char* DEMPath);
	~RefDem();
	void getData(int x0, int y0, int ww, int hh, short* demBuffer);
	void getData(double lat_min, double lat_max, double lon_min, double lon_max, double extralat,
		double extralon, short*&demBuffer, int& Lines, int& Pixels);
	void getData(double lat_min, double lat_max, double lon_min, double lon_max, double extralat,
		double extralon, int*&demBuffer, int& Lines, int& Pixels);
	void getIndex(double lat, double lon, double Res[2]);
	
};

//6 parameters to describe the translation relationship between master and slave
class TransFormCoef
{
public:
	int burst0;
	int Nbust;
	double *CpmAz=NULL;
	double *CpmRg = NULL;

	void Init(int BusrtBeign, int BusrtEnd);
	void setBurstCoeff(int BurstId, double CoeffAz[6], double CoeffRg[6]);
	void getBurstCoeff(int BurstId, double CoeffAz[6], double CoeffRg[6]);
	double* getAzCoeff(int BurstId);
	double* getRgCoeff(int BurstId);
	void clear();
	
	
};
/*************************************************************
*                     Class: TiffRead                       *
*************************************************************/

class TiffRead
{
private:
	GDALDatasetH pData_In;
public:
	void Init(const char* TiffIn);
	void ReadFloat(int x0, int y0, int Lines, int Pixels, float *Buffer);
	void ReadDouble(int x0, int y0, int Lines, int Pixels, double *Buffer);
	void ReadCpxShort(int x0, int y0, int Lines, int Pixels, complex<short> *Buffer);
	void ReadCpxShort(complex<short>*dataFor, complex<short>*dataBack,
		const int Numoverlap, int* overlapSizeArray, int linesPerBurst,
		int x0, int y0, int width);
	void ReadCpxFloat(int x0, int y0, int Lines, int Pixels, complex<float> *Buffer);
	void ReadCpxFloat(complex<float>*dataFor, complex<float>*dataBack,
		const int Numoverlap, int* overlapSizeArray, int linesPerBurst,
		int x0, int y0, int width);
	void Close();
};

/*************************************************************
*                     Class: TiffWrite                       *
*************************************************************/
class TiffWrite
{
private:
	GDALDatasetH pData_Out;
	int dataTypes;
public:
	void Init(const char* TiffOut, int types, int Pixels, int Lines);
	void WriteFloat(int x0, int y0, int Lines, int Pixels, float *Buffer);
	void WriteDouble(int x0, int y0, int Lines, int Pixels, double *Buffer);
	void WriteCpxFloat(int x0, int y0, int Lines, int Pixels, complex<float> *Buffer);
	void Close();
};
#endif // DATATYPES_H
