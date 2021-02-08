#ifndef FUNC_H
#define FUNC_H
#include "DataTypes.h"
#include "io.h"
#include "tinyxml2.h"
#include "Constants.h"
#include "cuComplex.h"
#include <ctime>

using namespace tinyxml2;

//Read ConfigFile
bool ReadInput(ConfigSet & ConfigInput, const char* ConfigPath);
bool CheckDir(string Dir);
//Check Image and Xml files
int GetFiles(string path, string * files, string polar);
//Initialize S1 product
bool S1_Initialize(SentinelTOPS& S1, const char* Path, string polar);

//Read Orbit
double utc2seconds(const char* utc);
void getIntArray(string str, char gap, int *Array);
void getDoubleArray(string str, char gap, double *Array);
int DiffUTCTime(string firstUTC, string lastUTC);
bool CheckOribtFile(string PreciseOrbit, string UTCtime, string MissionID);
bool S1_OrbitInitialize(S1PreciseOrbit &Orbit,
	const  char * Path, const char *firstUTC, const char *lastUTC, string MissionID);

//ATA
void matTxmat(double *A, double *B, double *Res,
	int ALines, int ResLines, int ResPixels);

// Cholesky factorization 
void choles(double * A, int N);

//Solving  Cholesky factorization 
void solvechol(
	double *A, double* B, int ALines, int APixels, int BLines, int BPixels);

void invertchol(
	double* A, int N);

bool Polyfit_Orbit(S1PreciseOrbit &Orbit);

double getLatitude(SentinelTOPS& S1, double aztime, double rgtime, int SubSwathID);
double getLontitude(SentinelTOPS& S1, double aztime, double rgtime, int SubSwathID);
void Locating(SentinelTOPS& S1, double aztime, double rgtime, int SubSwathID,
	int*Index, double * coef);
void Geodetic2Geographic(ellipsoid_WGS84 &ell, double lat, double lon, double height,
	double Res[3]);
inline  void getXyz(double t, double* possat, double* coef_x,
	double *coef_y, double* coef_z, int n);
inline  void getVelocity(double t, double* velsat, double* coef_x,
	double *coef_y, double* coef_z, int n);
inline void getAcceleration(double t, double*accsat, double* coef_x,
	double *coef_y, double* coef_z, int n);
double xyz2aztime(double *xyz, int NumOfCoef,
	double *coef_x, double* coef_y, double* coef_z,
	double Initialaztime, int midAztime);

bool ComputerBurstOverlap(int InputBurst0, int InputBurstN, int &M_burst0, int &M_burstN,
	int &S_burst0, int &S_burstN, SentinelTOPS &M_S1, SentinelTOPS &S_S1,
	int SubSwathIndex, S1PreciseOrbit& Morbit, S1PreciseOrbit& Sorbit);

bool CheckSpecificDem(const char *dempath, double lat_min, double lat_max,
	double lon_min, double lon_max);
bool ComputeDerampPara(SentinelTOPS &S1, S1PreciseOrbit &PreciseOrbit);
//Init resample look up table
bool InKernelInit(ResampleTable &ReTable, int Npoints, double prf, double abw, double rsr2x);
//Geometrical Coregistration
int GeometricCoreg(SubSwathInfo& M_Swath, SubSwathInfo& S_Swath, S1PreciseOrbit &Morbit, S1PreciseOrbit &Sorbit,
	RefDem& TopoDem, TransFormCoef TFCoef, int Mburst0, int MburstN, int Sburst0, int SburstN, int SwathId);

//deramping and resampling
int DerampAndResample(SubSwathInfo& M_Swath, SubSwathInfo& S_Swath, TransFormCoef& TFCoef,
	int Mburst0, int MburstN, int Sburst0, int SburstN, int SwathId, ResampleTable& ReTable,
	const char* OutPath, int RePoints);

//ESD and coherence Estimation
int ESDAndCoh(SubSwathInfo& M_Swath, SubSwathInfo& S_Swath, TransFormCoef& TFCoef,
	int Mburst0, int MburstN, int Sburst0, int SburstN, int SwathId, ResampleTable& ReTable,
	const char *ReSlaveFile, const char* ReSlaveESDFile, const char*Cohout);

int Stitch_OneSubswath(SubSwathInfo& M_Swath, const char* GappedImg, const char* StitchFile, int m_burst0, int m_burstN,
	int data_type);





#endif // DATATYPES_H