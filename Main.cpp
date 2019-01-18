#include <iostream>
#include<stdio.h>
#include<fstream>
#include "DataTypes.h"
#include "cuComplex.h"
#include "Func.h"



int main()
{
	ConfigSet S1Config;
	string ConfigPath="config.txt";
	if (!ReadInput(S1Config, ConfigPath.c_str()))
	{
		cout << "Failed in Openning Config File!" << endl;
		system("pause");
		return 1;
	}
		
	CheckDir(S1Config.process_dir);

	cout << "Current Working Catalogue:" << S1Config.process_dir << endl;


	SentinelTOPS M_S1;
	SentinelTOPS S_S1;
	//Reading InFo of Master and Slave SLCs
	S1_Initialize(M_S1, S1Config.masterpath.c_str(), S1Config.polarisation);
	S1_Initialize(S_S1, S1Config.slavepath.c_str(), S1Config.polarisation);

	S1PreciseOrbit Morbit;
	S1PreciseOrbit Sorbit;

	//Reading Precise Orbit Info 
	S1_OrbitInitialize(Morbit, S1Config.preciseOrbitMaster.c_str(), M_S1.SubSwath[0].firstLineUTC.c_str(),
		M_S1.SubSwath[0].lastLineUTC.c_str(),M_S1.MissionID);
	S1_OrbitInitialize(Sorbit, S1Config.preciseOrbitSlave.c_str(), S_S1.SubSwath[0].firstLineUTC.c_str(),
		S_S1.SubSwath[0].lastLineUTC.c_str(),S_S1.MissionID);



	if (!CheckSpecificDem(S1Config.SpecificDemPath.c_str(), M_S1.lat_min, M_S1.lat_max,
		M_S1.lon_min, M_S1.lon_max))
	{
		cout << "Something wrong with DEM!" << endl;
		system("pause");
		exit(0);
	}

	/********************************************************************/
	//Estimate the polynomials based on the orbit information
	/********************************************************************/
	Polyfit_Orbit(Morbit);
	Polyfit_Orbit(Sorbit);

	/********************************************************************/
	/*Compute Deramping Parameters*/
	/********************************************************************/
	ComputeDerapPara(M_S1, Morbit);
	ComputeDerapPara(S_S1, Sorbit);


	/********************************************************************/
	/*InitializeDemFile and ResampleTable*/
	/********************************************************************/
	RefDem TopoDEM;
	TopoDEM.Init(S1Config.SpecificDemPath.c_str());

	ResampleTable S1ReTable;
	InKernelInit(S1ReTable, 12, M_S1.getPRF(), M_S1.getABW(),
		M_S1.getRSR2X());

	TransFormCoef TransC;
	char buff[5];
	/********************************************************************/
	/*InSAR processing for S1 TOPS data. Coregistration, resampling, ESD, 
	Coherence Estimation modules are serially performed.*/
	/********************************************************************/
	for (int NSubS = S1Config.firstsubswath; NSubS <= S1Config.lastsubswath; NSubS++)
	{
		sprintf(buff, "-%01d", NSubS);
		string ReSlaveFile = S1Config.process_dir + "\\ReSlave" + string(buff) + ".tif";
		string ReSlaveESDFile = S1Config.process_dir + "\\ReSlaveESD" + string(buff) + ".tif";
		string CohFile = S1Config.process_dir + "\\Coh" + string(buff) + ".tif";

		/********************************************************************/
		/*Check and estimate the burst overlap between Master and Slave*/
		/********************************************************************/
		int M_burst0, S_burst0;
		int M_burstN, S_burstN;

		ComputerBurstOverlap(S1Config.burst0, S1Config.burstN, M_burst0, M_burstN,
			S_burst0, S_burstN, M_S1, S_S1, NSubS, Morbit, Sorbit);
		
		
		TransC.Init(M_burst0 , M_burstN );

		printf("Start geometric coregistration on subswath %d \n", NSubS);

		 //Geomteric Coregistration
		GeometricCoreg(M_S1.SubSwath[NSubS-1], S_S1.SubSwath[NSubS-1], Morbit, Sorbit, TopoDEM, TransC, M_burst0,
			M_burstN, S_burst0, S_burstN, NSubS);

		printf("Start deramping and resampling on subswath %d \n", NSubS);
		////The first alligment: 1.Deramping and Demodulation; 2.Resampling 3.Reramping
		DerampAndResample(M_S1.SubSwath[NSubS - 1], S_S1.SubSwath[NSubS - 1], TransC, M_burst0, M_burstN,
			S_burst0, S_burstN, NSubS, S1ReTable,ReSlaveFile.c_str(), 12);

//		system("pause");

		printf("Start ESD and Coherence estimation on subswath %d \n", NSubS);
		//Perform ESD estimation, resampling, and Coherence estimation
		ESDAndCoh(M_S1.SubSwath[NSubS - 1], S_S1.SubSwath[NSubS - 1], TransC,
			M_burst0, M_burstN, S_burst0, S_burstN, NSubS, S1ReTable,
			ReSlaveFile.c_str() , ReSlaveESDFile.c_str(), CohFile.c_str());

	}

	//clearing
	M_S1.clear();
	S_S1.clear();
	Morbit.clear();
	Sorbit.clear();
	TransC.clear();
	S1ReTable.clear();






	cout << "Finished!" << endl;




	return 1;
}





































 













 

 