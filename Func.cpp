
#include "Func.h"

/*****************************************************************/
/*Check config file                                              */
/*****************************************************************/
bool ReadInput(ConfigSet &ConfigInput, const char* ConfigPath)
{
	ifstream ConfigFile(ConfigPath);
	if (!ConfigFile)
	{
		cout << "The path dose not existed! Please check it!" << endl;
		return false;
	}
	bool GoOnReading = 1;

	while (GoOnReading)
	{

		if (ConfigFile.getline(eachLine, 4 * 512, '\n'));
		const int maxwords = 2;
		char *word[maxwords];

		char* c = eachLine;
		for (int i = 0; i < maxwords; i++)
		{
			while (*c == ' ' || *c == '\t' || *c == '=') // get rid of leadin white space
			{
				c++;                            // next address
			}
			word[i] = c;                    // pass first char address of a word

			while (*c != '\t' && *c != '\0'&&*c != '=')
			{
				c++;
			}

			if (*c != '\0') // at last char
			{
				*c = '\0';
				c++;
			}
		}

		char*keyword = word[0];
		if (!strcmp(keyword, "masterpath"))
		{

			ConfigInput.masterpath = word[1];

		}
		else if (!strcmp(keyword, "slavepath"))
		{

			ConfigInput.slavepath = word[1];

		}
		else if (!strcmp(keyword, "firstsubswath"))
		{

			ConfigInput.firstsubswath = atoi(word[1]);

		}
		else if (!strcmp(keyword, "lastsubswath"))
		{

			ConfigInput.lastsubswath = atoi(word[1]);

		}
		else if (!strcmp(keyword, "lastsubswath"))
		{

			ConfigInput.lastsubswath = atoi(word[1]);

		}
		else if (!strcmp(keyword, "polarisation"))
		{

			ConfigInput.polarisation = word[1];

		}
		else if (!strcmp(keyword, "preciseOrbitMaster"))
		{

			ConfigInput.preciseOrbitMaster = word[1];

		}
		else if (!strcmp(keyword, "preciseOrbitSlave"))
		{
			ConfigInput.preciseOrbitSlave = word[1];
		}
		else if (!strcmp(keyword, "SpecificDemPath"))
		{
			ConfigInput.SpecificDemPath = word[1];
		}
		else if (!strcmp(keyword, "process_dir"))
		{
			ConfigInput.process_dir = word[1];
		}

		else if (!strcmp(keyword, "multi_az"))
		{
			ConfigInput.multi_az = atoi(word[1]);
		}

		else if (!strcmp(keyword, "multi_rg"))
		{
			ConfigInput.multi_rg = atoi(word[1]);
		}

		else if (!strcmp(keyword, "burst0"))
		{
			ConfigInput.burst0 = atoi(word[1]);
		}
		else if (!strcmp(keyword, "burstN"))
		{
			ConfigInput.burstN = atoi(word[1]);
		}


		if (ConfigFile.eof())
		{
			GoOnReading = false;
		}

	}

	return true;

}

/*****************************************************************/
/*Initialize the S1 product                                      */
/*****************************************************************/
bool S1_Initialize(SentinelTOPS& S1, const char* Path, string polar)
{
	string polarFlag = "-slc-" + polar;
	string imgpath = string(Path) + "\\measurement";

	string ImgFiles[12];
	int NumImg;

	NumImg = GetFiles(imgpath, ImgFiles, polar);
	S1.NumSubSwath = NumImg;
	string XmlFiles[12];
	int NumXml;

	string xmlpath = string(Path) + "\\annotation";
	NumXml = GetFiles(xmlpath, XmlFiles, polar);

	//Set the lat/lon box of the frame
	S1.lat_max = S1.lon_max = -999;
	S1.lat_min = S1.lon_min = 999;



	if (NumImg <= 0)
	{
		cout << "number of subswath = " << NumImg << endl;
		system("pause");
		exit(0);
	}


	if (NumXml < NumImg)
	{
		cout << "annotation xml files not equal numOfSubSwath \n";
		system("pause");
		exit(0);
	}
	S1.SubSwath = new SubSwathInfo[NumImg];

	for (int i = 0; i < NumImg; i++)
	{
		S1.SubSwath[i].ImgPath = ImgFiles[i];
		S1.SubSwath[i].XmlPath = XmlFiles[i];
	}


	if (strcmp(S1.SubSwath[0].ImgPath.c_str(), "S1A") != NULL)
	{
		S1.MissionID = "S1A";
	}
	else if (strcmp(S1.SubSwath[0].ImgPath.c_str(), "S1B") != NULL)
	{
		S1.MissionID = "S1B";
	}
	else
	{
		cout << "Unkown Input SLC , It's assumed to be S1A or S1B!" << endl;
		system("pause");
		exit(0);
	}

	double SumRangeSpacing = 0.0;
	double SumAzimuthSpacing = 0.0;
	for (int i = 0; i < NumImg; i++)
	{


		XMLDocument xmldoc;

		bool LoadOK = xmldoc.LoadFile(XmlFiles[i].c_str());
		if (LoadOK)
		{
			cout << "The path is not valid! Please check it!" << endl;
			cout << XmlFiles[i] << endl;
			system("pause");
			exit(0);
		}

		XMLElement* root = xmldoc.RootElement();

		//Polarisation
		XMLElement* elem = root->FirstChildElement("adsHeader")->FirstChildElement("polarisation");
		XMLNode* node = elem->FirstChild();
		S1.SubSwath[i].polarisation = node->Value();

		//Number of Lines
		elem = root->FirstChildElement("imageAnnotation")->FirstChildElement("imageInformation")->FirstChildElement("numberOfLines");
		node = elem->FirstChild();
		S1.SubSwath[i].numOfLines = atoi(node->Value());

		//Number of Samples
		elem = root->FirstChildElement("imageAnnotation")->FirstChildElement("imageInformation")->FirstChildElement("numberOfSamples");
		node = elem->FirstChild();
		S1.SubSwath[i].numOfSamples = atoi(node->Value());

		//Time of First Line
		elem = root->FirstChildElement("imageAnnotation")->FirstChildElement("imageInformation")->FirstChildElement("productFirstLineUtcTime");
		node = elem->FirstChild();
		S1.SubSwath[i].firstLineUTC = node->Value();
		S1.SubSwath[i].firstLineTime = utc2seconds(node->Value());


		//Time of Last Line
		elem = root->FirstChildElement("imageAnnotation")->FirstChildElement("imageInformation")->FirstChildElement("productLastLineUtcTime");
		node = elem->FirstChild();
		S1.SubSwath[i].lastLineUTC = node->Value();
		S1.SubSwath[i].lastLineTime = utc2seconds(node->Value());

		//AzimithTime Interval
		elem = root->FirstChildElement("imageAnnotation")->FirstChildElement("imageInformation")->FirstChildElement("azimuthTimeInterval");
		node = elem->FirstChild();
		S1.SubSwath[i].azimuthTimeInterval = atof(node->Value());

		//Range Pixel Spacing 
		elem = root->FirstChildElement("imageAnnotation")->FirstChildElement("imageInformation")->FirstChildElement("rangePixelSpacing");
		node = elem->FirstChild();
		S1.SubSwath[i].rangePixelSpacing = atof(node->Value());
		SumRangeSpacing += S1.SubSwath[i].rangePixelSpacing;

		//Azimuth Pixel Spacing
		elem = root->FirstChildElement("imageAnnotation")->FirstChildElement("imageInformation")->FirstChildElement("azimuthPixelSpacing");
		node = elem->FirstChild();
		S1.SubSwath[i].azimuthPixelSpacing = atof(node->Value());
		SumAzimuthSpacing += S1.SubSwath[i].azimuthPixelSpacing;

		//Radar Carrier Frequency
		elem = root->FirstChildElement("generalAnnotation")->FirstChildElement("productInformation")->FirstChildElement("radarFrequency");
		node = elem->FirstChild();
		S1.SubSwath[i].radarFrequency = atof(node->Value());
		S1.SubSwath[i].wavelength = (double)SOL / S1.SubSwath[i].radarFrequency;

		//Range Sampling Rate
		elem = root->FirstChildElement("generalAnnotation")->FirstChildElement("productInformation")->FirstChildElement("rangeSamplingRate");
		node = elem->FirstChild();
		S1.SubSwath[i].rangeSamplingRate = atof(node->Value());

		//Azimuth Steering Rate
		elem = root->FirstChildElement("generalAnnotation")->FirstChildElement("productInformation")->FirstChildElement("azimuthSteeringRate");
		node = elem->FirstChild();
		S1.SubSwath[i].azimuthSteeringRate = atof(node->Value());

		//Satellite Heading
		elem = root->FirstChildElement("generalAnnotation")->FirstChildElement("productInformation")->FirstChildElement("platformHeading");
		node = elem->FirstChild();
		S1.SubSwath[i].headingAngle = atof(node->Value());

		//ascendingNodeTime
		elem = root->FirstChildElement("imageAnnotation")->FirstChildElement("imageInformation")->FirstChildElement("ascendingNodeTime");
		node = elem->FirstChild();
		S1.SubSwath[i].ascendingNodeTime = utc2seconds(node->Value());

		//prf= 1/AzimuthInterval
		S1.SubSwath[i].prf = double(1.0) / S1.SubSwath[i].azimuthTimeInterval;

		//azimuth BandWidth
		elem = root->FirstChildElement("imageAnnotation")->FirstChildElement("processingInformation")
			->FirstChildElement("swathProcParamsList")->FirstChildElement("swathProcParams")->
			FirstChildElement("azimuthProcessing")->FirstChildElement("processingBandwidth");


		node = elem->FirstChild();
		S1.SubSwath[i].azimuth_bandwidth = atof(node->Value());

		//range BandWidth
		elem = root->FirstChildElement("imageAnnotation")->FirstChildElement("processingInformation")
			->FirstChildElement("swathProcParamsList")->FirstChildElement("swathProcParams")->
			FirstChildElement("rangeProcessing")->FirstChildElement("processingBandwidth");
		node = elem->FirstChild();
		S1.SubSwath[i].range_bandwidth = atof(node->Value());


		//Fast time of first Pixels
		elem = root->FirstChildElement("imageAnnotation")->
			FirstChildElement("imageInformation")->FirstChildElement("slantRangeTime");
		node = elem->FirstChild();
		S1.SubSwath[i].slrTimeToFirstPixel = atof(node->Value()) / 2.0;

		S1.SubSwath[i].slrTimeToLastPixel = S1.SubSwath[i].slrTimeToFirstPixel +
			(S1.SubSwath[i].numOfSamples - 1)*S1.SubSwath[i].rangePixelSpacing / SOL;

		//Number of bursts
		elem = root->FirstChildElement("swathTiming")->
			FirstChildElement("burstList");
		S1.SubSwath[i].numOfBursts = atoi(elem->FirstAttribute()->Value());

		//Lines Per Burst
		elem = root->FirstChildElement("swathTiming")->
			FirstChildElement("linesPerBurst");
		S1.SubSwath[i].linesPerBurst = atoi(elem->FirstChild()->Value());

		//Pixels Per Burst
		elem = root->FirstChildElement("swathTiming")->
			FirstChildElement("samplesPerBurst");
		S1.SubSwath[i].samplesPerBurst = atoi(elem->FirstChild()->Value());


		//Read Burst Info
		int NumBursts = S1.SubSwath[i].numOfBursts;
		int LinesPerBurst = S1.SubSwath[i].linesPerBurst;

		S1.SubSwath[i].firstValidSample = new int[NumBursts*LinesPerBurst];
		S1.SubSwath[i].lastValidSample = new int[NumBursts*LinesPerBurst];

		S1.SubSwath[i].firstValidPixel = 0;
		S1.SubSwath[i].lastValidPixel = S1.SubSwath[i].numOfSamples;

		int firstValidPixel_temp = S1.SubSwath[i].numOfSamples;
		int lastValidPixel_temp = 0;

		elem = root->FirstChildElement("swathTiming")->FirstChildElement("burstList");
		XMLElement* elem1, *elem2;
		elem1 = elem->FirstChildElement();
		string TmpStr;
		S1.SubSwath[i].firstValidSample = new int[S1.SubSwath[i].numOfBursts*S1.SubSwath[i].linesPerBurst];
		S1.SubSwath[i].lastValidSample = new int[S1.SubSwath[i].numOfBursts*S1.SubSwath[i].linesPerBurst];
		for (int k = 0; k < S1.SubSwath[i].numOfBursts; k++)
		{
			//Burst First Azimuth Time
			elem2 = elem1->FirstChildElement("azimuthTime");;
			S1.SubSwath[i].burstFirstLineTime[k] = utc2seconds(elem2->GetText());

			//Burst Last Azimuth Time
			S1.SubSwath[i].burstLastLineTime[k] = S1.SubSwath[i].burstFirstLineTime[k] +
				(S1.SubSwath[i].linesPerBurst - 1)*S1.SubSwath[i].azimuthTimeInterval;

			//Array of First Valid Samples for each burst
			TmpStr = elem1->FirstChildElement("firstValidSample")->GetText();
			getIntArray(TmpStr, ' ', &S1.SubSwath[i].firstValidSample[k*LinesPerBurst]);


			//Array of last valid samples for each burst
			TmpStr = elem1->FirstChildElement("lastValidSample")->GetText();
			getIntArray(TmpStr, ' ', &S1.SubSwath[i].lastValidSample[k*LinesPerBurst]);



			int firstValidLineIdx = -1;
			int lastValidLineIdx = -1;

			for (int lineIdx = 0; lineIdx < LinesPerBurst; lineIdx++)
			{
				if (S1.SubSwath[i].firstValidSample[k*LinesPerBurst + lineIdx] != -1)
				{
					if (S1.SubSwath[i].firstValidSample[k*LinesPerBurst + lineIdx] < firstValidPixel_temp)
					{
						firstValidPixel_temp = S1.SubSwath[i].firstValidSample[k*LinesPerBurst + lineIdx];
					}
					if (firstValidLineIdx == -1)
					{
						firstValidLineIdx = lineIdx;
						lastValidLineIdx = lineIdx;
					}
					else
					{
						lastValidLineIdx++;
					}

				}

			}

			for (int lineIdx = 0; lineIdx < LinesPerBurst; lineIdx++)
			{
				if (S1.SubSwath[i].lastValidSample[k*LinesPerBurst + lineIdx] != -1 && S1.SubSwath[i].lastValidSample[k*LinesPerBurst + lineIdx]>lastValidPixel_temp)
					lastValidPixel_temp = S1.SubSwath[i].lastValidSample[k*LinesPerBurst + lineIdx];
			}

			S1.SubSwath[i].burstFirstValidLineTime[k] = S1.SubSwath[i].burstFirstLineTime[k] + firstValidLineIdx*S1.SubSwath[i].azimuthTimeInterval;

			S1.SubSwath[i].burstLastValidLineTime[k] = S1.SubSwath[i].burstFirstLineTime[k] + lastValidLineIdx*S1.SubSwath[i].azimuthTimeInterval;

			S1.SubSwath[i].firstValidLine[k] = firstValidLineIdx;
			S1.SubSwath[i].lastValidLine[k] = lastValidLineIdx;




			elem1 = elem1->NextSiblingElement();
		}

		//Valid Lines
		S1.SubSwath[i].firstValidLineTime = S1.SubSwath[i].burstFirstValidLineTime[0];
		S1.SubSwath[i].firstValidLineTime = S1.SubSwath[i].burstFirstValidLineTime[NumBursts - 1];

		//Valid Pixels
		S1.SubSwath[i].firstValidPixel = firstValidPixel_temp;
		S1.SubSwath[i].lastValidPixel = lastValidPixel_temp;
		//fast time
		S1.SubSwath[i].slrTimeToFirstValidPixel = S1.SubSwath[i].slrTimeToFirstPixel +
			S1.SubSwath[i].firstValidPixel*S1.SubSwath[i].rangePixelSpacing / SOL;
		S1.SubSwath[i].slrTimeToLastValidPixel = S1.SubSwath[i].slrTimeToFirstPixel +
			S1.SubSwath[i].lastValidPixel*S1.SubSwath[i].rangePixelSpacing / SOL;



		/************************************************************/
		//Reading GeoLocationGridPoints
		/************************************************************/
		elem = root->FirstChildElement("geolocationGrid");
		elem = elem->FirstChildElement("geolocationGridPointList");
		int numOfLocationPoints = atoi(elem->FirstAttribute()->Value());
		if (numOfLocationPoints <= 0)
		{
			cout << "numOfGeoLocationGridPoints < 0" << endl;
			system("pause");
			exit(0);
		}


		int numOfGeoPointsPerLine = 0;
		int line = 0;
		elem1 = elem->FirstChildElement();
		for (int k = 0; k < numOfLocationPoints; k++)
		{
			string str = elem1->FirstChildElement("line")->GetText();

			if (numOfLocationPoints == 0)
			{
				line = atoi(str.c_str());
				numOfGeoPointsPerLine++;
			}
			else if (line == atoi(str.c_str()))
			{
				numOfGeoPointsPerLine++;
			}
			else
			{
				break;
			}
			elem1 = elem1->NextSiblingElement();
		}



		S1.SubSwath[i].numOfGeoLines = numOfLocationPoints / numOfGeoPointsPerLine;
		S1.SubSwath[i].numOfGeoPointsPerLine = numOfGeoPointsPerLine;

		//Allocate memory space
		S1.SubSwath[i].azimuthTime = new double[numOfLocationPoints];
		S1.SubSwath[i].slantRangeTime = new double[numOfLocationPoints];
		S1.SubSwath[i].latitude = new double[numOfLocationPoints];
		S1.SubSwath[i].longitude = new double[numOfLocationPoints];
		S1.SubSwath[i].incidenceAngle = new double[numOfLocationPoints];
		S1.SubSwath[i].latlonBurst = new LatLonMaxMin[NumBursts];



		elem1 = elem->FirstChildElement();
		int Pixels = numOfGeoPointsPerLine;
		for (int k = 0; k < numOfLocationPoints; k++)
		{
			int row = k / Pixels;
			int col = k % Pixels;

			string str = elem1->FirstChildElement("azimuthTime")->GetText();
			S1.SubSwath[i].azimuthTime[row*Pixels + col] = utc2seconds(str.c_str());

			str = elem1->FirstChildElement("slantRangeTime")->GetText();
			S1.SubSwath[i].slantRangeTime[row*Pixels + col] = atof(str.c_str()) / 2.0;

			str = elem1->FirstChildElement("latitude")->GetText();
			S1.SubSwath[i].latitude[row*Pixels + col] = atof(str.c_str());

			str = elem1->FirstChildElement("longitude")->GetText();
			S1.SubSwath[i].longitude[row*Pixels + col] = atof(str.c_str());

			str = elem1->FirstChildElement("incidenceAngle")->GetText();
			S1.SubSwath[i].incidenceAngle[row*Pixels + col] = atof(str.c_str());







			elem1 = elem1->NextSiblingElement();
		}
		S1.SubSwath[i].nearRangeOnLeft = (S1.SubSwath[i].incidenceAngle[0] < S1.SubSwath[i].incidenceAngle[1]);
		/************************************************************/

		/************************************************************/
		//Finding the range box of each busrt
		/************************************************************/
		for (int b = 0; b < NumBursts; b++)
		{
			double lat1, lat2, lat3, lat4;
			double lon1, lon2, lon3, lon4;

			//Read latitudes of four corners
			lat1 = S1.SubSwath[i].latitude[b*Pixels + 0];//Left Top
			lat2 = S1.SubSwath[i].latitude[b*Pixels + Pixels - 1];//Right Top
			lat3 = S1.SubSwath[i].latitude[(b + 1)*Pixels + 0];//Left Bottom
			lat4 = S1.SubSwath[i].latitude[(b + 1)*Pixels + Pixels - 1];//Left Bottom

			//Read longitudes of four corners
			lon1 = S1.SubSwath[i].longitude[b*Pixels + 0];//Left Top
			lon2 = S1.SubSwath[i].longitude[b*Pixels + Pixels - 1];//Right Top
			lon3 = S1.SubSwath[i].longitude[(b + 1)*Pixels + 0];//Left Bottom
			lon4 = S1.SubSwath[i].longitude[(b + 1)*Pixels + Pixels - 1];//Left Bottom

			//Find the biggest and smallest lat/lon values
			double lat_max = max4(lat1, lat2, lat3, lat4);
			double lat_min = min4(lat1, lat2, lat3, lat4);
			double lon_max = max4(lon1, lon2, lon3, lon4);
			double lon_min = min4(lon1, lon2, lon3, lon4);

			if (S1.lat_max < lat_max)
				S1.lat_max = lat_max;

			if (S1.lat_min>lat_min)
				S1.lat_min = lat_min;

			if (S1.lon_max < lon_max)
				S1.lon_max = lon_max;

			if (S1.lon_min>lon_min)
				S1.lon_min = lon_min;

			S1.SubSwath[i].latlonBurst[b].lat_max = lat_max;
			S1.SubSwath[i].latlonBurst[b].lat_min = lat_min;
			S1.SubSwath[i].latlonBurst[b].lon_max = lon_max;
			S1.SubSwath[i].latlonBurst[b].lon_min = lon_min;

		}
		/************************************************************/





		xmldoc.Clear();

	}



	S1.azimuthSpacing = SumAzimuthSpacing / NumImg;
	S1.rangeSpacing = SumRangeSpacing / NumImg;


}

/*****************************************************************/
/*find S1 images and xml files                                   */
/*****************************************************************/
int GetFiles(string path, string * files, string polar)
{
	//Handle  
	intptr_t   hFile = 0;
	//文件信息  
	struct _finddata_t fileinfo;

	string p;
	int i = 0;
	if ((hFile = _findfirst(p.assign(path).append("\\*").c_str(), &fileinfo)) != -1)
	{




		do
		{
			//如果是目录,迭代之  
			//如果不是,加入列表  
			if ((fileinfo.attrib &  _A_SUBDIR))
			{

				continue;
				//if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
				//	getFiles(p.assign(path).append("\\").append(fileinfo.name), files);
			}
			else
			{
				string Temp;
				Temp.assign(path).append("\\").append(fileinfo.name);
				if (Temp.find(polar) != string::npos)
				{
					files[i] = Temp;
					i++;
				}


			}

		} while (_findnext(hFile, &fileinfo) == 0);

		_findclose(hFile);
	}
	return i;
}


/*****************************************************************/
/*Read orbit and fit orbit with polynomial                       */
/*****************************************************************/

double utc2seconds(const char* utc)
{
	int year, month, day;
	double hour, min, sec;
	sscanf(utc, "%d-%d-%dT%lf:%lf:%lf", &year, &month, &day, &hour, &min, &sec);
	return hour * 60 * 60 + min * 60 + sec;
}

void getIntArray(string str, char gap, int *Array)
{
	int lengths = str.length();
	string temp;
	int countt = 0;
	int i = 0;
	for (i = 0; i < lengths; i++)
	{
		temp += str[i];
		if (str[i] == gap)
		{
			Array[countt] = (atoi(temp.c_str()));
			countt++;
			temp.clear();
			continue;
		}
	}
	Array[countt] = (atoi(temp.c_str()));


}

int DiffUTCTime(string firstUTC, string lastUTC)
{
	int years1, months1, days1, hours1, mins1, secs1;
	sscanf_s(firstUTC.c_str(), "%d-%d-%dT%d:%d:%d", &years1, &months1,
		&days1, &hours1, &mins1, &secs1);

	int years11, months11, days11, hours11, mins11, secs11;
	sscanf_s(lastUTC.c_str(), "%d-%d-%dT%d:%d:%d", &years11, &months11,
		&days11, &hours11, &mins11, &secs11);

	int secInday = 60 * 60 * 24;
	int tt;
	if (years1 == years11&&months1 == months11&&days1 == days11)
		tt = hours11 * 3600 + mins11 * 60 + secs11 - (hours1 * 3600 + mins1 * 60 + secs1);
	else
		tt = hours11 * 3600 + mins11 * 60 + secs11 + secInday - (hours1 * 3600 + mins1 * 60 + secs1);

	return tt;
}

void getDoubleArray(string str, char gap, double *Array)
{
	int lengths = str.length();
	string temp;
	int countt = 0;
	int i = 0;
	for (i = 0; i < lengths; i++)
	{
		temp += str[i];
		if (str[i] == gap)
		{
			Array[countt] = (atof(temp.c_str()));
			countt++;
			temp.clear();
			continue;
		}
	}
	Array[countt] = (atof(temp.c_str()));


}



bool S1_OrbitInitialize(S1PreciseOrbit &Orbit,
	const  char * Path, const char *firstUTC, const char *lastUTC, string MissionID)
{

	CheckOribtFile(Path, firstUTC, MissionID);

	int years, months, days, hours, mins, secs;
	char buff[100];
	sscanf_s(firstUTC, "%d-%d-%dT%d:%d:%d", &years, &months,
		&days, &hours, &mins, &secs);

	sprintf(buff, "%04d-%02d-%02dT%02d:%02d:%01d", years, months, days, hours, mins, secs / 10);
	string firstUTCtime = buff;

	int tt = DiffUTCTime(firstUTC, lastUTC);

	int numOfOrbit = tt / 10 + 10;
	Orbit.NumPoints = numOfOrbit;

	XMLDocument xmldoc;

	bool LoadOK = xmldoc.LoadFile(Path);
	if (LoadOK)
	{
		cout << "The path is not valid! Please check it!" << endl;
		cout << Path << endl;
		system("pause");
		exit(0);
	}

	XMLElement* root = xmldoc.RootElement();

	XMLElement* elem = root->FirstChildElement("Data_Block");
	elem = elem->FirstChildElement("List_of_OSVs");

	int index = 0;
	XMLElement* elem1;

	elem1 = elem->FirstChildElement();


	while (true)
	{
		string str = elem1->FirstChildElement("UTC")->GetText();
		if (str.find(firstUTCtime.c_str()) != string::npos)
			break;

		elem1 = elem1->NextSiblingElement();
		index++;

		if (index>9361)
			break;
	}



	//Allocate Memory Space
	Orbit.orbitAzTime = new double[numOfOrbit];
	Orbit.x_pos = new double[numOfOrbit];
	Orbit.y_pos = new double[numOfOrbit];
	Orbit.z_pos = new double[numOfOrbit];
	Orbit.x_vel = new double[numOfOrbit];
	Orbit.y_vel = new double[numOfOrbit];
	Orbit.z_vel = new double[numOfOrbit];


	elem1 = elem->FirstChildElement();
	int num = 0;

	while (true)
	{
		if (num == index - 4)
		{
			for (int i = 0; i < numOfOrbit; i++)
			{
				string str = elem1->FirstChildElement("UTC")->GetText();
				str = str.substr(4, 30);
				Orbit.orbitAzTime[i] = (utc2seconds(str.c_str()));

				str = elem1->FirstChildElement("X")->GetText();
				Orbit.x_pos[i] = (atof(str.c_str()));

				str = elem1->FirstChildElement("Y")->GetText();
				Orbit.y_pos[i] = (atof(str.c_str()));

				str = elem1->FirstChildElement("Z")->GetText();
				Orbit.z_pos[i] = (atof(str.c_str()));

				str = elem1->FirstChildElement("VX")->GetText();
				Orbit.x_vel[i] = (atof(str.c_str()));

				str = elem1->FirstChildElement("VY")->GetText();
				Orbit.y_vel[i] = (atof(str.c_str()));

				str = elem1->FirstChildElement("VZ")->GetText();
				Orbit.z_vel[i] = (atof(str.c_str()));

				elem1 = elem1->NextSiblingElement();
			}

			break;
		}
		num++;
		elem1 = elem1->NextSiblingElement();

	}
	xmldoc.Clear();

	Orbit.TimeInerval = (Orbit.orbitAzTime[numOfOrbit - 1] - Orbit.orbitAzTime[0]) / (numOfOrbit - 1);



	return true;

	
}

time_t strTime2unix(const char*timeStamp)
{
	struct tm tm;
	memset(&tm, 0, sizeof(tm));

	sscanf(timeStamp, "%d-%d-%d %d:%d:%d",
		&tm.tm_year, &tm.tm_mon, &tm.tm_mday,
		&tm.tm_hour, &tm.tm_min, &tm.tm_sec);

	tm.tm_year -= 1900;
	tm.tm_mon--;

	return mktime(&tm);
}
time_t  addDay(time_t time1, int days)
{
	return (time1 + days * 60 * 60 * 24 + 28800);

}
time_t minusDay(time_t time1, int days)
{
	return (time1 - (days * 60 * 60 * 24 - 28800));

}

bool CheckOribtFile(string PreciseOrbit, string UTCtime, string MissionID)
{
	int monthDays[12] = { 31,28,31,30,31,30,31,31,30,31,30,31 };

	//Check Date
	char buff[100];
	char buff1[100];

	int years, months, days, hours, mins, secs;

	sscanf_s(UTCtime.c_str(), "%d-%d-%dT%d:%d:%d", &years, &months,
		&days, &hours, &mins, &secs);

	
	time_t t_unix = strTime2unix(UTCtime.c_str());
	
	
	time_t t_unix_before = minusDay(t_unix, 1);
	

	struct tm *p_before = gmtime(&t_unix_before);


	strftime(buff, sizeof(buff), "%Y%m%d", p_before);

	

	time_t t_unix_after = addDay(t_unix, 1);
	struct tm *p_after = gmtime(&t_unix_after);

	


	strftime(buff1, sizeof(buff1), "%Y%m%d", p_after);
	




	
		//sprintf(buff, "%04d%02d%02d", years, months, days - 1);
	
	
		//sprintf(buff, "%04d%02d%02d", years, months-1, monthDays[months - 1-1]);
	

	//sprintf(buff1, "%04d%02d%02d", years, months, days + 1);

	string str_start, str_end;
	int length = PreciseOrbit.length();

	str_start = PreciseOrbit.substr(length - 35, 8);
	str_end = PreciseOrbit.substr(length - 19, 8);

	


	


	if ((str_start != string(buff) || str_end != string(buff1))
		)
	{
		cout << "The date of file is not valid. Please check it." << endl;
		
		return false;
	}
	


	if (strstr(PreciseOrbit.c_str(), MissionID.c_str()) == NULL)
	{
		cout << "Mission is not mathched! Please check!." << endl;
		return false;
	}

	FILE *fp = fopen(PreciseOrbit.c_str(), "r");
	if (!fp) return -1;
	fseek(fp, 0L, SEEK_END);
	int size = ftell(fp);
	fclose(fp);

	if (size < 4096000)//should be bigger than 4000KBytes
	{
		cout << "The volume of data file is not correct! Please check!." << endl;
		return false;

	}


	return 1;
}


void matTxmat(double *A, double *B, double *Res,
	int ALines, int ResLines, int ResPixels)
{
	double sum = 0.0;
	for (int i = 0; i<ResLines; i++)
	{
		for (int j = 0; j<ResPixels; j++)
		{
			for (int k = 0; k<ALines; k++)
			{
				sum += A[k*ResLines + i] * B[k*ResPixels + j];
			}
			Res[i*ResPixels + j] = sum;
			sum = 0.0;
		}
	}

}

void choles(double * A, int N)
{

	double sum;
	for (int i = 0; i<N; ++i)
	{
		for (int j = i; j<N; ++j)
		{
			sum = A[i*N + j];
			for (int k = i - 1; k >= 0; --k)
			{
				sum -= A[i*N + k] * A[j*N + k];
			}
			if (i == j)
			{
				if (sum <= 0.) { cout << "choles: internal: A not pos. def." << endl; }
				A[i*N + i] = sqrt(sum);
			}
			else
			{
				A[j*N + i] = sum / A[i*N + i];
			}
		}
	}
}

void solvechol(
	double *A, double* B, int ALines, int APixels, int BLines, int BPixels)
{


	const int N = ALines;
	double sum;
	int i, j;

	for (i = 0; i<N; ++i)
	{
		sum = B[i*BPixels + 0];
		for (j = i - 1; j >= 0; --j)
		{
			sum -= A[i*APixels + j] * B[j*BPixels + 0];
		}
		B[i*BPixels + 0] = sum / A[i*APixels + i];
	}

	for (i = N - 1; i >= 0; --i)
	{
		sum = B[i*BPixels + 0];
		for (j = i + 1; j<N; ++j)
		{
			sum -= A[j*APixels + i] * B[j*BPixels + 0];
		}
		B[i*BPixels + 0] = sum / A[i*APixels + i];
	}
}

void invertchol(
	double* A, int N)
{

	double sum;
	int i, j, k;

	for (i = 0; i<N; ++i)
	{
		A[i*N + i] = 1. / A[i*N + i];
		for (j = i + 1; j<N; ++j)
		{
			sum = 0.;
			for (k = i; k<j; ++k)
			{
				sum -= A[j*N + k] * A[k*N + i];
			}
			A[j*N + i] = sum / A[j*N + j];
		}
	}

	for (i = 0; i<N; ++i)
	{
		for (j = i; j<N; ++j)
		{
			sum = 0.;
			for (k = j; k<N; ++k)
			{
				sum += A[k*N + i] * A[k*N + j];			// transpose
			}
			A[j*N + i] = sum;
		}
	}
}

bool Polyfit_Orbit(S1PreciseOrbit &Orbit)
{

	const int Npoints = Orbit.NumPoints;

	// Normalize tha azimuth time

	double * normal_t = new double[Npoints];
	for (int i = 0; i < Npoints; i++)
	{
		normal_t[i] = (Orbit.orbitAzTime[i] - Orbit.orbitAzTime[Npoints / 2]) / 10.0;
	}

	int Degree = (Npoints > 6) ? 5 : Npoints - 2;

	int Nunkown = Degree + 1;
	Orbit.NumCoeff = Nunkown;

	Orbit.coef_x = new double[Nunkown];
	Orbit.coef_y = new double[Nunkown];
	Orbit.coef_z = new double[Nunkown];

	double *A = new double[Npoints*Nunkown];
	double *N = new double[Nunkown*Nunkown];
	double *Qx_hat = new double[Nunkown*Nunkown];
	double *rhs_x = new double[Nunkown];
	double *rhs_y = new double[Nunkown];
	double *rhs_z = new double[Nunkown];


	for (int i = 0; i < Npoints; i++)
	{
		for (int j = 0; j < Nunkown; j++)
		{
			A[i*Nunkown + j] = pow(normal_t[i], j);
		}
	}


	matTxmat(A, A, N, Npoints, Nunkown, Nunkown);
	matTxmat(A, Orbit.x_pos, rhs_x, Npoints, Nunkown, 1);
	matTxmat(A, Orbit.y_pos, rhs_y, Npoints, Nunkown, 1);
	matTxmat(A, Orbit.z_pos, rhs_z, Npoints, Nunkown, 1);

	for (int i = 0; i < Nunkown*Nunkown; i++)
	{
		Qx_hat[i] = N[i];
	}

	choles(Qx_hat, Nunkown);
	solvechol(Qx_hat, rhs_x, Nunkown, Nunkown, Nunkown, 1);
	solvechol(Qx_hat, rhs_y, Nunkown, Nunkown, Nunkown, 1);
	solvechol(Qx_hat, rhs_z, Nunkown, Nunkown, Nunkown, 1);
	for (int i = 0; i < Nunkown; i++)
	{
		Orbit.coef_x[i] = rhs_x[i];
		Orbit.coef_y[i] = rhs_y[i];
		Orbit.coef_z[i] = rhs_z[i];

	}



	invertchol(Qx_hat, Nunkown);




	delete[] N;
	delete[] normal_t;
	delete[] A;
	delete[] rhs_x;
	delete[] rhs_y;
	delete[] rhs_z;
	delete[] Qx_hat;
	return true;
}


/*****************************************************************/
/*Compute the overlap area between master and slave images       */
/*****************************************************************/

double getLatitude(SentinelTOPS& S1, double aztime, double rgtime, int SubSwathID)
{

	int Index[4];
	double InterpolatedIndex[2];

	int SubSwathIndex = SubSwathID;
	int NumGeoPixels = S1.SubSwath[SubSwathIndex].numOfGeoPointsPerLine;
	Locating(S1, aztime, rgtime, SubSwathIndex, Index, InterpolatedIndex);

	double lat00 = S1.SubSwath[SubSwathIndex].latitude[Index[0] * NumGeoPixels + Index[2]];
	double lat01 = S1.SubSwath[SubSwathIndex].latitude[Index[0] * NumGeoPixels + Index[3]];
	double lat10 = S1.SubSwath[SubSwathIndex].latitude[Index[1] * NumGeoPixels + Index[2]];
	double lat11 = S1.SubSwath[SubSwathIndex].latitude[Index[1] * NumGeoPixels + Index[3]];

	int muX = InterpolatedIndex[0];
	int muY = InterpolatedIndex[1];

	return (1 - muY) * ((1 - muX) * lat00 + muX * lat01) +
		muY * ((1 - muX) * lat10 + muX * lat11);

}

double getLontitude(SentinelTOPS& S1, double aztime, double rgtime, int SubSwathID)
{

	int Index[4];
	double InterpolatedIndex[2];

	int SubSwathIndex = SubSwathID;
	int NumGeoPixels = S1.SubSwath[SubSwathIndex].numOfGeoPointsPerLine;
	Locating(S1, aztime, rgtime, SubSwathIndex, Index, InterpolatedIndex);

	double lon00 = S1.SubSwath[SubSwathIndex].longitude[Index[0] * NumGeoPixels + Index[2]];
	double lon01 = S1.SubSwath[SubSwathIndex].longitude[Index[0] * NumGeoPixels + Index[3]];
	double lon10 = S1.SubSwath[SubSwathIndex].longitude[Index[1] * NumGeoPixels + Index[2]];
	double lon11 = S1.SubSwath[SubSwathIndex].longitude[Index[1] * NumGeoPixels + Index[3]];

	int muX = InterpolatedIndex[0];
	int muY = InterpolatedIndex[1];

	return (1 - muY) * ((1 - muX) * lon00 + muX * lon01) +
		muY * ((1 - muX) * lon10 + muX * lon11);

}

void Locating(SentinelTOPS& S1, double aztime, double rgtime, int SubSwathID,
	int*Index, double * coef)
{

	int numOfGeoLines = S1.SubSwath[SubSwathID].numOfGeoLines;
	int numOfGeoPointsPerLine = S1.SubSwath[SubSwathID].numOfGeoPointsPerLine;
	int j0 = -1, j1 = -1;
	double muX = 0;
	if (rgtime < S1.SubSwath[SubSwathID].slantRangeTime[0])
	{
		j0 = 0;
		j1 = 1;
	}
	else if (rgtime >
		S1.SubSwath[SubSwathID].slantRangeTime[numOfGeoPointsPerLine - 1])
	{
		j0 = numOfGeoPointsPerLine - 2;
		j1 = numOfGeoPointsPerLine - 1;
	}
	else
	{
		for (int j = 0; j < numOfGeoPointsPerLine - 1; j++)
		{
			if (S1.SubSwath[SubSwathID].slantRangeTime[j] <= rgtime &&
				S1.SubSwath[SubSwathID].slantRangeTime[j + 1]> rgtime)
			{
				j0 = j;
				j1 = j + 1;
				break;
			}
		}
	}

	muX = (rgtime - S1.SubSwath[SubSwathID].slantRangeTime[j0]) /
		(S1.SubSwath[SubSwathID].slantRangeTime[j1] -
		S1.SubSwath[SubSwathID].slantRangeTime[j0]);

	int i0 = -1, i1 = -1;
	double muY = 0;
	for (int i = 0; i < S1.SubSwath[SubSwathID].numOfGeoLines - 1; i++) {
		double i0AzTime = (1 - muX) * S1.SubSwath[SubSwathID].azimuthTime[i*numOfGeoPointsPerLine + j0] +
			muX * S1.SubSwath[SubSwathID].azimuthTime[i*numOfGeoPointsPerLine + j1];

		double i1AzTime = (1 - muX) * S1.SubSwath[SubSwathID].azimuthTime[(i + 1)*numOfGeoPointsPerLine + j0] +
			muX * S1.SubSwath[SubSwathID].azimuthTime[(i + 1)*numOfGeoPointsPerLine + j1];

		if ((i == 0 && aztime < i0AzTime) ||
			(i == numOfGeoLines - 2 && aztime >= i1AzTime) ||
			(i0AzTime <= aztime && i1AzTime > aztime)) {

			i0 = i;
			i1 = i + 1;
			muY = (aztime - i0AzTime) / (i1AzTime - i0AzTime);
			break;
		}
	}

	Index[0] = i0;
	Index[1] = i1;
	Index[2] = j0;
	Index[3] = j1;
	coef[0] = muX;
	coef[1] = muY;



}

void Geodetic2Geographic(ellipsoid_WGS84 &ell, double lat, double lon, double height,
	double Res[3])
{
	double a = ell.a;
	double e2 = ell.e2;

	double N = a / sqrt(1.0 - e2*(sin(lat)*sin(lat)));
	double Nph = N + height;

	Res[0] = Nph * cos(lat) * cos(lon),
		Res[1] = Nph *  cos(lat) * sin(lon),
		Res[2] = (Nph - e2*N)    * sin(lat);

}

inline  void getXyz(double t, double* possat, double* coef_x,
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

inline  void getVelocity(double t, double* velsat, double* coef_x,
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

	velsat[0] /= (double)10.0;
	velsat[1] /= (double)10.0;
	velsat[2] /= (double)10.0;


}

inline void getAcceleration(double t, double*accsat, double* coef_x,
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

	accsat[0] /= (double) 100.0;
	accsat[1] /= (double)100.0;
	accsat[2] /= (double)100.0;



}

double xyz2aztime(double *xyz, int NumOfCoef,
	double *coef_x, double* coef_y, double* coef_z,
	double Initialaztime, int midAztime)
{
	double aztime = Initialaztime;
	double sol = 0.0;
	int iter;
	double possat[3];
	double velsat[3];
	double accsat[3];
	double dsat_P[3];
	for (iter = 0; iter <= 10; ++iter)           // break at convergence
	{
		//  Update equations 
		double normalized_Aztime = (aztime - midAztime) / 10.0;
		getXyz(normalized_Aztime, possat, coef_x, coef_y, coef_z, NumOfCoef);
		getVelocity(normalized_Aztime, velsat, coef_x, coef_y, coef_z, NumOfCoef);
		getAcceleration(normalized_Aztime, accsat, coef_x, coef_y, coef_z, NumOfCoef);
		dsat_P[0] = xyz[0] - possat[0];
		dsat_P[1] = xyz[1] - possat[1];
		dsat_P[2] = xyz[2] - possat[2];

		//  Update solution 

		sol = -(velsat[0] * dsat_P[0] + velsat[1] * dsat_P[1] + velsat[2] * dsat_P[2]) /
			((accsat[0] * dsat_P[0] + accsat[1] * dsat_P[1] + accsat[2] * dsat_P[2])
			- velsat[0] * velsat[0] - velsat[1] * velsat[1] - velsat[2] * velsat[2]);

		aztime += sol;

		// ______ Check convergence ______
		if (abs(sol) < 1e-010)                   // dta
			break;
	}

	return aztime;

}

/* 
The following codes are modified from Doris
add by ysli,201905 @ Newcastle Uni.
*/

double polyval1d(double x, double  *coeff, int num)
{
	double sum = 0.0;
	if (num > 100000)
	{
		return sum;
	}

	for (register long d = num - 1; d >= 0; --d)
	{
		sum *= x;
		sum += coeff[d];
	}
	return sum;
} // END polyval1d

int getxyz(double *satpos, S1PreciseOrbit  &orb, double t) // const // not, klo,khi
{


	int interp_method = 0;

	interp_method = (orb.NumPoints>6) ? 5 : orb.NumPoints - 2;
	 
	//const double t_tmp = (t - time[numberofpoints / 2]) / real8(10.0);     
	double t_tmp = (t - orb.orbitAzTime[orb.NumPoints / 2]) / 10.0;
 
		
		satpos[0] = polyval1d(t_tmp, orb.coef_x, interp_method + 1);
		satpos[1] = polyval1d(t_tmp, orb.coef_y, interp_method + 1);
		satpos[2] = polyval1d(t_tmp, orb.coef_z, interp_method + 1);
	 
	return 0;
} // END getxyz


int getxyzdot(double *velsat, S1PreciseOrbit  &orb, double t)
{
	
	 int interp_method = 0;

	  interp_method = (orb.NumPoints>6) ? 5 : orb.NumPoints - 2;


	  velsat[0] = orb.coef_x[1];// a_1*t^0 + 2a_2*t^1 + 3a_3*t^2 + ...
	  velsat[1] = orb.coef_y[1];
	  velsat[2] = orb.coef_z[1];
	   double t_tmp = (t - orb.orbitAzTime[orb.NumPoints / 2]) / 10.0;
		for (long i = 2; i <= interp_method; ++i)
		{
			double powt = double(i)*pow(t_tmp, double(i - 1));
			velsat[0] += orb.coef_x[i] * powt;
			velsat[1] += orb.coef_y[i] * powt;
			velsat[2] += orb.coef_z[i] * powt;
		}
		velsat[0] /= double(10.0);// normalization
		velsat[1] /= double(10.0);// normalization
		velsat[2] /= double(10.0);// normalization
	 
	return 0;
} // END getxyzdot


double  eq_in(double *inpos, double *P)            // scalar product: r=P.in(Q)
{
         // scalar product: r=P.in(Q)
	 return inpos[0] * P[0] + inpos[1] * P[1] + inpos[2] * P[2];
}

void solve33(
	double      RESULT[3],
	double rhs[3],
	double A[9])
{


	// ______  real8 L10, L20, L21: used lower matrix elements
	// ______  real8 U11, U12, U22: used upper matrix elements
	// ______  real8 b0,  b1,  b2:  used Ux=b
	const double L10 = A[3] / A[0];
	const double L20 = A[6] / A[0];
	const double U11 = A[4] - L10*A[1];
	const double L21 = (A[7] - (A[1] * L20)) / U11;
	const double U12 = A[5] - L10*A[2];
	const double U22 = A[8] - L20*A[2] - L21*U12;

	// ______ Solution: forward substitution ______
	const double b0 = rhs[0];
	const double b1 = rhs[1] - b0*L10;
	const double b2 = rhs[2] - b0*L20 - b1*L21;

	// ______ Solution: backwards substitution ______
	RESULT[2] = b2 / U22;
	RESULT[1] = (b1 - U12*RESULT[2]) / U11;
	RESULT[0] = (b0 - A[1] * RESULT[1] - A[2] * RESULT[2]) / A[0];
} // END solve33



long lp2xyz(
	double            line,
	double            pixel,
	const ellipsoid_WGS84 &ell,
	const SubSwathInfo  &image,
	S1PreciseOrbit       &orb,
	long  burstNum,
	double         *returnpos,
	long            MAXITER,
	double            CRITERPOS)
{
	// TRACE_FUNCTION("lp2xyz (BK 04-Jan-1999)")
	
	// ______ Convert lp to time ______
	const double aztime = image.burstFirstLineTime[burstNum] + (line - 1.0) / image.prf;       //image.line2ta(line);

	const double rangetime = image.slrTimeToFirstPixel + (pixel - 1.0) / (image.rangeSamplingRate*2);

	double possat[3];

	getxyz(&possat[0], orb, aztime); // const // not, klo,khi

	double velsat[3];

	getxyzdot(&velsat[0], orb, aztime); // const // not, klo,khi
 

	// ______ Set up system of equations and iterate ______
	//returnpos[0] = 0;  // iteration 0
	//returnpos[1] = 0;  // iteration 0
	//returnpos[2] = 0;  // iteration 0


	// ______ Save some alloc, init and dealloc time by declaring static (15%?) ______
	double solxyz[3];             // solution of system
	double equationset[3];        // observations
	double partialsxyz[3 * 3];        // partials to position due to nonlinear

	register long iter;
	for (iter = 0; iter <= MAXITER; ++iter)
	{


		double dsat_P[3];
		dsat_P[0] = returnpos[0] - possat[0];
		dsat_P[1] = returnpos[1] - possat[1];
		dsat_P[2] = returnpos[2] - possat[2];

		//equationset[0] = -eq1_doppler(velsat, dsat_P);
		equationset[0] = -eq_in(velsat, dsat_P);
 

		//equationset[1] = -eq2_range(dsat_P, ratime);
		equationset[1] = -1 *(eq_in(dsat_P, dsat_P) - sqr(SOL*rangetime));

       //equationset[2] = -eq3_ellipsoid(returnpos, ell.a, ell.b);

		equationset[2] = -1 * 	((sqr(returnpos[0]) + sqr(returnpos[1])) / sqr(ell.a) + sqr(returnpos[2] / ell.b) - 1.0);


		partialsxyz[0] = velsat[0];
		partialsxyz[1] = velsat[1];
		partialsxyz[2] = velsat[2];
		partialsxyz[3] = 2 * dsat_P[0];
		partialsxyz[4] = 2 * dsat_P[1];
		partialsxyz[5] = 2 * dsat_P[2];
		partialsxyz[6] = (2 * returnpos[0]) / sqr(ell.a);
		partialsxyz[7] = (2 * returnpos[1]) / sqr(ell.a);
		partialsxyz[8] = (2 * returnpos[2]) / sqr(ell.b);
	

		// ______ Solve system ______
		solve33(solxyz, equationset, partialsxyz);
	 
		// ______Update solution______
		returnpos[0] += solxyz[0];                         // update approx. value
		returnpos[1] += solxyz[1];                         // update approx. value
		returnpos[2] += solxyz[2];                         // update approx. value

		// ______Check convergence______
		if (abs(solxyz[0]) < CRITERPOS &&                         // dx
			abs(solxyz[1]) < CRITERPOS &&                         // dy
			abs(solxyz[2]) < CRITERPOS)                        // dz
			break; // converged

	
	}

	 
	return iter;
	
	return 0;
} // END lp2xyz

long lp2ell(
	double            line,
	double            pixel,
	const ellipsoid_WGS84 &ellips,
	const SubSwathInfo  &image,
	long  burstNum,
	S1PreciseOrbit &orb,
	double *returnpos,
	long            MAXITER = 10,
	double            CRITERPOS = 1e-6)
{

	long iter = lp2xyz(line, pixel, ellips, image, orb, burstNum,&returnpos[0], MAXITER, CRITERPOS);
	 
	return iter;
} // END lp2ell


// end of additional code//
 

 


bool ComputerBurstOverlap(int InputBurst0, int InputBurstN, int &M_burst0, int &M_burstN,
	int &S_burst0, int &S_burstN, SentinelTOPS &M_S1, SentinelTOPS &S_S1,
	int SubSwathIndex, S1PreciseOrbit& Morbit, S1PreciseOrbit& Sorbit)
{

	int Burst0 = InputBurst0 - 1;
	int BurstN = InputBurstN - 1;
	int SubSwathId = SubSwathIndex - 1;
	ellipsoid_WGS84 ell;

	if (Burst0 < 0)
	{
		cout << "Started number of burst is not valid!" << endl;
		system("pause");
		exit(0);
	}

	double aztime_burst0 = M_S1.SubSwath[SubSwathId].burstFirstLineTime[Burst0] +
		M_S1.SubSwath[SubSwathId].linesPerBurst / 2.0*M_S1.SubSwath[SubSwathId].azimuthTimeInterval;

	double rgtime_burst0 = (M_S1.SubSwath[SubSwathId].slrTimeToFirstPixel +
		M_S1.SubSwath[SubSwathId].slrTimeToFirstPixel) / 2.0;
	double lat, lon;
	lat = getLatitude(M_S1, aztime_burst0, rgtime_burst0, SubSwathId);
	lon = getLontitude(M_S1, aztime_burst0, rgtime_burst0, SubSwathId);

     double xyz[3];
 
	 Geodetic2Geographic(ell, deg2rad(lat), deg2rad(lon), 0.0, xyz); //Geo-location points is not precise
	
	 // use RD model to calculate the XYZ, add by ysli 20190517
	lp2ell(M_S1.SubSwath[SubSwathId].linesPerBurst / 2, M_S1.SubSwath[SubSwathId].samplesPerBurst/2, 
		ell, M_S1.SubSwath[SubSwathId], Burst0, Morbit, &xyz[0]);


	double zeroDopplerTime = xyz2aztime(xyz, Sorbit.NumCoeff, Sorbit.coef_x, Sorbit.coef_y, Sorbit.coef_z,
		aztime_burst0, Sorbit.orbitAzTime[Sorbit.NumPoints / 2]);



	if (zeroDopplerTime>S_S1.SubSwath[SubSwathId].lastLineTime)
	{
		cout << "No Overlap! Please check the coverage of the image pair, or reset the burst parameters!" << endl;
		system("pause");
		return false;
	}
	else if (zeroDopplerTime <S_S1.SubSwath[SubSwathId].lastLineTime
		&&zeroDopplerTime >S_S1.SubSwath[SubSwathId].firstLineTime)
	{
		//Locate the start burst of slave image
		M_burst0 = Burst0;

		int cc = 0;
		while (zeroDopplerTime>S_S1.SubSwath[SubSwathId].burstFirstLineTime[cc])
		{
			cc++;
		}
		S_burst0 = cc -1;




		if (BurstN - M_burst0 <= S_S1.SubSwath[SubSwathId].numOfBursts - 1 - S_burst0)
		{
			M_burstN = BurstN;
			S_burstN = S_burst0 + (M_burstN - M_burst0);
		}
		else
		{
			S_burstN = S_S1.SubSwath[SubSwathId].numOfBursts - 1;
			M_burstN = M_burst0 + (S_burstN - S_burst0);
		}

	}
	// Locate the start burst of master image from start burst of slave image

	else{


		aztime_burst0 = S_S1.SubSwath[SubSwathId].burstFirstLineTime[0] + S_S1.SubSwath[SubSwathId].linesPerBurst /
			2 * S_S1.SubSwath[SubSwathId].azimuthTimeInterval;
		rgtime_burst0 = (S_S1.SubSwath[SubSwathId].slrTimeToFirstPixel + S_S1.SubSwath[SubSwathId].slrTimeToLastPixel) / 2.0;

		lat = getLatitude(S_S1, aztime_burst0, rgtime_burst0, SubSwathId);
		lon = getLontitude(S_S1, aztime_burst0, rgtime_burst0, SubSwathId);

		double xyz[3];

		Geodetic2Geographic(ell, deg2rad(lat), deg2rad(lon), 0.0, xyz);

		lp2ell(S_S1.SubSwath[SubSwathId].linesPerBurst / 2, S_S1.SubSwath[SubSwathId].samplesPerBurst / 2,
			ell, S_S1.SubSwath[SubSwathId], 0, Sorbit, &xyz[0]);

		double zeroDopplerTime = xyz2aztime(xyz, Morbit.NumCoeff, Morbit.coef_x, Morbit.coef_y, Morbit.coef_z,
			aztime_burst0, Morbit.orbitAzTime[Morbit.NumPoints / 2]);



		if (zeroDopplerTime < M_S1.SubSwath[SubSwathId].firstLineTime)
		{
			cout << "No Overlap! Please check the coverage of the image pair, or reset the burst parameters!" << endl;
			system("pause");
			return false;
		}
		else if (zeroDopplerTime <M_S1.SubSwath[SubSwathId].lastLineTime
			&&zeroDopplerTime >M_S1.SubSwath[SubSwathId].firstLineTime)
		{
			//Locate the start burst of slave image
			S_burst0 = 0;

			int cc = 0;
			while (zeroDopplerTime > M_S1.SubSwath[SubSwathId].burstFirstLineTime[cc])
			{
				cc++;
			}
			M_burst0 = cc - 1;

			if (BurstN - M_burst0 <= S_S1.SubSwath[SubSwathId].numOfBursts - 1)
			{
				M_burstN = BurstN;
				S_burstN = S_burst0 + (M_burstN - M_burst0);
			}
			else
			{
				S_burstN = S_S1.SubSwath[SubSwathId].numOfBursts - 1;
				M_burstN = M_burst0 + (S_burstN - S_burst0);
			}


		}
	}





}


/*****************************************************************/
/*Check Input DEM                                                */
/*****************************************************************/

bool CheckSpecificDem(const char *dempath, double lat_min, double lat_max,
	double lon_min, double lon_max)
{
	bool res = 1;
	GDALAllRegister();
	//打开数据集
	GDALDatasetH hdem1;
	hdem1 = GDALOpen(dempath, GA_ReadOnly);

	if (hdem1 == NULL)
	{

		cout << " Unvalid SRTM file " << endl;
		cout << "It should be the format of Tiff!" << endl;
		GDALClose(hdem1);
		return 0;
	}


	//获取地理变换
	double adfGeoTransform[6];
	double deltaLat, deltaLon, demlonLeft, demlatUpper, demlonRight, demlatLow;
	if (GDALGetGeoTransform(hdem1, adfGeoTransform) == CE_None)
	{
		demlonLeft = adfGeoTransform[0];//pixel_left_up
		deltaLon = fabs(adfGeoTransform[1]);//pixel_interval
		//adfGeoTransform[2] = 0.0;
		demlatUpper = adfGeoTransform[3];//line_left_up
		//adfGeoTransform[4] = 0;
		deltaLat = fabs(adfGeoTransform[5]);//line_interval
	}

	int Ysize = GDALGetRasterYSize(hdem1);
	int Xsize = GDALGetRasterXSize(hdem1);

	demlatLow = demlatUpper - (Ysize - 1)*deltaLat;
	demlonRight = demlonLeft + (Xsize - 1)*deltaLon;

	if (demlatUpper <= lat_max)
	{
		cout << "Scene Outside DEM: most North latitude: " << demlatUpper << "[deg];Scene requires:" << lat_max <<
			"[deg]" << endl;

		res = 0;
	}

	if (demlatLow >= lat_min)
	{
		cout << "Scene Outside DEM: most South latitude: " << demlatLow << "[deg];Scene requires:" << lat_min <<
			"[deg]" << endl;
		res = 0;
	}

	if (demlonLeft >= lon_min)
	{
		cout << "Scene Outside DEM: most West longtitude: " << demlonLeft << "[deg];Scene requires:" << lon_min <<
			"[deg]" << endl;
		res = 0;
	}
	if (demlonRight <= lon_max)
	{
		cout << "Scene Outside DEM: most East longtitude: " << demlonRight << "[deg];Scene requires:" << lon_max <<
			"[deg]" << endl;
		res = 0;
	}

	GDALClose(hdem1);
	return res;

}


/*****************************************************************/
/*Compute deramping parameters                                   */
/*****************************************************************/
void computeDC(double AzTime, Polynomial* Input, int N, Polynomial &Res)
{
	int i0 = 0, i1 = 0;
	if (AzTime < Input[0].aztime)
	{
		i0 = 0;
		i1 = 1;
	}
	else if (AzTime > Input[N - 1].aztime)
	{
		i0 = N - 2;
		i1 = N;
	}
	else
	{
		for (int i = 0; i < N - 1; i++)
		{
			if (AzTime >= Input[i].aztime && AzTime < Input[i + 1].aztime)
			{
				i0 = i;
				i1 = i + 1;
				break;
			}
		}
	}



	Res.aztime = AzTime;
	Res.t0 = Input[i0].t0;
	double mu = (AzTime - Input[i0].aztime) / (Input[i1].aztime - Input[i0].aztime);

	Res.c0 = (1 - mu)* Input[i0].c0 + mu*Input[i1].c0;
	Res.c1 = (1 - mu)* Input[i0].c1 + mu*Input[i1].c1;
	Res.c2 = (1 - mu)* Input[i0].c2 + mu*Input[i1].c2;


}


bool ComputeDerampPara(SentinelTOPS &S1, S1PreciseOrbit &PreciseOrbit)
{
	Polynomial* DCEstimateList;
	for (int i = 0; i < S1.NumSubSwath; i++)
	{

		XMLDocument xmldoc;

		bool LoadOK = xmldoc.LoadFile(S1.SubSwath[i].XmlPath.c_str());
		if (LoadOK)
		{
			cout << "The path is not valid! Please check it!" << endl;
			cout << S1.SubSwath[i].XmlPath << endl;
			system("pause");
			exit(0);
		}

		XMLElement* root = xmldoc.RootElement();
		XMLElement* elem = root->FirstChildElement("imageAnnotation");
		elem = elem->FirstChildElement("processingInformation");
		elem = elem->FirstChildElement("dcMethod");
		string DCmethod = elem->GetText();

		elem = root->FirstChildElement("dopplerCentroid");
		elem = elem->FirstChildElement("dcEstimateList");
		string str = elem->FirstAttribute()->Value();
		int numbers = atoi(str.c_str());
		DCEstimateList = new Polynomial[numbers];

		int NDcPolynomial = atoi(elem->FirstChildElement()->FirstChildElement("dataDcPolynomial")->FirstAttribute()->Value());
		double *data = new double[NDcPolynomial + 5];//5 more space in case of emergencies

		XMLElement* elem1;
		elem1 = elem->FirstChildElement();
		for (int k = 0; k < numbers; k++)
		{
			str = elem1->FirstChildElement("azimuthTime")->GetText();
			DCEstimateList[k].aztime = utc2seconds(str.c_str());

			DCEstimateList[k].t0 = atof(elem1->FirstChildElement("t0")->GetText());

			if (DCmethod == "Data Analysis")
			{
				str = elem1->FirstChildElement("dataDcPolynomial")->GetText();
				int N = atoi(elem1->FirstChildElement("dataDcPolynomial")->FirstAttribute()->Value());

				getDoubleArray(str, ' ', data);
				DCEstimateList[k].c0 = data[0];
				DCEstimateList[k].c1 = data[1];
				DCEstimateList[k].c2 = data[2];

			}
			else
			{
				str = elem1->FirstChildElement("geometryDcPolynomial")->GetText();
				int N = atoi(elem1->FirstChildElement("dataDcPolynomial")->FirstAttribute()->Value());

				getDoubleArray(str, ' ', data);
				DCEstimateList[k].c0 = data[0];
				DCEstimateList[k].c1 = data[1];
				DCEstimateList[k].c2 = data[2];

			}
			elem1 = elem1->NextSiblingElement();
		}
		if (data)
		{
			delete[] data;
			data = NULL;
		}
		Polynomial *dcBurstList = new Polynomial[S1.SubSwath[i].numOfBursts];
		if (numbers > S1.SubSwath[i].numOfBursts)
		{
			for (int b = 0; b < S1.SubSwath[i].numOfBursts; b++)
			{
				dcBurstList[b] = DCEstimateList[b];
			}
		}
		else
		{
			for (int b = 0; b < S1.SubSwath[i].numOfBursts; b++)
			{
				double centreTime = 0.5*(S1.SubSwath[i].burstFirstLineTime[b] +
					S1.SubSwath[i].burstLastLineTime[b]);
				computeDC(centreTime, DCEstimateList, numbers, dcBurstList[b]);
			}

		}
		/******************Compute Doppler Centroid******************/

		S1.SubSwath[i].dopplerCentroid = new double[S1.SubSwath[i].numOfBursts*S1.SubSwath[i].samplesPerBurst];

		for (int b = 0; b < S1.SubSwath[i].numOfBursts; b++)
		{
			for (int x = 0; x < S1.SubSwath[i].samplesPerBurst; x++)
			{
				double slrt = (S1.SubSwath[i].slrTimeToFirstPixel +
					x*S1.SubSwath[i].rangePixelSpacing / SOL) * 2;//1-way to 2-way
				double dt = slrt - dcBurstList[b].t0;
				double dcvalue = dcBurstList[b].c0 + dcBurstList[b].c1*dt + dcBurstList[b].c2*sqr(dt);
				S1.SubSwath[i].dopplerCentroid[b*S1.SubSwath[i].samplesPerBurst + x] = dcvalue;
			}
		}
		/******************Compute Doppler Centroid******************/

		elem = root->FirstChildElement("generalAnnotation");
		elem = elem->FirstChildElement("azimuthFmRateList");
		str = elem->FirstAttribute()->Value();

		numbers = atoi(str.c_str());
		Polynomial* azFmRateList = new Polynomial[numbers];


		elem1 = elem->FirstChildElement();
		double *datas = NULL;
		if (elem1->FirstChildElement("azimuthFmRatePolynomial") != NULL)
		{
			int N = atoi(elem1->FirstChildElement("azimuthFmRatePolynomial")->FirstAttribute()->Value());

			datas = new double[N + 5];// 5 more space in case of emergencies


		}

		string str0;
		for (int k = 0; k < numbers; k++)
		{
			str = elem1->FirstChildElement("azimuthTime")->GetText();
			azFmRateList[k].aztime = utc2seconds(str.c_str());
			azFmRateList[k].t0 = atof(elem1->FirstChildElement("t0")->GetText());

			if (elem1->FirstChildElement("azimuthFmRatePolynomial") == NULL)
			{
				azFmRateList[k].c0 = atof(elem1->FirstChildElement("c0")->GetText());
				azFmRateList[k].c1 = atof(elem1->FirstChildElement("c1")->GetText());
				azFmRateList[k].c2 = atof(elem1->FirstChildElement("c2")->GetText());
			}
			else
			{
				str0 = elem1->FirstChildElement("azimuthFmRatePolynomial")->GetText();
				getDoubleArray(str0, ' ', datas);
				azFmRateList[k].c0 = datas[0];
				azFmRateList[k].c1 = datas[1];
				azFmRateList[k].c2 = datas[2];
			}
			elem1 = elem1->NextSiblingElement();
		}

		if (datas != NULL)delete[] datas;

		S1.SubSwath[i].rangeDependDopplerRate = new double[S1.SubSwath[i].numOfBursts*S1.SubSwath[i].samplesPerBurst];

		/******************Compute Range Depended Doppler Rate******************/
		for (int b = 0; b < S1.SubSwath[i].numOfBursts; b++)
		{

			for (int x = 0; x < S1.SubSwath[i].samplesPerBurst; x++)
			{
				double slrt = (S1.SubSwath[i].slrTimeToFirstPixel + x*S1.SubSwath[i].rangePixelSpacing / SOL) * 2;//1-way to 2-way
				double dt = slrt - azFmRateList[b].t0;
				double azdcvalue = azFmRateList[b].c0 + azFmRateList[b].c1*dt + azFmRateList[b].c2*dt*dt;
				S1.SubSwath[i].rangeDependDopplerRate[b*S1.SubSwath[i].samplesPerBurst + x] = azdcvalue;
			}
		}
		/******************Compute Range Depended Doppler Rate******************/

		/******************Compute Reference Time******************/
		S1.SubSwath[i].referenceTime = new double[S1.SubSwath[i].numOfBursts*S1.SubSwath[i].samplesPerBurst];

		double temp1 = S1.SubSwath[i].linesPerBurst*S1.SubSwath[i].azimuthTimeInterval / 2.0;

		int SamplesPerBurst = S1.SubSwath[i].samplesPerBurst;
		int NBursts = S1.SubSwath[i].numOfBursts;

		for (int b = 0; b < NBursts; b++)
		{


			double temp2 = temp1 + S1.SubSwath[i].dopplerCentroid[b*SamplesPerBurst + S1.SubSwath[i].firstValidPixel] /
				S1.SubSwath[i].rangeDependDopplerRate[b*SamplesPerBurst + S1.SubSwath[i].firstValidPixel];

			for (int x = 0; x < SamplesPerBurst; x++)
			{
				double reftime = temp2 - S1.SubSwath[i].dopplerCentroid[b*SamplesPerBurst + x] /
					S1.SubSwath[i].rangeDependDopplerRate[b*SamplesPerBurst + x];

				S1.SubSwath[i].referenceTime[b*SamplesPerBurst + x] = reftime;
			}
		}

		if (azFmRateList)
		{
			delete[] azFmRateList;
			azFmRateList = NULL;
		}
		if (dcBurstList)
		{
			delete[] dcBurstList;
			dcBurstList = NULL;
		}

		/*************************Compute Reference Time*************************/


		/****************** Compute Doppler Rate Kt(r) for each burst******************/
		double azTime = (S1.SubSwath[i].firstLineTime + S1.SubSwath[i].lastLineTime) / 2.0;
		azTime = (azTime - PreciseOrbit.orbitAzTime[PreciseOrbit.NumPoints / 2]) / 10.0;
		S1.SubSwath[i].dopplerRate = new double[NBursts*SamplesPerBurst];
		for (int b = 0; b < NBursts; b++)
		{
			double Velo[3];
			getVelocity(azTime, Velo, PreciseOrbit.coef_x, PreciseOrbit.coef_y, PreciseOrbit.coef_z,
				PreciseOrbit.NumCoeff);
			double v = abs3D(Velo);
			double steeringRate = deg2rad(S1.SubSwath[i].azimuthSteeringRate);
			double krot = 2 * v*steeringRate / S1.SubSwath[i].wavelength; // doppler rate by antenna steering

			for (int x = 0; x < SamplesPerBurst; x++)
			{
				double rate = S1.SubSwath[i].rangeDependDopplerRate[b*SamplesPerBurst + x] * krot
					/ (S1.SubSwath[i].rangeDependDopplerRate[b*SamplesPerBurst + x] - krot);
				S1.SubSwath[i].dopplerRate[b*SamplesPerBurst + x] = rate;
			}
		}

		/****************** Compute Doppler Rate Kt(r) for each burst******************/

		xmldoc.Clear();
	}

	return true;
}


/*****************************************************************/
/*Init a resample look-up table                                  */
/*****************************************************************/
bool InKernelInit(ResampleTable &ReTable, int Npoints, double prf, double abw, double rsr2x)
{



	const double CHI_az = prf / abw;
	const double CHI_rg = (rsr2x / 2.0) / abw;
	const double dx = 1.0 / Interval;
	const int Npointsd2 = Npoints / 2;
	const int Npointsd2m1 = Npointsd2 - 1;
	const int NInterval = Interval + 1;


	ReTable.KernelAz = new float[Npoints*NInterval];
	ReTable.KernelRg = new float[Npoints*NInterval];
	double Coeff_Az;
	double Coeff_Rg;
	//Using Raised Cosine Interpolation Kernels

	for (int i = 0; i < NInterval; i++)
	{

		double Sum_Az = 0.0;
		double Sum_Rg = 0.0;

		for (int j = 0; j < Npoints; j++)
		{
			double x_axis = 1.0 - Npointsd2 + j - dx*i;
			double v_az = 1. - 1. / CHI_az;
			double v_rg = 1. - 1. / CHI_rg;
			Coeff_Az = sinc(x_axis) * rect(x_axis / double(Npoints))*
				cos(v_az*PI*x_axis) / (1.0 - sqr(2.0*v_az*x_axis));
			Sum_Az += Coeff_Az;
			ReTable.KernelAz[i*Npoints + j] = Coeff_Az;

			Coeff_Rg = sinc(x_axis) * rect(x_axis / double(Npoints))*
				cos(v_rg*PI*x_axis) / (1.0 - sqr(2.0*v_rg*x_axis));
			Sum_Rg += Coeff_Rg;

			ReTable.KernelRg[i*Npoints + j] = Coeff_Rg;
		}

		for (int p = 0; p < Npoints; p++)
		{
			ReTable.KernelAz[i*Npoints + p] /= Sum_Az;
			ReTable.KernelRg[i*Npoints + p] /= Sum_Rg;
		}



	}

	return true;
}





/*************************************************************
*                Geometric Coregistration                    *
*************************************************************/
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

int GeometricCoreg(SubSwathInfo& M_Swath, SubSwathInfo& S_Swath, S1PreciseOrbit &Morbit, S1PreciseOrbit &Sorbit,
	RefDem& TopoDem, TransFormCoef TFCoef, int Mburst0, int MburstN, int Sburst0, int SburstN, int SwathId)
{

	int NumBurstProcessed = MburstN - Mburst0 + 1;
	int BurstGap = Sburst0 - Mburst0;
	int BurstPixels = M_Swath.samplesPerBurst;
	int BurstLines = M_Swath.linesPerBurst;

	//burst index need to be checked 
	for (int burstId = Mburst0; burstId <= MburstN; burstId++)
	{

		int xmin = 0;
		int ymin = burstId*BurstLines;
		int xmax = BurstPixels - 1;
		int ymax = (burstId + 1)*BurstLines - 1;


		double m_InitAztime = (M_Swath.burstFirstLineTime[burstId] + M_Swath.burstLastLineTime[burstId]) / 2.0;
		double s_InitAztime = (S_Swath.burstFirstLineTime[burstId] + S_Swath.burstLastLineTime[burstId]) / 2.0;


		double lat_max, lat_min, lon_max, lon_min;

		lat_min = M_Swath.latlonBurst[burstId].lat_min;
		lat_max = M_Swath.latlonBurst[burstId].lat_max;
		lon_min = M_Swath.latlonBurst[burstId].lon_min;
		lon_max = M_Swath.latlonBurst[burstId].lon_max;



		double extralat = 24 * TopoDem.deltaLat;
		double extralon = 24 * TopoDem.deltaLon;
		int DemLines, DemPixels;
		int* DemArray = NULL;
		TopoDem.getData(lat_min, lat_max, lon_min, lon_max, extralat, extralon,
			DemArray, DemLines, DemPixels);


		int InvalidCounters = 0;
		for (int i = 0; i < DemLines*DemPixels; i++)
		{

			if (DemArray[i] == demNoData)InvalidCounters++;

		}
		if (InvalidCounters > DemLines*DemPixels - 12)
		{
			cout << "Geometric Coregistration failed for this burst!" << endl;
		}



		

		GeoCoreg_warp_cuda(burstId, burstId + BurstGap, xmin, xmax, ymin, ymax, TopoDem.deltaLat,
			TopoDem.deltaLon, DemLines, DemPixels, TFCoef.getAzCoeff(burstId), TFCoef.getRgCoeff(burstId),
			lat_max, lon_min, demNoData, M_Swath.azimuthTimeInterval, S_Swath.azimuthTimeInterval, M_Swath.wavelength,
			M_Swath.linesPerBurst, S_Swath.linesPerBurst, M_Swath.burstFirstLineTime[burstId], S_Swath.burstFirstLineTime[burstId],
			M_Swath.samplesPerBurst, S_Swath.samplesPerBurst, Morbit.orbitAzTime, Sorbit.orbitAzTime, Morbit.NumPoints, Morbit.NumCoeff,
			Sorbit.NumPoints, Sorbit.NumCoeff, M_Swath.slrTimeToFirstPixel, S_Swath.slrTimeToFirstPixel, M_Swath.rangePixelSpacing,
			S_Swath.rangePixelSpacing, M_Swath.nearRangeOnLeft, S_Swath.nearRangeOnLeft, DemArray, Morbit.coef_x,
			Morbit.coef_y, Morbit.coef_z, Sorbit.coef_x, Sorbit.coef_y, Sorbit.coef_z, m_InitAztime, s_InitAztime);





		if (DemArray != NULL)
		{
			delete[] DemArray;
			DemArray = NULL;
		}
		cout << "Finish the Geometrical Coregistration on burst No." << burstId + 1 << " of SubSwath No."
			<< SwathId << "!" << endl;
	}




	return 0;
}



/*************************************************************
*                Deramping and Resmapling                    *
*************************************************************/
extern "C"   void DerampDemodResample(
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

int DerampAndResample(SubSwathInfo& M_Swath, SubSwathInfo& S_Swath, TransFormCoef& TFCoef,
	int Mburst0, int MburstN, int Sburst0, int SburstN, int SwathId, ResampleTable& ReTable,
	const char* OutPath, int RePoints)
{

	int MburstPixels = M_Swath.samplesPerBurst;
	int MburstLines = M_Swath.linesPerBurst;

	int SburstPixels = S_Swath.samplesPerBurst;
	int SburstLines = S_Swath.linesPerBurst;

	int burstOffset = Sburst0 - Mburst0;

	int OutLines = (MburstN - Mburst0 + 1)*MburstLines;
	TiffRead SlaveIn;
	TiffWrite ReSlaveOut;

	SlaveIn.Init(S_Swath.ImgPath.c_str());
	ReSlaveOut.Init(OutPath, GDT_CFloat32, MburstPixels, OutLines);

	complex<short>* Slave = new complex<short>[SburstLines*SburstPixels];//Allocate memory for slave input image
	complex<float>* ReSlave = new complex<float>[MburstLines*MburstPixels];// Allocate memory for resampled output image

	int StartLines = Mburst0*MburstLines;

	for (int burstId = Mburst0; burstId <= MburstN; burstId++)
	{

		int MasterBox[4];
		int SlaveBox[4];
		//Master range box
		MasterBox[0] = 0;
		MasterBox[1] = MburstPixels - 1;
		MasterBox[2] = burstId*MburstLines;
		MasterBox[3] = (burstId + 1)*MburstLines - 1;
		int mPixels = MasterBox[1] - MasterBox[0] + 1;
		int mLines = MasterBox[3] - MasterBox[2] + 1;

		//Slave range box
		SlaveBox[0] = 0;
		SlaveBox[1] = SburstPixels - 1;
		SlaveBox[2] = (burstId + burstOffset)*SburstLines;
		SlaveBox[3] = (burstId + burstOffset + 1)*SburstLines - 1;
		int sPixels = SlaveBox[1] - SlaveBox[0] + 1;
		int sLines = SlaveBox[3] - SlaveBox[2] + 1;
		int SburstId = burstId + burstOffset;

		//get transformation coefficients
		double CpmAz[6];
		double CpmRg[6];
		TFCoef.getBurstCoeff(burstId, CpmAz, CpmRg);


		SlaveIn.ReadCpxShort(SlaveBox[0], SlaveBox[2], sLines, sPixels, Slave);


		DerampDemodResample(Slave, CpmAz, CpmRg, 0.0, ReSlave, ReTable.KernelAz, ReTable.KernelRg, SburstId, sPixels,
			sLines, MasterBox, SlaveBox, SburstLines, SburstPixels, S_Swath.azimuthTimeInterval, S_Swath.dopplerRate,
			S_Swath.referenceTime, S_Swath.dopplerCentroid, RePoints);



		ReSlaveOut.WriteCpxFloat(MasterBox[0], MasterBox[2] - StartLines, mLines, mPixels, ReSlave);

		cout << "Finish the Deramping and Resampling on burst No." << burstId + 1 << " of SubSwath No."
			<< SwathId << "!" << endl;

	}
	SlaveIn.Close();
	ReSlaveOut.Close();
	if (Slave)
	{
		delete[] Slave;
		Slave = NULL;
	}
	if (ReSlave)
	{
		delete[] ReSlave;
	}
	return 0;
}



/*************************************************************
*1.Estimate azimuth coregistration residual.                 *
*2.Resmaple the slave image using the refined shift.         *
*3.Estimate coherence.                                       *
*************************************************************/
double GetSpectralSep(SubSwathInfo& M_Swath);

void cpOverlapSize(SubSwathInfo& M_Swath, int* OverlapSize, int Mburst0, int MburstN);

void EstOverlapSize(SubSwathInfo& M_Swath, int* OverlapSize, int Mburst0, int MburstN, int& TotalLines,
	int& MaximumLines);

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
);
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
);

int ESDAndCoh(SubSwathInfo& M_Swath, SubSwathInfo& S_Swath, TransFormCoef& TFCoef,
	int Mburst0, int MburstN, int Sburst0, int SburstN, int SwathId, ResampleTable& ReTable,
	const char *ReSlaveFile, const char* ReSlaveESDFile, const char*Cohout)
{

	int MburstPixels = M_Swath.samplesPerBurst;
	int MburstLines = M_Swath.linesPerBurst;

	int SburstPixels = S_Swath.samplesPerBurst;
	int SburstLines = S_Swath.linesPerBurst;
	int numBursts = MburstN - Mburst0 + 1;

	int burstOffset = Sburst0 - Mburst0;

	//Get overlap size of every burst
	int numOverlap = MburstN - Mburst0;
	int * OverlapSize = NULL;
	OverlapSize = new int[numOverlap];
	int TotalOverlapLines;
	int MaximumLines;
	EstOverlapSize(M_Swath, OverlapSize, Mburst0, MburstN, TotalOverlapLines, MaximumLines);



	//For Output File
	TiffWrite ReSlaveESD;
	TiffWrite CohRes;
	int hh = (MburstN - Mburst0 + 1)*MburstLines;
	int ww = MburstPixels;
	ReSlaveESD.Init(ReSlaveESDFile, GDT_CFloat32, ww, hh);
	CohRes.Init(Cohout, GDT_Float32, ww, hh);


	TiffRead MasterOverlap;
	TiffRead ReSlaveOverlap;
	MasterOverlap.Init(M_Swath.ImgPath.c_str());
	ReSlaveOverlap.Init(ReSlaveFile);




	//Output Azimuth Coregistration residual
	

	
	/*string logfiles = string(ReSlaveESDFile) + ".log";
	fstream fout;
	fout.open(logfiles.c_str(), ios::out);
	if (!fout.is_open())
	{
		cout << "output azoffset logfile error !\n";
		return 1;
	}*/
	char buff[100];

	double* BurstShifts = new double[numOverlap];



	double spectralSep = GetSpectralSep(M_Swath);




	//Read overlap area from master and resampled slave images
	int samplesPerBurst = M_Swath.samplesPerBurst;
	int linesPerBurst = M_Swath.linesPerBurst;

	complex<short>*masterFor = new complex<short>[TotalOverlapLines*samplesPerBurst];
	complex<short>*masterBack = new complex<short>[TotalOverlapLines*samplesPerBurst];

	complex<float>*ReslaveFor = new complex<float>[TotalOverlapLines*samplesPerBurst];
	complex<float>*ReslaveBack = new complex<float>[TotalOverlapLines*samplesPerBurst];


	MasterOverlap.ReadCpxShort(masterFor, masterBack, numOverlap, OverlapSize, linesPerBurst,
		0, Mburst0*linesPerBurst, samplesPerBurst);



	ReSlaveOverlap.ReadCpxFloat(ReslaveFor, ReslaveBack, numOverlap, OverlapSize, linesPerBurst,
		0, 0, samplesPerBurst);//y0 is alwasy zero.




	//Estimate ESD shifts
	int coh_rg = 24;
	int coh_az = 6;
	float cohThreshold = 0.3; //Coherence Threshold

	EstESDShifts(masterFor, masterBack, ReslaveFor, ReslaveBack, OverlapSize, BurstShifts,
		numOverlap, TotalOverlapLines, samplesPerBurst, coh_az, coh_rg, MaximumLines, cohThreshold,
		spectralSep, M_Swath.azimuthTimeInterval);

	if (masterFor)
	{
		delete[]masterFor;
		masterFor = NULL;
	}
	if (masterBack)
	{
		delete[]masterBack;
		masterBack = NULL;
	}
	if (ReslaveFor)
	{
		delete[]ReslaveFor;
		ReslaveFor = NULL;
	}
	if (ReslaveBack)
	{
		delete[]ReslaveBack;
		ReslaveBack = NULL;
	}

	for (int k = 0; k < numOverlap; k++)
	{
		cout << "Shift Residual:" << BurstShifts[k];
		if (fabs(BurstShifts[k])>100)
		{
			cout << "----Invalid Value, will be reset as zero!";
			BurstShifts[k] = 0;
		}
		cout << endl;
	}


	//remove the invalid shifts.





	//1.Resample the slave image using the refined shifts
	//2.Estimate Coherence

	TiffRead SlaveRaw;
	SlaveRaw.Init(S_Swath.ImgPath.c_str());


	cuComplex* d_ReSlave;
	complex<short>*SlaveIn = new complex<short>[S_Swath.samplesPerBurst*S_Swath.linesPerBurst];
	complex<float>*ReSlaveOut = new complex<float>[samplesPerBurst*linesPerBurst];
	complex<short>* MasterIn = new complex<short>[samplesPerBurst*linesPerBurst];
	float* Coh = new float[samplesPerBurst*linesPerBurst];

	int dx = (coh_rg - 1) / 2;
	int dy = (coh_az - 1) / 2;
	int startLines = Mburst0*linesPerBurst;
	for (int b = Mburst0; b <= MburstN; b++)
	{
		printf("Progress Info :  burst %d  \n", b + 1 );

		int y0 = b*linesPerBurst;

		//For first burst. Just loading it into GPU memory
		if (b - Mburst0 == 0)
		{

			ReSlaveOverlap.ReadCpxFloat(0, 0, linesPerBurst, samplesPerBurst,
				ReSlaveOut);




			d_ReSlave = ResampleFirstBurst(ReSlaveOut, samplesPerBurst, linesPerBurst);



			ReSlaveESD.WriteCpxFloat(0, 0, linesPerBurst, samplesPerBurst, ReSlaveOut);//y0 is always zero

			MasterOverlap.ReadCpxShort(0, startLines, linesPerBurst, samplesPerBurst, MasterIn);

			CoherenceEst_ESD(0, startLines, samplesPerBurst, linesPerBurst, dx, dy, coh_rg, coh_az, MasterIn,
				d_ReSlave, Coh);


			CohRes.WriteFloat(0, 0, linesPerBurst, samplesPerBurst, Coh);


		}
		//For other bursts
		else
		{

			int sPixels = S_Swath.samplesPerBurst;
			int sLines = S_Swath.linesPerBurst;

			int MasterBox[4] = { 0, samplesPerBurst - 1,
				b*linesPerBurst, (b + 1)*linesPerBurst - 1 };


			int SlaveBox[4] = { 0, sPixels - 1,
				(b + burstOffset)*sLines, (b + burstOffset + 1)*sLines - 1 };
			SlaveRaw.ReadCpxShort(SlaveBox[0], SlaveBox[2], sLines, sPixels,
				SlaveIn);

			//Read coregistration translation coefficients
			double CpmAz[6];
			double CpmRg[6];
			TFCoef.getBurstCoeff(b, CpmAz, CpmRg);



			//Resample
			d_ReSlave = DerampDemodResample_ESD(SlaveIn, CpmAz, CpmRg, BurstShifts[b - Mburst0 - 1], ReSlaveOut, ReTable.KernelAz,
				ReTable.KernelRg, b + burstOffset, sPixels, sLines, MasterBox, SlaveBox, linesPerBurst, samplesPerBurst,
				M_Swath.azimuthTimeInterval, S_Swath.dopplerRate, S_Swath.referenceTime, S_Swath.dopplerCentroid,
				12);


			//Read the resampled slave image after ESD correction
			ReSlaveESD.WriteCpxFloat(MasterBox[0], MasterBox[2] - startLines, linesPerBurst, samplesPerBurst, ReSlaveOut);


			//Read master image
			MasterOverlap.ReadCpxShort(MasterBox[0], MasterBox[2], linesPerBurst, samplesPerBurst, MasterIn);


			//Estimate Coherence
			CoherenceEst(0, y0, samplesPerBurst, linesPerBurst, dx, dy, coh_rg, coh_az, MasterIn,
				d_ReSlave, Coh);

			

			//Output Coherence image
			CohRes.WriteFloat(MasterBox[0], MasterBox[2] - startLines, linesPerBurst, samplesPerBurst, Coh);

		}


	}
	if (Coh)
	{
		delete[] Coh;
		Coh = NULL;
	}
	if (MasterIn)
	{
		delete[] MasterIn;
		MasterIn = NULL;
	}
	if (SlaveIn)
	{
		delete[] SlaveIn;
		SlaveIn = NULL;
	}
	if (ReSlaveOut)
	{
		delete[] ReSlaveOut;
		ReSlaveOut = NULL;
	}



	SlaveRaw.Close();
	//fout.close();
	MasterOverlap.Close();
	ReSlaveOverlap.Close();
	ReSlaveESD.Close();
	CohRes.Close();
	delete[] OverlapSize;
	delete[] BurstShifts;


	return 0;
}


int getLineIdxInMetaData(SubSwathInfo &mSubswath, double targetLineTime,int m_burst0, int m_burstN)
{
	int sy0 = -1;
	int sy1 = -1;
	int burstNum0 = -1;
	int burstNum1 = -1;
	double midTime = 0;
	int k = 0;
	for (int i = 0; i < mSubswath.numOfBursts; i++)
	{
		if (targetLineTime >= mSubswath.burstFirstLineTime[i] && targetLineTime < mSubswath.burstLastLineTime[i])
		{
			int sy = i * mSubswath.linesPerBurst +
				(int)(((targetLineTime - mSubswath.burstFirstLineTime[i]) / mSubswath.azimuthTimeInterval) + 0.5);

			if (k == 0)
			{
				sy0 = sy;
				burstNum0 = i;
			}
			else
			{
				sy1 = sy;
				burstNum1 = i;
				break;
			}
			++k;
		}
	}

	if (sy0 != -1 && sy1 != -1)
	{
		if (burstNum1 == m_burst0)
			return sy1;

		if (burstNum1 == m_burstN + 1)
			return sy0;

		// find time between bursts midTime
		// use first burst if targetLineTime is before midTime
		midTime = (mSubswath.burstLastLineTime[burstNum0] +
			mSubswath.burstFirstLineTime[burstNum1]) / 2.0;

		if (targetLineTime > midTime)
			return sy1;
	}

	if (targetLineTime < mSubswath.burstLastLineTime[mSubswath.numOfBursts - 1] + mSubswath.azimuthTimeInterval&&
		targetLineTime >= mSubswath.burstLastLineTime[mSubswath.numOfBursts - 1])
	{
		sy0 = mSubswath.numOfLines - 1;
	}

	return sy0;
}

int Stitch_OneSubswath(SubSwathInfo& M_Swath, const char* GappedImg, const char* StitchFile,int m_burst0, int m_burstN,
	int data_type)
{
	double targetFirstLineTime = M_Swath.burstFirstValidLineTime[m_burst0];
	double targetLastLineTime = M_Swath.burstLastValidLineTime[m_burstN];
	double targetLineTimeInterval = M_Swath.azimuthTimeInterval;

	double targetSlantRangeTimeToFirstPixel = M_Swath.slrTimeToFirstPixel;
	double targetSlantRangeTimeToLastPixel = M_Swath.slrTimeToLastPixel;
	double targetDeltaSlantRangeTime = M_Swath.rangePixelSpacing / SOL;


	int yMin = (targetFirstLineTime - M_Swath.firstLineTime) / targetLineTimeInterval;
	int yMax = (targetLastLineTime - M_Swath.firstLineTime) / targetLineTimeInterval;

	int xMin = 0;
	int xMax = M_Swath.samplesPerBurst - 1;

	int width = xMax - xMin + 1;
	int height = yMax - yMin + 1;

	TiffRead GappedIn;
	TiffWrite ReStitchOut;

	GappedIn.Init(GappedImg);
	
	ReStitchOut.Init(StitchFile, data_type, width, height);

	//only support three kinds of data
	complex<short> * data_CInt16 = new complex<short> [width];
	float* data_Float32 = new float[width];
	complex<float>* data_Cpxfloat32 = new complex<float>[width];


	for (int y = yMin; y <= yMax; y++)
	{
		int yy = y - yMin;
		printf("Stitching Progress Info : %.2f%%\r", yy*100.0 / height);


		double LineTime = targetFirstLineTime + yy * targetLineTimeInterval;
	
		int lineIndex = getLineIdxInMetaData(M_Swath, LineTime, m_burst0,m_burstN);

		if (lineIndex == -1)
		{
			printf("error occured in the function (Stitch_OneSubswath) when \
					trying to get line index in meta data \n");
			
			system("pause");
			exit(0);
			
		}

		
		if (data_type == GDT_CInt16)
		{
			lineIndex -= m_burst0 * M_Swath.linesPerBurst;
			
			memset(data_CInt16, 0, sizeof(complex<short>)*width);
			GappedIn.ReadCpxShort( xMin, lineIndex, 1,width, data_CInt16);

			ReStitchOut.WriteCpxShort(0, yy,1,width,data_CInt16);
		}

		else if (data_type == GDT_Float32)
		{
			lineIndex -= m_burst0 * M_Swath.linesPerBurst;

			memset(data_Float32, 0, sizeof(float)*width);
			GappedIn.ReadFloat(xMin, lineIndex, 1, width, data_Float32);

			ReStitchOut.WriteFloat(0, yy, 1, width, data_Float32);
		
		}

		else if (data_type == GDT_CFloat32)
		{
			lineIndex -= m_burst0 * M_Swath.linesPerBurst;

			memset(data_Cpxfloat32, 0, sizeof(complex<float>)*width);

			GappedIn.ReadCpxFloat(xMin, lineIndex, 1, width, data_Cpxfloat32);

			ReStitchOut.WriteCpxFloat(0, yy, 1, width, data_Cpxfloat32);
		}

		else
		{
			cout << "unsupported image type for the stitching step!\n";
			system("pause");
			exit(0);
		}

	}
	cout << endl;







	delete[] data_CInt16; data_CInt16 = NULL;
	delete[] data_Float32; data_Float32 = NULL;
	delete[] data_Cpxfloat32; data_Cpxfloat32 = NULL;

	GappedIn.Close();
	ReStitchOut.Close();

	return 1;
}

void EstOverlapSize(SubSwathInfo& M_Swath, int* OverlapSize, int Mburst0, int MburstN, int& TotalLines,
	int& MaximumLines)
{

	double dt = M_Swath.azimuthTimeInterval;
	if (OverlapSize == NULL)
	{
		cout << "The pointer has not been initialiazed! Please Intialize it!" << endl;
		system("pause");
		exit(0);
	}
	TotalLines = 0;
	MaximumLines = -9999;

	for (int i = Mburst0; i < MburstN; i++)
	{
		double endtime = M_Swath.burstLastLineTime[i];
		double starttime = M_Swath.burstFirstLineTime[i + 1];
		int OverlapLines = floor((endtime - starttime) / dt);
		OverlapSize[i - Mburst0] = OverlapLines;
		TotalLines += OverlapLines;
		if (OverlapLines > MaximumLines)MaximumLines = OverlapLines;

	}

}

double GetSpectralSep(SubSwathInfo& M_Swath)
{
	double tCycle = (M_Swath.linesPerBurst)*M_Swath.azimuthTimeInterval;

	double sumSpectralSep = 0.0;

	for (int i = 0; i < M_Swath.numOfBursts; i++)
	{
		for (int p = 0; p < M_Swath.samplesPerBurst; p++)
		{
			sumSpectralSep += M_Swath.dopplerRate[i*M_Swath.samplesPerBurst + p] * tCycle;
		}
	}
	return sumSpectralSep / M_Swath.samplesPerBurst / M_Swath.numOfBursts;
}



//Check if the directory exists. If not, create it.
bool CheckDir(string Dir)
{
	if (0 != access(Dir.c_str(), 0))
	{
		// if this folder not exist, create a new one.
		if (mkdir(Dir.c_str()) != 0)return false;
		//换成 ::_mkdir  ::_access 也行，不知道什么意思
	}

	return true;

}