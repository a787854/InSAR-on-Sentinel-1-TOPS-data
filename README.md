# InSAR-on-Sentinel-1-TOPS-data
This code is developed under Visual Studio 2013 and CUDA toolkit 8.0.
The software makes use of Graphic Processing Unit to perform Geometric coregistration, resampling, ESD and Coherence Estimation on Sentinel-1 TOPS data.

This software dependents on two libraries:

1.TinyXml2:

Please find it at:https://github.com/leethomason/tinyxml2

2.GDAL.


Platform:

It is suggested to use the GPUs with powerful double floating computing.
Otherwise, the accuracy of coregistration is not guaranteed. 
To achieve a better performance, the driver mode of GPU should be set as TCC mode. 


Our testing platform:
GPU: Nvidia GTX Titan black  CPU:Intel I7-6700K


Config File:

process_dir= Current Working Directory

masterpath=The Folder of Master Image

slavepath=The Folder of Slave Image

preciseOrbitMaster=The Path to the Master Precise Orbit File

preciseOrbitSlave=The Path to the Slave Precise Orbit File

SpecificDemPath=The Path to DEM File (Tiff format only)

burst0=The Start Burst

burstN=The Stop Burst

firstsubswath=The Start Proccesing Subswath

lastsubswath=The End Proccesing Subswath

polarisation=Specific polarisation types

Research Paper:
"GPU accelerated interferometric SAR processing for Sentinel-1 TOPS data" (2019), Computers and Geosciences, Doi: https://doi.org/10.1016/j.cageo.2019.04.010.




How to compile this file:

Regarding my environment, I used the VS2017 as the IDE and CUDA 10.1 as the GPU toolkit to compile the all codes (VS2013 and CUDA 8.0 also works). 

Here are the steps which should be noted in compilation:

1. Create a new empty win console application project (without pre-compiled header).

2. add all the source code files and two additional files tinyxml2.cpp and tinyxml2.h (from library TinyXml2)

3. Link the CUDA toolkit to the project: right click the project-> build customization->click CUDA toolkit.

4. Specify all the .cu files as the CUDA C/C++ files, so these files can be included into the compilation.

5. Set the include directories (include TinyXml2 and GDAL libraries), and lib directories (include GDAL libs).

6. link to the GDAL libraries and also to CUDA libraries (cudart.lib; cusolver.lib; cublas.lib).

7. Project -> Properties -> Configuration Properties -> CUDA C/C++ -> Device -> Code Generation -> compute_35, sm_35 (higher than 35, and choose an appropriate option according to your card’s compute capability);

8. Note : "This function or variable may be unsafe. Consider using sprintf_s instead", please include this two lines in the Pre-processor Definitions : _CRT_SECURE_NO_WARNINGS , _CRT_NONSTDC_NO_DEPRECATE.

If you meet some specific problems, please contact me!

Contact:
filippoyu0717@gmail.com

