# InSAR-on-Sentinel-1-TOPS-data

The software makes use of Graphic Processing Unit (GPUs) to carry out Geometric coregistration, resampling, ESD and Coherence Estimation on Sentinel-1 TOPS data.
Relevant algorithm details and speed test you can refer to the paper:
"GPU accelerated interferometric SAR processing for Sentinel-1 TOPS data" (2019), Computers and Geosciences, Doi: https://doi.org/10.1016/j.cageo.2019.04.010.
There is also a speed test reported at [here ](https://github.com/glemoine62/InSAR-on-Sentinel-1-TOPS-data-POSIX).


# Platform:
The code is developed under windows operating system .  
For a Linux user,  you may go to [the project by Guido Lemoine ](https://github.com/glemoine62/InSAR-on-Sentinel-1-TOPS-data-POSIX)  to check the linux version of this code.


# GPU hardware:
A GPU (with >=3.0 compute capability) is sufficient.
Double precision computing capability is appreciated to improve the accuracy of coregistration. 


# Basic Configuration File:

process_dir= Current Working Directory

masterpath=The Path to the Master Image

slavepath=The Path to the Slave Image

preciseOrbitMaster=The Path to the Master Precise Orbit File

preciseOrbitSlave=The Path to the Slave Precise Orbit File

SpecificDemPath=The Path to DEM File (Tiff format only)

burst0=The Start Burst

burstN=The Stop Burst

firstsubswath=The Start Proccesing Subswath

lastsubswath=The End Proccesing Subswath

polarisation=Specific polarisation types



# How to compile 

library dependency:

1.TinyXml2:

Please go to [here](https://github.com/leethomason/tinyxml2).


2.GDAL.


Compilation using VS and cuda toolkit:

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

# Contact:
filippoyu0717@gmail.com

