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

Contact:
filippoyu0717@gmail.com
