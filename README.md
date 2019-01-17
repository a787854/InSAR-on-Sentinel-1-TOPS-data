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

Our testing platform:
GPU: Nvidia GTX Titan black  CPU:Intel I7-6700K

