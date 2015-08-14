## NeedMat, a Matlab Package for Spherical Needlets
### Description
NeedMat provides a Matlab Package that implements fast spherical needlet transforms and fast spherical needlet evaluations. 
### Dependencies
NeedMat is developed on Matlab R2013a. It depends on the following two packages.
* MEALPix 3.0 (The original link is broken), which is a Matlab implementation of HEALPix based on HEALPix-F90 original source code.
* Spherical-Harmonic-Transform (https://github.com/polarch/Spherical-Harmonic-Transform), which is a collection of MATLAB routines for the Spherical Harmonic Transform and related manipulations in the spherical harmonic spectrum. Note that NeedMat only uses this package to compute Voronoi diagrams on the sphere.

### Installation
1. Download the codes by cloning the repository `git clone https://github.com/minjay/NeedMat`.
2. Add the repository (including subfolders) to the search path of Matlab.
