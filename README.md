## NeedMat, a Matlab Package for Spherical Needlets
### Description
NeedMat provides a Matlab Package that implements fast spherical needlet transforms and fast spherical needlet evaluations. 
### Dependencies
NeedMat is developed on Matlab R2013a. It depends on the following two packages.
* MEALPix 3.0 (the original link is broken), which is a Matlab implementation of HEALPix based on HEALPix-F90 original source code.
* Spherical-Harmonic-Transform (https://github.com/polarch/Spherical-Harmonic-Transform), which is a collection of MATLAB routines for the Spherical Harmonic Transform and related manipulations in the spherical harmonic spectrum. Note that NeedMat only uses this package to compute Voronoi diagrams on the sphere.

To ensure compatibility, these two packages are already included in the repository.

### Installation
1. Download the codes by cloning the repository `git clone https://github.com/minjay/NeedMat`.
2. Add the repository (including subfolders) to the search path of Matlab.

### Functions (in the folder named core)
1. `A = get_A(B, j_min, j_max, theta, phi, n_dist)`
Compute the design matrix A, where A is an N-by-M matrix, N is the number of observations, and M is the number of spherical needlets. (theta, phi) gives the locations of these observations, and the spherical needlets are from frequency level j_min to j_max (inclusively). For fast computation, this function first evaluate the spherical needlets at a very fine grid and then interpolate the values to the query points.
2. 
