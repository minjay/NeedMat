## NeedMat, a Matlab Package for Spherical Needlets
## Description
NeedMat provides a Matlab Package that implements fast spherical needlet transforms and fast spherical needlet evaluations. For more details on the spherical needlets and the package, please refer to the technical report:
**A Note on Spherical Needlets**
http://arxiv.org/abs/1508.05406.
## Dependencies
NeedMat is developed on Matlab R2013a. It depends on the following three packages.
* RBFSPHERE Matlab package (http://math.boisestate.edu/~wright/montestigliano/index.html).
According to the author of RBFSPHERE, Grady Wright, the functions are integrated into the package SpherePts (https://github.com/gradywright/spherepts).
* MEALPix 3.0 (the original link is broken), which is a Matlab implementation of HEALPix based on HEALPix-F90 original source code.
* Spherical-Harmonic-Transform (https://github.com/polarch/Spherical-Harmonic-Transform), which is a collection of MATLAB routines for the Spherical Harmonic Transform and related manipulations in the spherical harmonic spectrum. Note that NeedMat only uses this package to compute Voronoi diagrams on the sphere.

To ensure compatibility, the second package is already included in the repository.

## Installation
1. Download the codes by cloning the repository `git clone https://github.com/minjay/NeedMat`.
2. Add the repository (including subfolders) to the search path of Matlab.

## Main Functions (in the folder named core)
1. `A = get_A(B, j_min, j_max, theta, phi, n_dist)`
Compute the design matrix A, where A is an N-by-M matrix, N is the number of observations, and M is the number of spherical needlets. `(theta, phi)` gives the locations of these observations, and the spherical needlets are from frequency level `j_min` to `j_max` (inclusively). For fast computation, this function first evaluates the spherical needlets on a very fine grid and then interpolates the values for the query points.
2. `Nside = get_Nside(B, j)` 
Compute `Nside`.
3. `j_max = get_j_max(B, l_max)` 
Compute the maximal j, `j_max`.
4. `[dist, psi] = get_psi(B, j, k, theta, phi)` 
A wrapper of the function `spneedlet_eval_fast`. It evaluates the spherical needlet with subscripts `j` and `k` at locations `(theta, phi)`.
5. `map = inv_spharmonic_tran_naive(alm, theta, phi, l_max)`
A *naive* implementation of the inverse spherical harmonic transform.
6. `plot_needlets(B, j, k, res)`
Plot the spherical needlet with subscripts `j` and `k`.
7. `Y = spharmonic_eval(l, m, theta, phi)`
Evaluate the spherical harmonic with subscripts `l` and `m` at locations `(theta, phi)`.
8. `alm = spharmonic_tran_irr(theta, phi, f, l_max)`
Spherical harmonic transform for irregularly spaced observations of function `f` at locations `(theta, phi)`. It estimates the spherical harmonic coefficients empirically using the weights determined by the Voronoi diagram of these locations on the sphere.
9. `beta = spneedlet_tran(alm, l_max, B)`
Fast spherical needlet transform. It computes the needlet coefficients based on the spherical harmonic coefficients `alm`.

## Contact
Please report any bugs to 
mjfan@ucdavis.edu

My research: [Google Scholar](https://scholar.google.com/citations?user=UlavjOkAAAAJ&hl=en)
