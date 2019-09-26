# GLORIA
Matlab code for "Hyperspectral Super-Resolution via Global-Local Low-Rank Matrix Estimation"

# Hybrid inexact Block Coordinate Descent (HiBCD)
Matlab code for "Hybrid Inexact BCD for Coupled Structured Matrix Factorization in Hyperspectral Super-Resolution", submitted to IEEE Transaction on Signal Processing, 2019.

## Usage
1. Run "demo_synthetic.m" to 

## Reference
Ruiyuan Wu, Wing-Kin Ma, Xiao Fu, and Qiang Li, "Hyperspectral Super-Resolution via Global-Local Low-Rank Matrix Estimation" [[pdf]](https://arxiv.org/pdf/1907.01149.pdf)

### Abstract
Hyperspectral super-resolution (HSR) is a problem that aims to estimate an image of high spectral and spatial resolutions from a pair of co-registered multispectral (MS) and hyperspectral (HS) images, which have coarser spectral and spatial resolutions, respectively. In this paper we pursue a low-rank matrix estimation approach for HSR. We assume that the spectral-spatial matrices associated with the whole image and the local areas of the image have low rank structures. The local low-rank assumption, in particular, has the aim of providing a more flexible model for accounting for local variation effects due to endmember variability. We formulate the HSR problem as a global-local rank-regularized least-squares problem. By leveraging on the recent advances in non-convex large-scale optimization, namely, the smooth Schatten-p approximation and the accelerated majorization-minimization method, we developed an efficient algorithm for the global-local low-rank problem. Numerical experiments on synthetic and semi-real data show that the proposed algorithm outperforms a number of benchmark algorithms in terms of recovery performance.

## Sources of the tested algorithms
1. Hyperspectral and multispectral data fusion toolbox/CNMF: http://naotoyokoya.com/Download.html

2. CNMF: http://naotoyokoya.com/Download.html


### Please advise to remove immediately if any infringement caused.
