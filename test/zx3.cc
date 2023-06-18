//
// Created by louis on 2021/8/29.
//


#include <iostream>
#include <mpi.h>
#include "target/hip/hipTarget.h"
#include <vector> 

// hip header file
#include "hip/hip_runtime.h"
#include "hip/hip_complex.h"
#include "operator.h"


__global__ void transfer_x_f(double *src, double *dest_b, double *U,
                             const int v_x, const int v_y, const int v_z, const int v_t,
                             const int s_x, const int s_y, const int s_z, const int s_t,
                             const int rank, const int cb, const int flag) {
/*
    int y = blockDim.x * blockIdx.x + threadIdx.x;
    int z = blockDim.y * blockIdx.y + threadIdx.y;
    int tp = blockDim.z * blockIdx.z + threadIdx.z;
*/

//    if (y >= s_y || z >= s_z) { return; }

    int subgrid_vol = (s_x * s_y * s_z * s_t);
    int subgrid_vol_cb = (subgrid_vol) >> 1;

    int N_sub[4] = {
            v_x / s_x,
            v_y / s_y,
            v_z / s_z,
            v_t / s_t
    };

    const double half = 0.5;
//    const hipDoubleComplex I(0, 1);

    const int x_p = ((rank / N_sub[0]) % N_sub[1]) * s_y +
                    ((rank / (N_sub[1] * N_sub[0])) % N_sub[2]) * s_z +
                    (rank / (N_sub[2] * N_sub[1] * N_sub[0])) * s_t;

    const int s_x_cb = s_x >> 1;

    int xyz = blockDim.x * blockIdx.x + threadIdx.x;
    
    int tp = xyz / (s_z * s_y);
    int z = (xyz / s_y) % s_z;
    int y = xyz % s_y;	

    int t = (y + z + 2 * tp + x_p) % 2 == cb ? 2 * tp + 1 : 2 * tp;
    if (t >= s_t) { return; }

    int cont = s_y * s_z * tp + s_y * z + y;

    int x = 0;

    double tmp[2];

    double *srcO = src + (s_x_cb * s_y * s_z * t +
                          s_x_cb * s_y * z +
                          s_x_cb * y +
                          x + (1 - cb) * subgrid_vol_cb) * 12 * 2;

//    hipDoubleComplex tmp;
/*    
    hipDoubleComplex *srcO = (hipDoubleComplex *) (src + (s_x_cb * s_y * s_z * t +
                                                          s_x_cb * s_y * z +
                                                          s_x_cb * y + x +
                                                          (1 - cb) * subgrid_vol_cb) * 12 * 2);
*/
    for (int c2 = 0; c2 < 3; c2++) {
/*
        tmp = -(srcO[0 * 3 + c2] - flag * I * srcO[3 * 3 + c2]) * half;
        dest_b[cont * 6 * 2 + (0 * 3 + c2) * 2 + 0] = tmp.x;
        dest_b[cont * 6 * 2 + (0 * 3 + c2) * 2 + 1] = tmp.y;
        tmp = -(srcO[1 * 3 + c2] - flag * I * srcO[2 * 3 + c2]) * half;
        dest_b[cont * 6 * 2 + (1 * 3 + c2) * 2 + 0] = tmp.x;
        dest_b[cont * 6 * 2 + (1 * 3 + c2) * 2 + 1] = tmp.y;
*/
	tmp[0] = -(srcO[(0 * 3 + c2) * 2 + 0] + flag * srcO[(3 * 3 + c2) * 2 + 1]) * half;
	tmp[1] = -(srcO[(0 * 3 + c2) * 2 + 1] - flag * srcO[(3 * 3 + c2) * 2 + 0]) * half;	
	dest_b[cont * 6 * 2 + (0 * 3 + c2) * 2 + 0] = tmp[0];
	dest_b[cont * 6 * 2 + (0 * 3 + c2) * 2 + 1] = tmp[1];	
	tmp[0] = -(srcO[(1 * 3 + c2) * 2 + 0] + flag * srcO[(2 * 3 + c2) * 2 + 1]) * half;
	tmp[1] = -(srcO[(1 * 3 + c2) * 2 + 1] - flag * srcO[(2 * 3 + c2) * 2 + 0]) * half;
	dest_b[cont * 6 * 2 + (1 * 3 + c2) * 2 + 0] = tmp[0];
	dest_b[cont * 6 * 2 + (1 * 3 + c2) * 2 + 1] = tmp[1];

    }
/*
    tmp = -(srcO[0 * 3 + 0] - flag * I * srcO[3 * 3 + 0]) * half;
    dest_b[cont * 6 * 2 + (0 * 3 + 0) * 2 + 0] = tmp.x;
    dest_b[cont * 6 * 2 + (0 * 3 + 0) * 2 + 1] = tmp.y;
    tmp = -(srcO[1 * 3 + 0] - flag * I * srcO[2 * 3 + 0]) * half;
    dest_b[cont * 6 * 2 + (1 * 3 + 0) * 2 + 0] = tmp.x;
    dest_b[cont * 6 * 2 + (1 * 3 + 0) * 2 + 1] = tmp.y;

    tmp = -(srcO[0 * 3 + 1] - flag * I * srcO[3 * 3 + 1]) * half;
    dest_b[cont * 6 * 2 + (0 * 3 + 1) * 2 + 0] = tmp.x;
    dest_b[cont * 6 * 2 + (0 * 3 + 1) * 2 + 1] = tmp.y;
    tmp = -(srcO[1 * 3 + 1] - flag * I * srcO[2 * 3 + 1]) * half;
    dest_b[cont * 6 * 2 + (1 * 3 + 1) * 2 + 0] = tmp.x;
    dest_b[cont * 6 * 2 + (1 * 3 + 1) * 2 + 1] = tmp.y;

    tmp = -(srcO[0 * 3 + 2] - flag * I * srcO[3 * 3 + 2]) * half;
    dest_b[cont * 6 * 2 + (0 * 3 + 2) * 2 + 0] = tmp.x;
    dest_b[cont * 6 * 2 + (0 * 3 + 2) * 2 + 1] = tmp.y;
    tmp = -(srcO[1 * 3 + 2] - flag * I * srcO[2 * 3 + 2]) * half;
    dest_b[cont * 6 * 2 + (1 * 3 + 2) * 2 + 0] = tmp.x;
    dest_b[cont * 6 * 2 + (1 * 3 + 2) * 2 + 1] = tmp.y;
*/
}

__global__ void transfer_x_b(double *src, double *dest_f, double *U,
                             const int v_x, const int v_y, const int v_z, const int v_t,
                             const int s_x, const int s_y, const int s_z, const int s_t,
                             const int rank, const int cb, const int flag) {
/*
    int y = blockDim.x * blockIdx.x + threadIdx.x;
    int z = blockDim.y * blockIdx.y + threadIdx.y;
    int tp = blockDim.z * blockIdx.z + threadIdx.z;
*/
  
//  if (y >= s_y || z >= s_z) { return; }

    int subgrid_vol = (s_x * s_y * s_z * s_t);
    int subgrid_vol_cb = (subgrid_vol) >> 1;

    int N_sub[4] = {v_x / s_x, v_y / s_y, v_z / s_z, v_t / s_t};

    const double half = 0.5;
//    const hipDoubleComplex I(0, 1);

    const int x_p = ((rank / N_sub[0]) % N_sub[1]) * s_y +
                    ((rank / (N_sub[1] * N_sub[0])) % N_sub[2]) * s_z +
                    (rank / (N_sub[2] * N_sub[1] * N_sub[0])) * s_t;

    const int s_x_cb = s_x >> 1;

    int xyz = blockDim.x * blockIdx.x + threadIdx.x;	

    int tp = xyz / (s_z * s_y);
    int z = (xyz / s_y) % s_z;
    int y = xyz % s_y;

    int t = (y + z + 2 * tp + x_p) % 2 != cb ? 2 * tp + 1 : 2 * tp;
    if (t >= s_t) { return; }

    int cont = s_y * s_z * tp + s_y * z + y;

    int x = s_x_cb - 1;

/*
    hipDoubleComplex tmp;
    hipDoubleComplex *AO = (hipDoubleComplex *) (U + (s_x_cb * s_y * s_z * t +
                                                      s_x_cb * s_y * z +
                                                      s_x_cb * y +
                                                      x + (1 - cb) * subgrid_vol_cb) * 9 * 2);

    hipDoubleComplex *srcO = (hipDoubleComplex *) (src + (s_x_cb * s_y * s_z * t +
                                                          s_x_cb * s_y * z +
                                                          s_x_cb * y + x +
                                                          (1 - cb) * subgrid_vol_cb) * 12 * 2);
*/

    double tmp[2];
    double destE[12];
    double srcO[24];
    double AO[18];

    for (int i = 0; i < 24; i++) {

        srcO[i] = src[(s_x_cb * s_y * s_z * t +
                       s_x_cb * s_y * z +
                       s_x_cb * y +
                       x + (1 - cb) * subgrid_vol_cb) * 12 * 2 + i];

    }

    for (int i = 0; i < 12; i++) {

        destE[i] = 0;

    }

    for (int i = 0; i < 18; i++) {

        AO[i] = U[(s_x_cb * s_y * s_z * t +
                   s_x_cb * s_y * z +
                   s_x_cb * y +
                   x + (1 - cb) * subgrid_vol_cb) * 9 * 2 + i];

    }	

    for (int c1 = 0; c1 < 3; c1++) {
        for (int c2 = 0; c2 < 3; c2++) {
/*
            tmp = -(srcO[0 * 3 + c2] + flag * I * srcO[3 * 3 + c2]) * half * hipConj(AO[c2 * 3 + c1]);

            dest_f[cont * 6 * 2 + (0 * 3 + c1) * 2 + 0] += tmp.x;
            dest_f[cont * 6 * 2 + (0 * 3 + c1) * 2 + 1] += tmp.y;

            tmp = -(srcO[1 * 3 + c2] + flag * I * srcO[2 * 3 + c2]) * half * hipConj(AO[c2 * 3 + c1]);

            dest_f[cont * 6 * 2 + (1 * 3 + c1) * 2 + 0] += tmp.x;
            dest_f[cont * 6 * 2 + (1 * 3 + c1) * 2 + 1] += tmp.y;
*/
	    tmp[0] = -(srcO[(0 * 3 + c2) * 2 + 0] - flag * srcO[(3 * 3 + c2) * 2 + 1]) * half * AO[(c2 * 3 + c1) * 2 + 0]
	 	     -(srcO[(0 * 3 + c2) * 2 + 1] + flag * srcO[(3 * 3 + c2) * 2 + 0]) * half * AO[(c2 * 3 + c1) * 2 + 1];

	    tmp[1] = +(srcO[(0 * 3 + c2) * 2 + 0] - flag * srcO[(3 * 3 + c2) * 2 + 1]) * half * AO[(c2 * 3 + c1) * 2 + 1]
                     -(srcO[(0 * 3 + c2) * 2 + 1] + flag * srcO[(3 * 3 + c2) * 2 + 0]) * half * AO[(c2 * 3 + c1) * 2 + 0];
	    	
            destE[(0 * 3 + c1) * 2 + 0] += tmp[0];
            destE[(0 * 3 + c1) * 2 + 1] += tmp[1];

	    tmp[0] = -(srcO[(1 * 3 + c2) * 2 + 0] - flag * srcO[(2 * 3 + c2) * 2 + 1]) * half * AO[(c2 * 3 + c1) * 2 + 0]
                     -(srcO[(1 * 3 + c2) * 2 + 1] + flag * srcO[(2 * 3 + c2) * 2 + 0]) * half * AO[(c2 * 3 + c1) * 2 + 1];		
	   
	    tmp[1] = +(srcO[(1 * 3 + c2) * 2 + 0] - flag * srcO[(2 * 3 + c2) * 2 + 1]) * half * AO[(c2 * 3 + c1) * 2 + 1]
                     -(srcO[(1 * 3 + c2) * 2 + 1] + flag * srcO[(2 * 3 + c2) * 2 + 0]) * half * AO[(c2 * 3 + c1) * 2 + 0];

            destE[(1 * 3 + c1) * 2 + 0] += tmp[0];
            destE[(1 * 3 + c1) * 2 + 1] += tmp[1];

        }
    }

    for (int i = 0; i < 12; i++) {
            dest_f[cont * 6 * 2 + i] = destE[i];
    }		
}

__global__ void transfer_y_f(double *src, double *dest_b, double *U,
                             const int v_x, const int v_y, const int v_z, const int v_t,
                             const int s_x, const int s_y, const int s_z, const int s_t,
                             const int rank, const int cb, const int flag) {
/*
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int z = blockDim.y * blockIdx.y + threadIdx.y;
    int t = blockDim.z * blockIdx.z + threadIdx.z;
*/
    const int s_x_cb = s_x >> 1;

//    if (x >= s_x_cb || z >= s_z || t >= s_t) { return; }

    const double half = 0.5;
//    const hipDoubleComplex I(0, 1);

    int subgrid_vol = (s_x * s_y * s_z * s_t);
    int subgrid_vol_cb = (subgrid_vol) >> 1;

    int y = 0;

    int xyz = blockDim.x * blockIdx.x + threadIdx.x;	

    int t = xyz / (s_z * s_x_cb);
    int z = (xyz / s_x_cb) % s_z;
    int x = xyz % s_x_cb;

    double tmp[2];	

/*
    hipDoubleComplex tmp;
    hipDoubleComplex *srcO = (hipDoubleComplex *) (src + (s_x_cb * s_y * s_z * t +
                                                          s_x_cb * s_y * z +
                                                          s_x_cb * y +
                                                          x + (1 - cb) * subgrid_vol_cb) * 12 * 2);
*/

    double *srcO = src + (s_x_cb * s_y * s_z * t +
                          s_x_cb * s_y * z +
                          s_x_cb * y +
                          x + (1 - cb) * subgrid_vol_cb) * 12 * 2;	

    int cont = s_x_cb * s_z * t + s_x_cb * z + x;

    for (int c2 = 0; c2 < 3; c2++) {
/*
        tmp = -(srcO[0 * 3 + c2] + flag * srcO[3 * 3 + c2]) * half;
        dest_b[cont * 6 * 2 + (0 * 3 + c2) * 2 + 0] = tmp.x;
        dest_b[cont * 6 * 2 + (0 * 3 + c2) * 2 + 1] = tmp.y;
        tmp = -(srcO[1 * 3 + c2] - flag * srcO[2 * 3 + c2]) * half;
        dest_b[cont * 6 * 2 + (1 * 3 + c2) * 2 + 0] = tmp.x;
        dest_b[cont * 6 * 2 + (1 * 3 + c2) * 2 + 1] = tmp.y;
*/

        tmp[0] = -(srcO[(0 * 3 + c2) * 2 + 0] + flag * srcO[(3 * 3 + c2) * 2 + 0]) * half;
        tmp[1] = -(srcO[(0 * 3 + c2) * 2 + 1] + flag * srcO[(3 * 3 + c2) * 2 + 1]) * half;

        dest_b[cont * 6 * 2 + (0 * 3 + c2) * 2 + 0] = tmp[0];
        dest_b[cont * 6 * 2 + (0 * 3 + c2) * 2 + 1] = tmp[1];

        tmp[0] = -(srcO[(1 * 3 + c2) * 2 + 0] - flag * srcO[(2 * 3 + c2) * 2 + 0]) * half;
        tmp[1] = -(srcO[(1 * 3 + c2) * 2 + 1] - flag * srcO[(2 * 3 + c2) * 2 + 1]) * half;

        dest_b[cont * 6 * 2 + (1 * 3 + c2) * 2 + 0] = tmp[0];
        dest_b[cont * 6 * 2 + (1 * 3 + c2) * 2 + 1] = tmp[1];

    }
}

__global__ void transfer_y_b(double *src, double *dest_f, double *U,
                             const int v_x, const int v_y, const int v_z, const int v_t,
                             const int s_x, const int s_y, const int s_z, const int s_t,
                             const int rank, const int cb, const int flag) {
/*
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int z = blockDim.y * blockIdx.y + threadIdx.y;
    int t = blockDim.z * blockIdx.z + threadIdx.z;
*/

    const int s_x_cb = s_x >> 1;

//    if (x >= s_x_cb || z >= s_z || t >= s_t) { return; }

    const double half = 0.5;
//    const hipDoubleComplex I(0, 1);

    int subgrid_vol = (s_x * s_y * s_z * s_t);
    int subgrid_vol_cb = (subgrid_vol) >> 1;

//    hipDoubleComplex tmp;

    int xyz = blockDim.x * blockIdx.x + threadIdx.x;	

    int y = s_y - 1;

    int t = xyz / (s_z * s_x_cb);
    int z = (xyz / s_x_cb) % s_z;
    int x = xyz % s_x_cb;	
/*
    hipDoubleComplex *srcO = (hipDoubleComplex *) (src + (s_x_cb * s_y * s_z * t +
                                                          s_x_cb * s_y * z +
                                                          s_x_cb * y +
                                                          x + (1 - cb) * subgrid_vol_cb) * 12 * 2);

    hipDoubleComplex *AO = (hipDoubleComplex *) (U + (s_x_cb * s_y * s_z * t +
                                                      s_x_cb * s_y * z +
                                                      s_x_cb * y +
                                                      x + (1 - cb) * subgrid_vol_cb) * 9 * 2);
*/
    double tmp[2];
    double destE[12];
    double srcO[24];
    double AO[18];

    for (int i = 0; i < 24; i++) {

        srcO[i] = src[(s_x_cb * s_y * s_z * t +
                       s_x_cb * s_y * z +
                       s_x_cb * y +
                       x + (1 - cb) * subgrid_vol_cb) * 12 * 2 + i];

    }

    for (int i = 0; i < 12; i++) {

        destE[i] = 0;

    }

    for (int i = 0; i < 18; i++) {

        AO[i] = U[(s_x_cb * s_y * s_z * t +
                   s_x_cb * s_y * z +
                   s_x_cb * y +
                   x + (1 - cb) * subgrid_vol_cb) * 9 * 2 + i];

    }

    int cont = s_x_cb * s_z * t + s_x_cb * z + x;
    for (int c1 = 0; c1 < 3; c1++) {
        for (int c2 = 0; c2 < 3; c2++) {
/*
            tmp = -(srcO[0 * 3 + c2] - flag * srcO[3 * 3 + c2]) * half * hipConj(AO[c2 * 3 + c1]);
            dest_f[cont * 6 * 2 + (0 * 3 + c1) * 2 + 0] += tmp.x;
            dest_f[cont * 6 * 2 + (0 * 3 + c1) * 2 + 1] += tmp.y;
            tmp = -(srcO[1 * 3 + c2] + flag * srcO[2 * 3 + c2]) * half * hipConj(AO[c2 * 3 + c1]);
            dest_f[cont * 6 * 2 + (1 * 3 + c1) * 2 + 0] += tmp.x;
            dest_f[cont * 6 * 2 + (1 * 3 + c1) * 2 + 1] += tmp.y;
*/

            tmp[0] = -(srcO[(0 * 3 + c2) * 2 + 0] - flag * srcO[(3 * 3 + c2) * 2 + 0]) * half * AO[(c2 * 3 + c1) * 2 + 0]
                     -(srcO[(0 * 3 + c2) * 2 + 1] - flag * srcO[(3 * 3 + c2) * 2 + 1]) * half * AO[(c2 * 3 + c1) * 2 + 1];
            tmp[1] = +(srcO[(0 * 3 + c2) * 2 + 0] - flag * srcO[(3 * 3 + c2) * 2 + 0]) * half * AO[(c2 * 3 + c1) * 2 + 1]
                     -(srcO[(0 * 3 + c2) * 2 + 1] - flag * srcO[(3 * 3 + c2) * 2 + 1]) * half * AO[(c2 * 3 + c1) * 2 + 0];

            destE[(0 * 3 + c1) * 2 + 0] += tmp[0];
            destE[(0 * 3 + c1) * 2 + 1] += tmp[1];

            tmp[0] = -(srcO[(1 * 3 + c2) * 2 + 0] + flag * srcO[(2 * 3 + c2) * 2 + 0]) * half * AO[(c2 * 3 + c1) * 2 + 0]
                     -(srcO[(1 * 3 + c2) * 2 + 1] + flag * srcO[(2 * 3 + c2) * 2 + 1]) * half * AO[(c2 * 3 + c1) * 2 + 1];
            tmp[1] = +(srcO[(1 * 3 + c2) * 2 + 0] + flag * srcO[(2 * 3 + c2) * 2 + 0]) * half * AO[(c2 * 3 + c1) * 2 + 1]
                     -(srcO[(1 * 3 + c2) * 2 + 1] + flag * srcO[(2 * 3 + c2) * 2 + 1]) * half * AO[(c2 * 3 + c1) * 2 + 0];

            destE[(1 * 3 + c1) * 2 + 0] += tmp[0];
            destE[(1 * 3 + c1) * 2 + 1] += tmp[1];

        }
    }

    for (int i = 0; i < 12; i++) {
            dest_f[cont * 6 * 2 + i] = destE[i];
    }	
}

__global__ void transfer_z_f(double *src, double *dest_b, double *U,
                             const int v_x, const int v_y, const int v_z, const int v_t,
                             const int s_x, const int s_y, const int s_z, const int s_t,
                             const int rank, const int cb, const int flag) {

//    int x = blockDim.x * blockIdx.x + threadIdx.x;
//    int y = blockDim.y * blockIdx.y + threadIdx.y;
//    int t = blockDim.z * blockIdx.z + threadIdx.z;

    const int s_x_cb = s_x >> 1;

    const double half = 0.5;
//    const hipDoubleComplex I(0, 1);

    int subgrid_vol = (s_x * s_y * s_z * s_t);
    int subgrid_vol_cb = (subgrid_vol) >> 1;

    int xyz = blockDim.x * blockIdx.x + threadIdx.x;

    int t = xyz / (s_y * s_x_cb);
    int y = (xyz / s_x_cb) % s_y;
    int x = xyz % s_x_cb;	

    int z = 0;
//    hipDoubleComplex tmp;

    double *srcO = src + (s_x_cb * s_y * s_z * t +
                          s_x_cb * s_y * z +
                          s_x_cb * y +
                          x + (1 - cb) * subgrid_vol_cb) * 12 * 2;

    double tmp[2];

    int cont = s_x_cb * s_y * t + s_x_cb * y + x;

    for (int c2 = 0; c2 < 3; c2++) {
/*
        tmp = -(srcO[0 * 3 + c2] - flag * I * srcO[2 * 3 + c2]) * half;
        dest_b[cont * 6 * 2 + (0 * 3 + c2) * 2 + 0] += tmp.x;
        dest_b[cont * 6 * 2 + (0 * 3 + c2) * 2 + 1] += tmp.y;
        tmp = -(srcO[1 * 3 + c2] + flag * I * srcO[3 * 3 + c2]) * half;
        dest_b[cont * 6 * 2 + (1 * 3 + c2) * 2 + 0] += tmp.x;
        dest_b[cont * 6 * 2 + (1 * 3 + c2) * 2 + 1] += tmp.y;
*/
	tmp[0] = -(srcO[(0 * 3 + c2) * 2 + 0] + flag * srcO[(2 * 3 + c2) * 2 + 1]) * half;
	tmp[1] = -(srcO[(0 * 3 + c2) * 2 + 1] - flag * srcO[(2 * 3 + c2) * 2 + 0]) * half;

	dest_b[cont * 6 * 2 + (0 * 3 + c2) * 2 + 0] = tmp[0];
	dest_b[cont * 6 * 2 + (0 * 3 + c2) * 2 + 1] = tmp[1];

	tmp[0] = -(srcO[(1 * 3 + c2) * 2 + 0] - flag * srcO[(3 * 3 + c2) * 2 + 1]) * half;
	tmp[1] = -(srcO[(1 * 3 + c2) * 2 + 1] + flag * srcO[(3 * 3 + c2) * 2 + 0]) * half;

        dest_b[cont * 6 * 2 + (1 * 3 + c2) * 2 + 0] = tmp[0];
        dest_b[cont * 6 * 2 + (1 * 3 + c2) * 2 + 1] = tmp[1];
    }
}

__global__ void transfer_z_b(double *src, double *dest_f, double *U,
                             const int v_x, const int v_y, const int v_z, const int v_t,
                             const int s_x, const int s_y, const int s_z, const int s_t,
                             const int rank, const int cb, const int flag) {

//    int x = blockDim.x * blockIdx.x + threadIdx.x;
//    int y = blockDim.y * blockIdx.y + threadIdx.y;
//    int t = blockDim.z * blockIdx.z + threadIdx.z;

    const int s_x_cb = s_x >> 1;

    const double half = 0.5;
//    const hipDoubleComplex I(0, 1);

    int subgrid_vol = (s_x * s_y * s_z * s_t);
    int subgrid_vol_cb = (subgrid_vol) >> 1;

    int xyz = blockDim.x * blockIdx.x + threadIdx.x;

    int t = xyz / (s_y * s_x_cb);
    int y = (xyz / s_x_cb) % s_y;
    int x = xyz % s_x_cb;

//    hipDoubleComplex tmp;

    int z = s_z - 1;
/*
    hipDoubleComplex *srcO = (hipDoubleComplex *) (src + (s_x_cb * s_y * s_z * t +
                                                          s_x_cb * s_y * z +
                                                          s_x_cb * y +
                                                          x + (1 - cb) * subgrid_vol_cb) * 12 * 2);

    hipDoubleComplex *AO = (hipDoubleComplex *) (U + (s_x_cb * s_y * s_z * t +
                                                      s_x_cb * s_y * z +
                                                      s_x_cb * y +
                                                      x + (1 - cb) * subgrid_vol_cb) * 9 * 2);
*/

    double tmp[2];
    double destE[12];
    double srcO[24];
    double AO[18];

    for (int i = 0; i < 24; i++) {

        srcO[i] = src[(s_x_cb * s_y * s_z * t +
                       s_x_cb * s_y * z +
                       s_x_cb * y +
                       x + (1 - cb) * subgrid_vol_cb) * 12 * 2 + i];

    }

    for (int i = 0; i < 12; i++) {

        destE[i] = 0;

    }

    for (int i = 0; i < 18; i++) {

        AO[i] = U[(s_x_cb * s_y * s_z * t +
                   s_x_cb * s_y * z +
                   s_x_cb * y +
                   x + (1 - cb) * subgrid_vol_cb) * 9 * 2 + i];

    }	

    int cont = s_x_cb * s_y * t + s_x_cb * y + x;
    for (int c1 = 0; c1 < 3; c1++) {
        for (int c2 = 0; c2 < 3; c2++) {
/*
            tmp = -(srcO[0 * 3 + c2] + flag * I * srcO[2 * 3 + c2]) * half * hipConj(AO[c2 * 3 + c1]);
            dest_f[cont * 6 * 2 + (0 * 3 + c1) * 2 + 0] += tmp.x;
            dest_f[cont * 6 * 2 + (0 * 3 + c1) * 2 + 1] += tmp.y;
            tmp = -(srcO[1 * 3 + c2] - flag * I * srcO[3 * 3 + c2]) * half * hipConj(AO[c2 * 3 + c1]);
            dest_f[cont * 6 * 2 + (1 * 3 + c1) * 2 + 0] += tmp.x;
            dest_f[cont * 6 * 2 + (1 * 3 + c1) * 2 + 1] += tmp.y;
*/

	    tmp[0] = -(srcO[(0 * 3 + c2) * 2 + 0] - flag * srcO[(2 * 3 + c2) * 2 + 1]) * half * AO[(c2 * 3 + c1) * 2 + 0]
		     -(srcO[(0 * 3 + c2) * 2 + 1] + flag * srcO[(2 * 3 + c2) * 2 + 0]) * half * AO[(c2 * 3 + c1) * 2 + 1];
	    tmp[1] = +(srcO[(0 * 3 + c2) * 2 + 0] - flag * srcO[(2 * 3 + c2) * 2 + 1]) * half * AO[(c2 * 3 + c1) * 2 + 1]
	             -(srcO[(0 * 3 + c2) * 2 + 1] + flag * srcO[(2 * 3 + c2) * 2 + 0]) * half * AO[(c2 * 3 + c1) * 2 + 0]; 	

	    destE[(0 * 3 + c1) * 2 + 0] += tmp[0];
            destE[(0 * 3 + c1) * 2 + 1] += tmp[1];

	    tmp[0] = -(srcO[(1 * 3 + c2) * 2 + 0] + flag * srcO[(3 * 3 + c2) * 2 + 1]) * half * AO[(c2 * 3 + c1) * 2 + 0]
		     -(srcO[(1 * 3 + c2) * 2 + 1] - flag * srcO[(3 * 3 + c2) * 2 + 0]) * half * AO[(c2 * 3 + c1) * 2 + 1];
	    tmp[1] = +(srcO[(1 * 3 + c2) * 2 + 0] + flag * srcO[(3 * 3 + c2) * 2 + 1]) * half * AO[(c2 * 3 + c1) * 2 + 1]
		     -(srcO[(1 * 3 + c2) * 2 + 1] - flag * srcO[(3 * 3 + c2) * 2 + 0]) * half * AO[(c2 * 3 + c1) * 2 + 0];	

	    destE[(1 * 3 + c1) * 2 + 0] += tmp[0];
            destE[(1 * 3 + c1) * 2 + 1] += tmp[1];
        }
    }

    for (int i = 0; i < 12; i++) {
            dest_f[cont * 6 * 2 + i] = destE[i];
    }	
}

__global__ void transfer_t_f(double *src, double *dest_b, double *U,
                             const int v_x, const int v_y, const int v_z, const int v_t,
                             const int s_x, const int s_y, const int s_z, const int s_t,
                             const int rank, const int cb, const int flag) {
/*
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;
    int z = blockDim.z * blockIdx.z + threadIdx.z;
*/


    const int s_x_cb = s_x >> 1;

    const double half = 0.5;
//    const hipDoubleComplex I(0, 1);

    int subgrid_vol = (s_x * s_y * s_z * s_t);
    int subgrid_vol_cb = (subgrid_vol) >> 1;

    int xyz = blockDim.x * blockIdx.x + threadIdx.x;

    int z = xyz / (s_y * s_x_cb);
    int y = (xyz / s_x_cb) % s_y;
    int x = xyz % s_x_cb;


    int t = 0;
/*
    hipDoubleComplex tmp;
    hipDoubleComplex *srcO = (hipDoubleComplex * )(src + (s_x_cb * s_y * s_z * t +
                                                          s_x_cb * s_y * z +
                                                          s_x_cb * y +
                                                          x + (1 - cb) * subgrid_vol_cb) * 12 * 2);
*/

    double *srcO = src + (s_x_cb * s_y * s_z * t +
                          s_x_cb * s_y * z +
                          s_x_cb * y +
                          x + (1 - cb) * subgrid_vol_cb) * 12 * 2;

    double tmp[2];


    int cont = s_x_cb * s_y * z + s_x_cb * y + x;

    for (int c2 = 0; c2 < 3; c2++) {
//        tmp = -(srcO[0 * 3 + c2] - flag * srcO[2 * 3 + c2]) * half;
        tmp[0] = -(srcO[(0 * 3 + c2) * 2 + 0] - flag * srcO[(2 * 3 + c2) * 2 + 0]) * half;
        tmp[1] = -(srcO[(0 * 3 + c2) * 2 + 1] - flag * srcO[(2 * 3 + c2) * 2 + 1]) * half;

//        dest_b[cont * 6 * 2 + (0 * 3 + c2) * 2 + 0] += tmp.x;
//        dest_b[cont * 6 * 2 + (0 * 3 + c2) * 2 + 1] += tmp.y;

        dest_b[cont * 6 * 2 + (0 * 3 + c2) * 2 + 0] = tmp[0];
        dest_b[cont * 6 * 2 + (0 * 3 + c2) * 2 + 1] = tmp[1];

//        tmp = -(srcO[1 * 3 + c2] - flag * srcO[3 * 3 + c2]) * half;

        tmp[0] = -(srcO[(1 * 3 + c2) * 2 + 0] - flag * srcO[(3 * 3 + c2) * 2 + 0]) * half;
        tmp[1] = -(srcO[(1 * 3 + c2) * 2 + 1] - flag * srcO[(3 * 3 + c2) * 2 + 1]) * half;

        dest_b[cont * 6 * 2 + (1 * 3 + c2) * 2 + 0] = tmp[0];
        dest_b[cont * 6 * 2 + (1 * 3 + c2) * 2 + 1] = tmp[1];
    }
}
         
__global__ void transfer_t_b(double *src, double *dest_f, double *U,
                             const int v_x, const int v_y, const int v_z, const int v_t,
                             const int s_x, const int s_y, const int s_z, const int s_t,
                             const int rank, const int cb, const int flag) {

//    int x = blockDim.x * blockIdx.x + threadIdx.x;
//    int y = blockDim.y * blockIdx.y + threadIdx.y;
//    int z = blockDim.z * blockIdx.z + threadIdx.z;

    const int s_x_cb = s_x >> 1;

    const double half = 0.5;
//    const hipDoubleComplex I(0, 1);

    int subgrid_vol = (s_x * s_y * s_z * s_t);
    int subgrid_vol_cb = (subgrid_vol) >> 1;

//    hipDoubleComplex tmp;

    int xyz = blockDim.x * blockIdx.x + threadIdx.x;

    int z = xyz / (s_y * s_x_cb);
    int y = (xyz / s_x_cb) % s_y;
    int x = xyz % s_x_cb;

    int t = s_t - 1;

/*
    hipDoubleComplex *srcO = (hipDoubleComplex * )(src + (s_x_cb * s_y * s_z * t +
                                                          s_x_cb * s_y * z +
                                                          s_x_cb * y +
                                                          x + (1 - cb) * subgrid_vol_cb) * 12 * 2);

    hipDoubleComplex *AO = (hipDoubleComplex * )(U + (s_x_cb * s_y * s_z * t +
                                                      s_x_cb * s_y * z +
                                                      s_x_cb * y +
                                                      x + (1 - cb) * subgrid_vol_cb) * 9 * 2);
*/

    double tmp[2];
    double destE[12];
    double srcO[24];
    double AO[18];

    for (int i = 0; i < 24; i++) {

        srcO[i] = src[(s_x_cb * s_y * s_z * t +
                       s_x_cb * s_y * z +
                       s_x_cb * y +
                       x + (1 - cb) * subgrid_vol_cb) * 12 * 2 + i];

    }

    for (int i = 0; i < 12; i++) {

        destE[i] = 0;

    }

    for (int i = 0; i < 18; i++) {

        AO[i] = U[(s_x_cb * s_y * s_z * t +
                   s_x_cb * s_y * z +
                   s_x_cb * y +
                   x + (1 - cb) * subgrid_vol_cb) * 9 * 2 + i];

    }

    int cont = s_x_cb * s_y * z + s_x_cb * y + x;
    for (int c1 = 0; c1 < 3; c1++) {
        for (int c2 = 0; c2 < 3; c2++) {
/*
            tmp = -(srcO[0 * 3 + c2] + flag * srcO[2 * 3 + c2]) * half * hipConj(AO[c2 * 3 + c1]);
            dest_f[cont * 6 * 2 + (0 * 3 + c1) * 2 + 0] += tmp.x;
            dest_f[cont * 6 * 2 + (0 * 3 + c1) * 2 + 1] += tmp.y;
            tmp = -(srcO[1 * 3 + c2] + flag * srcO[3 * 3 + c2]) * half * hipConj(AO[c2 * 3 + c1]);
            dest_f[cont * 6 * 2 + (1 * 3 + c1) * 2 + 0] += tmp.x;
            dest_f[cont * 6 * 2 + (1 * 3 + c1) * 2 + 1] += tmp.y;
*/
            tmp[0] = (-(srcO[(0 * 3 + c2) * 2 + 0] + flag * srcO[(2 * 3 + c2) * 2 + 0]) *
                      half *
                      AO[(c2 * 3 + c1) * 2 + 0]) -
                     ((srcO[(0 * 3 + c2) * 2 + 1] + flag * srcO[(2 * 3 + c2) * 2 + 1]) *
                      half *
                      AO[(c2 * 3 + c1) * 2 + 1]);

            tmp[1] = (+(srcO[(0 * 3 + c2) * 2 + 0] + flag * srcO[(2 * 3 + c2) * 2 + 0]) *
                      half *
                      AO[(c2 * 3 + c1) * 2 + 1]) -
                     ((srcO[(0 * 3 + c2) * 2 + 1] + flag * srcO[(2 * 3 + c2) * 2 + 1]) *
                      half *
                      AO[(c2 * 3 + c1) * 2 + 0]);

            destE[(0 * 3 + c1) * 2 + 0] += tmp[0];
            destE[(0 * 3 + c1) * 2 + 1] += tmp[1];

            tmp[0] = (-(srcO[(1 * 3 + c2) * 2 + 0] + flag * srcO[(3 * 3 + c2) * 2 + 0]) *
                      half *
                      AO[(c2 * 3 + c1) * 2 + 0]) -
                     ((srcO[(1 * 3 + c2) * 2 + 1] + flag * srcO[(3 * 3 + c2) * 2 + 1]) *
                      half *
                      AO[(c2 * 3 + c1) * 2 + 1]);

            tmp[1] = (+(srcO[(1 * 3 + c2) * 2 + 0] + flag * srcO[(3 * 3 + c2) * 2 + 0]) *
                      half *
                      AO[(c2 * 3 + c1) * 2 + 1])
                     - ((srcO[(1 * 3 + c2) * 2 + 1] + flag * srcO[(3 * 3 + c2) * 2 + 1]) *
                        half *
                        AO[(c2 * 3 + c1) * 2 + 0]);

            destE[(1 * 3 + c1) * 2 + 0] += tmp[0];
            destE[(1 * 3 + c1) * 2 + 1] += tmp[1];

        }
    }

    for (int i = 0; i < 12; i++) {
        dest_f[cont * 6 * 2 + i] = destE[i];
    }

}

__global__ void main_xyzt_abp5(double *src, double *dest,
                               double *U_x, double *U_y, double *U_z, double *U_t,
                               const int v_x, const int v_y, const int v_z, const int v_t,
                               const int s_x, const int s_y, const int s_z, const int s_t,
                               const int rank, const int cb, const int flag, const double a, const double b) {

    int subgrid_vol = (s_x * s_y * s_z * s_t);
    int subgrid_vol_cb = (subgrid_vol) >> 1;

    int N_sub[4] = {v_x / s_x, v_y / s_y, v_z / s_z, v_t / s_t};

    const double half = 0.5;
    const hipDoubleComplex I(0, 1);

    const int x_p = ((rank / N_sub[0]) % N_sub[1]) * s_y +
                    ((rank / (N_sub[1] * N_sub[0])) % N_sub[2]) * s_z +
                    (rank / (N_sub[2] * N_sub[1] * N_sub[0])) * s_t;

    const int s_x_cb = s_x >> 1;

    int xyzt = blockDim.x * blockIdx.x + threadIdx.x;

    int t = xyzt / (s_z * s_y * s_x_cb);
    int z = (xyzt / (s_y * s_x_cb)) % s_z;
    int y = (xyzt / s_x_cb) % s_y;
    int x = xyzt % s_x_cb;

    double tmp[2];
    double destE[24];
    double srcO[24];
    double AE[18];

    for (int i = 0; i < 24; i++) {
        destE[i] = 0;
    }

    int x_u = ((y + z + t + x_p) % 2 == cb || N_sub[0] == 1) ? s_x_cb : s_x_cb - 1;
    if (x < x_u) {
        int f_x = ((y + z + t + x_p) % 2 == cb) ? x : (x + 1) % s_x_cb;

        for (int i = 0; i < 24; i++) {

            srcO[i] = src[(s_x_cb * s_y * s_z * t +
                           s_x_cb * s_y * z +
                           s_x_cb * y +
                           f_x + (1 - cb) * subgrid_vol_cb) * 12 * 2 + i];

        }

        for (int i = 0; i < 18; i++) {

            AE[i] = U_x[(s_x_cb * s_y * s_z * t +
                         s_x_cb * s_y * z +
                         s_x_cb * y +
                         x + cb * subgrid_vol_cb) * 9 * 2 + i];

        }

        for (int c1 = 0; c1 < 3; c1++) {
            for (int c2 = 0; c2 < 3; c2++) {

                tmp[0] = -(srcO[(0 * 3 + c2) * 2 + 0] + flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                         AE[(c1 * 3 + c2) * 2 + 0]
                         + (srcO[(0 * 3 + c2) * 2 + 1] - flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                           AE[(c1 * 3 + c2) * 2 + 1];
                tmp[1] = -(srcO[(0 * 3 + c2) * 2 + 0] + flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                         AE[(c1 * 3 + c2) * 2 + 1]
                         - (srcO[(0 * 3 + c2) * 2 + 1] - flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                           AE[(c1 * 3 + c2) * 2 + 0];

                destE[(0 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(0 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(3 * 3 + c1) * 2 + 0] += -flag * tmp[1];
                destE[(3 * 3 + c1) * 2 + 1] += flag * tmp[0];

                tmp[0] = -(srcO[(1 * 3 + c2) * 2 + 0] + flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                         AE[(c1 * 3 + c2) * 2 + 0]
                         + (srcO[(1 * 3 + c2) * 2 + 1] - flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                           AE[(c1 * 3 + c2) * 2 + 1];
                tmp[1] = -(srcO[(1 * 3 + c2) * 2 + 0] + flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                         AE[(c1 * 3 + c2) * 2 + 1]
                         - (srcO[(1 * 3 + c2) * 2 + 1] - flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                           AE[(c1 * 3 + c2) * 2 + 0];

                destE[(1 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(1 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(2 * 3 + c1) * 2 + 0] += -flag * tmp[1];
                destE[(2 * 3 + c1) * 2 + 1] += flag * tmp[0];

            }
        }
    }

    int x_d = (((y + z + t + x_p) % 2) != cb || N_sub[0] == 1) ? 0 : 1;
    if (x >= x_d) {
        int b_x = ((t + z + y + x_p) % 2 == cb) ? (x - 1 + s_x_cb) % s_x_cb : x;

        for (int i = 0; i < 24; i++) {

            srcO[i] = src[(s_x_cb * s_y * s_z * t +
                           s_x_cb * s_y * z +
                           s_x_cb * y +
                           b_x + (1 - cb) * subgrid_vol_cb) * 12 * 2 + i];
        }

        for (int i = 0; i < 18; i++) {

            AE[i] = U_x[(s_x_cb * s_y * s_z * t +
                         s_x_cb * s_y * z +
                         s_x_cb * y +
                         b_x + (1 - cb) * subgrid_vol_cb) * 9 * 2 + i];
        }


        for (int c1 = 0; c1 < 3; c1++) {
            for (int c2 = 0; c2 < 3; c2++) {
                tmp[0] = -(srcO[(0 * 3 + c2) * 2 + 0] - flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                         AE[(c2 * 3 + c1) * 2 + 0]
                         - (srcO[(0 * 3 + c2) * 2 + 1] + flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                           AE[(c2 * 3 + c1) * 2 + 1];
                tmp[1] = +(srcO[(0 * 3 + c2) * 2 + 0] - flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                         AE[(c2 * 3 + c1) * 2 + 1]
                         - (srcO[(0 * 3 + c2) * 2 + 1] + flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                           AE[(c2 * 3 + c1) * 2 + 0];

                destE[(0 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(0 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(3 * 3 + c1) * 2 + 0] += flag * tmp[1];
                destE[(3 * 3 + c1) * 2 + 1] += -flag * tmp[0];

                tmp[0] = -(srcO[(1 * 3 + c2) * 2 + 0] - flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                         AE[(c2 * 3 + c1) * 2 + 0]
                         - (srcO[(1 * 3 + c2) * 2 + 1] + flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                           AE[(c2 * 3 + c1) * 2 + 1];
                tmp[1] = +(srcO[(1 * 3 + c2) * 2 + 0] - flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                         AE[(c2 * 3 + c1) * 2 + 1]
                         - (srcO[(1 * 3 + c2) * 2 + 1] + flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                           AE[(c2 * 3 + c1) * 2 + 0];

                destE[(1 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(1 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(2 * 3 + c1) * 2 + 0] += flag * tmp[1];
                destE[(2 * 3 + c1) * 2 + 1] += -flag * tmp[0];
            }
        }
    }

    int y_u = (N_sub[1] == 1) ? s_y : s_y - 1;
    if (y < y_u) {

        int f_y = (y + 1) % s_y;


        for (int i = 0; i < 24; i++) {

            srcO[i] = src[(s_x_cb * s_y * s_z * t +
                           s_x_cb * s_y * z +
                           s_x_cb * f_y +
                           x + (1 - cb) * subgrid_vol_cb) * 12 * 2 + i];
        }

        for (int i = 0; i < 18; i++) {


            AE[i] = U_y[(s_x_cb * s_y * s_z * t +
                         s_x_cb * s_y * z +
                         s_x_cb * y +
                         x + cb * subgrid_vol_cb) * 9 * 2 + i];

        }

        for (int c1 = 0; c1 < 3; c1++) {
            for (int c2 = 0; c2 < 3; c2++) {
                tmp[0] = -(srcO[(0 * 3 + c2) * 2 + 0] + flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                         AE[(c1 * 3 + c2) * 2 + 0]
                         + (srcO[(0 * 3 + c2) * 2 + 1] + flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                           AE[(c1 * 3 + c2) * 2 + 1];
                tmp[1] = -(srcO[(0 * 3 + c2) * 2 + 0] + flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                         AE[(c1 * 3 + c2) * 2 + 1]
                         - (srcO[(0 * 3 + c2) * 2 + 1] + flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                           AE[(c1 * 3 + c2) * 2 + 0];

                destE[(0 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(0 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(3 * 3 + c1) * 2 + 0] += flag * tmp[0];
                destE[(3 * 3 + c1) * 2 + 1] += flag * tmp[1];

                tmp[0] = -(srcO[(1 * 3 + c2) * 2 + 0] - flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                         AE[(c1 * 3 + c2) * 2 + 0]
                         + (srcO[(1 * 3 + c2) * 2 + 1] - flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                           AE[(c1 * 3 + c2) * 2 + 1];
                tmp[1] = -(srcO[(1 * 3 + c2) * 2 + 0] - flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                         AE[(c1 * 3 + c2) * 2 + 1]
                         - (srcO[(1 * 3 + c2) * 2 + 1] - flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                           AE[(c1 * 3 + c2) * 2 + 0];

                destE[(1 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(1 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(2 * 3 + c1) * 2 + 0] -= flag * tmp[0];
                destE[(2 * 3 + c1) * 2 + 1] -= flag * tmp[1];


            }
        }
    }

    int y_d = (N_sub[1] == 1) ? 0 : 1;
    if (y >= y_d) {
        int b_y = (y - 1 + s_y) % s_y;


        for (int i = 0; i < 24; i++) {

            srcO[i] = src[(s_x_cb * s_y * s_z * t +
                           s_x_cb * s_y * z +
                           s_x_cb * b_y +
                           x + (1 - cb) * subgrid_vol_cb) * 12 * 2 + i];
        }

        for (int i = 0; i < 18; i++) {

            AE[i] = U_y[(s_x_cb * s_y * s_z * t +
                         s_x_cb * s_y * z +
                         s_x_cb * b_y +
                         x + (1 - cb) * subgrid_vol_cb) * 9 * 2 + i];

        }

        for (int c1 = 0; c1 < 3; c1++) {
            for (int c2 = 0; c2 < 3; c2++) {
                tmp[0] = -(srcO[(0 * 3 + c2) * 2 + 0] - flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                         AE[(c2 * 3 + c1) * 2 + 0]
                         - (srcO[(0 * 3 + c2) * 2 + 1] - flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                           AE[(c2 * 3 + c1) * 2 + 1];
                tmp[1] = +(srcO[(0 * 3 + c2) * 2 + 0] - flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                         AE[(c2 * 3 + c1) * 2 + 1]
                         - (srcO[(0 * 3 + c2) * 2 + 1] - flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                           AE[(c2 * 3 + c1) * 2 + 0];

                destE[(0 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(0 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(3 * 3 + c1) * 2 + 0] -= flag * tmp[0];
                destE[(3 * 3 + c1) * 2 + 1] -= flag * tmp[1];

                tmp[0] = -(srcO[(1 * 3 + c2) * 2 + 0] + flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                         AE[(c2 * 3 + c1) * 2 + 0]
                         - (srcO[(1 * 3 + c2) * 2 + 1] + flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                           AE[(c2 * 3 + c1) * 2 + 1];
                tmp[1] = +(srcO[(1 * 3 + c2) * 2 + 0] + flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                         AE[(c2 * 3 + c1) * 2 + 1]
                         - (srcO[(1 * 3 + c2) * 2 + 1] + flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                           AE[(c2 * 3 + c1) * 2 + 0];

                destE[(1 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(1 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(2 * 3 + c1) * 2 + 0] += flag * tmp[0];
                destE[(2 * 3 + c1) * 2 + 1] += flag * tmp[1];

            }
        }
    }


    int z_u = (N_sub[2] == 1) ? s_z : s_z - 1;
    if (z < z_u) {
        int f_z = (z + 1) % s_z;

        for (int i = 0; i < 24; i++) {

            srcO[i] = src[(s_x_cb * s_y * s_z * t +
                           s_x_cb * s_y * f_z +
                           s_x_cb * y +
                           x + (1 - cb) * subgrid_vol_cb) * 12 * 2 + i];
        }

        for (int i = 0; i < 18; i++) {

            AE[i] = U_z[(s_x_cb * s_y * s_z * t +
                         s_x_cb * s_y * z +
                         s_x_cb * y +
                         x + cb * subgrid_vol_cb) * 9 * 2 + i];

        }

        for (int c1 = 0; c1 < 3; c1++) {
            for (int c2 = 0; c2 < 3; c2++) {
                tmp[0] = -(srcO[(0 * 3 + c2) * 2 + 0] + flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                         AE[(c1 * 3 + c2) * 2 + 0]
                         + (srcO[(0 * 3 + c2) * 2 + 1] - flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                           AE[(c1 * 3 + c2) * 2 + 1];
                tmp[1] = -(srcO[(0 * 3 + c2) * 2 + 0] + flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                         AE[(c1 * 3 + c2) * 2 + 1]
                         - (srcO[(0 * 3 + c2) * 2 + 1] - flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                           AE[(c1 * 3 + c2) * 2 + 0];

                destE[(0 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(0 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(2 * 3 + c1) * 2 + 0] += -flag * tmp[1];
                destE[(2 * 3 + c1) * 2 + 1] += flag * tmp[0];

                tmp[0] = -(srcO[(1 * 3 + c2) * 2 + 0] - flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                         AE[(c1 * 3 + c2) * 2 + 0]
                         + (srcO[(1 * 3 + c2) * 2 + 1] + flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                           AE[(c1 * 3 + c2) * 2 + 1];
                tmp[1] = -(srcO[(1 * 3 + c2) * 2 + 0] - flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                         AE[(c1 * 3 + c2) * 2 + 1]
                         - (srcO[(1 * 3 + c2) * 2 + 1] + flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                           AE[(c1 * 3 + c2) * 2 + 0];

                destE[(1 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(1 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(3 * 3 + c1) * 2 + 0] += flag * tmp[1];
                destE[(3 * 3 + c1) * 2 + 1] += -flag * tmp[0];
            }
        }
    }

    int z_d = (N_sub[2] == 1) ? 0 : 1;
    if (z >= z_d) {

        int b_z = (z - 1 + s_z) % s_z;


        for (int i = 0; i < 24; i++) {

            srcO[i] = src[(s_x_cb * s_y * s_z * t +
                           s_x_cb * s_y * b_z +
                           s_x_cb * y +
                           x + (1 - cb) * subgrid_vol_cb) * 12 * 2 + i];
        }

        for (int i = 0; i < 18; i++) {

            AE[i] = U_z[(s_x_cb * s_y * s_z * t +
                         s_x_cb * s_y * b_z +
                         s_x_cb * y +
                         x + (1 - cb) * subgrid_vol_cb) * 9 * 2 + i];
        }


        for (int c1 = 0; c1 < 3; c1++) {
            for (int c2 = 0; c2 < 3; c2++) {
                tmp[0] = -(srcO[(0 * 3 + c2) * 2 + 0] - flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                         AE[(c2 * 3 + c1) * 2 + 0]
                         - (srcO[(0 * 3 + c2) * 2 + 1] + flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                           AE[(c2 * 3 + c1) * 2 + 1];
                tmp[1] = +(srcO[(0 * 3 + c2) * 2 + 0] - flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                         AE[(c2 * 3 + c1) * 2 + 1]
                         - (srcO[(0 * 3 + c2) * 2 + 1] + flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                           AE[(c2 * 3 + c1) * 2 + 0];

                destE[(0 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(0 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(2 * 3 + c1) * 2 + 0] += flag * tmp[1];
                destE[(2 * 3 + c1) * 2 + 1] += -flag * tmp[0];

                tmp[0] = -(srcO[(1 * 3 + c2) * 2 + 0] + flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                         AE[(c2 * 3 + c1) * 2 + 0]
                         - (srcO[(1 * 3 + c2) * 2 + 1] - flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                           AE[(c2 * 3 + c1) * 2 + 1];
                tmp[1] = +(srcO[(1 * 3 + c2) * 2 + 0] + flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                         AE[(c2 * 3 + c1) * 2 + 1]
                         - (srcO[(1 * 3 + c2) * 2 + 1] - flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                           AE[(c2 * 3 + c1) * 2 + 0];

                destE[(1 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(1 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(3 * 3 + c1) * 2 + 0] += -flag * tmp[1];
                destE[(3 * 3 + c1) * 2 + 1] += flag * tmp[0];
            }
        }
    }

    int t_u = (N_sub[3] == 1) ? s_t : s_t - 1;
    if (t < t_u) {

        int f_t = (t + 1) % s_t;


        for (int i = 0; i < 24; i++) {

            srcO[i] = src[(s_x_cb * s_y * s_z * f_t +
                           s_x_cb * s_y * z +
                           s_x_cb * y +
                           x + (1 - cb) * subgrid_vol_cb) * 12 * 2 + i];
        }

        for (int i = 0; i < 18; i++) {

            AE[i] = U_t[(s_x_cb * s_y * s_z * t +
                         s_x_cb * s_y * z +
                         s_x_cb * y +
                         x + cb * subgrid_vol_cb) * 9 * 2 + i];
        }

        for (int c1 = 0; c1 < 3; c1++) {
            for (int c2 = 0; c2 < 3; c2++) {
                tmp[0] = -(srcO[(0 * 3 + c2) * 2 + 0] - flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                         AE[(c1 * 3 + c2) * 2 + 0]
                         + (srcO[(0 * 3 + c2) * 2 + 1] - flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                           AE[(c1 * 3 + c2) * 2 + 1];
                tmp[1] = -(srcO[(0 * 3 + c2) * 2 + 0] - flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                         AE[(c1 * 3 + c2) * 2 + 1]
                         - (srcO[(0 * 3 + c2) * 2 + 1] - flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                           AE[(c1 * 3 + c2) * 2 + 0];

                destE[(0 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(0 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(2 * 3 + c1) * 2 + 0] -= flag * tmp[0];
                destE[(2 * 3 + c1) * 2 + 1] -= flag * tmp[1];

                tmp[0] = -(srcO[(1 * 3 + c2) * 2 + 0] - flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                         AE[(c1 * 3 + c2) * 2 + 0]
                         + (srcO[(1 * 3 + c2) * 2 + 1] - flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                           AE[(c1 * 3 + c2) * 2 + 1];
                tmp[1] = -(srcO[(1 * 3 + c2) * 2 + 0] - flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                         AE[(c1 * 3 + c2) * 2 + 1]
                         - (srcO[(1 * 3 + c2) * 2 + 1] - flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                           AE[(c1 * 3 + c2) * 2 + 0];

                destE[(1 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(1 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(3 * 3 + c1) * 2 + 0] -= flag * tmp[0];
                destE[(3 * 3 + c1) * 2 + 1] -= flag * tmp[1];
            }
        }
    }

    int t_d = (N_sub[3] == 1) ? 0 : 1;
    if (t >= t_d) {
        int b_t = (t - 1 + s_t) % s_t;


        for (int i = 0; i < 24; i++) {

            srcO[i] = src[(s_x_cb * s_y * s_z * b_t +
                           s_x_cb * s_y * z +
                           s_x_cb * y +
                           x + (1 - cb) * subgrid_vol_cb) * 12 * 2 + i];
        }

        for (int i = 0; i < 18; i++) {

            AE[i] = U_t[(s_x_cb * s_y * s_z * b_t +
                         s_x_cb * s_y * z +
                         s_x_cb * y +
                         x + (1 - cb) * subgrid_vol_cb) * 9 * 2 + i];

        }

        for (int c1 = 0; c1 < 3; c1++) {
            for (int c2 = 0; c2 < 3; c2++) {
                tmp[0] = -(srcO[(0 * 3 + c2) * 2 + 0] + flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                         AE[(c2 * 3 + c1) * 2 + 0]
                         - (srcO[(0 * 3 + c2) * 2 + 1] + flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                           AE[(c2 * 3 + c1) * 2 + 1];
                tmp[1] = +(srcO[(0 * 3 + c2) * 2 + 0] + flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                         AE[(c2 * 3 + c1) * 2 + 1]
                         - (srcO[(0 * 3 + c2) * 2 + 1] + flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                           AE[(c2 * 3 + c1) * 2 + 0];

                destE[(0 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(0 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(2 * 3 + c1) * 2 + 0] += flag * tmp[0];
                destE[(2 * 3 + c1) * 2 + 1] += flag * tmp[1];

                tmp[0] = -(srcO[(1 * 3 + c2) * 2 + 0] + flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                         AE[(c2 * 3 + c1) * 2 + 0]
                         - (srcO[(1 * 3 + c2) * 2 + 1] + flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                           AE[(c2 * 3 + c1) * 2 + 1];
                tmp[1] = +(srcO[(1 * 3 + c2) * 2 + 0] + flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                         AE[(c2 * 3 + c1) * 2 + 1]
                         - (srcO[(1 * 3 + c2) * 2 + 1] + flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                           AE[(c2 * 3 + c1) * 2 + 0];

                destE[(1 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(1 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(3 * 3 + c1) * 2 + 0] += flag * tmp[0];
                destE[(3 * 3 + c1) * 2 + 1] += flag * tmp[1];
            }
        }
    }

    for (int i = 0; i < 24; i++) {

        srcO[i] = src[(s_x_cb * s_y * s_z * t +
                       s_x_cb * s_y * z +
                       s_x_cb * y +
                       x + cb * subgrid_vol_cb) * 12 * 2 + i];
    }

    for (int i = 0; i < 24; i++) {
        destE[i] = a * srcO[i] + b * destE[i];
    }


    for (int i = 0; i < 3; i++) {
        destE[2 * 3 * 2 + i * 2 + 0] = -destE[2 * 3 * 2 + i * 2 + 0];
        destE[2 * 3 * 2 + i * 2 + 1] = -destE[2 * 3 * 2 + i * 2 + 1];
        destE[3 * 3 * 2 + i * 2 + 0] = -destE[3 * 3 * 2 + i * 2 + 0];
        destE[3 * 3 * 2 + i * 2 + 1] = -destE[3 * 3 * 2 + i * 2 + 1];
    }

    for (int i = 0; i < 24; i++) {

        dest[(s_x_cb * s_y * s_z * t +
              s_x_cb * s_y * z +
              s_x_cb * y +
              x + cb * subgrid_vol_cb) * 12 * 2 + i] = destE[i];
    }
}

__global__ void main_xyzt(double *src, double *dest,
                          double *U_x, double *U_y, double *U_z, double *U_t,
                          const int v_x, const int v_y, const int v_z, const int v_t,
                          const int s_x, const int s_y, const int s_z, const int s_t,  // s_x s_y s_z s_t  from 16 to 24   
                          const int rank, const int cb, const int flag) {


    int subgrid_vol = (s_x * s_y * s_z * s_t);
    int subgrid_vol_cb = (subgrid_vol) >> 1;

    int N_sub[4] = {v_x / s_x, v_y / s_y, v_z / s_z, v_t / s_t};

    const double half = 0.5;

    const int x_p = ((rank / N_sub[0]) % N_sub[1]) * s_y +
                    ((rank / (N_sub[1] * N_sub[0])) % N_sub[2]) * s_z +
                    (rank / (N_sub[2] * N_sub[1] * N_sub[0])) * s_t;

    const int s_x_cb = s_x >> 1;

    int xyzt = blockDim.x * blockIdx.x + threadIdx.x;

    int t = xyzt / (s_z * s_y * s_x_cb);
    int z = (xyzt / (s_y * s_x_cb)) % s_z;
    int y = (xyzt / s_x_cb) % s_y;
    int x = xyzt % s_x_cb;

    double tmp[2];

    double destE[24];
    double srcO[24];
    double AE[18];

    for (int i = 0; i < 24; i++) {
        destE[i] = 0;
    }

    int x_u = ((y + z + t + x_p) % 2 == cb || N_sub[0] == 1) ? s_x_cb : s_x_cb - 1;
    if (x < x_u) {

        int f_x = ((y + z + t + x_p) % 2 == cb) ? x : (x + 1) % s_x_cb;

        for (int i = 0; i < 24; i++) {

            srcO[i] = src[(s_x_cb * s_y * s_z * t +
                           s_x_cb * s_y * z +
                           s_x_cb * y +
                           f_x + (1 - cb) * subgrid_vol_cb) * 12 * 2 + i];
        }

        for (int i = 0; i < 18; i++) {

            AE[i] = U_x[(s_x_cb * s_y * s_z * t +
                         s_x_cb * s_y * z +
                         s_x_cb * y +
                         x + cb * subgrid_vol_cb) * 9 * 2 + i];

        }


        for (int c1 = 0; c1 < 3; c1++) {
            for (int c2 = 0; c2 < 3; c2++) {

                tmp[0] = -(srcO[(0 * 3 + c2) * 2 + 0] + flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                         AE[(c1 * 3 + c2) * 2 + 0]
                         + (srcO[(0 * 3 + c2) * 2 + 1] - flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                           AE[(c1 * 3 + c2) * 2 + 1];
                tmp[1] = -(srcO[(0 * 3 + c2) * 2 + 0] + flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                         AE[(c1 * 3 + c2) * 2 + 1]
                         - (srcO[(0 * 3 + c2) * 2 + 1] - flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                           AE[(c1 * 3 + c2) * 2 + 0];

                destE[(0 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(0 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(3 * 3 + c1) * 2 + 0] += -flag * tmp[1];
                destE[(3 * 3 + c1) * 2 + 1] += flag * tmp[0];

                tmp[0] = -(srcO[(1 * 3 + c2) * 2 + 0] + flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                         AE[(c1 * 3 + c2) * 2 + 0]
                         + (srcO[(1 * 3 + c2) * 2 + 1] - flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                           AE[(c1 * 3 + c2) * 2 + 1];
                tmp[1] = -(srcO[(1 * 3 + c2) * 2 + 0] + flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                         AE[(c1 * 3 + c2) * 2 + 1]
                         - (srcO[(1 * 3 + c2) * 2 + 1] - flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                           AE[(c1 * 3 + c2) * 2 + 0];

                destE[(1 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(1 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(2 * 3 + c1) * 2 + 0] += -flag * tmp[1];
                destE[(2 * 3 + c1) * 2 + 1] += flag * tmp[0];
            }
        }
    }

    int x_d = (((y + z + t + x_p) % 2) != cb || N_sub[0] == 1) ? 0 : 1;
    if (x >= x_d) {

        int b_x = ((t + z + y + x_p) % 2 == cb) ? (x - 1 + s_x_cb) % s_x_cb : x;

        for (int i = 0; i < 24; i++) {

            srcO[i] = src[(s_x_cb * s_y * s_z * t +
                           s_x_cb * s_y * z +
                           s_x_cb * y +
                           b_x + (1 - cb) * subgrid_vol_cb) * 12 * 2 + i];

        }

        for (int i = 0; i < 18; i++) {

            AE[i] = U_x[(s_x_cb * s_y * s_z * t +
                         s_x_cb * s_y * z +
                         s_x_cb * y +
                         b_x + (1 - cb) * subgrid_vol_cb) * 9 * 2 + i];

        }

        for (int c1 = 0; c1 < 3; c1++) {
            for (int c2 = 0; c2 < 3; c2++) {

                tmp[0] = -(srcO[(0 * 3 + c2) * 2 + 0] - flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                         AE[(c2 * 3 + c1) * 2 + 0]
                         - (srcO[(0 * 3 + c2) * 2 + 1] + flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                           AE[(c2 * 3 + c1) * 2 + 1];
                tmp[1] = +(srcO[(0 * 3 + c2) * 2 + 0] - flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                         AE[(c2 * 3 + c1) * 2 + 1]
                         - (srcO[(0 * 3 + c2) * 2 + 1] + flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                           AE[(c2 * 3 + c1) * 2 + 0];

                destE[(0 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(0 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(3 * 3 + c1) * 2 + 0] += flag * tmp[1];
                destE[(3 * 3 + c1) * 2 + 1] += -flag * tmp[0];

                tmp[0] = -(srcO[(1 * 3 + c2) * 2 + 0] - flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                         AE[(c2 * 3 + c1) * 2 + 0]
                         - (srcO[(1 * 3 + c2) * 2 + 1] + flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                           AE[(c2 * 3 + c1) * 2 + 1];
                tmp[1] = +(srcO[(1 * 3 + c2) * 2 + 0] - flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                         AE[(c2 * 3 + c1) * 2 + 1]
                         - (srcO[(1 * 3 + c2) * 2 + 1] + flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                           AE[(c2 * 3 + c1) * 2 + 0];

                destE[(1 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(1 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(2 * 3 + c1) * 2 + 0] += flag * tmp[1];
                destE[(2 * 3 + c1) * 2 + 1] += -flag * tmp[0];
            }
        }
    }

    int y_u = (N_sub[1] == 1) ? s_y : s_y - 1;
    if (y < y_u) {

        int f_y = (y + 1) % s_y;

        for (int i = 0; i < 24; i++) {

            srcO[i] = src[(s_x_cb * s_y * s_z * t +
                           s_x_cb * s_y * z +
                           s_x_cb * f_y +
                           x + (1 - cb) * subgrid_vol_cb) * 12 * 2 + i];

        }

        for (int i = 0; i < 18; i++) {

            AE[i] = U_y[(s_x_cb * s_y * s_z * t +
                         s_x_cb * s_y * z +
                         s_x_cb * y +
                         x + cb * subgrid_vol_cb) * 9 * 2 + i];

        }

        for (int c1 = 0; c1 < 3; c1++) {
            for (int c2 = 0; c2 < 3; c2++) {

                tmp[0] = -(srcO[(0 * 3 + c2) * 2 + 0] + flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                         AE[(c1 * 3 + c2) * 2 + 0]
                         + (srcO[(0 * 3 + c2) * 2 + 1] + flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                           AE[(c1 * 3 + c2) * 2 + 1];
                tmp[1] = -(srcO[(0 * 3 + c2) * 2 + 0] + flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                         AE[(c1 * 3 + c2) * 2 + 1]
                         - (srcO[(0 * 3 + c2) * 2 + 1] + flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                           AE[(c1 * 3 + c2) * 2 + 0];

                destE[(0 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(0 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(3 * 3 + c1) * 2 + 0] += flag * tmp[0];
                destE[(3 * 3 + c1) * 2 + 1] += flag * tmp[1];

                tmp[0] = -(srcO[(1 * 3 + c2) * 2 + 0] - flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                         AE[(c1 * 3 + c2) * 2 + 0]
                         + (srcO[(1 * 3 + c2) * 2 + 1] - flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                           AE[(c1 * 3 + c2) * 2 + 1];
                tmp[1] = -(srcO[(1 * 3 + c2) * 2 + 0] - flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                         AE[(c1 * 3 + c2) * 2 + 1]
                         - (srcO[(1 * 3 + c2) * 2 + 1] - flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                           AE[(c1 * 3 + c2) * 2 + 0];

                destE[(1 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(1 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(2 * 3 + c1) * 2 + 0] -= flag * tmp[0];
                destE[(2 * 3 + c1) * 2 + 1] -= flag * tmp[1];


            }
        }
    }

    int y_d = (N_sub[1] == 1) ? 0 : 1;

    if (y >= y_d) {

        int b_y = (y - 1 + s_y) % s_y;

        for (int i = 0; i < 24; i++) {

            srcO[i] = src[(s_x_cb * s_y * s_z * t +
                           s_x_cb * s_y * z +
                           s_x_cb * b_y +
                           x + (1 - cb) * subgrid_vol_cb) * 12 * 2 + i];

        }

        for (int i = 0; i < 18; i++) {

            AE[i] = U_y[(s_x_cb * s_y * s_z * t +
                         s_x_cb * s_y * z +
                         s_x_cb * b_y +
                         x + (1 - cb) * subgrid_vol_cb) * 9 * 2 + i];

        }

        for (int c1 = 0; c1 < 3; c1++) {
            for (int c2 = 0; c2 < 3; c2++) {

                tmp[0] = -(srcO[(0 * 3 + c2) * 2 + 0] - flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                         AE[(c2 * 3 + c1) * 2 + 0]
                         - (srcO[(0 * 3 + c2) * 2 + 1] - flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                           AE[(c2 * 3 + c1) * 2 + 1];
                tmp[1] = +(srcO[(0 * 3 + c2) * 2 + 0] - flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                         AE[(c2 * 3 + c1) * 2 + 1]
                         - (srcO[(0 * 3 + c2) * 2 + 1] - flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                           AE[(c2 * 3 + c1) * 2 + 0];

                destE[(0 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(0 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(3 * 3 + c1) * 2 + 0] -= flag * tmp[0];
                destE[(3 * 3 + c1) * 2 + 1] -= flag * tmp[1];

                tmp[0] = -(srcO[(1 * 3 + c2) * 2 + 0] + flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                         AE[(c2 * 3 + c1) * 2 + 0]
                         - (srcO[(1 * 3 + c2) * 2 + 1] + flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                           AE[(c2 * 3 + c1) * 2 + 1];
                tmp[1] = +(srcO[(1 * 3 + c2) * 2 + 0] + flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                         AE[(c2 * 3 + c1) * 2 + 1]
                         - (srcO[(1 * 3 + c2) * 2 + 1] + flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                           AE[(c2 * 3 + c1) * 2 + 0];

                destE[(1 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(1 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(2 * 3 + c1) * 2 + 0] += flag * tmp[0];
                destE[(2 * 3 + c1) * 2 + 1] += flag * tmp[1];

            }
        }
    }


    int z_u = (N_sub[2] == 1) ? s_z : s_z - 1;
    if (z < z_u) {
        int f_z = (z + 1) % s_z;

        for (int i = 0; i < 24; i++) {

            srcO[i] = src[(s_x_cb * s_y * s_z * t +
                           s_x_cb * s_y * f_z +
                           s_x_cb * y +
                           x + (1 - cb) * subgrid_vol_cb) * 12 * 2 + i];
        }

        for (int i = 0; i < 18; i++) {

            AE[i] = U_z[(s_x_cb * s_y * s_z * t +
                         s_x_cb * s_y * z +
                         s_x_cb * y +
                         x + cb * subgrid_vol_cb) * 9 * 2 + i];

        }

        for (int c1 = 0; c1 < 3; c1++) {
            for (int c2 = 0; c2 < 3; c2++) {

                tmp[0] = -(srcO[(0 * 3 + c2) * 2 + 0] + flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                         AE[(c1 * 3 + c2) * 2 + 0]
                         + (srcO[(0 * 3 + c2) * 2 + 1] - flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                           AE[(c1 * 3 + c2) * 2 + 1];
                tmp[1] = -(srcO[(0 * 3 + c2) * 2 + 0] + flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                         AE[(c1 * 3 + c2) * 2 + 1]
                         - (srcO[(0 * 3 + c2) * 2 + 1] - flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                           AE[(c1 * 3 + c2) * 2 + 0];

                destE[(0 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(0 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(2 * 3 + c1) * 2 + 0] += -flag * tmp[1];
                destE[(2 * 3 + c1) * 2 + 1] += flag * tmp[0];

                tmp[0] = -(srcO[(1 * 3 + c2) * 2 + 0] - flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                         AE[(c1 * 3 + c2) * 2 + 0]
                         + (srcO[(1 * 3 + c2) * 2 + 1] + flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                           AE[(c1 * 3 + c2) * 2 + 1];
                tmp[1] = -(srcO[(1 * 3 + c2) * 2 + 0] - flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                         AE[(c1 * 3 + c2) * 2 + 1]
                         - (srcO[(1 * 3 + c2) * 2 + 1] + flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                           AE[(c1 * 3 + c2) * 2 + 0];

                destE[(1 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(1 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(3 * 3 + c1) * 2 + 0] += flag * tmp[1];
                destE[(3 * 3 + c1) * 2 + 1] += -flag * tmp[0];
            }
        }
    }


    int z_d = (N_sub[2] == 1) ? 0 : 1;
    if (z >= z_d) {

        int b_z = (z - 1 + s_z) % s_z;

        for (int i = 0; i < 24; i++) {

            srcO[i] = src[(s_x_cb * s_y * s_z * t +
                           s_x_cb * s_y * b_z +
                           s_x_cb * y +
                           x + (1 - cb) * subgrid_vol_cb) * 12 * 2 + i];

        }

        for (int i = 0; i < 18; i++) {

            AE[i] = U_z[(s_x_cb * s_y * s_z * t +
                         s_x_cb * s_y * b_z +
                         s_x_cb * y +
                         x + (1 - cb) * subgrid_vol_cb) * 9 * 2 + i];
        }


        for (int c1 = 0; c1 < 3; c1++) {
            for (int c2 = 0; c2 < 3; c2++) {

                tmp[0] = -(srcO[(0 * 3 + c2) * 2 + 0] - flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                         AE[(c2 * 3 + c1) * 2 + 0]
                         - (srcO[(0 * 3 + c2) * 2 + 1] + flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                           AE[(c2 * 3 + c1) * 2 + 1];
                tmp[1] = +(srcO[(0 * 3 + c2) * 2 + 0] - flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                         AE[(c2 * 3 + c1) * 2 + 1]
                         - (srcO[(0 * 3 + c2) * 2 + 1] + flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                           AE[(c2 * 3 + c1) * 2 + 0];

                destE[(0 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(0 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(2 * 3 + c1) * 2 + 0] += flag * tmp[1];
                destE[(2 * 3 + c1) * 2 + 1] += -flag * tmp[0];

                tmp[0] = -(srcO[(1 * 3 + c2) * 2 + 0] + flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                         AE[(c2 * 3 + c1) * 2 + 0]
                         - (srcO[(1 * 3 + c2) * 2 + 1] - flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                           AE[(c2 * 3 + c1) * 2 + 1];
                tmp[1] = +(srcO[(1 * 3 + c2) * 2 + 0] + flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                         AE[(c2 * 3 + c1) * 2 + 1]
                         - (srcO[(1 * 3 + c2) * 2 + 1] - flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                           AE[(c2 * 3 + c1) * 2 + 0];

                destE[(1 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(1 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(3 * 3 + c1) * 2 + 0] += -flag * tmp[1];
                destE[(3 * 3 + c1) * 2 + 1] += flag * tmp[0];
            }
        }
    }

    int t_u = (N_sub[3] == 1) ? s_t : s_t - 1;
    if (t < t_u) {

        int f_t = (t + 1) % s_t;

        for (int i = 0; i < 24; i++) {

            srcO[i] = src[(s_x_cb * s_y * s_z * f_t +
                           s_x_cb * s_y * z +
                           s_x_cb * y +
                           x + (1 - cb) * subgrid_vol_cb) * 12 * 2 + i];

        }

        for (int i = 0; i < 18; i++) {

            AE[i] = U_t[(s_x_cb * s_y * s_z * t +
                         s_x_cb * s_y * z +
                         s_x_cb * y +
                         x + cb * subgrid_vol_cb) * 9 * 2 + i];

        }

        for (int c1 = 0; c1 < 3; c1++) {
            for (int c2 = 0; c2 < 3; c2++) {

                tmp[0] = -(srcO[(0 * 3 + c2) * 2 + 0] - flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                         AE[(c1 * 3 + c2) * 2 + 0]
                         + (srcO[(0 * 3 + c2) * 2 + 1] - flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                           AE[(c1 * 3 + c2) * 2 + 1];
                tmp[1] = -(srcO[(0 * 3 + c2) * 2 + 0] - flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                         AE[(c1 * 3 + c2) * 2 + 1]
                         - (srcO[(0 * 3 + c2) * 2 + 1] - flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                           AE[(c1 * 3 + c2) * 2 + 0];

                destE[(0 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(0 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(2 * 3 + c1) * 2 + 0] -= flag * tmp[0];
                destE[(2 * 3 + c1) * 2 + 1] -= flag * tmp[1];

                tmp[0] = -(srcO[(1 * 3 + c2) * 2 + 0] - flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                         AE[(c1 * 3 + c2) * 2 + 0]
                         + (srcO[(1 * 3 + c2) * 2 + 1] - flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                           AE[(c1 * 3 + c2) * 2 + 1];
                tmp[1] = -(srcO[(1 * 3 + c2) * 2 + 0] - flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                         AE[(c1 * 3 + c2) * 2 + 1]
                         - (srcO[(1 * 3 + c2) * 2 + 1] - flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                           AE[(c1 * 3 + c2) * 2 + 0];

                destE[(1 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(1 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(3 * 3 + c1) * 2 + 0] -= flag * tmp[0];
                destE[(3 * 3 + c1) * 2 + 1] -= flag * tmp[1];
            }
        }
    }

    int t_d = (N_sub[3] == 1) ? 0 : 1;
    if (t >= t_d) {

        int b_t = (t - 1 + s_t) % s_t;

        for (int i = 0; i < 24; i++) {

            srcO[i] = src[(s_x_cb * s_y * s_z * b_t +
                           s_x_cb * s_y * z +
                           s_x_cb * y +
                           x + (1 - cb) * subgrid_vol_cb) * 12 * 2 + i];

        }

        for (int i = 0; i < 18; i++) {

            AE[i] = U_t[(s_x_cb * s_y * s_z * b_t +
                         s_x_cb * s_y * z +
                         s_x_cb * y +
                         x + (1 - cb) * subgrid_vol_cb) * 9 * 2 + i];

        }

        for (int c1 = 0; c1 < 3; c1++) {
            for (int c2 = 0; c2 < 3; c2++) {

                tmp[0] = -(srcO[(0 * 3 + c2) * 2 + 0] + flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                         AE[(c2 * 3 + c1) * 2 + 0]
                         - (srcO[(0 * 3 + c2) * 2 + 1] + flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                           AE[(c2 * 3 + c1) * 2 + 1];
                tmp[1] = +(srcO[(0 * 3 + c2) * 2 + 0] + flag * srcO[(2 * 3 + c2) * 2 + 0]) * half *
                         AE[(c2 * 3 + c1) * 2 + 1]
                         - (srcO[(0 * 3 + c2) * 2 + 1] + flag * srcO[(2 * 3 + c2) * 2 + 1]) * half *
                           AE[(c2 * 3 + c1) * 2 + 0];

                destE[(0 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(0 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(2 * 3 + c1) * 2 + 0] += flag * tmp[0];
                destE[(2 * 3 + c1) * 2 + 1] += flag * tmp[1];

                tmp[0] = -(srcO[(1 * 3 + c2) * 2 + 0] + flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                         AE[(c2 * 3 + c1) * 2 + 0]
                         - (srcO[(1 * 3 + c2) * 2 + 1] + flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                           AE[(c2 * 3 + c1) * 2 + 1];
                tmp[1] = +(srcO[(1 * 3 + c2) * 2 + 0] + flag * srcO[(3 * 3 + c2) * 2 + 0]) * half *
                         AE[(c2 * 3 + c1) * 2 + 1]
                         - (srcO[(1 * 3 + c2) * 2 + 1] + flag * srcO[(3 * 3 + c2) * 2 + 1]) * half *
                           AE[(c2 * 3 + c1) * 2 + 0];

                destE[(1 * 3 + c1) * 2 + 0] += tmp[0];
                destE[(1 * 3 + c1) * 2 + 1] += tmp[1];
                destE[(3 * 3 + c1) * 2 + 0] += flag * tmp[0];
                destE[(3 * 3 + c1) * 2 + 1] += flag * tmp[1];

            }
        }
    }

    for (int i = 0; i < 24; i++) {

        dest[(s_x_cb * s_y * s_z * t +
              s_x_cb * s_y * z +
              s_x_cb * y +
              x + cb * subgrid_vol_cb) * 12 * 2 + i] = destE[i];
    }
}

__global__ void ghost_x_f(double *src_f, double *dest, double *U,
                          const int v_x, const int v_y, const int v_z, const int v_t,
                          const int s_x, const int s_y, const int s_z, const int s_t,
                          const int rank, const int cb, const int flag) {

    int subgrid_vol = (s_x * s_y * s_z * s_t);
    int subgrid_vol_cb = (subgrid_vol) >> 1;

    int N_sub[4] = {v_x / s_x, v_y / s_y, v_z / s_z, v_t / s_t};

    const double half = 0.5;

    const int x_p = ((rank / N_sub[0]) % N_sub[1]) * s_y +
                    ((rank / (N_sub[1] * N_sub[0])) % N_sub[2]) * s_z +
                    (rank / (N_sub[2] * N_sub[1] * N_sub[0])) * s_t;

    int xyz = blockDim.x * blockIdx.x + threadIdx.x;

    int tp =  xyz / ( s_z * s_y  ) ;
    int z =  (xyz / s_y ) % s_z ;
    int y =  xyz % s_y ;	

    int t = (y + z + 2 * tp + x_p) % 2 == cb ? 2 * tp + 1 : 2 * tp;
    if (t >= s_t) { return; }

    const int s_x_cb = s_x >> 1;

    int cont = s_y * s_z * tp + s_y * z + y;

    int x = s_x_cb - 1;

    double tmp[2];
    double destE[24];
    double srcO[12];
    double AE[18];

    for (int i = 0; i < 24; i++) {
        destE[i] = 0;
    }

    for (int i = 0; i < 12; i++) {
        srcO[i] = src_f[cont * 6 * 2 + i];
    }

    for (int i = 0; i < 18; i++) {
        AE[i] = U[(s_x_cb * s_y * s_z * t +
                   s_x_cb * s_y * z +
                   s_x_cb * y +
                   x + cb * subgrid_vol_cb) * 9 * 2 + i];
    }	

    for (int c1 = 0; c1 < 3; c1++) {
        for (int c2 = 0; c2 < 3; c2++) {

	    tmp[0] = srcO[(0 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 0]
	    	    -srcO[(0 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 1];
	    tmp[1] = srcO[(0 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 1]
                    +srcO[(0 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 0];

            destE[(0 * 3 + c1) * 2 + 0] += tmp[0];
            destE[(0 * 3 + c1) * 2 + 1] += tmp[1];
            destE[(3 * 3 + c1) * 2 + 0] +=-flag * tmp[1];
            destE[(3 * 3 + c1) * 2 + 1] += flag * tmp[0];	

            tmp[0] = srcO[(1 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 0]
                    -srcO[(1 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 1];
            tmp[1] = srcO[(1 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 1]
                    +srcO[(1 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 0];

            destE[(1 * 3 + c1) * 2 + 0] += tmp[0];
            destE[(1 * 3 + c1) * 2 + 1] += tmp[1];
            destE[(2 * 3 + c1) * 2 + 0] +=-flag * tmp[1];
            destE[(2 * 3 + c1) * 2 + 1] += flag * tmp[0];    
  
         }
    }

    for (int i = 0; i < 24; i++) {
        dest[(s_x_cb * s_y * s_z * t +
              s_x_cb * s_y * z +
              s_x_cb * y +
              x + cb * subgrid_vol_cb) * 12 * 2 + i] += destE[i];
    }	

}

__global__ void ghost_x_b(double *src_b, double *dest, double *U,
                          const int v_x, const int v_y, const int v_z, const int v_t,
                          const int s_x, const int s_y, const int s_z, const int s_t,
                          const int rank, const int cb, const int flag) {

    int subgrid_vol = (s_x * s_y * s_z * s_t);
    int subgrid_vol_cb = (subgrid_vol) >> 1;

    int N_sub[4] = {v_x / s_x, v_y / s_y, v_z / s_z, v_t / s_t};

    const double half = 0.5;

    const int x_p = ((rank / N_sub[0]) % N_sub[1]) * s_y +
                    ((rank / (N_sub[1] * N_sub[0])) % N_sub[2]) * s_z +
                    (rank / (N_sub[2] * N_sub[1] * N_sub[0])) * s_t;

    int xyz = blockDim.x * blockIdx.x + threadIdx.x;

    int tp =  xyz / ( s_z * s_y  ) ;
    int z =  (xyz / s_y ) % s_z ;
    int y =  xyz % s_y ;	

    int t = (y + z + 2 * tp + x_p) % 2 != cb ? 2 * tp + 1 : 2 * tp;
    if (t >= s_t) { return; }

    const int s_x_cb = s_x >> 1;

    int cont = s_y * s_z * tp + s_y * z + y;

    int x = 0;

    double srcO[12];

    for ( int i=0 ; i < 12 ; i++ ){
        srcO[i] = src_b[cont * 6 * 2 + i];
    }

    double *destE = dest + (s_x_cb * s_y * s_z * t +
                            s_x_cb * s_y * z +
                            s_x_cb * y +
                            x + cb * subgrid_vol_cb) * 12 * 2;

    for (int c1 = 0; c1 < 3; c1++) {

        destE[(0 * 3 + c1) * 2 + 0] += srcO[(0 * 3 + c1) * 2 + 0];
        destE[(0 * 3 + c1) * 2 + 1] += srcO[(0 * 3 + c1) * 2 + 1];
        destE[(3 * 3 + c1) * 2 + 0] += flag * srcO[(0 * 3 + c1) * 2 + 1];
        destE[(3 * 3 + c1) * 2 + 1] +=-flag * srcO[(0 * 3 + c1) * 2 + 0];
        destE[(1 * 3 + c1) * 2 + 0] += srcO[(1 * 3 + c1) * 2 + 0];
        destE[(1 * 3 + c1) * 2 + 1] += srcO[(1 * 3 + c1) * 2 + 1];
        destE[(2 * 3 + c1) * 2 + 0] += flag * srcO[(1 * 3 + c1) * 2 + 1];
        destE[(2 * 3 + c1) * 2 + 1] +=-flag * srcO[(1 * 3 + c1) * 2 + 0];
    }
}

__global__ void ghost_y_f(double *src_f, double *dest, double *U,
                          const int v_x, const int v_y, const int v_z, const int v_t,
                          const int s_x, const int s_y, const int s_z, const int s_t,
                          const int rank, const int cb, const int flag) {

    const int s_x_cb = s_x >> 1;

    int subgrid_vol = (s_x * s_y * s_z * s_t);
    int subgrid_vol_cb = (subgrid_vol) >> 1;

    const double half = 0.5;

    int y = s_y - 1;

    int xyz = blockDim.x * blockIdx.x + threadIdx.x;	

    int t =  xyz / ( s_z * s_x_cb  ) ;
    int z =  (xyz / s_x_cb ) % s_z ;
    int x =  xyz % s_x_cb ;

    double tmp[2];
    double destE[24];
    double srcO[12];
    double AE[18];	

    int cont = s_x_cb * s_z * t + s_x_cb * z + x;

    for (int i = 0; i < 24; i++) {
        destE[i] = 0;
    }

    for (int i = 0; i < 12; i++) {
        srcO[i] = src_f[cont * 6 * 2 + i];
    }

    for (int i = 0; i < 18; i++) {
        AE[i] = U[(s_x_cb * s_y * s_z * t +
                   s_x_cb * s_y * z +
                   s_x_cb * y +
                   x + cb * subgrid_vol_cb) * 9 * 2 + i];
    }

    for (int c1 = 0; c1 < 3; c1++) {
        for (int c2 = 0; c2 < 3; c2++) {

            tmp[0] = + srcO[(0 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 0]
                     - srcO[(0 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 1];

            tmp[1] = + srcO[(0 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 1]
                     + srcO[(0 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 0];

            destE[(0 * 3 + c1) * 2 + 0] += tmp[0];
            destE[(0 * 3 + c1) * 2 + 1] += tmp[1];
            destE[(3 * 3 + c1) * 2 + 0] += flag * tmp[0];
            destE[(3 * 3 + c1) * 2 + 1] += flag * tmp[1];

            tmp[0] = + srcO[(1 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 0]
                     - srcO[(1 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 1];

            tmp[1] = + srcO[(1 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 1]
                     + srcO[(1 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 0];

            destE[(1 * 3 + c1) * 2 + 0] += tmp[0];
            destE[(1 * 3 + c1) * 2 + 1] += tmp[1];
            destE[(2 * 3 + c1) * 2 + 0] -= flag * tmp[0];
            destE[(2 * 3 + c1) * 2 + 1] -= flag * tmp[1];

        }
    }

    for (int i = 0; i < 24; i++) {
        dest[(s_x_cb * s_y * s_z * t +
              s_x_cb * s_y * z +
              s_x_cb * y +
              x + cb * subgrid_vol_cb) * 12 * 2 + i] += destE[i];
    }	

}

__global__ void ghost_y_b(double *src_b, double *dest, double *U,
                          const int v_x, const int v_y, const int v_z, const int v_t,
                          const int s_x, const int s_y, const int s_z, const int s_t,
                          const int rank, const int cb, const int flag) {

    const int s_x_cb = s_x >> 1;

    int subgrid_vol = (s_x * s_y * s_z * s_t);
    int subgrid_vol_cb = (subgrid_vol) >> 1;

    const double half = 0.5;

    int y = 0;

    int xyz = blockDim.x * blockIdx.x + threadIdx.x;	

    int t =  xyz / ( s_z * s_x_cb  ) ;
    int z =  (xyz / s_x_cb ) % s_z ;
    int x =  xyz % s_x_cb ;	

    int cont = s_x_cb * s_z * t + s_x_cb * z + x;

    double srcO[12];

    for ( int i=0 ; i < 12 ; i++ ){
        srcO[i] = src_b[cont * 6 * 2 + i];
    }	

    double *destE = dest + (s_x_cb * s_y * s_z * t +
                            s_x_cb * s_y * z +
                            s_x_cb * y +
                            x + cb * subgrid_vol_cb) * 12 * 2;

    for (int c1 = 0; c1 < 3; c1++) {
        destE[(0 * 3 + c1) * 2 + 0] += srcO[(0 * 3 + c1) * 2 + 0];
        destE[(0 * 3 + c1) * 2 + 1] += srcO[(0 * 3 + c1) * 2 + 1];
        destE[(3 * 3 + c1) * 2 + 0] -= flag * srcO[(0 * 3 + c1) * 2 + 0];
        destE[(3 * 3 + c1) * 2 + 1] -= flag * srcO[(0 * 3 + c1) * 2 + 1];
        destE[(1 * 3 + c1) * 2 + 0] += srcO[(1 * 3 + c1) * 2 + 0];
        destE[(1 * 3 + c1) * 2 + 1] += srcO[(1 * 3 + c1) * 2 + 1];
        destE[(2 * 3 + c1) * 2 + 0] += flag * srcO[(1 * 3 + c1) * 2 + 0];
        destE[(2 * 3 + c1) * 2 + 1] += flag * srcO[(1 * 3 + c1) * 2 + 1];
    }
}

__global__ void ghost_z_f(double *src_f, double *dest, double *U,
                          const int v_x, const int v_y, const int v_z, const int v_t,
                          const int s_x, const int s_y, const int s_z, const int s_t,
                          const int rank, const int cb, const int flag) {

    int xyz = blockDim.x * blockIdx.x + threadIdx.x;

    const int s_x_cb = s_x >> 1;

    int subgrid_vol = (s_x * s_y * s_z * s_t);
    int subgrid_vol_cb = (subgrid_vol) >> 1;

    const double half = 0.5;

    int z = s_z - 1;

    int t =  xyz / ( s_y * s_x_cb  ) ;
    int y =  (xyz / s_x_cb ) % s_y ;
    int x =  xyz % s_x_cb ;	

    int cont = s_x_cb * s_y * t + s_x_cb * y + x;

    double tmp[2];
    double destE[24];
    double srcO[12];
    double AE[18];

    for (int i = 0; i < 24; i++) {
        destE[i] = 0;
    }

    for (int i = 0; i < 12; i++) {
        srcO[i] = src_f[cont * 6 * 2 + i];
    }

    for (int i = 0; i < 18; i++) {
        AE[i] = U[(s_x_cb * s_y * s_z * t +
                   s_x_cb * s_y * z +
                   s_x_cb * y +
                   x + cb * subgrid_vol_cb) * 9 * 2 + i];
    }

    for (int c1 = 0; c1 < 3; c1++) {
        for (int c2 = 0; c2 < 3; c2++) {

            tmp[0] = + srcO[(0 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 0]
                     - srcO[(0 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 1];

            tmp[1] = + srcO[(0 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 1]
                     + srcO[(0 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 0];

            destE[(0 * 3 + c1) * 2 + 0] += tmp[0];
            destE[(0 * 3 + c1) * 2 + 1] += tmp[1];
            destE[(2 * 3 + c1) * 2 + 0] +=-flag * tmp[1];
            destE[(2 * 3 + c1) * 2 + 1] += flag * tmp[0];

            tmp[0] = + srcO[(1 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 0]
                     - srcO[(1 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 1];

            tmp[1] = + srcO[(1 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 1]
                     + srcO[(1 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 0];

            destE[(1 * 3 + c1) * 2 + 0] += tmp[0];
            destE[(1 * 3 + c1) * 2 + 1] += tmp[1];
            destE[(3 * 3 + c1) * 2 + 0] += flag * tmp[1];
            destE[(3 * 3 + c1) * 2 + 1] +=-flag * tmp[0];	

        }
    }

    for (int i = 0; i < 24; i++) {
        dest[(s_x_cb * s_y * s_z * t +
              s_x_cb * s_y * z +
              s_x_cb * y +
              x + cb * subgrid_vol_cb) * 12 * 2 + i] += destE[i];
    }	

}

__global__ void ghost_z_b(double *src_b, double *dest, double *U,
                          const int v_x, const int v_y, const int v_z, const int v_t,
                          const int s_x, const int s_y, const int s_z, const int s_t,
                          const int rank, const int cb, const int flag) {

    int xyz = blockDim.x * blockIdx.x + threadIdx.x;
   	
    const int s_x_cb = s_x >> 1;

    int subgrid_vol = (s_x * s_y * s_z * s_t);
    int subgrid_vol_cb = (subgrid_vol) >> 1;

    const double half = 0.5;

    int z = 0;

    int t =  xyz / ( s_y * s_x_cb  ) ;
    int y =  (xyz / s_x_cb ) % s_y ;
    int x =  xyz % s_x_cb ;	

    int cont = s_x_cb * s_y * t + s_x_cb * y + x;

    double srcO[12];

    for ( int i=0 ; i < 12 ; i++ ){
        srcO[i] = src_b[cont * 6 * 2 + i];
    }	

    double *destE = dest + (s_x_cb * s_y * s_z * t +
                            s_x_cb * s_y * z +
                            s_x_cb * y +
                            x + cb * subgrid_vol_cb) * 12 * 2;

    for (int c1 = 0; c1 < 3; c1++) {
	destE[(0 * 3 + c1) * 2 + 0] += srcO[(0 * 3 + c1) * 2 + 0];
	destE[(0 * 3 + c1) * 2 + 1] += srcO[(0 * 3 + c1) * 2 + 1];

	destE[(2 * 3 + c1) * 2 + 0] += flag * srcO[(0 * 3 + c1) * 2 + 1];
	destE[(2 * 3 + c1) * 2 + 1] +=-flag * srcO[(0 * 3 + c1) * 2 + 0];

	destE[(1 * 3 + c1) * 2 + 0] += srcO[(1 * 3 + c1) * 2 + 0];
	destE[(1 * 3 + c1) * 2 + 1] += srcO[(1 * 3 + c1) * 2 + 1];
	
	destE[(3 * 3 + c1) * 2 + 0] +=-flag * srcO[(1 * 3 + c1) * 2 + 1];
	destE[(3 * 3 + c1) * 2 + 1] += flag * srcO[(1 * 3 + c1) * 2 + 0];

    }
}

__global__ void ghost_t_f(double *src_f, double *dest, double *U,
                          const int v_x, const int v_y, const int v_z, const int v_t,
                          const int s_x, const int s_y, const int s_z, const int s_t,
                          const int rank, const int cb, const int flag) {

    int xyz = blockDim.x * blockIdx.x + threadIdx.x;

    const int s_x_cb = s_x >> 1;

    int subgrid_vol = (s_x * s_y * s_z * s_t);
    int subgrid_vol_cb = (subgrid_vol) >> 1;

    const double half = 0.5;

    int t = s_t - 1;

    int z =  xyz / ( s_y * s_x_cb  ) ;
    int y =  (xyz / s_x_cb ) % s_y ;
    int x =  xyz % s_x_cb ;	

    double tmp[2];
    double destE[24];
    double srcO[12];
    double AE[18];

    int cont = s_x_cb * s_y * z + s_x_cb * y + x;

    for (int i = 0; i < 24; i++) {
        destE[i] = 0;
    }

    for (int i = 0; i < 12; i++) {
        srcO[i] = src_f[cont * 6 * 2 + i];
    }

    for (int i = 0; i < 18; i++) {
        AE[i] = U[(s_x_cb * s_y * s_z * t +
                   s_x_cb * s_y * z +
                   s_x_cb * y +
                   x + cb * subgrid_vol_cb) * 9 * 2 + i];
    }

    for (int c1 = 0; c1 < 3; c1++) {
        for (int c2 = 0; c2 < 3; c2++) {

            tmp[0] = +srcO[(0 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 0]
                     - srcO[(0 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 1];

            tmp[1] = +srcO[(0 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 1]
                     + srcO[(0 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 0];

            destE[(0 * 3 + c1) * 2 + 0] += tmp[0];
            destE[(0 * 3 + c1) * 2 + 1] += tmp[1];
            destE[(2 * 3 + c1) * 2 + 0] -= flag * tmp[0];
            destE[(2 * 3 + c1) * 2 + 1] -= flag * tmp[1];

            tmp[0] = +srcO[(1 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 0]
                     - srcO[(1 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 1];

            tmp[1] = +srcO[(1 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 1]
                     + srcO[(1 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 0];

            destE[(1 * 3 + c1) * 2 + 0] += tmp[0];
            destE[(1 * 3 + c1) * 2 + 1] += tmp[1];
            destE[(3 * 3 + c1) * 2 + 0] -= flag * tmp[0];
            destE[(3 * 3 + c1) * 2 + 1] -= flag * tmp[1];

        }
    }

    for (int i = 0; i < 24; i++) {
        dest[(s_x_cb * s_y * s_z * t +
              s_x_cb * s_y * z +
              s_x_cb * y +
              x + cb * subgrid_vol_cb) * 12 * 2 + i] += destE[i];
    }
}

__global__ void ghost_t_b(double *src_b, double *dest, double *U,
                          const int v_x, const int v_y, const int v_z, const int v_t,
                          const int s_x, const int s_y, const int s_z, const int s_t,
                          const int rank, const int cb, const int flag) {
    const int s_x_cb = s_x >> 1;

    int subgrid_vol = (s_x * s_y * s_z * s_t);
    int subgrid_vol_cb = (subgrid_vol) >> 1;

    const double half = 0.5;

    int t = 0;

    int xyz = blockDim.x * blockIdx.x + threadIdx.x;

    int z =  xyz / ( s_y * s_x_cb  ) ;
    int y =  (xyz / s_x_cb ) % s_y ;
    int x =  xyz % s_x_cb ;	

    int cont = s_x_cb * s_y * z + s_x_cb * y + x;

    double srcO[12];

    for (int i = 0; i < 12; i++) {
        srcO[i] = src_b[cont * 6 * 2 + i];
    }	

    double *destE = dest + (s_x_cb * s_y * s_z * t +
                            s_x_cb * s_y * z +
                            s_x_cb * y +
                            x + cb * subgrid_vol_cb) * 12 * 2;

    for (int c1 = 0; c1 < 3; c1++) {
        destE[(0 * 3 + c1) * 2 + 0] += srcO[(0 * 3 + c1) * 2 + 0];
        destE[(0 * 3 + c1) * 2 + 1] += srcO[(0 * 3 + c1) * 2 + 1];
        destE[(2 * 3 + c1) * 2 + 0] += flag * srcO[(0 * 3 + c1) * 2 + 0];
        destE[(2 * 3 + c1) * 2 + 1] += flag * srcO[(0 * 3 + c1) * 2 + 1];
        destE[(1 * 3 + c1) * 2 + 0] += srcO[(1 * 3 + c1) * 2 + 0];
        destE[(1 * 3 + c1) * 2 + 1] += srcO[(1 * 3 + c1) * 2 + 1];
        destE[(3 * 3 + c1) * 2 + 0] += flag * srcO[(1 * 3 + c1) * 2 + 0];
        destE[(3 * 3 + c1) * 2 + 1] += flag * srcO[(1 * 3 + c1) * 2 + 1];
    }
}

///////////////////////////////////////////////////////////// for overlap  /////////////////////////////////////////////////////////////////////////////////////////////
//
//                sink * gamma5 * b
//
//
/////////////////////////////////////////////////////////////////

__global__ void ghost_x_f_abp5(double *src_f, double *dest, double *U,
                               const int v_x, const int v_y, const int v_z, const int v_t,
                               const int s_x, const int s_y, const int s_z, const int s_t,
                               const int rank, const int cb, const int flag, const double b) {

    int subgrid_vol = (s_x * s_y * s_z * s_t);
    int subgrid_vol_cb = (subgrid_vol) >> 1;

    int N_sub[4] = {v_x / s_x, v_y / s_y, v_z / s_z, v_t / s_t};

    const double half = 0.5;

    const int x_p = ((rank / N_sub[0]) % N_sub[1]) * s_y +
                    ((rank / (N_sub[1] * N_sub[0])) % N_sub[2]) * s_z +
                    (rank / (N_sub[2] * N_sub[1] * N_sub[0])) * s_t;

    int xyz = blockDim.x * blockIdx.x + threadIdx.x;

    int tp =  xyz / ( s_z * s_y  ) ;
    int z =  (xyz / s_y ) % s_z ;
    int y =  xyz % s_y ;

    int t = (y + z + 2 * tp + x_p) % 2 == cb ? 2 * tp + 1 : 2 * tp;
    if (t >= s_t) { return; }

    const int s_x_cb = s_x >> 1;

    int cont = s_y * s_z * tp + s_y * z + y;

    int x = s_x_cb - 1;

    double tmp[2];
    double destE[24];
    double srcO[12];
    double AE[18];

    for (int i = 0; i < 24; i++) {
        destE[i] = 0;
    }

    for (int i = 0; i < 12; i++) {
        srcO[i] = src_f[cont * 6 * 2 + i];
    }

    for (int i = 0; i < 18; i++) {
        AE[i] = U[(s_x_cb * s_y * s_z * t +
                   s_x_cb * s_y * z +
                   s_x_cb * y +
                   x + cb * subgrid_vol_cb) * 9 * 2 + i];
    }

    for (int c1 = 0; c1 < 3; c1++) {
        for (int c2 = 0; c2 < 3; c2++) {

            tmp[0] = srcO[(0 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 0]
                    -srcO[(0 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 1];
            tmp[1] = srcO[(0 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 1]
                    +srcO[(0 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 0];

            destE[(0 * 3 + c1) * 2 + 0] += tmp[0] * b;
            destE[(0 * 3 + c1) * 2 + 1] += tmp[1] * b;
            destE[(3 * 3 + c1) * 2 + 0] += flag * tmp[1] * b;
            destE[(3 * 3 + c1) * 2 + 1] +=-flag * tmp[0] * b;

            tmp[0] = srcO[(1 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 0]
                    -srcO[(1 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 1];
            tmp[1] = srcO[(1 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 1]
                    +srcO[(1 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 0];

            destE[(1 * 3 + c1) * 2 + 0] += tmp[0] * b;
            destE[(1 * 3 + c1) * 2 + 1] += tmp[1] * b;
            destE[(2 * 3 + c1) * 2 + 0] += flag * tmp[1] * b;
            destE[(2 * 3 + c1) * 2 + 1] +=-flag * tmp[0] * b;

         }
    }

    for (int i = 0; i < 24; i++) {
        dest[(s_x_cb * s_y * s_z * t +
              s_x_cb * s_y * z +
              s_x_cb * y +
              x + cb * subgrid_vol_cb) * 12 * 2 + i] += destE[i];
    }

}
	

__global__ void ghost_x_b_abp5(double *src_b, double *dest, double *U,
                               const int v_x, const int v_y, const int v_z, const int v_t,
                               const int s_x, const int s_y, const int s_z, const int s_t,
                               const int rank, const int cb, const int flag, const double b) {

    int subgrid_vol = (s_x * s_y * s_z * s_t);
    int subgrid_vol_cb = (subgrid_vol) >> 1;

    int N_sub[4] = {v_x / s_x, v_y / s_y, v_z / s_z, v_t / s_t};

    const double half = 0.5;

    const int x_p = ((rank / N_sub[0]) % N_sub[1]) * s_y +
                    ((rank / (N_sub[1] * N_sub[0])) % N_sub[2]) * s_z +
                    (rank / (N_sub[2] * N_sub[1] * N_sub[0])) * s_t;

    int xyz = blockDim.x * blockIdx.x + threadIdx.x;

    int tp =  xyz / ( s_z * s_y  ) ;
    int z =  (xyz / s_y ) % s_z ;
    int y =  xyz % s_y ;

    int t = (y + z + 2 * tp + x_p) % 2 != cb ? 2 * tp + 1 : 2 * tp;
    if (t >= s_t) { return; }

    const int s_x_cb = s_x >> 1;

    int cont = s_y * s_z * tp + s_y * z + y;

    int x = 0;

    double srcO[12];

    for ( int i=0 ; i < 12 ; i++ ){
        srcO[i] = src_b[cont * 6 * 2 + i];
    }

    double *destE = dest + (s_x_cb * s_y * s_z * t +
                            s_x_cb * s_y * z +
                            s_x_cb * y +
                            x + cb * subgrid_vol_cb) * 12 * 2;

    for (int c1 = 0; c1 < 3; c1++) {

        destE[(0 * 3 + c1) * 2 + 0] += srcO[(0 * 3 + c1) * 2 + 0] * b;
        destE[(0 * 3 + c1) * 2 + 1] += srcO[(0 * 3 + c1) * 2 + 1] * b;
        destE[(3 * 3 + c1) * 2 + 0] +=-flag * srcO[(0 * 3 + c1) * 2 + 1] * b;
        destE[(3 * 3 + c1) * 2 + 1] += flag * srcO[(0 * 3 + c1) * 2 + 0] * b;
        destE[(1 * 3 + c1) * 2 + 0] += srcO[(1 * 3 + c1) * 2 + 0] * b;
        destE[(1 * 3 + c1) * 2 + 1] += srcO[(1 * 3 + c1) * 2 + 1] * b;
        destE[(2 * 3 + c1) * 2 + 0] +=-flag * srcO[(1 * 3 + c1) * 2 + 1] * b;
        destE[(2 * 3 + c1) * 2 + 1] += flag * srcO[(1 * 3 + c1) * 2 + 0] * b;
    }
}

__global__ void ghost_y_f_abp5(double *src_f, double *dest, double *U,
                               const int v_x, const int v_y, const int v_z, const int v_t,
                               const int s_x, const int s_y, const int s_z, const int s_t,
                               const int rank, const int cb, const int flag, const double b) {

    const int s_x_cb = s_x >> 1;

    int subgrid_vol = (s_x * s_y * s_z * s_t);
    int subgrid_vol_cb = (subgrid_vol) >> 1;

    const double half = 0.5;

    int y = s_y - 1;

    int xyz = blockDim.x * blockIdx.x + threadIdx.x;

    int t =  xyz / ( s_z * s_x_cb  ) ;
    int z =  (xyz / s_x_cb ) % s_z ;
    int x =  xyz % s_x_cb ;

    double tmp[2];
    double destE[24];
    double srcO[12];
    double AE[18];

    int cont = s_x_cb * s_z * t + s_x_cb * z + x;

    for (int i = 0; i < 24; i++) {
        destE[i] = 0;
    }

    for (int i = 0; i < 12; i++) {
        srcO[i] = src_f[cont * 6 * 2 + i];
    }

    for (int i = 0; i < 18; i++) {
        AE[i] = U[(s_x_cb * s_y * s_z * t +
                   s_x_cb * s_y * z +
                   s_x_cb * y +
                   x + cb * subgrid_vol_cb) * 9 * 2 + i];
    }

    for (int c1 = 0; c1 < 3; c1++) {
        for (int c2 = 0; c2 < 3; c2++) {

            tmp[0] = + srcO[(0 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 0]
                     - srcO[(0 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 1];

            tmp[1] = + srcO[(0 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 1]
                     + srcO[(0 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 0];

            destE[(0 * 3 + c1) * 2 + 0] += tmp[0] * b;
            destE[(0 * 3 + c1) * 2 + 1] += tmp[1] * b;
            destE[(3 * 3 + c1) * 2 + 0] +=-flag * tmp[0] * b;
            destE[(3 * 3 + c1) * 2 + 1] +=-flag * tmp[1] * b;

            tmp[0] = + srcO[(1 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 0]
                     - srcO[(1 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 1];

            tmp[1] = + srcO[(1 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 1]
                     + srcO[(1 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 0];

            destE[(1 * 3 + c1) * 2 + 0] += tmp[0] * b;
            destE[(1 * 3 + c1) * 2 + 1] += tmp[1] * b;
            destE[(2 * 3 + c1) * 2 + 0] += flag * tmp[0] * b;
            destE[(2 * 3 + c1) * 2 + 1] += flag * tmp[1] * b;

        }
    }

    for (int i = 0; i < 24; i++) {
        dest[(s_x_cb * s_y * s_z * t +
              s_x_cb * s_y * z +
              s_x_cb * y +
              x + cb * subgrid_vol_cb) * 12 * 2 + i] += destE[i];
    }

}	

__global__ void ghost_y_b_abp5(double *src_b, double *dest, double *U,
                               const int v_x, const int v_y, const int v_z, const int v_t,
                      	       const int s_x, const int s_y, const int s_z, const int s_t,
                     	       const int rank, const int cb, const int flag, const double b) {

    const int s_x_cb = s_x >> 1;

    int subgrid_vol = (s_x * s_y * s_z * s_t);
    int subgrid_vol_cb = (subgrid_vol) >> 1;

    const double half = 0.5;

    int y = 0;

    int xyz = blockDim.x * blockIdx.x + threadIdx.x;

    int t =  xyz / ( s_z * s_x_cb  ) ;
    int z =  (xyz / s_x_cb ) % s_z ;
    int x =  xyz % s_x_cb ;

    int cont = s_x_cb * s_z * t + s_x_cb * z + x;

    double srcO[12];

    for ( int i=0 ; i < 12 ; i++ ){
        srcO[i] = src_b[cont * 6 * 2 + i];
    }

    double *destE = dest + (s_x_cb * s_y * s_z * t +
                            s_x_cb * s_y * z +
                            s_x_cb * y +
                            x + cb * subgrid_vol_cb) * 12 * 2;

    for (int c1 = 0; c1 < 3; c1++) {
	
        destE[(0 * 3 + c1) * 2 + 0] += srcO[(0 * 3 + c1) * 2 + 0] * b;
        destE[(0 * 3 + c1) * 2 + 1] += srcO[(0 * 3 + c1) * 2 + 1] * b;
        destE[(3 * 3 + c1) * 2 + 0] += flag * srcO[(0 * 3 + c1) * 2 + 0] * b;
        destE[(3 * 3 + c1) * 2 + 1] += flag * srcO[(0 * 3 + c1) * 2 + 1] * b;
        destE[(1 * 3 + c1) * 2 + 0] += srcO[(1 * 3 + c1) * 2 + 0] * b;
        destE[(1 * 3 + c1) * 2 + 1] += srcO[(1 * 3 + c1) * 2 + 1] * b;
        destE[(2 * 3 + c1) * 2 + 0] +=-flag * srcO[(1 * 3 + c1) * 2 + 0] * b;
        destE[(2 * 3 + c1) * 2 + 1] +=-flag * srcO[(1 * 3 + c1) * 2 + 1] * b;
    }
}	

__global__ void ghost_z_f_abp5 (double *src_f, double *dest, double *U,
                         	const int v_x, const int v_y, const int v_z, const int v_t,
                        	const int s_x, const int s_y, const int s_z, const int s_t,
                       		const int rank, const int cb, const int flag, const double b) {

    int xyz = blockDim.x * blockIdx.x + threadIdx.x;

    const int s_x_cb = s_x >> 1;

    int subgrid_vol = (s_x * s_y * s_z * s_t);
    int subgrid_vol_cb = (subgrid_vol) >> 1;

    const double half = 0.5;

    int z = s_z - 1;

    int t =  xyz / ( s_y * s_x_cb  ) ;
    int y =  (xyz / s_x_cb ) % s_y ;
    int x =  xyz % s_x_cb ;

    int cont = s_x_cb * s_y * t + s_x_cb * y + x;

    double tmp[2];
    double destE[24];
    double srcO[12];
    double AE[18];

    for (int i = 0; i < 24; i++) {
        destE[i] = 0;
    }

    for (int i = 0; i < 12; i++) {
        srcO[i] = src_f[cont * 6 * 2 + i];
    }

    for (int i = 0; i < 18; i++) {
        AE[i] = U[(s_x_cb * s_y * s_z * t +
                   s_x_cb * s_y * z +
                   s_x_cb * y +
                   x + cb * subgrid_vol_cb) * 9 * 2 + i];
    }

    for (int c1 = 0; c1 < 3; c1++) {
        for (int c2 = 0; c2 < 3; c2++) {

            tmp[0] = + srcO[(0 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 0]
                     - srcO[(0 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 1];

            tmp[1] = + srcO[(0 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 1]
                     + srcO[(0 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 0];

            destE[(0 * 3 + c1) * 2 + 0] += tmp[0] * b;
            destE[(0 * 3 + c1) * 2 + 1] += tmp[1] * b;
            destE[(2 * 3 + c1) * 2 + 0] += flag * tmp[1] * b;
            destE[(2 * 3 + c1) * 2 + 1] +=-flag * tmp[0] * b;

            tmp[0] = + srcO[(1 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 0]
                     - srcO[(1 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 1];

            tmp[1] = + srcO[(1 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 1]
                     + srcO[(1 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 0];

            destE[(1 * 3 + c1) * 2 + 0] += tmp[0] * b;
            destE[(1 * 3 + c1) * 2 + 1] += tmp[1] * b;
            destE[(3 * 3 + c1) * 2 + 0] +=-flag * tmp[1] * b;
            destE[(3 * 3 + c1) * 2 + 1] += flag * tmp[0] * b;

        }
    }

    for (int i = 0; i < 24; i++) {
        dest[(s_x_cb * s_y * s_z * t +
              s_x_cb * s_y * z +
              s_x_cb * y +
              x + cb * subgrid_vol_cb) * 12 * 2 + i] += destE[i];
    }

}

__global__ void ghost_z_b_abp5 (double *src_b, double *dest, double *U,
                          const int v_x, const int v_y, const int v_z, const int v_t,
                          const int s_x, const int s_y, const int s_z, const int s_t,
                          const int rank, const int cb, const int flag, const double b) {

    int xyz = blockDim.x * blockIdx.x + threadIdx.x;

    const int s_x_cb = s_x >> 1;

    int subgrid_vol = (s_x * s_y * s_z * s_t);
    int subgrid_vol_cb = (subgrid_vol) >> 1;

    const double half = 0.5;

    int z = 0;

    int t =  xyz / ( s_y * s_x_cb  ) ;
    int y =  (xyz / s_x_cb ) % s_y ;
    int x =  xyz % s_x_cb ;

    int cont = s_x_cb * s_y * t + s_x_cb * y + x;

    double srcO[12];

    for ( int i=0 ; i < 12 ; i++ ){
        srcO[i] = src_b[cont * 6 * 2 + i];
    }

    double *destE = dest + (s_x_cb * s_y * s_z * t +
                            s_x_cb * s_y * z +
                            s_x_cb * y +
                            x + cb * subgrid_vol_cb) * 12 * 2;

    for (int c1 = 0; c1 < 3; c1++) {
        destE[(0 * 3 + c1) * 2 + 0] += srcO[(0 * 3 + c1) * 2 + 0] * b;
        destE[(0 * 3 + c1) * 2 + 1] += srcO[(0 * 3 + c1) * 2 + 1] * b;

        destE[(2 * 3 + c1) * 2 + 0] +=-flag * srcO[(0 * 3 + c1) * 2 + 1] * b;
        destE[(2 * 3 + c1) * 2 + 1] += flag * srcO[(0 * 3 + c1) * 2 + 0] * b;

        destE[(1 * 3 + c1) * 2 + 0] += srcO[(1 * 3 + c1) * 2 + 0] * b;
        destE[(1 * 3 + c1) * 2 + 1] += srcO[(1 * 3 + c1) * 2 + 1] * b;

        destE[(3 * 3 + c1) * 2 + 0] += flag * srcO[(1 * 3 + c1) * 2 + 1] * b;
        destE[(3 * 3 + c1) * 2 + 1] +=-flag * srcO[(1 * 3 + c1) * 2 + 0] * b;

    }
}

__global__ void ghost_t_f_abp5(double *src_f, double *dest, double *U,
                          const int v_x, const int v_y, const int v_z, const int v_t,
                          const int s_x, const int s_y, const int s_z, const int s_t,
                          const int rank, const int cb, const int flag, const double b) {

    int xyz = blockDim.x * blockIdx.x + threadIdx.x;

    const int s_x_cb = s_x >> 1;

    int subgrid_vol = (s_x * s_y * s_z * s_t);
    int subgrid_vol_cb = (subgrid_vol) >> 1;

    int t = s_t - 1;

    int z = xyz / (s_y * s_x_cb);
    int y = (xyz / s_x_cb) % s_y;
    int x = xyz % s_x_cb;

    double tmp[2];
    double destE[24];
    double srcO[12];
    double AE[18];

    int cont = s_x_cb * s_y * z + s_x_cb * y + x;

    for (int i = 0; i < 24; i++) {
        destE[i] = 0;
    }

    for (int i = 0; i < 12; i++) {
        srcO[i] = src_f[cont * 6 * 2 + i];
    }

    for (int i = 0; i < 18; i++) {
        AE[i] = U[(s_x_cb * s_y * s_z * t +
                   s_x_cb * s_y * z +
                   s_x_cb * y +
                   x + cb * subgrid_vol_cb) * 9 * 2 + i];
    }

    for (int c1 = 0; c1 < 3; c1++) {
        for (int c2 = 0; c2 < 3; c2++) {

            tmp[0] = +srcO[(0 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 0]
                     - srcO[(0 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 1];

            tmp[1] = +srcO[(0 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 1]
                     + srcO[(0 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 0];

            destE[(0 * 3 + c1) * 2 + 0] += tmp[0] * b;
            destE[(0 * 3 + c1) * 2 + 1] += tmp[1] * b;
            destE[(2 * 3 + c1) * 2 + 0] += flag * tmp[0] * b;
            destE[(2 * 3 + c1) * 2 + 1] += flag * tmp[1] * b;

            tmp[0] = +srcO[(1 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 0]
                     - srcO[(1 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 1];

            tmp[1] = +srcO[(1 * 3 + c2) * 2 + 0] * AE[(c1 * 3 + c2) * 2 + 1]
                     + srcO[(1 * 3 + c2) * 2 + 1] * AE[(c1 * 3 + c2) * 2 + 0];

            destE[(1 * 3 + c1) * 2 + 0] += tmp[0] * b;
            destE[(1 * 3 + c1) * 2 + 1] += tmp[1] * b;
            destE[(3 * 3 + c1) * 2 + 0] += flag * tmp[0] * b;
            destE[(3 * 3 + c1) * 2 + 1] += flag * tmp[1] * b;

        }
    }
/*
    for (int i = 0; i < 3; i++) {
        destE[0 * 3 * 2 + i * 2 + 0] = destE[0 * 3 * 2 + i * 2 + 0] * b;
        destE[0 * 3 * 2 + i * 2 + 1] = destE[0 * 3 * 2 + i * 2 + 1] * b;
        destE[1 * 3 * 2 + i * 2 + 0] = destE[1 * 3 * 2 + i * 2 + 0] * b;
        destE[1 * 3 * 2 + i * 2 + 1] = destE[1 * 3 * 2 + i * 2 + 1] * b;

        destE[2 * 3 * 2 + i * 2 + 0] = -destE[2 * 3 * 2 + i * 2 + 0] * b;
        destE[2 * 3 * 2 + i * 2 + 1] = -destE[2 * 3 * 2 + i * 2 + 1] * b;
        destE[3 * 3 * 2 + i * 2 + 0] = -destE[3 * 3 * 2 + i * 2 + 0] * b;
        destE[3 * 3 * 2 + i * 2 + 1] = -destE[3 * 3 * 2 + i * 2 + 1] * b;
    }
*/
    for (int i = 0; i < 24; i++) {
        dest[(s_x_cb * s_y * s_z * t +
              s_x_cb * s_y * z +
              s_x_cb * y +
              x + cb * subgrid_vol_cb) * 12 * 2 + i] += destE[i];
    }
}

__global__ void ghost_t_b_abp5(double *src_b, double *dest, double *U,
                          const int v_x, const int v_y, const int v_z, const int v_t,
                          const int s_x, const int s_y, const int s_z, const int s_t,
                          const int rank, const int cb, const int flag, const double b) {

    const int s_x_cb = s_x >> 1;

    int subgrid_vol = (s_x * s_y * s_z * s_t);
    int subgrid_vol_cb = (subgrid_vol) >> 1;

    int t = 0;

    int xyz = blockDim.x * blockIdx.x + threadIdx.x;

    int z =  xyz / ( s_y * s_x_cb  ) ;
    int y =  (xyz / s_x_cb ) % s_y ;
    int x =  xyz % s_x_cb ;	

    int cont = s_x_cb * s_y * z + s_x_cb * y + x;

    double srcO[12];

    for ( int i=0 ; i < 12 ; i++ ){
        srcO[i] = src_b[cont * 6 * 2 + i];
    }

    double *destE = dest + (s_x_cb * s_y * s_z * t +
                            s_x_cb * s_y * z +
                            s_x_cb * y +
                            x + cb * subgrid_vol_cb) * 12 * 2;

    for (int c1 = 0; c1 < 3; c1++) {
        destE[(0 * 3 + c1) * 2 + 0] += srcO[(0 * 3 + c1) * 2 + 0] * b;
        destE[(0 * 3 + c1) * 2 + 1] += srcO[(0 * 3 + c1) * 2 + 1] * b;
        destE[(2 * 3 + c1) * 2 + 0] += -flag * srcO[(0 * 3 + c1) * 2 + 0] * b;
        destE[(2 * 3 + c1) * 2 + 1] += -flag * srcO[(0 * 3 + c1) * 2 + 1] * b;
        destE[(1 * 3 + c1) * 2 + 0] += srcO[(1 * 3 + c1) * 2 + 0] * b;
        destE[(1 * 3 + c1) * 2 + 1] += srcO[(1 * 3 + c1) * 2 + 1] * b;
        destE[(3 * 3 + c1) * 2 + 0] += -flag * srcO[(1 * 3 + c1) * 2 + 0] * b;
        destE[(3 * 3 + c1) * 2 + 1] += -flag * srcO[(1 * 3 + c1) * 2 + 1] * b;
    }
}

__global__ void Dslash_d(double *src, double *dest,
                         const int s_x, const int s_y, const int s_z, const int s_t,
                         const double mass, const int cb) {

    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;
    int z = blockDim.z * blockIdx.z + threadIdx.z;

    const double a = 4.0;
    int subgrid_vol = (s_x * s_y * s_z * s_t);
    int subgrid_vol_cb = (subgrid_vol) >> 1;

    int s_x_cb = s_x >> 1;

    for (int t = 0; t < s_t; t++) {

        double *src_d = src + (s_x_cb * s_y * s_z * t +
                               s_x_cb * s_y * z +
                               s_x_cb * y +
                               x + cb * subgrid_vol_cb) * 12 * 2;

        double *dest_d = dest + (s_x_cb * s_y * s_z * t +
                                 s_x_cb * s_y * z +
                                 s_x_cb * y +
                                 x + cb * subgrid_vol_cb) * 12 * 2;

        for (int i = 0; i < 24; i++) {
            dest_d[i] = (a + mass) * src_d[i];
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const int Dx=2;
const int Dy=8;
const int Dz=8;
const int Dt=4;

inline int test_2(double *src_g, double *dest_g, double *U_x_g, double *U_y_g, double *U_z_g, double *U_t_g,
           const int v_x, const int v_y, const int v_z, const int v_t,
           const int s_x, const int s_y, const int s_z, const int s_t,
           const int cb, const int flag) {

//hipDeviceProp_t devProp;
//hipGetDeviceProperties(&devProp, 0);

//std::cout << "Device name " << devProp.name << std::endl;

//    int rank =0;
//    int s_x_cb = s_x >> 1;	


    int N_sub[4] = {v_x / s_x, v_y / s_y, v_z / s_z, v_t / s_t};

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Request req[8 * size];
    MPI_Request reqr[8 * size];
    MPI_Status sta[8 * size];


    int site_x_f[4] = {(rank + 1) % N_sub[0],
                       (rank / N_sub[0]) % N_sub[1],
                       (rank / (N_sub[1] * N_sub[0])) % N_sub[2],
                       rank / (N_sub[2] * N_sub[1] * N_sub[0])};

    int site_x_b[4] = {(rank - 1 + N_sub[0]) % N_sub[0],
                       (rank / N_sub[0]) % N_sub[1],
                       (rank / (N_sub[1] * N_sub[0])) % N_sub[2],
                       rank / (N_sub[2] * N_sub[1] * N_sub[0])};

    const int nodenum_x_b = get_nodenum(site_x_b, N_sub, 4);
    const int nodenum_x_f = get_nodenum(site_x_f, N_sub, 4);

    int site_y_f[4] = {(rank % N_sub[0]),
                       (rank / N_sub[0] + 1) % N_sub[1],
                       (rank / (N_sub[1] * N_sub[0])) % N_sub[2],
                       rank / (N_sub[2] * N_sub[1] * N_sub[0])};

    int site_y_b[4] = {(rank % N_sub[0]),
                       (rank / N_sub[0] - 1 + N_sub[1]) % N_sub[1],
                       (rank / (N_sub[1] * N_sub[0])) % N_sub[2],
                       rank / (N_sub[2] * N_sub[1] * N_sub[0])};

    const int nodenum_y_b = get_nodenum(site_y_b, N_sub, 4);
    const int nodenum_y_f = get_nodenum(site_y_f, N_sub, 4);

    int site_z_f[4] = {(rank % N_sub[0]),
                       (rank / N_sub[0]) % N_sub[1],
                       (rank / (N_sub[1] * N_sub[0]) + 1) % N_sub[2],
                       rank / (N_sub[2] * N_sub[1] * N_sub[0])};

    int site_z_b[4] = {(rank % N_sub[0]),
                       (rank / N_sub[0]) % N_sub[1],
                       (rank / (N_sub[1] * N_sub[0]) - 1 + N_sub[2]) % N_sub[2],
                       rank / (N_sub[2] * N_sub[1] * N_sub[0])};

    const int nodenum_z_b = get_nodenum(site_z_b, N_sub, 4);
    const int nodenum_z_f = get_nodenum(site_z_f, N_sub, 4);

    int site_t_f[4] = {(rank % N_sub[0]),
                       (rank / N_sub[0]) % N_sub[1],
                       (rank / (N_sub[1] * N_sub[0])) % N_sub[2],
                       (rank / (N_sub[2] * N_sub[1] * N_sub[0]) + 1) % N_sub[3]};

    int site_t_b[4] = {(rank % N_sub[0]),
                       (rank / N_sub[0]) % N_sub[1],
                       (rank / (N_sub[1] * N_sub[0])) % N_sub[2],
                       (rank / (N_sub[2] * N_sub[1] * N_sub[0]) - 1 + N_sub[3]) %  N_sub[3]};

    const int nodenum_t_b = get_nodenum(site_t_b, N_sub, 4);
    const int nodenum_t_f = get_nodenum(site_t_f, N_sub, 4);	

    int s_x_cb = s_x >> 1;

    int len_x = (s_y * s_z * s_t + cb) >> 1;

    double *resv_x_f ;
    double *send_x_b ;
    double *resv_x_b ;
    double *send_x_f ;


//    resv_x_f = new double[len_x * 6 * 2];
//    send_x_b = new double[len_x * 6 * 2];
//    resv_x_b = new double[len_x * 6 * 2];
//    send_x_f = new double[len_x * 6 * 2];



    if (N_sub[0] != 1) {

        resv_x_f = new double[len_x * 6 * 2];
        send_x_b = new double[len_x * 6 * 2];
        resv_x_b = new double[len_x * 6 * 2];
        send_x_f = new double[len_x * 6 * 2];

        double *tran_f;
        double *tran_b;
        int size_T = len_x * 6 * 2;

        hipMalloc((void **) &tran_f, size_T * sizeof(double));
        hipMemset(tran_f, 0, size_T * sizeof(double));

        hipMalloc((void **) &tran_b, size_T * sizeof(double));
        hipMemset(tran_b, 0, size_T * sizeof(double));

        transfer_x_f<<<dim3(s_y / Dy, s_z / Dz, s_t / Dt / 2), dim3(Dy, Dz, Dt) >>>
                                                        (src_g, tran_b, U_x_g, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);
        transfer_x_b<<<dim3(s_y / Dy, s_z / Dz, s_t / Dt / 2), dim3(Dy, Dz, Dt) >>>
                                                        (src_g, tran_f, U_x_g, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);

        hipMemcpy(send_x_b, tran_b, size_T * sizeof(double), hipMemcpyDeviceToHost);
        hipMemcpy(send_x_f, tran_f, size_T * sizeof(double), hipMemcpyDeviceToHost);

        hipFree(tran_f);
        hipFree(tran_b);

        MPI_Isend(send_x_b, len_x * 6 * 2, MPI_DOUBLE, nodenum_x_b, 8 * rank, MPI_COMM_WORLD, &req[8 * rank]);
        MPI_Irecv(resv_x_f, len_x * 6 * 2, MPI_DOUBLE, nodenum_x_f, 8 * nodenum_x_f, MPI_COMM_WORLD,
                  &reqr[8 * nodenum_x_f]);

        MPI_Isend(send_x_f, len_x * 6 * 2, MPI_DOUBLE, nodenum_x_f, 8 * rank + 1, MPI_COMM_WORLD, &req[8 * rank + 1]);
        MPI_Irecv(resv_x_b, len_x * 6 * 2, MPI_DOUBLE, nodenum_x_b, 8 * nodenum_x_b + 1, MPI_COMM_WORLD,
                  &reqr[8 * nodenum_x_b + 1]);

    }

    int len_y = s_x_cb * s_z * s_t;



    double *resv_y_f;
    double *send_y_b;
    double *resv_y_b;
    double *send_y_f;

//    resv_y_f = new double[len_y * 6 * 2];
//    send_y_b = new double[len_y * 6 * 2];
//    resv_y_b = new double[len_y * 6 * 2];
//    send_y_f = new double[len_y * 6 * 2];

    if (N_sub[1] != 1) {

        resv_y_f = new double[len_y * 6 * 2];
        send_y_b = new double[len_y * 6 * 2];
        resv_y_b = new double[len_y * 6 * 2];
        send_y_f = new double[len_y * 6 * 2];

        double *tran_f;
        double *tran_b;
        int size_T = len_y * 6 * 2;

        hipMalloc((void **) &tran_f, size_T * sizeof(double));
        hipMemset(tran_f, 0, size_T * sizeof(double));

        hipMalloc((void **) &tran_b, size_T * sizeof(double));
        hipMemset(tran_b, 0, size_T * sizeof(double));

        transfer_y_f<<<dim3(s_x_cb / Dx, s_z / Dz, s_t / Dt), dim3(Dx, Dz, Dt) >>>(src_g, tran_b, U_y_g, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);
        transfer_y_b<<<dim3(s_x_cb / Dx, s_z / Dz, s_t / Dt), dim3(Dx, Dz, Dt) >>>(src_g, tran_f, U_y_g, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);

        hipMemcpy(send_y_b, tran_b, size_T * sizeof(double), hipMemcpyDeviceToHost);
        hipMemcpy(send_y_f, tran_f, size_T * sizeof(double), hipMemcpyDeviceToHost);

        hipFree(tran_f);
        hipFree(tran_b);

        MPI_Isend(send_y_b, len_y * 6 * 2, MPI_DOUBLE, nodenum_y_b, 8 * rank + 2, MPI_COMM_WORLD, &req[8 * rank + 2]);
        MPI_Irecv(resv_y_f, len_y * 6 * 2, MPI_DOUBLE, nodenum_y_f, 8 * nodenum_y_f + 2, MPI_COMM_WORLD,
                  &reqr[8 * nodenum_y_f + 2]);

        MPI_Isend(send_y_f, len_y * 6 * 2, MPI_DOUBLE, nodenum_y_f, 8 * rank + 3, MPI_COMM_WORLD, &req[8 * rank + 3]);
        MPI_Irecv(resv_y_b, len_y * 6 * 2, MPI_DOUBLE, nodenum_y_b, 8 * nodenum_y_b + 3, MPI_COMM_WORLD,
                  &reqr[8 * nodenum_y_b + 3]);

    }

    int len_z = s_x_cb * s_y * s_t;

    double *resv_z_f;
    double *send_z_b;
    double *resv_z_b;
    double *send_z_f;

//    resv_z_f = new double[len_z * 6 * 2];
//    send_z_b = new double[len_z * 6 * 2];
//    resv_z_b = new double[len_z * 6 * 2];
//    send_z_f = new double[len_z * 6 * 2];

    if (N_sub[2] != 1) {

        resv_z_f = new double[len_z * 6 * 2];
        send_z_b = new double[len_z * 6 * 2];
        resv_z_b = new double[len_z * 6 * 2];
        send_z_f = new double[len_z * 6 * 2];

        double *tran_f;
        double *tran_b;
        int size_T = len_z * 6 * 2;

        hipMalloc((void **) &tran_f, size_T * sizeof(double));
        hipMemset(tran_f, 0, size_T * sizeof(double));

        hipMalloc((void **) &tran_b, size_T * sizeof(double));
        hipMemset(tran_b, 0, size_T * sizeof(double));

        transfer_z_f<<<dim3(s_x_cb / Dx, s_y / Dy, s_t / Dt), dim3(Dx, Dy, Dt) >>>
                                                       (src_g, tran_b, U_z_g, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);
        transfer_z_b<<<dim3(s_x_cb / Dx, s_y / Dy, s_t / Dt), dim3(Dx, Dy, Dt) >>>
                                                       (src_g, tran_f, U_z_g, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);

        hipMemcpy(send_z_b, tran_b, size_T * sizeof(double), hipMemcpyDeviceToHost);
        hipMemcpy(send_z_f, tran_f, size_T * sizeof(double), hipMemcpyDeviceToHost);

        hipFree(tran_f);
        hipFree(tran_b);

        MPI_Isend(send_z_b, len_z * 6 * 2, MPI_DOUBLE, nodenum_z_b, 8 * rank + 4, MPI_COMM_WORLD, &req[8 * rank + 4]);
        MPI_Irecv(resv_z_f, len_z * 6 * 2, MPI_DOUBLE, nodenum_z_f, 8 * nodenum_z_f + 4, MPI_COMM_WORLD,
                  &reqr[8 * nodenum_z_f + 4]);

        MPI_Isend(send_z_f, len_z * 6 * 2, MPI_DOUBLE, nodenum_z_f, 8 * rank + 5, MPI_COMM_WORLD, &req[8 * rank + 5]);
        MPI_Irecv(resv_z_b, len_z * 6 * 2, MPI_DOUBLE, nodenum_z_b, 8 * nodenum_z_b + 5, MPI_COMM_WORLD,
                  &reqr[8 * nodenum_z_b + 5]);

    }	

    int len_t = s_x_cb * s_y * s_z;

    double *resv_t_f;
    double *send_t_b;
    double *resv_t_b;
    double *send_t_f;	

//    resv_t_f = new double[len_t * 6 * 2];
//    send_t_b = new double[len_t * 6 * 2];
//    resv_t_b = new double[len_t * 6 * 2];
//    send_t_f = new double[len_t * 6 * 2];

    if (N_sub[3] != 1) {

        resv_t_f = new double[len_t * 6 * 2];
        send_t_b = new double[len_t * 6 * 2];
        resv_t_b = new double[len_t * 6 * 2];
        send_t_f = new double[len_t * 6 * 2];

        double *tran_f;
        double *tran_b;
        int size_T = len_t * 6 * 2;

        hipMalloc((void **) &tran_f, size_T * sizeof(double));
//        hipMemset(tran_f, 0, size_T * sizeof(double));

        hipMalloc((void **) &tran_b, size_T * sizeof(double));
//        hipMemset(tran_b, 0, size_T * sizeof(double));

//        transfer_t_f<<<dim3(s_x_cb / Dx, s_y / Dy, s_z / Dz), dim3(Dx, Dy, Dz) >>>
//                                                       (src_g, tran_b, U_t_g, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);


        transfer_t_f<<<s_x_cb * s_y * s_z / Dx / Dy / Dz, Dx * Dy * Dz >>>
                                                          (src_g, tran_b, U_t_g, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);


//        transfer_t_b<<<dim3(s_x_cb / Dx, s_y / Dy, s_z / Dz), dim3(Dx, Dy, Dz) >>>
//                                                       (src_g, tran_f, U_t_g, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);


        transfer_t_b<<<s_x_cb * s_y * s_z / Dx / Dy / Dz, Dx * Dy * Dz >>>
                                                          (src_g, tran_f, U_t_g, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);

        hipMemcpy(send_t_b, tran_b, size_T * sizeof(double), hipMemcpyDeviceToHost);
        hipMemcpy(send_t_f, tran_f, size_T * sizeof(double), hipMemcpyDeviceToHost);

        hipFree(tran_f);
        hipFree(tran_b);

        MPI_Isend(send_t_b, len_t * 6 * 2, MPI_DOUBLE, nodenum_t_b, 8 * rank + 6, MPI_COMM_WORLD, &req[8 * rank + 6]);
        MPI_Irecv(resv_t_f, len_t * 6 * 2, MPI_DOUBLE, nodenum_t_f, 8 * nodenum_t_f + 6, MPI_COMM_WORLD,
                  &reqr[8 * nodenum_t_f + 6]);

        MPI_Isend(send_t_f, len_t * 6 * 2, MPI_DOUBLE, nodenum_t_f, 8 * rank + 7, MPI_COMM_WORLD, &req[8 * rank + 7]);
        MPI_Irecv(resv_t_b, len_t * 6 * 2, MPI_DOUBLE, nodenum_t_b, 8 * nodenum_t_b + 7, MPI_COMM_WORLD,
                  &reqr[8 * nodenum_t_b + 7]);

    }

/* 	
    float eventMs = 1.0f;
    hipEvent_t start, stop;
    hipEventCreate(&start);
    hipEventCreate(&stop);
    hipEventRecord(start, 0);
*/


    main_xyzt<<< s_x_cb * s_y * s_z * s_t / Dx / Dy / Dz , Dx * Dy * Dz >>>
                                             (src_g, dest_g, U_x_g, U_y_g, U_z_g, U_t_g, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);	

/*
    hipEventRecord(stop, 0);
    hipEventSynchronize(stop);
    hipEventElapsedTime(&eventMs, start, stop);
    printf("main_xyzt time taken  = %6.3fms\n", eventMs);    	
*/

    if (N_sub[0] != 1) {

        MPI_Wait(&reqr[8 * nodenum_x_f], &sta[8 * nodenum_x_f]);
        MPI_Wait(&reqr[8 * nodenum_x_b + 1], &sta[8 * nodenum_x_b + 1]);

        double *tran_f;
        double *tran_b;
        int size_T = len_x * 6 * 2;
        hipMalloc((void **) &tran_f, size_T * sizeof(double));
        hipMalloc((void **) &tran_b, size_T * sizeof(double));

        hipMemcpy(tran_f, resv_x_f, size_T * sizeof(double), hipMemcpyHostToDevice);
        hipMemcpy(tran_b, resv_x_b, size_T * sizeof(double), hipMemcpyHostToDevice);

        ghost_x_f<<<dim3(s_y / Dy, s_z / Dz, s_t / Dt / 2), dim3(Dy, Dz, Dt) >>>
                                                     (tran_f, dest_g, U_x_g, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);
        ghost_x_b<<<dim3(s_y / Dy, s_z / Dz, s_t / Dt / 2), dim3(Dy, Dz, Dt) >>>
                                                     (tran_b, dest_g, U_x_g, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);

        hipFree(tran_f);
        hipFree(tran_b);
    }


    if (N_sub[1] != 1) {

        MPI_Wait(&reqr[8 * nodenum_y_f + 2], &sta[8 * nodenum_y_f + 2]);
        MPI_Wait(&reqr[8 * nodenum_y_b + 3], &sta[8 * nodenum_y_b + 3]);

        double *tran_f;
        double *tran_b;
        int size_T = len_y * 6 * 2;
        hipMalloc((void **) &tran_f, size_T * sizeof(double));
        hipMalloc((void **) &tran_b, size_T * sizeof(double));

        hipMemcpy(tran_f, resv_y_f, size_T * sizeof(double), hipMemcpyHostToDevice);
        hipMemcpy(tran_b, resv_y_b, size_T * sizeof(double), hipMemcpyHostToDevice);

        ghost_y_f<<<dim3(s_x_cb / Dx, s_z / Dz, s_t / Dt), dim3(Dx, Dz, Dt) >>>
                                                    (tran_f, dest_g, U_y_g, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);
        ghost_y_b<<<dim3(s_x_cb / Dx, s_z / Dz, s_t / Dt), dim3(Dx, Dz, Dt) >>>
                                                    (tran_b, dest_g, U_y_g, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);

        hipFree(tran_f);
        hipFree(tran_b);
    }

    if (N_sub[2] != 1) {

        MPI_Wait(&reqr[8 * nodenum_z_f + 4], &sta[8 * nodenum_z_f + 4]);
        MPI_Wait(&reqr[8 * nodenum_z_b + 5], &sta[8 * nodenum_z_b + 5]);

        double *tran_f;
        double *tran_b;
        int size_T = len_z * 6 * 2;
        hipMalloc((void **) &tran_f, size_T * sizeof(double));
        hipMalloc((void **) &tran_b, size_T * sizeof(double));

        hipMemcpy(tran_f, resv_z_f, size_T * sizeof(double), hipMemcpyHostToDevice);
        hipMemcpy(tran_b, resv_z_b, size_T * sizeof(double), hipMemcpyHostToDevice);

        ghost_z_f<<<dim3(s_x_cb / Dx, s_y / Dy, s_t / Dt), dim3(Dx, Dy, Dt) >>>
                                                    (tran_f, dest_g, U_z_g, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);
        ghost_z_b<<<dim3(s_x_cb / Dx, s_y / Dy, s_t / Dt), dim3(Dx, Dy, Dt) >>>
                                                    (tran_b, dest_g, U_z_g, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);

        hipFree(tran_f);
        hipFree(tran_b);
    }

    if (N_sub[3] != 1) {

        MPI_Wait(&reqr[8 * nodenum_t_f + 6], &sta[8 * nodenum_t_f + 6]);
        MPI_Wait(&reqr[8 * nodenum_t_b + 7], &sta[8 * nodenum_t_b + 7]);

        double *tran_f;
        double *tran_b;
        int size_T = len_t * 6 * 2;
        hipMalloc((void **) &tran_f, size_T * sizeof(double));
        hipMalloc((void **) &tran_b, size_T * sizeof(double));

        hipMemcpy(tran_f, resv_t_f, size_T * sizeof(double), hipMemcpyHostToDevice);
        hipMemcpy(tran_b, resv_t_b, size_T * sizeof(double), hipMemcpyHostToDevice);

//        ghost_t_f<<<dim3(s_x_cb / Dx, s_y / Dy, s_z / Dz), dim3(Dx, Dy, Dz) >>>
//                                                    (tran_f, dest_g, U_t_g, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);



        ghost_t_f<<<s_x_cb * s_y * s_z / Dx / Dy / Dz, Dx * Dy * Dz >>>
                                                       (tran_f, dest_g, U_t_g, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);

//        ghost_t_b<<<dim3(s_x_cb / Dx, s_y / Dy, s_z / Dz), dim3(Dx, Dy, Dz) >>>
//                                                    (tran_b, dest_g, U_t_g, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);


        ghost_t_b<<<s_x_cb * s_y * s_z / Dx / Dy / Dz, Dx * Dy * Dz >>>
                                                       (tran_b, dest_g, U_t_g, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);


        hipFree(tran_f);
        hipFree(tran_b);

    }


    MPI_Barrier(MPI_COMM_WORLD);

    if (N_sub[0] != 1) {
        delete[] resv_x_f;
        delete[] send_x_b;
        delete[] resv_x_b;
        delete[] send_x_f;
    }

    if (N_sub[1] != 1) {
        delete[] resv_y_f;
        delete[] send_y_b;
        delete[] resv_y_b;
        delete[] send_y_f;
    }

    if (N_sub[2] != 1) {
        delete[] resv_z_f;
        delete[] send_z_b;
        delete[] resv_z_b;
        delete[] send_z_f;
    }

    if (N_sub[3] != 1) {
        delete[] resv_t_f;
        delete[] send_t_b;
        delete[] resv_t_b;
        delete[] send_t_f;
    }
//    printf("PASSED!\n");

    return 0;
}

int test(double *src, double *dest, double *U_x, double *U_y, double *U_z, double *U_t,
         const int v_x, const int v_y, const int v_z, const int v_t,
         const int s_x, const int s_y, const int s_z, const int s_t,
         const int cb, const int flag) {

    hipDeviceProp_t devProp;
    hipGetDeviceProperties(&devProp, 0);

//    std::cout << "Device name " << devProp.name << std::endl;

    double *src_g;
    double *dest_g;
    double *U_x_g;
    double *U_y_g;
    double *U_z_g;
    double *U_t_g;

    int size_f = s_x * s_y * s_z * s_t * 12 * 2;
    int size_u = s_x * s_y * s_z * s_t * 9 * 2;
//    allocate the memory on the device side
    hipMalloc((void **) &src_g, size_f * sizeof(double));
    hipMalloc((void **) &dest_g, size_f * sizeof(double));
    hipMalloc((void **) &U_x_g, size_u * sizeof(double));
    hipMalloc((void **) &U_y_g, size_u * sizeof(double));
    hipMalloc((void **) &U_z_g, size_u * sizeof(double));
    hipMalloc((void **) &U_t_g, size_u * sizeof(double));

    hipMemcpy(src_g, src, size_f * sizeof(double), hipMemcpyHostToDevice);
    hipMemcpy(U_x_g, U_x, size_u * sizeof(double), hipMemcpyHostToDevice);
    hipMemcpy(U_y_g, U_y, size_u * sizeof(double), hipMemcpyHostToDevice);
    hipMemcpy(U_z_g, U_z, size_u * sizeof(double), hipMemcpyHostToDevice);
    hipMemcpy(U_t_g, U_t, size_u * sizeof(double), hipMemcpyHostToDevice);
    hipMemset(dest_g, 0, size_f * sizeof(double));

    test_2(src_g, dest_g, U_x_g, U_y_g, U_z_g, U_t_g, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, cb, flag);

//    hipDeviceSynchronize();

    hipMemcpy(dest, dest_g, size_f * sizeof(double), hipMemcpyDeviceToHost);

    hipFree(src_g);
    hipFree(dest_g);
    hipFree(U_x_g);
    hipFree(U_y_g);
    hipFree(U_z_g);
    hipFree(U_t_g);

//    printf("test PASSED!\n");
    return 0;
}

__global__ void psi_g5(double *src, const int s_x, const int s_y, const int s_z, const int s_t) {

    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;
    int z = blockDim.z * blockIdx.z + threadIdx.z;

    for (int t = 0; t < s_t; t++) {
        double *src_d = src + (s_x * s_y * s_z * t +
                               s_x * s_y * z +
                               s_x * y +
                               x) * 12 * 2;


        for (int i = 0; i < 3; i++) {
//            src_d[0 * 3 * 2 + i * 2 + 0] = src_d[0 * 3 * 2 + i * 2 + 0];
//            src_d[0 * 3 * 2 + i * 2 + 1] = src_d[0 * 3 * 2 + i * 2 + 1];
//            src_d[1 * 3 * 2 + i * 2 + 0] = src_d[1 * 3 * 2 + i * 2 + 0];
//            src_d[1 * 3 * 2 + i * 2 + 1] = src_d[1 * 3 * 2 + i * 2 + 1];
            src_d[2 * 3 * 2 + i * 2 + 0] = -src_d[2 * 3 * 2 + i * 2 + 0];
            src_d[2 * 3 * 2 + i * 2 + 1] = -src_d[2 * 3 * 2 + i * 2 + 1];
            src_d[3 * 3 * 2 + i * 2 + 0] = -src_d[3 * 3 * 2 + i * 2 + 0];
            src_d[3 * 3 * 2 + i * 2 + 1] = -src_d[3 * 3 * 2 + i * 2 + 1];
        }
    }
}

__device__ void psi_g5_d(double *src) {

    for (int i = 0; i < 3; i++) {
        src[2 * 3 * 2 + i * 2 + 0] = -src[2 * 3 * 2 + i * 2 + 0];
        src[2 * 3 * 2 + i * 2 + 1] = -src[2 * 3 * 2 + i * 2 + 1];
        src[3 * 3 * 2 + i * 2 + 0] = -src[3 * 3 * 2 + i * 2 + 0];
        src[3 * 3 * 2 + i * 2 + 1] = -src[3 * 3 * 2 + i * 2 + 1];
    }
}



__global__ void overlapLinop(double *out, double *in, const double k0, const double k1, const double k2) {

    int i = blockDim.x * blockIdx.x + threadIdx.x;

    int j = i % 24;

    if (j < 12) {
        out[i] = k0 * in[i] + (k2 + k1) * out[i];
    } else {
        out[i] = k0 * in[i] + (k2 - k1) * out[i];
    }

}

__global__ void axpbyz_g(const double a, double *x_i, const double b, double *y_i, double *z_i,
                         const int s_x, const int s_y, const int s_z, const int s_t) {

    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;
    int z = blockDim.z * blockIdx.z + threadIdx.z;

    for (int t = 0; t < s_t; t++) {
        int pos = s_x * s_y * s_z * t + s_x * s_y * z + s_x * y + x;
        double *x_t = x_i + pos * 12 * 2;
        double *y_t = y_i + pos * 12 * 2;
        double *z_t = z_i + pos * 12 * 2;
        for (int i = 0; i < 12; i++) {
            z_t[i * 2 + 0] = a * x_t[i * 2 + 0] + b * y_t[i * 2 + 0];
            z_t[i * 2 + 1] = a * x_t[i * 2 + 1] + b * y_t[i * 2 + 1];
        }
    }
}

__global__ void axpbyz_g2(const double a, double *x, const double b, double *y, double *z) {

    int i = blockDim.x * blockIdx.x + threadIdx.x;

    z[i] = a * x[i] + b * y[i];
}

__global__ void axpbyczw_g2(const double a, double *x, const double b, double *y, const double c, double *z, double *w) {

    int i = blockDim.x * blockIdx.x + threadIdx.x;

    w[i] = a * x[i] + b * y[i] + c * z[i];
}


__device__ void axpbyz_g_d(const double a, double *x, const double b, double *y, double *z) {

    for (int i = 0; i < 12; i++) {
        z[i * 2 + 0] = a * x[i * 2 + 0] + b * y[i * 2 + 0];
        z[i * 2 + 1] = a * x[i * 2 + 1] + b * y[i * 2 + 1];
    }
}

__global__ void caxpbyz_g(double *a, double *x_i, double *b, double *y_i, double *z_i,
                          const int s_x, const int s_y, const int s_z, const int s_t) {

    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;
    int z = blockDim.z * blockIdx.z + threadIdx.z;

    hipDoubleComplex *a_t = (hipDoubleComplex *) a;
    hipDoubleComplex *b_t = (hipDoubleComplex *) b;

    for (int t = 0; t < s_t; t++) {

        int pos = s_x * s_y * s_z * t + s_x * s_y * z + s_x * y + x;

        hipDoubleComplex *x_t = (hipDoubleComplex *) (x_i + pos * 12 * 2);
        hipDoubleComplex *y_t = (hipDoubleComplex *) (y_i + pos * 12 * 2);
        hipDoubleComplex *z_t = (hipDoubleComplex *) (z_i + pos * 12 * 2);

        for (int i = 0; i < 12; i++) {
            z_t[i] = (*a_t) * x_t[i] + (*b_t) * y_t[i];
        }
    }
}

__global__ void caxmbyz_g(double *a, double *x_i, double *b, double *y_i, double *z_i,
                          const int s_x, const int s_y, const int s_z, const int s_t) {

    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;
    int z = blockDim.z * blockIdx.z + threadIdx.z;

    hipDoubleComplex *a_t = (hipDoubleComplex *) a;
    hipDoubleComplex *b_t = (hipDoubleComplex *) b;

    for (int t = 0; t < s_t; t++) {

        int pos = s_x * s_y * s_z * t + s_x * s_y * z + s_x * y + x;

        hipDoubleComplex *x_t = (hipDoubleComplex *) (x_i + pos * 12 * 2);
        hipDoubleComplex *y_t = (hipDoubleComplex *) (y_i + pos * 12 * 2);
        hipDoubleComplex *z_t = (hipDoubleComplex *) (z_i + pos * 12 * 2);

        for (int i = 0; i < 12; i++) {
            z_t[i] = (*a_t) * x_t[i] - (*b_t) * y_t[i];
        }
    }
}

__global__ void cxpbyz_g(double *x_i, double *b, double *y_i, double *z_i,
                         const int s_x, const int s_y, const int s_z, const int s_t) {

    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;
    int z = blockDim.z * blockIdx.z + threadIdx.z;

    hipDoubleComplex *b_t = (hipDoubleComplex *) b;

    for (int t = 0; t < s_t; t++) {

        int pos = s_x * s_y * s_z * t + s_x * s_y * z + s_x * y + x;

        hipDoubleComplex *x_t = (hipDoubleComplex *) (x_i + pos * 12 * 2);
        hipDoubleComplex *y_t = (hipDoubleComplex *) (y_i + pos * 12 * 2);
        hipDoubleComplex *z_t = (hipDoubleComplex *) (z_i + pos * 12 * 2);

        for (int i = 0; i < 12; i++) {
            z_t[i] = x_t[i] + (*b_t) * y_t[i];
        }
    }
}

__global__ void cxpbyz_g2(double *x_i, double *b, double *y_i, double *z_i) {

    int x = blockDim.x * blockIdx.x + threadIdx.x;
	    	
            z_i[2 * x + 0] = x_i[2 * x + 0] + b[0] * y_i[2 * x + 0] - b[1] * y_i[2 * x + 1];
	    z_i[2 * x + 1] = x_i[2 * x + 1] + b[0] * y_i[2 * x + 1] + b[1] * y_i[2 * x + 0];

}

__global__ void cxpbyx_v1_g2(double *x_i, double *b, double *y_i, const int size, const int N) {

    int x = blockDim.x * blockIdx.x + threadIdx.x;

            double x_t_r = x_i[2 * x + 0];
            double x_t_i = x_i[2 * x + 1];


            for ( int i=0; i < N; i++ ){

            double b_t_r = b[ i * 2 + 0 ];
            double b_t_i = b[ i * 2 + 1 ];
            double y_t_r = y_i[ i * size * 2 + x * 2 + 0];
            double y_t_i = y_i[ i * size * 2 + x * 2 + 1];

                    x_t_r +=  + b_t_r * y_t_r - b_t_i * y_t_i;
                    x_t_i +=  + b_t_r * y_t_i + b_t_i * y_t_r;
            }

            x_i[2 * x + 0] = x_t_r;
            x_i[2 * x + 1] = x_t_i;

}


__global__ void cxmbyz_g(double *x_i, double *b, double *y_i, double *z_i,
                         const int s_x, const int s_y, const int s_z, const int s_t) {

    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;
    int z = blockDim.z * blockIdx.z + threadIdx.z;


    hipDoubleComplex *b_t = (hipDoubleComplex *) b;

    for (int t = 0; t < s_t; t++) {

        int pos = s_x * s_y * s_z * t + s_x * s_y * z + s_x * y + x;

        hipDoubleComplex *x_t = (hipDoubleComplex *) (x_i + pos * 12 * 2);
        hipDoubleComplex *y_t = (hipDoubleComplex *) (y_i + pos * 12 * 2);
        hipDoubleComplex *z_t = (hipDoubleComplex *) (z_i + pos * 12 * 2);

        for (int i = 0; i < 12; i++) {
            z_t[i] = x_t[i] - (*b_t) * y_t[i];
        }
    }
}

__global__ void cxmbyx_v1_g2(double *x_i, double *b, double *y_i, const int size, const int N) {

    int x = blockDim.x * blockIdx.x + threadIdx.x;

	    double x_t_r = x_i[2 * x + 0];
	    double x_t_i = x_i[2 * x + 1];
	    	

	    for ( int i=0; i < N; i++ ){

	    double b_t_r = b[ i * 2 + 0 ];
	    double b_t_i = b[ i * 2 + 1 ];
	    double y_t_r = y_i[ i * size * 2 + x * 2 + 0];
	    double y_t_i = y_i[ i * size * 2 + x * 2 + 1];			
	
	            x_t_r +=  - b_t_r * y_t_r + b_t_i * y_t_i;
	            x_t_i +=  - b_t_r * y_t_i - b_t_i * y_t_r;
            } 

	    x_i[2 * x + 0] = x_t_r;
	    x_i[2 * x + 1] = x_t_i; 	
	    
}

__global__ void cxmbyz_g2(double *x_i, double *b, double *y_i, double *z_i) {

    int x = blockDim.x * blockIdx.x + threadIdx.x;

            z_i[2 * x + 0] = x_i[2 * x + 0] - b[0] * y_i[2 * x + 0] + b[1] * y_i[2 * x + 1];
            z_i[2 * x + 1] = x_i[2 * x + 1] - b[0] * y_i[2 * x + 1] - b[1] * y_i[2 * x + 0];

}


void caxpbyz(std::complex<double> a, std::complex<double> *x_i,
             std::complex<double> b, std::complex<double> *y_i,
             std::complex<double> *z_i,
             const int s_x, const int s_y, const int s_z, const int s_t) {

    for (int t = 0; t < s_t; t++) {
        for (int z = 0; z < s_z; z++) {
            for (int y = 0; y < s_y; y++) {
                for (int x = 0; x < s_x; x++) {

                    int pos = s_x * s_y * s_z * t + s_x * s_y * z + s_x * y + x;

                    std::complex<double> *x_t = x_i + pos * 12;
                    std::complex<double> *y_t = y_i + pos * 12;
                    std::complex<double> *z_t = z_i + pos * 12;

                    for (int i = 0; i < 12; i++) {
                        z_t[i] = a * x_t[i] + b * y_t[i];
                    }
                }
            }
        }
    }
}

void caxz(std::complex<double> a, std::complex<double> *x_i, std::complex<double> *z_i,
          const int s_x, const int s_y, const int s_z, const int s_t) {

    for (int t = 0; t < s_t; t++) {
        for (int z = 0; z < s_z; z++) {
            for (int y = 0; y < s_y; y++) {
                for (int x = 0; x < s_x; x++) {

                    int pos = s_x * s_y * s_z * t + s_x * s_y * z + s_x * y + x;
                    std::complex<double> *x_t = x_i + pos * 12;
                    std::complex<double> *z_t = z_i + pos * 12;
                    for (int i = 0; i < 12; i++) {
                        z_t[i] = a * x_t[i];
                    }
                }
            }
        }
    }
}

__global__ void cDotProduct_g(double *result, double *x_i, double *y_i,
                              const int s_x, const int s_y, const int s_z, const int s_t) {

    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;
    int z = blockDim.z * blockIdx.z + threadIdx.z;

    __shared__ double r_t[8000];

    const int pos_loc_xyz = blockDim.x * blockDim.y * threadIdx.z + blockDim.x * threadIdx.y + threadIdx.x;
    const int vol_loc = blockDim.x * blockDim.y * blockDim.z;

    r_t[pos_loc_xyz * 2 + 0] = 0;
    r_t[pos_loc_xyz * 2 + 1] = 0;

    for (int t = 0; t < s_t; t++) {
        int pos = s_x * s_y * s_z * t + s_x * s_y * z + s_x * y + x;
        hipDoubleComplex *x_t = (hipDoubleComplex *) (x_i + pos * 24);
        hipDoubleComplex *y_t = (hipDoubleComplex *) (y_i + pos * 24);
        for (int j = 0; j < 12; j++) {

            r_t[pos_loc_xyz * 2 + 0] += (y_t[j] * hipConj(x_t[j])).x;
            r_t[pos_loc_xyz * 2 + 1] += (y_t[j] * hipConj(x_t[j])).y;

        }
    }

    __syncthreads();

    for (int stride = 1; stride < vol_loc; stride *= 2) {
	
        if (pos_loc_xyz % (2 * stride) == 0 && pos_loc_xyz + stride < vol_loc) {

            r_t[pos_loc_xyz * 2 + 0] += r_t[(pos_loc_xyz + stride) * 2 + 0];
            r_t[pos_loc_xyz * 2 + 1] += r_t[(pos_loc_xyz + stride) * 2 + 1];
        }

        __syncthreads();

    }

    double real = r_t[0];
    double imag = r_t[1];


    if (threadIdx.x + threadIdx.y + threadIdx.z == 0) {
        result[0] += real;
        result[1] += imag;
    }
    __syncthreads();
}

__global__ void cDotProduct_g2(double *result, double *x_i, double *y_i,
                              const int s_x, const int s_y, const int s_z, const int s_t) {

    int x = blockDim.x * blockIdx.x + threadIdx.x;

    int thread = blockDim.x;	

    __shared__ double r_t[8000];

    int tot = s_x * s_y * s_z * s_t * 12; 	   
    
    int sub = tot / thread;	


    r_t[x * 2 + 0] = 0;
    r_t[x * 2 + 1] = 0;

    for (int i = 0; i < sub ; i++) {

	double * x_t = x_i + ( i * thread + x ) * 2 ;
	double * y_t = y_i + ( i * thread + x ) * 2 ;

	r_t[x * 2 + 0] += y_t[0] * x_t[0] + y_t[1] * x_t[1]; 
     	r_t[x * 2 + 1] += y_t[1] * x_t[0] - y_t[0] * x_t[1];
    }
    __syncthreads();

    for (int stride = 1; stride < thread; stride *= 2) {

        if (x % (2 * stride) == 0 && x + stride < thread) {

            r_t[x * 2 + 0] += r_t[(x + stride) * 2 + 0];
            r_t[x * 2 + 1] += r_t[(x + stride) * 2 + 1];
        }
        __syncthreads();
    }

    if (threadIdx.x  == 0) {
        result[0] = r_t[0];
        result[1] = r_t[1];
    }
}

__global__ void cDotProduct_v1_g2(double *result, double *x_i, double *y_i,
                              const int size, const int N, const int block) {

//    int x = blockDim.x * blockIdx.x + threadIdx.x;

    int x = threadIdx.x;

    int x_t = blockDim.x * blockIdx.x + threadIdx.x;

    int thread = blockDim.x;

    __shared__ double r_t[8000];
    

//    const int sub = size / thread;

//    for (int j = 0; j < N ; j++){	

//	    r_t[x * 2 + 0] = 0;
//	    r_t[x * 2 + 1] = 0;

//     for (int i = 0; i < sub; i++){

	

        double y_t_r = y_i[x_t * 2 + 0];
  	double y_t_i = y_i[x_t * 2 + 1];

	for (int j = 0; j < N; j++){

        double x_t_r = x_i[(j * size  + x_t ) * 2 + 0];
        double x_t_i = x_i[(j * size  + x_t ) * 2 + 1];

        r_t[x * 2 + 0] = y_t_r * x_t_r + y_t_i * x_t_i;
        r_t[x * 2 + 1] = y_t_i * x_t_r - y_t_r * x_t_i;
  
//     }	
	
    __syncthreads();


    for (int stride = 1; stride < thread; stride *= 2) {

        if (x % (2 * stride) == 0 && x + stride < thread) {

            r_t[x * 2 + 0] += r_t[(x + stride) * 2 + 0];
            r_t[x * 2 + 1] += r_t[(x + stride) * 2 + 1];
	    	
        }
        __syncthreads();
    }

      if (threadIdx.x  == 0) {
	
	  result[block * j * 2 + blockIdx.x *  2  + 0] = r_t[0];
	  result[block * j * 2 + blockIdx.x *  2  + 1] = r_t[1];
      }	
  }
}

__global__ void cDotProduct_v2_g2(double *result, double *x_i,
                              const int size, const int N, const int block, const int thread) {

	int x_t = blockDim.x * blockIdx.x + threadIdx.x;

	int x = threadIdx.x;

//	int thread = blockDim.x;

	__shared__ double r_t[8000];

//	for (int j = 0; j < N ; j++){

//		r_t[x * 2 + 0] = x_i[j * size * 2 + x_t * 2 + 0];
//		r_t[x * 2 + 1] = x_i[j * size * 2 + x_t * 2 + 1];

		int j = ( x_t / size );
		int i = ( x_t / thread) % block;

		r_t[x * 2 + 0] = x_i[x_t * 2 + 0];
		r_t[x * 2 + 1] = x_i[x_t * 2 + 1];		

		 __syncthreads();

		int x_j = x % thread;

		for (int stride = 1; stride < thread; stride *= 2) {

			if (x_j % (2 * stride) == 0 && x_j + stride < thread) {

				r_t[x * 2 + 0] += r_t[(x + stride) * 2 + 0];
				r_t[x * 2 + 1] += r_t[(x + stride) * 2 + 1];

			}
			__syncthreads();
		}		

		if (x_j  == 0) {

			result[block * j * 2 + i *  2  + 0] = r_t[x * 2 + 0];
			result[block * j * 2 + i *  2  + 1] = r_t[x * 2 + 1];
		}
//	}
}

void cDotProduct(double *result, double *x_i, double *y_i,
                 const int s_x, const int s_y, const int s_z, const int s_t) {

    std::complex<double> *r_t = (std::complex<double> *) result;

    *r_t = 0;

    for (int t = 0; t < s_t; t++) {
        for (int z = 0; z < s_z; z++) {
            for (int y = 0; y < s_y; y++) {
                for (int x = 0; x < s_x; x++) {
                    int pos = s_x * s_y * s_z * t + s_x * s_y * z + s_x * y + x;
                    std::complex<double> *x_t = (std::complex<double> *) (x_i + pos * 24);
                    std::complex<double> *y_t = (std::complex<double> *) (y_i + pos * 24);
                    for (int j = 0; j < 12; j++) {
                        *r_t = *r_t + y_t[j] * std::conj(x_t[j]);
                    }
                }
            }
        }
    }
}

__global__ void equal(double *out, double *in,
                      const int s_x, const int s_y, const int s_z, const int s_t) {

    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;
    int z = blockDim.z * blockIdx.z + threadIdx.z;

    for (int t = 0; t < s_t; t++) {
        int pos = s_x * s_y * s_z * t + s_x * s_y * z + s_x * y + x;
        for (int i = 0; i < 24; i++) {
            out[pos * 24 + i] = in[pos * 24 + i];
        }
    }
}

__global__ void  sign( int * out, double * in){

         int x = blockDim.x * blockIdx.x + threadIdx.x;

         out[x]= ( in[x] > 0 )? 1: -1;
}

__global__ void  mult( double * out, int * in){

         int x = blockDim.x * blockIdx.x + threadIdx.x;

         out[x]*= in[x/2];
	 
}

class DiracSetup{

public:

    double *U_x;
    double *U_y;
    double *U_z;
    double *U_t;

    int v_x, v_y, v_z, v_t;
    int s_x, s_y, s_z, s_t;
    
    int N_sub[4];
    int rank;
    int size;	

    int nodenum_x_b;
    int nodenum_x_f;
    int nodenum_y_b;
    int nodenum_y_f;
    int nodenum_z_b;
    int nodenum_z_f;
    int nodenum_t_b;
    int nodenum_t_f;

    int s_x_cb;

    int len_x;
    int len_y;
    int len_z;
    int len_t;

    double *resv_x_f;
    double *send_x_b;
    double *resv_x_b;
    double *send_x_f;

    double *resv_y_f;
    double *send_y_b;
    double *resv_y_b;
    double *send_y_f;

    double *resv_z_f;
    double *send_z_b;
    double *resv_z_b;
    double *send_z_f;

    double *resv_t_f;
    double *send_t_b;
    double *resv_t_b;
    double *send_t_f;

    double *tran_x_f;
    double *tran_x_b;
    double *tran_y_f;
    double *tran_y_b;
    double *tran_z_f;
    double *tran_z_b;
    double *tran_t_f;
    double *tran_t_b;

    double *ghos_x_f;
    double *ghos_x_b;
    double *ghos_y_f;
    double *ghos_y_b;
    double *ghos_z_f;
    double *ghos_z_b;
    double *ghos_t_f;
    double *ghos_t_b;	

    hipStream_t stre_x_f;
    hipStream_t stre_x_b;
    hipStream_t stre_y_f;
    hipStream_t stre_y_b;	
    hipStream_t stre_z_f;	
    hipStream_t stre_z_b;
    hipStream_t stre_t_f;
    hipStream_t stre_t_b;

    hipEvent_t control;	

    int volume;	

    DiracSetup(double *U_x_in, double *U_y_in, double *U_z_in, double *U_t_in,
              const int v_x_in, const int v_y_in, const int v_z_in, const int v_t_in,
              const int s_x_in, const int s_y_in, const int s_z_in, const int s_t_in){
	
	volume = s_x_in * s_y_in * s_z_in * s_t_in;

        int size_u = s_x_in * s_y_in * s_z_in * s_t_in * 9 * 2;

        hipMalloc((void **) &U_x, size_u * sizeof(double));
        hipMalloc((void **) &U_y, size_u * sizeof(double));
        hipMalloc((void **) &U_z, size_u * sizeof(double));
        hipMalloc((void **) &U_t, size_u * sizeof(double));

        hipMemcpy(U_x, U_x_in, size_u * sizeof(double), hipMemcpyHostToDevice);
        hipMemcpy(U_y, U_y_in, size_u * sizeof(double), hipMemcpyHostToDevice);
        hipMemcpy(U_z, U_z_in, size_u * sizeof(double), hipMemcpyHostToDevice);
        hipMemcpy(U_t, U_t_in, size_u * sizeof(double), hipMemcpyHostToDevice);

        v_x = v_x_in;
        v_y = v_y_in;
        v_z = v_z_in;
        v_t = v_t_in;
        s_x = s_x_in;
        s_y = s_y_in;
        s_z = s_z_in;
        s_t = s_t_in;	

        N_sub[0] = v_x / s_x;
        N_sub[1] = v_y / s_y;
        N_sub[2] = v_z / s_z;
        N_sub[3] = v_t / s_t;

        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        int site_x_f[4] = {(rank + 1) % N_sub[0],
                           (rank / N_sub[0]) % N_sub[1],
                           (rank / (N_sub[1] * N_sub[0])) % N_sub[2],
                           rank / (N_sub[2] * N_sub[1] * N_sub[0])};

        int site_x_b[4] = {(rank - 1 + N_sub[0]) % N_sub[0],
                           (rank / N_sub[0]) % N_sub[1],
                           (rank / (N_sub[1] * N_sub[0])) % N_sub[2],
                           rank / (N_sub[2] * N_sub[1] * N_sub[0])};

        nodenum_x_b = get_nodenum(site_x_b, N_sub, 4);
        nodenum_x_f = get_nodenum(site_x_f, N_sub, 4);

        int site_y_f[4] = {(rank % N_sub[0]),
                           (rank / N_sub[0] + 1) % N_sub[1],
                           (rank / (N_sub[1] * N_sub[0])) % N_sub[2],
                           rank / (N_sub[2] * N_sub[1] * N_sub[0])};

        int site_y_b[4] = {(rank % N_sub[0]),
                           (rank / N_sub[0] - 1 + N_sub[1]) % N_sub[1],
                           (rank / (N_sub[1] * N_sub[0])) % N_sub[2],
                           rank / (N_sub[2] * N_sub[1] * N_sub[0])};

        nodenum_y_b = get_nodenum(site_y_b, N_sub, 4);
        nodenum_y_f = get_nodenum(site_y_f, N_sub, 4);

        int site_z_f[4] = {(rank % N_sub[0]),
                           (rank / N_sub[0]) % N_sub[1],
                           (rank / (N_sub[1] * N_sub[0]) + 1) % N_sub[2],
                           rank / (N_sub[2] * N_sub[1] * N_sub[0])};

        int site_z_b[4] = {(rank % N_sub[0]),
                           (rank / N_sub[0]) % N_sub[1],
                           (rank / (N_sub[1] * N_sub[0]) - 1 + N_sub[2]) % N_sub[2],
                           rank / (N_sub[2] * N_sub[1] * N_sub[0])};

        nodenum_z_b = get_nodenum(site_z_b, N_sub, 4);
        nodenum_z_f = get_nodenum(site_z_f, N_sub, 4);

        int site_t_f[4] = {(rank % N_sub[0]),
                           (rank / N_sub[0]) % N_sub[1],
                           (rank / (N_sub[1] * N_sub[0])) % N_sub[2],
                           (rank / (N_sub[2] * N_sub[1] * N_sub[0]) + 1) % N_sub[3]};

        int site_t_b[4] = {(rank % N_sub[0]),
                           (rank / N_sub[0]) % N_sub[1],
                           (rank / (N_sub[1] * N_sub[0])) % N_sub[2],
                           (rank / (N_sub[2] * N_sub[1] * N_sub[0]) - 1 + N_sub[3]) % N_sub[3]};

        nodenum_t_b = get_nodenum(site_t_b, N_sub, 4);
        nodenum_t_f = get_nodenum(site_t_f, N_sub, 4);

	s_x_cb = s_x >> 1;

        len_x = s_y * s_z * s_t >> 1;
        len_y = s_x_cb * s_z * s_t;
        len_z = s_x_cb * s_y * s_t;
        len_t = s_x_cb * s_y * s_z;

        if (N_sub[0] != 1) {
/*
            resv_x_f = new double[len_x * 6 * 2];
            send_x_b = new double[len_x * 6 * 2];
            resv_x_b = new double[len_x * 6 * 2];
            send_x_f = new double[len_x * 6 * 2];
*/
            int size_T = len_x * 6 * 2;

            hipMallocManaged((void **) &resv_x_f, size_T * sizeof(double));
            hipMallocManaged((void **) &send_x_b, size_T * sizeof(double));
            hipMallocManaged((void **) &resv_x_b, size_T * sizeof(double));
            hipMallocManaged((void **) &send_x_f, size_T * sizeof(double));


            hipMalloc((void **) &tran_x_f, size_T * sizeof(double));
            hipMalloc((void **) &tran_x_b, size_T * sizeof(double));

            hipMalloc((void **) &ghos_x_f, size_T * sizeof(double));
            hipMalloc((void **) &ghos_x_b, size_T * sizeof(double));

        }

	if (N_sub[1] != 1) {
/*
            resv_y_f = new double[len_y * 6 * 2];
            send_y_b = new double[len_y * 6 * 2];
            resv_y_b = new double[len_y * 6 * 2];
            send_y_f = new double[len_y * 6 * 2];
*/
            int size_T = len_y * 6 * 2;

	    hipMallocManaged((void **) &resv_y_f, size_T * sizeof(double));
            hipMallocManaged((void **) &send_y_b, size_T * sizeof(double));
            hipMallocManaged((void **) &resv_y_b, size_T * sizeof(double));
            hipMallocManaged((void **) &send_y_f, size_T * sizeof(double));	

            hipMalloc((void **) &tran_y_f, size_T * sizeof(double));
            hipMalloc((void **) &tran_y_b, size_T * sizeof(double));

            hipMalloc((void **) &ghos_y_f, size_T * sizeof(double));
            hipMalloc((void **) &ghos_y_b, size_T * sizeof(double));

        }

	if (N_sub[2] != 1) {
/*
            resv_z_f = new double[len_z * 6 * 2];
            send_z_b = new double[len_z * 6 * 2];
            resv_z_b = new double[len_z * 6 * 2];
            send_z_f = new double[len_z * 6 * 2];
*/
            int size_T = len_z * 6 * 2;

	    hipMallocManaged((void **) &resv_z_f, size_T * sizeof(double));
            hipMallocManaged((void **) &send_z_b, size_T * sizeof(double));
            hipMallocManaged((void **) &resv_z_b, size_T * sizeof(double));
            hipMallocManaged((void **) &send_z_f, size_T * sizeof(double));

            hipMalloc((void **) &tran_z_f, size_T * sizeof(double));
            hipMalloc((void **) &tran_z_b, size_T * sizeof(double));

            hipMalloc((void **) &ghos_z_f, size_T * sizeof(double));
            hipMalloc((void **) &ghos_z_b, size_T * sizeof(double));

        }

	if (N_sub[3] != 1) {
/*
            resv_t_f = new double[len_t * 6 * 2];
            send_t_b = new double[len_t * 6 * 2];
            resv_t_b = new double[len_t * 6 * 2];
            send_t_f = new double[len_t * 6 * 2];
*/
            int size_T = len_t * 6 * 2;

	    hipMallocManaged((void **) &resv_t_f, size_T * sizeof(double));
	    hipMallocManaged((void **) &send_t_b, size_T * sizeof(double));
	    hipMallocManaged((void **) &resv_t_b, size_T * sizeof(double));
            hipMallocManaged((void **) &send_t_f, size_T * sizeof(double));

/*
            hipHostMalloc((void **) &resv_t_f, size_T * sizeof(double));
            hipHostMalloc((void **) &send_t_b, size_T * sizeof(double));
            hipHostMalloc((void **) &resv_t_f, size_T * sizeof(double));
            hipHostMalloc((void **) &send_t_b, size_T * sizeof(double));	
*/
	    		
            hipMalloc((void **) &tran_t_f, size_T * sizeof(double));
            hipMalloc((void **) &tran_t_b, size_T * sizeof(double));

            hipMalloc((void **) &ghos_t_f, size_T * sizeof(double));
            hipMalloc((void **) &ghos_t_b, size_T * sizeof(double));

        }

        hipStreamCreate(&stre_x_f);
        hipStreamCreate(&stre_x_b);
        hipStreamCreate(&stre_y_f);
        hipStreamCreate(&stre_y_b);
        hipStreamCreate(&stre_t_f);
        hipStreamCreate(&stre_t_b);
	hipStreamCreate(&stre_z_f);
	hipStreamCreate(&stre_z_b);
        hipEventCreate(&control);
	
    }

	~DiracSetup(){

	        hipFree(U_x);
	        hipFree(U_y);
	        hipFree(U_z);
	        hipFree(U_t);

		if (N_sub[0] != 1) {
/*
        		delete[] resv_x_f;
                	delete[] send_x_b;
                	delete[] resv_x_b;
                	delete[] send_x_f;
*/
			hipHostFree(resv_x_f);
                        hipHostFree(send_x_b);
                        hipHostFree(resv_x_b);
                        hipHostFree(send_x_f);

            		hipFree(tran_x_f);
            		hipFree(tran_x_b);

            		hipFree(ghos_x_f);
            		hipFree(ghos_x_b);
        	}

	        if (N_sub[1] != 1) {
/*
        		delete[] resv_y_f;
                	delete[] send_y_b;
                	delete[] resv_y_b;
                	delete[] send_y_f;
*/
			hipHostFree(resv_y_f);
                        hipHostFree(send_y_b);
                        hipHostFree(resv_y_b);
                        hipHostFree(send_y_f);

                	hipFree(tran_y_f);
                	hipFree(tran_y_b);

                	hipFree(ghos_y_f);
                	hipFree(ghos_y_b);
        	}

		if (N_sub[2] != 1) {
/*
			delete[] resv_z_f;
			delete[] send_z_b;
			delete[] resv_z_b;
			delete[] send_z_f;
*/
			hipHostFree(resv_z_f);
                        hipHostFree(send_z_b);
                        hipHostFree(resv_z_b);
                        hipHostFree(send_z_f);

			hipFree(tran_z_f);
			hipFree(tran_z_b);

			hipFree(ghos_z_f);
			hipFree(ghos_z_b);
		}

		if (N_sub[3] != 1) {
/*
			delete[] resv_t_f;
			delete[] send_t_b;
			delete[] resv_t_b;
			delete[] send_t_f;
*/

			hipHostFree(resv_t_f);
			hipHostFree(send_t_b);
			hipHostFree(resv_t_b);
			hipHostFree(send_t_f);

			hipFree(tran_t_f);
			hipFree(tran_t_b);

			hipFree(ghos_t_f);
			hipFree(ghos_t_b);
        	}

                hipStreamDestroy(stre_x_f);
                hipStreamDestroy(stre_x_b);
                hipStreamDestroy(stre_y_f);
                hipStreamDestroy(stre_y_b);
		hipStreamDestroy(stre_t_f);
		hipStreamDestroy(stre_t_b);
		hipStreamDestroy(stre_z_f);
		hipStreamDestroy(stre_z_b);
		hipEventDestroy(control);
	}
};

class DiracWilson:public DiracSetup{

	public:

	DiracWilson(double *U_x_in, double *U_y_in, double *U_z_in, double *U_t_in,
              	    const int v_x_in, const int v_y_in, const int v_z_in, const int v_t_in,
                    const int s_x_in, const int s_y_in, const int s_z_in, const int s_t_in):
		    DiracSetup(U_x_in, U_y_in, U_z_in, U_t_in, v_x_in, v_y_in, v_z_in, v_t_in, s_x_in, s_y_in, s_z_in, s_t_in){
	}

	~DiracWilson(){}

	inline void dslash(double *src_g, double *dest_g, const int cb, const int flag) {

	        hipEventRecord(control, 0);

	        MPI_Request req[8 * size];
	        MPI_Request reqr[8 * size];
	        MPI_Status sta[8 * size];

		if (N_sub[0] != 1) {

                        hipStreamWaitEvent(stre_x_f, control, 0);

                        transfer_x_f<<<s_y * s_z * s_t / 64 / 2, 64, 0, stre_x_f>>>(src_g, tran_x_b, U_x, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);

                        hipStreamWaitEvent(stre_x_b, control, 0);

                        transfer_x_b<<<s_y * s_z * s_t / 64 / 2, 64, 0, stre_x_b>>>(src_g, tran_x_f, U_x, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);
                }

		if (N_sub[1] != 1) {

                        hipStreamWaitEvent(stre_y_f, control, 0);

                        transfer_y_f<<<s_x_cb * s_z * s_t / 64, 64, 0, stre_y_f>>>(src_g, tran_y_b, U_y, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);

                        hipStreamWaitEvent(stre_y_b, control, 0);

                        transfer_y_b<<<s_x_cb * s_z * s_t / 64, 64, 0, stre_y_b>>>(src_g, tran_y_f, U_y, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);
                }

		if (N_sub[2] != 1) {

			hipStreamWaitEvent(stre_z_f, control, 0);

			transfer_z_f<<<s_x_cb * s_y * s_t / 64, 64, 0, stre_z_f>>>(src_g, tran_z_b, U_z, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);

			hipStreamWaitEvent(stre_z_b, control, 0);

			transfer_z_b<<<s_x_cb * s_y * s_t / 64, 64, 0, stre_z_b>>>(src_g, tran_z_f, U_z, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);
		}

		if (N_sub[3] != 1) {

		  	hipStreamWaitEvent(stre_t_f, control, 0);

			transfer_t_f<<<s_x_cb * s_y * s_z / 64, 64, 0, stre_t_f>>>(src_g, tran_t_b, U_t, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);

			hipStreamWaitEvent(stre_t_b, control, 0);

			transfer_t_b<<<s_x_cb * s_y * s_z / 64, 64, 0, stre_t_b>>>(src_g, tran_t_f, U_t, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);
        	}

		main_xyzt<<< s_x_cb * s_y * s_z * s_t / 64, 64>>>(src_g, dest_g, U_x, U_y, U_z, U_t, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);

		if (N_sub[0] != 1) {

                        hipMemcpyAsync(send_x_b, tran_x_b, len_x * 6 * 2 * sizeof(double), hipMemcpyDeviceToHost, stre_x_f);

			hipStreamSynchronize(stre_x_f);

                        MPI_Isend(send_x_b, len_x * 6 * 2, MPI_DOUBLE, nodenum_x_b, 8 * rank + 0, MPI_COMM_WORLD, &req[8 * rank + 0]);
                        MPI_Irecv(resv_x_f, len_x * 6 * 2, MPI_DOUBLE, nodenum_x_f, 8 * nodenum_x_f + 0, MPI_COMM_WORLD, &reqr[8 * nodenum_x_f + 0]);

                        hipMemcpyAsync(send_x_f, tran_x_f, len_x * 6 * 2 * sizeof(double), hipMemcpyDeviceToHost, stre_x_b);

			hipStreamSynchronize(stre_x_b);

                        MPI_Isend(send_x_f, len_x * 6 * 2, MPI_DOUBLE, nodenum_x_f, 8 * rank + 1, MPI_COMM_WORLD, &req[8 * rank + 1]);
                        MPI_Irecv(resv_x_b, len_x * 6 * 2, MPI_DOUBLE, nodenum_x_b, 8 * nodenum_x_b + 1, MPI_COMM_WORLD, &reqr[8 * nodenum_x_b + 1]);

                        MPI_Wait(&reqr[8 * nodenum_x_f + 0], &sta[8 * nodenum_x_f + 0]);

                        hipMemcpyAsync(ghos_x_f, resv_x_f, len_x * 6 * 2 * sizeof(double), hipMemcpyHostToDevice, stre_x_f);

                        hipEventRecord(control, stre_x_f);
                        hipStreamWaitEvent(0, control, 0);

                        ghost_x_f<<<s_y * s_z * s_t / 64 / 2, 64>>>(ghos_x_f, dest_g, U_x, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);

                        MPI_Wait(&reqr[8 * nodenum_x_b + 1], &sta[8 * nodenum_x_b + 1]);

                        hipMemcpyAsync(ghos_x_b, resv_x_b, len_x * 6 * 2 * sizeof(double), hipMemcpyHostToDevice, stre_x_b);

                        hipEventRecord(control, stre_x_b);
                        hipStreamWaitEvent(0, control, 0);

                        ghost_x_b<<<s_y * s_z * s_t / 64 / 2, 64>>>(ghos_x_b, dest_g, U_x, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);
                }

		if (N_sub[1] != 1) {

                        hipMemcpyAsync(send_y_b, tran_y_b, len_y * 6 * 2 * sizeof(double), hipMemcpyDeviceToHost, stre_y_f);

			hipStreamSynchronize(stre_y_f);

                        MPI_Isend(send_y_b, len_y * 6 * 2, MPI_DOUBLE, nodenum_y_b, 8 * rank + 2, MPI_COMM_WORLD, &req[8 * rank + 2]);
                        MPI_Irecv(resv_y_f, len_y * 6 * 2, MPI_DOUBLE, nodenum_y_f, 8 * nodenum_y_f + 2, MPI_COMM_WORLD, &reqr[8 * nodenum_y_f + 2]);

                        hipMemcpyAsync(send_y_f, tran_y_f, len_y * 6 * 2 * sizeof(double), hipMemcpyDeviceToHost, stre_y_b);

			hipStreamSynchronize(stre_y_b);

                        MPI_Isend(send_y_f, len_y * 6 * 2, MPI_DOUBLE, nodenum_y_f, 8 * rank + 3, MPI_COMM_WORLD, &req[8 * rank + 3]);
                        MPI_Irecv(resv_y_b, len_y * 6 * 2, MPI_DOUBLE, nodenum_y_b, 8 * nodenum_y_b + 3, MPI_COMM_WORLD, &reqr[8 * nodenum_y_b + 3]);

                        MPI_Wait(&reqr[8 * nodenum_y_f + 2], &sta[8 * nodenum_y_f + 2]);

                        hipMemcpyAsync(ghos_y_f, resv_y_f, len_y * 6 * 2 * sizeof(double), hipMemcpyHostToDevice, stre_y_f);

                        hipEventRecord(control, stre_y_f);
                        hipStreamWaitEvent(0, control, 0);

                        ghost_y_f<<<s_x_cb * s_z * s_t / 64, 64>>>(ghos_y_f, dest_g, U_y, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);

                        MPI_Wait(&reqr[8 * nodenum_y_b + 3], &sta[8 * nodenum_y_b + 3]);

                        hipMemcpyAsync(ghos_y_b, resv_y_b, len_y * 6 * 2 * sizeof(double), hipMemcpyHostToDevice, stre_y_b);

                        hipEventRecord(control, stre_y_b);
                        hipStreamWaitEvent(0, control, 0);

                        ghost_y_b<<<s_x_cb * s_z * s_t / 64, 64>>>(ghos_y_b, dest_g, U_y, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);
                }

		if (N_sub[2] != 1) {

                        hipMemcpyAsync(send_z_b, tran_z_b, len_z * 6 * 2 * sizeof(double), hipMemcpyDeviceToHost, stre_z_f);

			hipStreamSynchronize(stre_z_f);

                        MPI_Isend(send_z_b, len_z * 6 * 2, MPI_DOUBLE, nodenum_z_b, 8 * rank + 4, MPI_COMM_WORLD, &req[8 * rank + 4]);
                        MPI_Irecv(resv_z_f, len_z * 6 * 2, MPI_DOUBLE, nodenum_z_f, 8 * nodenum_z_f + 4, MPI_COMM_WORLD, &reqr[8 * nodenum_z_f + 4]);

                        hipMemcpyAsync(send_z_f, tran_z_f, len_z * 6 * 2 * sizeof(double), hipMemcpyDeviceToHost, stre_z_b);

			hipStreamSynchronize(stre_z_b);

                        MPI_Isend(send_z_f, len_z * 6 * 2, MPI_DOUBLE, nodenum_z_f, 8 * rank + 5, MPI_COMM_WORLD, &req[8 * rank + 5]);
                        MPI_Irecv(resv_z_b, len_z * 6 * 2, MPI_DOUBLE, nodenum_z_b, 8 * nodenum_z_b + 5, MPI_COMM_WORLD, &reqr[8 * nodenum_z_b + 5]);

                        MPI_Wait(&reqr[8 * nodenum_z_f + 4], &sta[8 * nodenum_z_f + 4]);

                        hipMemcpyAsync(ghos_z_f, resv_z_f, len_z * 6 * 2 * sizeof(double), hipMemcpyHostToDevice, stre_z_f);

                        hipEventRecord(control, stre_z_f);
                        hipStreamWaitEvent(0, control, 0);

                        ghost_z_f<<<s_x_cb * s_y * s_t / 64, 64>>>(ghos_z_f, dest_g, U_z, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);

                        MPI_Wait(&reqr[8 * nodenum_z_b + 5], &sta[8 * nodenum_z_b + 5]);

                        hipMemcpyAsync(ghos_z_b, resv_z_b, len_z * 6 * 2 * sizeof(double), hipMemcpyHostToDevice, stre_z_b);

                        hipEventRecord(control, stre_z_b);
                        hipStreamWaitEvent(0, control, 0);

                        ghost_z_b<<<s_x_cb * s_y * s_t / 64, 64>>>(ghos_z_b, dest_g, U_z, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);
                }

	        if (N_sub[3] != 1) {

			hipMemcpyAsync(send_t_b, tran_t_b, len_t * 6 * 2 * sizeof(double), hipMemcpyDeviceToHost, stre_t_f);

			hipStreamSynchronize(stre_t_f);

			MPI_Isend(send_t_b, len_t * 6 * 2, MPI_DOUBLE, nodenum_t_b, 8 * rank + 6, MPI_COMM_WORLD, &req[8 * rank + 6]);
			MPI_Irecv(resv_t_f, len_t * 6 * 2, MPI_DOUBLE, nodenum_t_f, 8 * nodenum_t_f + 6, MPI_COMM_WORLD, &reqr[8 * nodenum_t_f + 6]);

			hipMemcpyAsync(send_t_f, tran_t_f, len_t * 6 * 2 * sizeof(double), hipMemcpyDeviceToHost, stre_t_b);

			hipStreamSynchronize(stre_t_b);

			MPI_Isend(send_t_f, len_t * 6 * 2, MPI_DOUBLE, nodenum_t_f, 8 * rank + 7, MPI_COMM_WORLD, &req[8 * rank + 7]);
			MPI_Irecv(resv_t_b, len_t * 6 * 2, MPI_DOUBLE, nodenum_t_b, 8 * nodenum_t_b + 7, MPI_COMM_WORLD, &reqr[8 * nodenum_t_b + 7]);

			MPI_Wait(&reqr[8 * nodenum_t_f + 6], &sta[8 * nodenum_t_f + 6]);

			hipMemcpyAsync(ghos_t_f, resv_t_f, len_t * 6 * 2 * sizeof(double), hipMemcpyHostToDevice, stre_t_f);

			hipEventRecord(control, stre_t_f);
			hipStreamWaitEvent(0, control, 0);

			ghost_t_f<<<s_x_cb * s_y * s_z / 64, 64>>>(ghos_t_f, dest_g, U_t, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);

			MPI_Wait(&reqr[8 * nodenum_t_b + 7], &sta[8 * nodenum_t_b + 7]);

			hipMemcpyAsync(ghos_t_b, resv_t_b, len_t * 6 * 2 * sizeof(double), hipMemcpyHostToDevice, stre_t_b);

			hipEventRecord(control, stre_t_b);
			hipStreamWaitEvent(0, control, 0);
	
			ghost_t_b<<<s_x_cb * s_y * s_z / 64, 64>>>(ghos_t_b, dest_g, U_t, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);
		}

	    MPI_Barrier(MPI_COMM_WORLD);	

	}

	void Kernal(double *out, double *in) {
/*
        float eventMs = 1.0f;
        hipEvent_t start, stop;
        hipEventCreate(&start);
        hipEventCreate(&stop);
        hipEventRecord(start, 0);
*/
        dslash(in, out, 0, 1);
        dslash(in, out, 1, 1);
/*
	hipEventRecord(stop, 0);
	hipEventSynchronize(stop);
	hipEventElapsedTime(&eventMs, start, stop);
	printf("dslash time taken  = %6.3fms\n", eventMs);
*/
	}  

};


class DiracOverlapWilson:public DiracSetup {

protected:

    double prec0;
    std::vector<double *> hw_evec;
    std::vector<double> hw_eval;
    std::vector<std::vector<double>> coef;
    std::vector<int> hw_size;
    bool build_hw;

    void calc_coef(double cut);

public:

    double rho;
    double kappa;

    double *tmp1;
    double *tmp2;
   
    double *hw_evec_g;
    int * sig_hw_eval_g;

    DiracOverlapWilson(std::vector<double *> &evecs,
                       std::vector<double> &evals,
                       std::vector<std::vector<double> > &coefs,
                       std::vector<int> &sizes,
                       double *U_x_in, double *U_y_in, double *U_z_in, double *U_t_in,
                       const int v_x_in, const int v_y_in, const int v_z_in, const int v_t_in,
                       const int s_x_in, const int s_y_in, const int s_z_in, const int s_t_in,
                       const double kappa_in) : prec0(1e-12), hw_evec(evecs), hw_eval(evals), coef(coefs), hw_size(sizes),
		       DiracSetup(U_x_in, U_y_in, U_z_in, U_t_in, v_x_in, v_y_in, v_z_in, v_t_in, s_x_in, s_y_in, s_z_in, s_t_in){

        build_hw = false;
        kappa = kappa_in;
        rho = 4 - 0.5 / kappa;

        int size_f = s_x * s_y * s_z * s_t * 12 * 2;
        hipMalloc((void **) &tmp1, size_f * sizeof(double));
        hipMalloc((void **) &tmp2, size_f * sizeof(double));

        hw_eval = evals;

 	hipMalloc((void **) &hw_evec_g, hw_evec.size() * size_f * sizeof(double));

	for ( int i=0; i < hw_evec.size(); i++  ){
	
	hipMemcpy(hw_evec_g + i * size_f , hw_evec[i], size_f * sizeof(double), hipMemcpyHostToDevice);

	}

	double * hw_eval_g;
	hipMalloc((void **) &hw_eval_g, hw_evec.size() * sizeof(double));
	hipMemcpy(hw_eval_g, &hw_eval[0], hw_evec.size() * sizeof(double), hipMemcpyHostToDevice);
	
	hipMalloc((void **) &sig_hw_eval_g, hw_evec.size() * sizeof(int));
	sign<<<1,hw_evec.size()>>>(sig_hw_eval_g,hw_eval_g);
	hipFree(hw_eval_g);
    }

    ~DiracOverlapWilson() {
      
	hipFree(tmp1);
        hipFree(tmp2);
	hipFree(hw_evec_g);
	hipFree(sig_hw_eval_g);

    }

    inline void dslash(double *src_g, double *dest_g, const int cb, const int flag, const double a, const double b) {

	hipEventRecord(control, 0);

        MPI_Request req[8 * size];
        MPI_Request reqr[8 * size];
        MPI_Status sta[8 * size];

	if (N_sub[0] != 1) {

            hipStreamWaitEvent(stre_x_f,control,0);

            transfer_x_f<<<s_y * s_z * s_t / 64 / 2, 64, 0, stre_x_f>>>(src_g, tran_x_b, U_x, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);

            hipStreamWaitEvent(stre_x_b,control,0);

            transfer_x_b<<<s_y * s_z * s_t / 64 / 2, 64, 0, stre_x_b>>>(src_g, tran_x_f, U_x, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);
        }

        if (N_sub[1] != 1) {

            hipStreamWaitEvent(stre_y_f,control,0);

            transfer_y_f<<<s_x_cb * s_z * s_t / 64, 64, 0, stre_y_f>>>(src_g, tran_y_b, U_y, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);

            hipStreamWaitEvent(stre_y_b,control,0);

            transfer_y_b<<<s_x_cb * s_z * s_t / 64, 64, 0, stre_y_b>>>(src_g, tran_y_f, U_y, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);
        }

        if (N_sub[2] != 1) {

            hipStreamWaitEvent(stre_z_f,control,0);

            transfer_z_f<<<s_x_cb * s_y * s_t / 64, 64, 0, stre_z_f>>>(src_g, tran_z_b, U_z, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);

            hipStreamWaitEvent(stre_z_b,control,0);

            transfer_z_b<<<s_x_cb * s_y * s_t / 64, 64, 0, stre_z_b>>>(src_g, tran_z_f, U_z, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);
        }
	
        if (N_sub[3] != 1) {

	    hipStreamWaitEvent(stre_t_f,control,0);

	    transfer_t_f<<<s_x_cb * s_y * s_z / 64, 64, 0, stre_t_f>>>(src_g, tran_t_b, U_t, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);

	    hipStreamWaitEvent(stre_t_b,control,0); 

            transfer_t_b<<<s_x_cb * s_y * s_z / 64, 64, 0, stre_t_b>>>(src_g, tran_t_f, U_t, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag);
	}

        main_xyzt_abp5<<< s_x_cb * s_y * s_z * s_t / 64, 64>>>(src_g, dest_g, U_x, U_y, U_z, U_t, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag, a, b);

	if (N_sub[0] != 1) {

            hipMemcpyAsync(send_x_b, tran_x_b, len_x * 6 * 2 * sizeof(double), hipMemcpyDeviceToHost, stre_x_f);

	    hipStreamSynchronize(stre_x_f);

            MPI_Isend(send_x_b, len_x * 6 * 2, MPI_DOUBLE, nodenum_x_b, 8 * rank + 0, MPI_COMM_WORLD,
                      &req[8 * rank + 0]);
            MPI_Irecv(resv_x_f, len_x * 6 * 2, MPI_DOUBLE, nodenum_x_f, 8 * nodenum_x_f + 0, MPI_COMM_WORLD,
                      &reqr[8 * nodenum_x_f + 0]);

            hipMemcpyAsync(send_x_f, tran_x_f, len_x * 6 * 2 * sizeof(double), hipMemcpyDeviceToHost, stre_x_b);

	    hipStreamSynchronize(stre_x_b);

            MPI_Isend(send_x_f, len_x * 6 * 2, MPI_DOUBLE, nodenum_x_f, 8 * rank + 1, MPI_COMM_WORLD,
                      &req[8 * rank + 1]);
            MPI_Irecv(resv_x_b, len_x * 6 * 2, MPI_DOUBLE, nodenum_x_b, 8 * nodenum_x_b + 1, MPI_COMM_WORLD,
                      &reqr[8 * nodenum_x_b + 1]);

            MPI_Wait(&reqr[8 * nodenum_x_f + 0], &sta[8 * nodenum_x_f + 0]);

            hipMemcpyAsync(ghos_x_f, resv_x_f, len_x * 6 * 2 * sizeof(double), hipMemcpyHostToDevice, stre_x_f);

            hipEventRecord(control, stre_x_f);
            hipStreamWaitEvent(0, control, 0);

            ghost_x_f_abp5<<<s_y * s_z * s_t / 64 / 2, 64>>>(ghos_x_f, dest_g, U_x, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag, b);

            MPI_Wait(&reqr[8 * nodenum_x_b + 1], &sta[8 * nodenum_x_b + 1]);

            hipMemcpyAsync(ghos_x_b, resv_x_b, len_x * 6 * 2 * sizeof(double), hipMemcpyHostToDevice, stre_x_b);

            hipEventRecord(control, stre_x_b);
            hipStreamWaitEvent(0, control, 0);

            ghost_x_b_abp5<<<s_y * s_z * s_t / 64 / 2, 64>>>(ghos_x_b, dest_g, U_x, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag, b);

        }

        if (N_sub[1] != 1) {

            hipMemcpyAsync(send_y_b, tran_y_b, len_y * 6 * 2 * sizeof(double), hipMemcpyDeviceToHost, stre_y_f);

	    hipStreamSynchronize(stre_y_f);

            MPI_Isend(send_y_b, len_y * 6 * 2, MPI_DOUBLE, nodenum_y_b, 8 * rank + 2, MPI_COMM_WORLD,
                      &req[8 * rank + 2]);
            MPI_Irecv(resv_y_f, len_y * 6 * 2, MPI_DOUBLE, nodenum_y_f, 8 * nodenum_y_f + 2, MPI_COMM_WORLD,
                      &reqr[8 * nodenum_y_f + 2]);

            hipMemcpyAsync(send_y_f, tran_y_f, len_y * 6 * 2 * sizeof(double), hipMemcpyDeviceToHost, stre_y_b);

	    hipStreamSynchronize(stre_y_b);

            MPI_Isend(send_y_f, len_y * 6 * 2, MPI_DOUBLE, nodenum_y_f, 8 * rank + 3, MPI_COMM_WORLD,
                      &req[8 * rank + 3]);
            MPI_Irecv(resv_y_b, len_y * 6 * 2, MPI_DOUBLE, nodenum_y_b, 8 * nodenum_y_b + 3, MPI_COMM_WORLD,
                      &reqr[8 * nodenum_y_b + 3]);

            MPI_Wait(&reqr[8 * nodenum_y_f + 2], &sta[8 * nodenum_y_f + 2]);

            hipMemcpyAsync(ghos_y_f, resv_y_f, len_y * 6 * 2 * sizeof(double), hipMemcpyHostToDevice, stre_y_f);

            hipEventRecord(control, stre_y_f);
            hipStreamWaitEvent(0, control, 0);

            ghost_y_f_abp5<<<s_x_cb * s_z * s_t / 64, 64>>>(ghos_y_f, dest_g, U_y, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag, b);

            MPI_Wait(&reqr[8 * nodenum_y_b + 3], &sta[8 * nodenum_y_b + 3]);

            hipMemcpyAsync(ghos_y_b, resv_y_b, len_y * 6 * 2 * sizeof(double), hipMemcpyHostToDevice, stre_y_b);

            hipEventRecord(control, stre_y_b);
            hipStreamWaitEvent(0, control, 0);

            ghost_y_b_abp5<<<s_x_cb * s_z * s_t / 64, 64>>>(ghos_y_b, dest_g, U_y, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag, b);

        }

	if (N_sub[2] != 1) {

            hipMemcpyAsync(send_z_b, tran_z_b, len_z * 6 * 2 * sizeof(double), hipMemcpyDeviceToHost, stre_z_f);

	    hipStreamSynchronize(stre_z_f);

            MPI_Isend(send_z_b, len_z * 6 * 2, MPI_DOUBLE, nodenum_z_b, 8 * rank + 4, MPI_COMM_WORLD,
                      &req[8 * rank + 4]);
            MPI_Irecv(resv_z_f, len_z * 6 * 2, MPI_DOUBLE, nodenum_z_f, 8 * nodenum_z_f + 4, MPI_COMM_WORLD,
                      &reqr[8 * nodenum_z_f + 4]);
   
            hipMemcpyAsync(send_z_f, tran_z_f, len_z * 6 * 2 * sizeof(double), hipMemcpyDeviceToHost, stre_z_b);

	    hipStreamSynchronize(stre_z_b);

            MPI_Isend(send_z_f, len_z * 6 * 2, MPI_DOUBLE, nodenum_z_f, 8 * rank + 5, MPI_COMM_WORLD,
                      &req[8 * rank + 5]);
            MPI_Irecv(resv_z_b, len_z * 6 * 2, MPI_DOUBLE, nodenum_z_b, 8 * nodenum_z_b + 5, MPI_COMM_WORLD,
                      &reqr[8 * nodenum_z_b + 5]);

            MPI_Wait(&reqr[8 * nodenum_z_f + 4], &sta[8 * nodenum_z_f + 4]);

            hipMemcpyAsync(ghos_z_f, resv_z_f, len_z * 6 * 2 * sizeof(double), hipMemcpyHostToDevice, stre_z_f);

            hipEventRecord(control, stre_z_f);
            hipStreamWaitEvent(0, control, 0);

            ghost_z_f_abp5<<<s_x_cb * s_y * s_t / 64, 64>>>(ghos_z_f, dest_g, U_z, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag, b);

            MPI_Wait(&reqr[8 * nodenum_z_b + 5], &sta[8 * nodenum_z_b + 5]);

            hipMemcpyAsync(ghos_z_b, resv_z_b, len_z * 6 * 2 * sizeof(double), hipMemcpyHostToDevice, stre_z_b);

            hipEventRecord(control, stre_z_b);
            hipStreamWaitEvent(0, control, 0);

            ghost_z_b_abp5<<<s_x_cb * s_y * s_t / 64, 64>>>(ghos_z_b, dest_g, U_z, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag, b);

        }

	if (N_sub[3] != 1) {

	    hipMemcpyAsync(send_t_b, tran_t_b, len_t * 6 * 2 * sizeof(double), hipMemcpyDeviceToHost, stre_t_f);            

 	    hipStreamSynchronize(stre_t_f);	    	

            MPI_Isend(send_t_b, len_t * 6 * 2, MPI_DOUBLE, nodenum_t_b, 8 * rank + 6, MPI_COMM_WORLD,
                      &req[8 * rank + 6]);
            MPI_Irecv(resv_t_f, len_t * 6 * 2, MPI_DOUBLE, nodenum_t_f, 8 * nodenum_t_f + 6, MPI_COMM_WORLD,
                      &reqr[8 * nodenum_t_f + 6]);
   
	    hipMemcpyAsync(send_t_f, tran_t_f, len_t * 6 * 2 * sizeof(double), hipMemcpyDeviceToHost, stre_t_b);

	    hipStreamSynchronize(stre_t_b);	

            MPI_Isend(send_t_f, len_t * 6 * 2, MPI_DOUBLE, nodenum_t_f, 8 * rank + 7, MPI_COMM_WORLD,
                      &req[8 * rank + 7]);
            MPI_Irecv(resv_t_b, len_t * 6 * 2, MPI_DOUBLE, nodenum_t_b, 8 * nodenum_t_b + 7, MPI_COMM_WORLD,
                      &reqr[8 * nodenum_t_b + 7]);

            MPI_Wait(&reqr[8 * nodenum_t_f + 6], &sta[8 * nodenum_t_f + 6]);

            hipMemcpyAsync(ghos_t_f, resv_t_f, len_t * 6 * 2 * sizeof(double), hipMemcpyHostToDevice, stre_t_f);

            hipEventRecord(control, stre_t_f);
            hipStreamWaitEvent(0, control, 0);

	    ghost_t_f_abp5<<<s_x_cb * s_y * s_z / 64, 64>>>(ghos_t_f, dest_g, U_t, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag, b);

            MPI_Wait(&reqr[8 * nodenum_t_b + 7], &sta[8 * nodenum_t_b + 7]);

            hipMemcpyAsync(ghos_t_b, resv_t_b, len_t * 6 * 2 * sizeof(double), hipMemcpyHostToDevice, stre_t_b);
	
            hipEventRecord(control, stre_t_b);
            hipStreamWaitEvent(0, control, 0);  

	    ghost_t_b_abp5<<<s_x_cb * s_y * s_z / 64, 64>>>(ghos_t_b, dest_g, U_t, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, rank, cb, flag, b);	 

        }

        MPI_Barrier(MPI_COMM_WORLD);

    }

    void Kernel(double *out, double *in) {
/*
        float eventMs = 1.0f;
        hipEvent_t start, stop;
        hipEventCreate(&start);
        hipEventCreate(&stop);
        hipEventRecord(start, 0);
*/
        const double a = -0.5 / (4 - rho);

        dslash(in, out, 0, 1, 1.0, -2 * a);
        dslash(in, out, 1, 1, 1.0, -2 * a);
/*
        hipEventRecord(stop, 0);
        hipEventSynchronize(stop);
        hipEventElapsedTime(&eventMs, start, stop);
        printf("dslash time taken  = %6.3fms\n", eventMs);
*/
    }

    void KernelSq_scaled(double *out, double *in, double cut) {

        const double sc1 = 2 / ((1 + 8 * kappa) * (1 + 8 * kappa) * (1 - cut));
        const double sc2 = (1 + cut) / (1 - cut);

        Kernel(tmp1, in);

        Kernel(tmp2, tmp1);

        axpbyz_g2<<<s_x * s_y * s_z * s_t * 24 / 1024, 1024 >>> (sc1, tmp2, -sc2, in, out);

    }

    void eps_l_g(double *out, double *in, int size) {

        const int size_f = s_x * s_y * s_z * s_t * 24;
        double *inner_g;
        hipMalloc((void **) &inner_g,   size * 2 * sizeof(double));

            double * &result_g = inner_g;
	
	    if(size_f == 24 * 24 * 24 * 16 * 24 ){	
	
		double * result_g1;

		hipMalloc((void **) &result_g1,  size * 2592 * 2 * sizeof(double));

		cDotProduct_v1_g2<<< size_f / 2 / 1024, 1024 >>> (result_g1, hw_evec_g, in, size_f / 2 , size, 2592);

		double * result_g2;

		hipMalloc((void **) &result_g2,  4 * size * 2 * sizeof(double));			

		cDotProduct_v2_g2<<< size * size_f / 2 / 1024 / 648, 648 >>> (result_g2, result_g1 , 2592  , size, 4, 648);						

		hipFree(result_g1);

		cDotProduct_v2_g2<<< size * size_f / 2 / 1024 / 648 / 4 / 100, 400 >>> (result_g, result_g2 , 4  , size, 1, 4);

		hipFree(result_g2);
	
	    }else if(size_f == 8 * 8 * 8 * 8 * 24){

		double * result_g1;

		hipMalloc((void **) &result_g1,  size * 48 * 2 * sizeof(double));

		cDotProduct_v1_g2<<< size_f / 2 / 1024, 1024 >>> (result_g1, hw_evec_g, in, size_f / 2 , size, 48);

		cDotProduct_v2_g2<<< size * size_f / 2 / 1024 / 48 / 10, 480 >>> (result_g, result_g1 , 48  , size, 1, 48);

		hipFree(result_g1);
	     	
	    } else if (size_f == 24 * 24 * 24 * 8 * 24 ){                                   

                double * result_g1;

                hipMalloc((void **) &result_g1,  size * 1296 * 2 * sizeof(double));

                cDotProduct_v1_g2<<< size_f / 2 / 1024, 1024 >>> (result_g1, hw_evec_g, in, size_f / 2 , size, 1296);

                double * result_g2;

                hipMalloc((void **) &result_g2, 2 * size * 2 * sizeof(double));

                cDotProduct_v2_g2<<< size * size_f / 2 / 1024 / 648, 648 >>> (result_g2, result_g1 , 1296  , size, 2, 648);

                hipFree(result_g1);

                cDotProduct_v2_g2<<< size * size_f / 2 / 1024 / 648 / 4 / 100, 400 >>> (result_g, result_g2 , 2  , size, 1, 2);

                hipFree(result_g2);

            }

            double result[2 * size];
            hipMemcpy(result, result_g, size * 2 * sizeof(double), hipMemcpyDeviceToHost);

            double inner[2 * size];

            MPI_Reduce(&result[0], &inner[0], 2 * size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            MPI_Bcast(&inner[0], 2 * size , MPI_DOUBLE, 0, MPI_COMM_WORLD);

            MPI_Barrier(MPI_COMM_WORLD);

            hipMemcpy(inner_g, inner, 2 * size * sizeof(double), hipMemcpyHostToDevice);

            cxmbyx_v1_g2<<<size_f / 2  / 768, 768 >>>(in, inner_g, hw_evec_g, size_f / 2, size);  	
	
	    mult<<<1, 2*size>>>(inner_g, sig_hw_eval_g);

	    cxpbyx_v1_g2<<<size_f / 2  / 768, 768 >>>(out, inner_g, hw_evec_g, size_f / 2, size);

            hipFree(inner_g);
    }

    void general_dov(double *out, double *in, double k0, double k1, double k2, double prec) {

        float eventMs = 1.0f;
        hipEvent_t start, stop;
        hipEventCreate(&start);
        hipEventCreate(&stop);
        hipEventRecord(start, 0);


        int is = -(int) log(prec / 0.2);
        if (is < 0)is = 0;
        if (is > coef.size() - 1)is = coef.size() - 1;

        int size_f = s_x * s_y * s_z * s_t * 24;

        double *src, *high;
        hipMalloc((void **) &src, size_f * sizeof(double));
        hipMemcpy(src, in, size_f * sizeof(double), hipMemcpyDeviceToDevice);

        hipMalloc((void **) &high, size_f * sizeof(double));
        hipMemset(high, 0, size_f * sizeof(double));

        hipMemset(out, 0, size_f * sizeof(double));

        eps_l_g(out, src, hw_size[is]);

        hipEventRecord(stop, 0);
        hipEventSynchronize(stop);
        hipEventElapsedTime(&eventMs, start, stop);
        if(rank==0) printf("Kernel 1 time taken  = %6.3fms\n", eventMs);

        hipEventCreate(&start);
        hipEventCreate(&stop);
        hipEventRecord(start, 0);


        double cut = pow(hw_eval[hw_size[is] - 1] / (1 + 8 * kappa), 2);

        double *pbn2, *pbn1, *pbn0, *ptmp;
        hipMalloc((void **) &pbn0, size_f * sizeof(double));
        hipMemset(pbn0, 0, size_f * sizeof(double));
        hipMalloc((void **) &pbn1, size_f * sizeof(double));
        hipMemset(pbn1, 0, size_f * sizeof(double));
        hipMalloc((void **) &pbn2, size_f * sizeof(double));
        hipMemset(pbn2, 0, size_f * sizeof(double));

        for (size_t i = coef[is].size() - 1; i >= 1; i--) {

            if (i < coef[is].size() - 1)KernelSq_scaled(high, pbn1, cut);

	    axpbyczw_g2<<<size_f / 1024 , 1024 >>>(2.0, high, -1.0, pbn0, coef[is][i], src, pbn2);  

            ptmp = pbn0;
            pbn0 = pbn1;
            pbn1 = pbn2;
            pbn2 = ptmp;
        }

        KernelSq_scaled(high, pbn1, cut);

	axpbyczw_g2<<<size_f / 1024 , 1024 >>>(1.0, high, -1.0, pbn0, coef[is][0], src, pbn2);

        Kernel(high, pbn2);

        axpbyz_g2<<< size_f / 1024, 1024  >>>(1.0 / (1 + 8 * kappa), high, 1.0, out, out);

        hipFree(pbn0);
        hipFree(pbn1);
        hipFree(pbn2);

        overlapLinop<<< size_f / 1024, 1024 >>> (out, in, k0, k1, k2);

        hipFree(src);
        hipFree(high);

        hipEventRecord(stop, 0);
        hipEventSynchronize(stop);
        hipEventElapsedTime(&eventMs, start, stop);
        if(rank==0) printf("Kernel 2 time taken  = %6.3fms\n", eventMs);
    }
};

void * newApplyOverlapQuda( std::vector<double *> &evecs, std::vector<double> &evals,
                      std::vector<std::vector<double> > &coefs, std::vector<int> &sizes,
                      double *U_x, double *U_y, double *U_z, double *U_t,
                      const int v_x, const int v_y, const int v_z, const int v_t,
                      const int s_x, const int s_y, const int s_z, const int s_t,
                      const double kappa  ){

	DiracOverlapWilson * ov_instance = new DiracOverlapWilson( evecs, evals, coefs, sizes, U_x, U_y, U_z, U_t, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, kappa);

	return  static_cast<void*>(ov_instance); 

}

void delApplyOverlapQuda(void *ov_instance){

	delete static_cast<DiracOverlapWilson*>(ov_instance);
}

void ApplyOverlapQuda(double *dest, double *src, double k0, double k1, double k2, double prec, void *ov_instance, int size){

DiracOverlapWilson *dirac=static_cast<DiracOverlapWilson*>(ov_instance);

double *src_g;
double *dest_g;

int size_f = size * 12 * 2;	

hipMalloc((void **) &src_g, size_f * sizeof(double));
hipMalloc((void **) &dest_g, size_f * sizeof(double));

hipMemcpy(src_g, src, size_f * sizeof(double), hipMemcpyHostToDevice);
hipMemset(dest_g, 0, size_f * sizeof(double));

dirac->general_dov(dest_g,src_g,k0,k1,k2,prec);

hipMemcpy(dest, dest_g, size_f * sizeof(double), hipMemcpyDeviceToHost);

hipFree(src_g);
hipFree(dest_g);

}

void * newApplyWilsonQuda(double *U_x, double *U_y, double *U_z, double *U_t,
             	          const int v_x, const int v_y, const int v_z, const int v_t,
                 	  const int s_x, const int s_y, const int s_z, const int s_t){

        DiracWilson * ov_instance = new DiracWilson(U_x, U_y, U_z, U_t, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t);

        return  static_cast<void*>(ov_instance);

}

void delApplyWilsonQuda(void *ov_instance){

        delete static_cast<DiracWilson*>(ov_instance);
}

void ApplyWilsonQuda(double *dest, double *src, void *ov_instance){

DiracWilson *dirac=static_cast<DiracWilson*>(ov_instance);

double *src_g;
double *dest_g;

int size_f = dirac->volume * 12 * 2;

hipMalloc((void **) &src_g, size_f * sizeof(double));
hipMalloc((void **) &dest_g, size_f * sizeof(double));

hipMemcpy(src_g, src, size_f * sizeof(double), hipMemcpyHostToDevice);
//hipMemset(dest_g, 0, size_f * sizeof(double));

dirac->Kernal(dest_g,src_g);

hipMemcpy(dest, dest_g, size_f * sizeof(double), hipMemcpyDeviceToHost);

hipFree(src_g);
hipFree(dest_g);

}

void ApplyOverlapQuda(double *dest, double *src, double k0, double k1, double k2, double prec,
                      std::vector<double *> &evecs, std::vector<double> &evals,
                      std::vector<std::vector<double> > &coefs, std::vector<int> &sizes,
                      double *U_x, double *U_y, double *U_z, double *U_t,
                      const int v_x, const int v_y, const int v_z, const int v_t,
                      const int s_x, const int s_y, const int s_z, const int s_t,
                      const double kappa) {

    double *src_g;
    double *dest_g;

    int size_f = s_x * s_y * s_z * s_t * 12 * 2;

    hipMalloc((void **) &src_g, size_f * sizeof(double));
    hipMalloc((void **) &dest_g, size_f * sizeof(double));

    hipMemcpy(src_g, src, size_f * sizeof(double), hipMemcpyHostToDevice);

    hipMemset(dest_g, 0, size_f * sizeof(double));

    DiracOverlapWilson OW(evecs, evals, coefs, sizes, U_x, U_y, U_z, U_t, v_x, v_y, v_z, v_t, s_x, s_y, s_z,
                          s_t, kappa);

    OW.general_dov(dest_g, src_g, k0, k1, k2, prec);

    hipMemcpy(dest, dest_g, size_f * sizeof(double), hipMemcpyDeviceToHost);

    hipFree(src_g);
    hipFree(dest_g);
}

void dslash_g5_kernal(double *src, double *dest, double *U_x, double *U_y, double *U_z, double *U_t,
                      const int v_x, const int v_y, const int v_z, const int v_t,
                      const int s_x, const int s_y, const int s_z, const int s_t,
                      const double a, const int flag) {

    test_2(src, dest, U_x, U_y, U_z, U_t, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, 0, flag);
    test_2(src, dest, U_x, U_y, U_z, U_t, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, 1, flag);

    axpbyz_g<<<dim3(s_x / 2, s_y / 2, s_z / 2), dim3(2, 2, 2) >>> (1.0, src, a, dest, dest, s_x, s_y, s_z, s_t);

    psi_g5<<<dim3(s_x / 2, s_y / 2, s_z / 2), dim3(2, 2, 2) >>> (dest, s_x, s_y, s_z, s_t);

}


int dslash_g5(double *src, double *dest, double *U_x, double *U_y, double *U_z, double *U_t,
              const int v_x, const int v_y, const int v_z, const int v_t,
              const int s_x, const int s_y, const int s_z, const int s_t,
              const double a, const int flag) {

    hipDeviceProp_t devProp;
    hipGetDeviceProperties(&devProp, 0);

//    std::cout << "Device name " << devProp.name << std::endl;

    double *src_g;
    double *dest_g;
    double *U_x_g;
    double *U_y_g;
    double *U_z_g;
    double *U_t_g;

    int size_f = s_x * s_y * s_z * s_t * 12 * 2;
    int size_u = s_x * s_y * s_z * s_t * 9 * 2;
//    allocate the memory on the device side
    hipMalloc((void **) &src_g, size_f * sizeof(double));
    hipMalloc((void **) &dest_g, size_f * sizeof(double));
    hipMalloc((void **) &U_x_g, size_u * sizeof(double));
    hipMalloc((void **) &U_y_g, size_u * sizeof(double));
    hipMalloc((void **) &U_z_g, size_u * sizeof(double));
    hipMalloc((void **) &U_t_g, size_u * sizeof(double));

    hipMemcpy(src_g, src, size_f * sizeof(double), hipMemcpyHostToDevice);
    hipMemcpy(U_x_g, U_x, size_u * sizeof(double), hipMemcpyHostToDevice);
    hipMemcpy(U_y_g, U_y, size_u * sizeof(double), hipMemcpyHostToDevice);
    hipMemcpy(U_z_g, U_z, size_u * sizeof(double), hipMemcpyHostToDevice);
    hipMemcpy(U_t_g, U_t, size_u * sizeof(double), hipMemcpyHostToDevice);
    hipMemset(dest_g, 0, size_f * sizeof(double));


    dslash_g5_kernal(src_g, dest_g, U_x_g, U_y_g, U_z_g, U_t_g, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, a, flag);

/*
    test_2(src_g, dest_g, U_x_g, U_y_g, U_z_g, U_t_g, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, 0, flag);
    test_2(src_g, dest_g, U_x_g, U_y_g, U_z_g, U_t_g, v_x, v_y, v_z, v_t, s_x, s_y, s_z, s_t, 1, flag);

    psi_g5<<<dim3(s_x/2, s_y/2, s_z/2), dim3(2, 2, 2) >>>(dest_g, s_x, s_y, s_z, s_t);	
*/

//    hipDeviceSynchronize();

    hipMemcpy(dest, dest_g, size_f * sizeof(double), hipMemcpyDeviceToHost);


    hipFree(src_g);
    hipFree(dest_g);
    hipFree(U_x_g);
    hipFree(U_y_g);
    hipFree(U_z_g);
    hipFree(U_t_g);

//    printf("test PASSED!\n");

    return 0;

}

