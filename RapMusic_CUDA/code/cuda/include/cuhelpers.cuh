//=============================================================================================================
/**
* @file     cusvd.cuh
* @author   Christoph Dinh <christoph.dinh@live.de>; Johannes Ruehle <johannes.ruehle@tu-ilmenau.de>
* @version  1.0
* @date     March, 2011
*
* @section  LICENSE
*
* Copyright (C) 2011 Christoph Dinh. All rights reserved.
*
* No part of this program may be photocopied, reproduced,
* or translated to another program language without the
* prior written consent of the author.
*
*
* @brief    Contains a svd gpu implementation.
*
*/


#ifndef CUHELPERS_CUH
#define CUHELPERS_CUH

//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include <math.h>


//*************************************************************************************************************
//=============================================================================================================
// DEFINE NAMESPACE HPCLib
//=============================================================================================================

namespace HPCLib
{


//*************************************************************************************************************
//=============================================================================================================
// SOME DEFINES
//=============================================================================================================

#define SIGN(a,b) ( (b) >= 0.0 ? fabs(a) : -fabs(a) )   /**< A macro which makes \a a equal signed as \a b. */

#define DMAX(a,b) ( (a) > (b) ? (a) : (b) ) /**< A macro that returns the maximum of \a a and \a b. */
#define IMIN(a,b) ( (a) < (b) ? (a) : (b) ) /**< A macro that returns the minimum of \a a and \a b. */

//static float sqrarg;
#define SQR(a) (a*a)    /**< A macro which calculates the square of \a a. */


//*************************************************************************************************************
//=============================================================================================================
// Kernels and Device Functions
//=============================================================================================================

//=============================================================================================================
/**
* Computes the pythagoras a2 + b2 = c2 without destructive underflow or overflow.
*
* CODE of 3rd EDITION OF Numerical Recipes in C (2007)
*
* @param[in] a  Coefficient one.
* @param[in] b  Coefficient two.
* @return c of a2 + b2 = c2.
*/
__device__ float pythag(float a, float b)
{
    float absa,absb;
    absa=fabs(a);
    absb=fabs(b);
    if (absa > absb) 
        return absa*sqrt(1.0+SQR(absb/absa));
    else 
        return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}


//*************************************************************************************************************
//for (s=0.0,k=i;k<m;k++) s += u[k+i*m]*u[k+j*m];
__device__ void cuScalarProductYZ_shared(   float *p_pVec_a,
                                            float *p_pVec_b,
                                            int p_iNumElements,
                                            float* p_sharedYZCache,
                                            float *result )
{
    int tid = threadIdx.y+blockDim.y*threadIdx.z;
    int cacheIndex = threadIdx.y+blockDim.y*threadIdx.z;

    float   temp = 0;
    while (tid < p_iNumElements) {
        temp += p_pVec_a[tid] * p_pVec_b[tid];
        tid += blockDim.y*blockDim.z;
    }
    
    // set the cache values
    p_sharedYZCache[cacheIndex] = temp;
    
    // synchronize threads in this block
    __syncthreads();

    // for reductions, threadsPerBlock must be a power of 2
    // because of the following code
    int i = (blockDim.y*blockDim.z)/2;
    while (i != 0) {
        if (cacheIndex < i)
            p_sharedYZCache[cacheIndex] += p_sharedYZCache[cacheIndex + i];
        __syncthreads();
        i /= 2;
    }

    if (cacheIndex == 0)
        *result = p_sharedYZCache[0];
}


//*************************************************************************************************************
//for (s=0.0,k=i;k<m;k++) s += u[k+i*m]*u[k+j*m];
__device__ void cuScalarProductY_shared( float *p_pVec_a,
                                        float *p_pVec_b,
                                        int p_iNumElements,
                                        float* p_sharedYCache,
                                        float *result )
{
    int tid = threadIdx.y;
    int cacheIndex = threadIdx.y;

    float   temp = 0;
    while (tid < p_iNumElements) {
        temp += p_pVec_a[tid] * p_pVec_b[tid];
        tid += blockDim.y;
    }
    
    // set the cache values
    p_sharedYCache[cacheIndex] = temp;
    
    // synchronize threads in this block
    __syncthreads();

    // for reductions, threadsPerBlock must be a power of 2
    // because of the following code
    int i = blockDim.y/2;
    while (i != 0) {
        if (cacheIndex < i)
            p_sharedYCache[cacheIndex] += p_sharedYCache[cacheIndex + i];
        __syncthreads();
        i /= 2;
    }

    if (cacheIndex == 0)
        *result = p_sharedYCache[0];
}


//*************************************************************************************************************
//for (s=0.0,k=i;k<m;k++) s += u[k+i*m]*u[k+j*m];
__device__ void cuScalarProductAddZ_incr_shared(   float *p_pVec_a,
                                                float *p_pVec_b,//in,out
                                                int p_iNumElements,
                                                float* p_sharedZCache,
                                                float *result,
                                                int incr,
                                                float* p_pRV1
                                            )
{
    int tid = threadIdx.z;
    int cacheIndex = threadIdx.z;

    float   temp = 0;
    while (tid < p_iNumElements) {
        temp += p_pVec_a[tid*incr] * p_pVec_b[tid*incr];
        tid += blockDim.z;
    }
    
    // set the cache values
    p_sharedZCache[cacheIndex] = temp;
    
    // synchronize threads in this block
    __syncthreads();

    // for reductions, threadsPerBlock must be a power of 2
    // because of the following code
    int i = blockDim.z/2;
    while (i != 0) {
        if (cacheIndex < i)
            p_sharedZCache[cacheIndex] += p_sharedZCache[cacheIndex + i];
        __syncthreads();
        i /= 2;
    }


    tid = threadIdx.z;

    //Increment
    while (tid < p_iNumElements) {
        p_pVec_a[tid*incr] += p_sharedZCache[0]*p_pRV1[tid];
        tid += blockDim.z;
    }
}


//*************************************************************************************************************
//                 for (k=i;k<m;k++) {
//                     u[k+i*m] /= *t_pScale;
//                     s += u[k+i*m]*u[k+i*m];
//                 }
__device__ void cuVecScaleSquareReduceYZ_shared( float* p_pScale,
                                                float *p_pVec,//in,out(scaled vec)
                                                int p_iNumElements,
                                                float* p_sharedYZCache,
                                                float *sum//out
                                                )
{
    int tid = threadIdx.y+blockDim.y*threadIdx.z;
    int cacheIndex = threadIdx.y+blockDim.y*threadIdx.z;

    float   temp = 0;
    while (tid < p_iNumElements) {
        p_pVec[tid] /= *p_pScale;
        temp += p_pVec[tid] * p_pVec[tid];
        tid += blockDim.y*blockDim.z;
    }
    
    // set the cache values
    p_sharedYZCache[cacheIndex] = temp;
    
    // synchronize threads in this block
    __syncthreads();

    // for reductions, threadsPerBlock must be a power of 2
    // because of the following code
    int i = (blockDim.y*blockDim.z)/2;
    while (i != 0) {
        if (cacheIndex < i)
            p_sharedYZCache[cacheIndex] += p_sharedYZCache[cacheIndex + i];
        __syncthreads();
        i /= 2;
    }

    if (cacheIndex == 0)
        *sum = p_sharedYZCache[0];
}


//*************************************************************************************************************
//                 for (k=i;k<m;k++) {
//                     u[k+i*m] /= *t_pScale;
//                     s += u[k+i*m]*u[k+i*m];
//                 }
__device__ void cuVecScaleSquareReduceYZ_incr_shared(    float* p_pScale,
                                                        float *p_pVec,//in,out(scaled vec)
                                                        int p_iNumElements,
                                                        float* p_sharedYZCache,
                                                        float *sum,//out
                                                        int incr
                                                    )
{
    int tid = threadIdx.y+blockDim.y*threadIdx.z;
    int cacheIndex = threadIdx.y+blockDim.y*threadIdx.z;

    float   temp = 0;
    while (tid < p_iNumElements) {
        p_pVec[tid*incr] /= *p_pScale;
        temp += p_pVec[tid*incr] * p_pVec[tid*incr];
        tid += blockDim.y*blockDim.z;
    }
    
    // set the cache values
    p_sharedYZCache[cacheIndex] = temp;
    
    // synchronize threads in this block
    __syncthreads();

    // for reductions, threadsPerBlock must be a power of 2
    // because of the following code
    int i = (blockDim.y*blockDim.z)/2;
    while (i != 0) {
        if (cacheIndex < i)
            p_sharedYZCache[cacheIndex] += p_sharedYZCache[cacheIndex + i];
        __syncthreads();
        i /= 2;
    }

    if (cacheIndex == 0)
        *sum = p_sharedYZCache[0];
}


//*************************************************************************************************************
// for (k=i;k<m;k++) *t_pScale += abs(u[k+i*m]);
__device__ void cuReduceAbsYZ_shared( float *p_pVec, int p_iNumElements, float* p_sharedYZCache, float *result )
{
    int tid = threadIdx.y+blockDim.y*threadIdx.z;
    int cacheIndex = threadIdx.y+blockDim.y*threadIdx.z;

    float   temp = 0;
    while (tid < p_iNumElements) {
        temp += abs(p_pVec[tid]);
        tid += blockDim.y*blockDim.z;
    }
    
    // set the cache values
    p_sharedYZCache[cacheIndex] = temp;

    // for reductions, threadsPerBlock must be a power of 2
    // because of the following code
    int i = (blockDim.y*blockDim.z)/2;
    while (i != 0) {
        if (cacheIndex < i)
            p_sharedYZCache[cacheIndex] += p_sharedYZCache[cacheIndex + i];
        __syncthreads();
        i /= 2;
    }

    if (cacheIndex == 0)
        *result = p_sharedYZCache[0];
}

}//Namespace

#endif //CUHELPERS_CUH
