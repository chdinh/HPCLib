//=============================================================================================================
/**
* @file     rapmusic_kernel.cu
* @author   Christoph Dinh <christoph.dinh@live.de>;
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
* @brief    ToDo Documentation.. 
*
*/

//*************************************************************************************************************
//=============================================================================================================
// CUDA INCLUDES
//=============================================================================================================

#include "../include/rapmusic_kernel.cuh"
#include "../include/cusvd.cuh"
#include "../include/cuhelpers.cuh"


//*************************************************************************************************************
//=============================================================================================================
// STL INCLUDES
//=============================================================================================================


#define SHDEBUG


//*************************************************************************************************************
//=============================================================================================================
// DEFINE NAMESPACE HPCLib
//=============================================================================================================

namespace HPCLib
{


//*************************************************************************************************************
//=============================================================================================================
// Kernels and Device Functions
//=============================================================================================================

__device__ void cuGetPointPair( const int p_iPoints,
                                const int p_iCurIdx,
                                int *p_pIdx1, int *p_pIdx2)
{
    int ii = p_iPoints*(p_iPoints+1)/2-1-p_iCurIdx;
    int K = (int)floor((sqrt((double)(8*ii+1))-1)/2);
    
    *p_pIdx1 = p_iPoints-1-K;

    *p_pIdx2 = (p_iCurIdx-p_iPoints*(p_iPoints+1)/2 + (K+1)*(K+2)/2)+(*p_pIdx1);
}


//*************************************************************************************************************

__global__ void cuCalcPairCombinations( int p_iNumPoints,
                                        int p_iNumCombinations,
                                        int* p_pPairIdxCombinations)
{
    int t_iCombIdx = threadIdx.x + blockIdx.x * blockDim.x;

    int *t_pIdx1 = new int;
    int *t_pIdx2 = new int;

    while (t_iCombIdx < p_iNumCombinations) {

        cuGetPointPair(p_iNumPoints, t_iCombIdx, t_pIdx1, t_pIdx2);

        p_pPairIdxCombinations[2*t_iCombIdx] = *t_pIdx1;
        p_pPairIdxCombinations[2*t_iCombIdx+1] = *t_pIdx2;

        t_iCombIdx += gridDim.x*blockDim.x;
    }

    delete t_pIdx2;
    delete t_pIdx1;
}


//*************************************************************************************************************

__global__ void cuCalcU_B(  float* p_pMatProj_Phi_s,
                            int p_iProj_Phi_sRows,
                            int p_iProj_Phi_sCols,
                            float* p_pMatU_B,
                            int* p_pRank)
{
    int t_iSizeProj_Phi_sMat = p_iProj_Phi_sRows * p_iProj_Phi_sCols;
    
    extern __shared__ float t_pSharedMem[];

    float* t_pU_B = t_pSharedMem;
    float* t_pW = t_pU_B + t_iSizeProj_Phi_sMat;
    float* t_pYZCache = t_pW + p_iProj_Phi_sCols;
    float* t_pSVDCache_all = t_pYZCache + blockDim.y*blockDim.z;

     if(threadIdx.z == 0)//To prevent other z threads performing memory access
     {
        //Copy first Lead Field point
        int i = threadIdx.y;//combination index
        while(i < t_iSizeProj_Phi_sMat)
        {
            t_pU_B[i] = p_pMatProj_Phi_s[i];
            i += blockDim.y;
        }
     }

    cuSVD_UW_shared(    t_pU_B,   /* [m x n ]*/
                        p_iProj_Phi_sRows,      /* rows */
                        p_iProj_Phi_sCols,      /* columns */
                        t_pW,   /* [nx1]*/
                        t_pSVDCache_all, /* [nx1] */
                        t_pYZCache );

    *p_pRank = 0;

    int t_iKey = 0;
    float t_vMax = 0;

    for(int n = 0; n < p_iProj_Phi_sCols; ++n)
    {
        if (t_pW[n] >= 0.00001f)
        {
            if(t_vMax < t_pW[n])
            {
                t_iKey = n;
                t_vMax = t_pW[n];
            }
            ++(*p_pRank);
        }
    }

    
    //order nonzero Singular values
    int* vecOrdKey = new int[*p_pRank];
    float* vecOrdVal = new float[*p_pRank];

    vecOrdKey[0] = t_iKey;
    vecOrdVal[0] = t_vMax;

    //very uneffective key-value-sorting
    for(int k = 1; k < *p_pRank; ++k)
    {
        vecOrdVal[k] = 0;
        for(int n = 0; n < p_iProj_Phi_sCols; ++n)
        {
            if( t_pW[n] <= vecOrdVal[k-1] && n != vecOrdKey[k-1] && t_pW[n] > vecOrdVal[k] )
            {
                vecOrdKey[k] = n;
                vecOrdVal[k] = t_pW[n];
            }
        }
    }


    int c = 0;
    for(int n = 0; n < *p_pRank; ++n)//ToDo Parallel
    {
        int m = threadIdx.y;
        while(m < p_iProj_Phi_sRows)
        {
            p_pMatU_B[c*p_iProj_Phi_sRows + m] = t_pU_B[vecOrdKey[n]*p_iProj_Phi_sRows + m]; 
            m += blockDim.y;
        }
        ++c;
    }

    __syncthreads();
}

//*************************************************************************************************************

__device__ void cuGetLeadFieldPair( float* p_pMatProjLeadField, //Input
                                    int p_iLeadFieldRows,
                                    float* p_pMatProj_G,
                                    int p_iIdx1, int p_iIdx2)
{
     if(threadIdx.z == 0)//To prevent other z threads performing memory access
     {
        int iidx1 = p_iIdx1 * p_iLeadFieldRows*3;//index with offset (idx1 * p_iLeadFieldRows)
        int iidx2 = p_iIdx2 * p_iLeadFieldRows*3;//index with offset (idx2 * p_iLeadFieldRows)

        int t_iSizePointMat = p_iLeadFieldRows*3;

        //Copy first Lead Field point
        int i = threadIdx.y;//combination index
        while(i < t_iSizePointMat)
        {
            p_pMatProj_G[i] = p_pMatProjLeadField[iidx1 + i];
            i += blockDim.y;
        }

        //Copy second Lead Field point
        i = threadIdx.y;//combination index
        while(i < t_iSizePointMat)
        {
            p_pMatProj_G[t_iSizePointMat+i] = p_pMatProjLeadField[iidx2 + i];
            i += blockDim.y;
        }
    }
}


//*************************************************************************************************************

__global__ void RapMusicSubcorr(    float* p_pMatProjLeadField, //Input
                                    int p_iLeadFieldRows,
                                    int p_iLeadFieldCols, 
                                    int* p_pPairIdxCombinations, //Combination
                                    int p_iNumOfCombinations,

                                    float* p_pMatU_B, //[rowsA x colsB] //from kernel part2
                                    int p_iColsB,
                                    float* p_pRoh )
{
    const int t_iSizePairMat = p_iLeadFieldRows * 6;

    const int t_iPairCols = 6;
    const int iColsA = 6;
    const int t_iSizeSVDCache = t_iPairCols+1+1;

    const int t_iSizeCorMat = iColsA * p_iColsB;

    //Create all Pair Mats in shared mem
    extern __shared__ float t_pSharedMem[];
    
    //Split shared Memory
    float* t_pMatProj_G_all = t_pSharedMem;//size = t_iSizePairMat*blockDim.x
    float* t_pW_all = t_pMatProj_G_all + t_iSizePairMat*blockDim.x;//size = t_iPairCols*blockDim.x
    float* t_pCor_all = t_pW_all + t_iPairCols*blockDim.x;//size = t_iSizeCorMat*blockDim.x
    float* t_pCacheYZ_all = t_pCor_all + t_iSizeCorMat*blockDim.x;//size = blockDim.y*blockDim.z*blockDim.x
    float* t_pSVDCache_all = t_pCacheYZ_all + blockDim.y*blockDim.z*blockDim.x;//size = t_iSizeSVDCache*blockDim.x

    //Split YZ Cache
    float* t_pCacheYZ = t_pCacheYZ_all +(threadIdx.x*blockDim.y*blockDim.z);
    float* t_pCacheY = t_pCacheYZ + threadIdx.z*blockDim.y;

    int t_iNumPairMatsPerBlock = blockDim.x;//Number of threads in x
    int t_iCurrentBlock = blockIdx.x;

    int t_iCombIdx = threadIdx.x + t_iCurrentBlock * t_iNumPairMatsPerBlock;


    while (t_iCombIdx < p_iNumOfCombinations) {

        int idx1 = p_pPairIdxCombinations[2*t_iCombIdx];//*3, 3x more cols -> x y z
        int idx2 = p_pPairIdxCombinations[2*t_iCombIdx+1];//*3, 3x more cols -> x y z

        float* t_pMatProj_G = t_pMatProj_G_all + t_iSizePairMat*threadIdx.x;

        cuGetLeadFieldPair( p_pMatProjLeadField, //Input
                            p_iLeadFieldRows,
                            t_pMatProj_G,
                            idx1, idx2);
        __syncthreads();

    //Part 1
        //führe svd auf paarmat aus
        float* w = t_pW_all + t_iPairCols*threadIdx.x;
        float* t_pSVDCache = t_pSVDCache_all + t_iSizeSVDCache*threadIdx.x;

        cuSVD_UW_shared( t_pMatProj_G, p_iLeadFieldRows, t_iPairCols, w, t_pSVDCache, t_pCacheYZ);
//        cuSVD_UW( t_pMatProj_G, p_iLeadFieldRows, t_iPairCols, w);
        __syncthreads();

        float* t_pMatU_A_full = t_pMatProj_G;

        //if once a singularvalue is smaller than epsilon = 10^-5 the following values are also smaller
        // -> because Singular values are ordered
        int t_iRank = threadIdx.z;
        while(t_iRank < t_iPairCols)
        {
            if (w[t_iRank] < 0.00001f)//set Eigenvectors with corresponding 0 eigenvalue to 0
            {
                int c = threadIdx.y;
                while(c < p_iLeadFieldRows)
                {
                    t_pMatU_A_full[t_iRank*p_iLeadFieldRows + c] = 0.0f; 
                    c += blockDim.y;
                }
            }
            t_iRank += blockDim.z;
        }
        __syncthreads();


    //Part 2
        float* U_B = p_pMatU_B;

        float* Cor = t_pCor_all + t_iSizeCorMat*threadIdx.x;//new float[t_iSizeCorMat];//p_pMatCor+(t_iCombIdx * t_iSizeCorMat);

        //lt. Mosher 1998: C = U_A^T * U_B

        //Cor.cols() >= Cor.rows() == U_B.cols > U_A.cols
        if(p_iColsB >= iColsA) //Bug ID 1 - fixed! changed from > to >=
        {
            //C = U_B^T * U_A
            for( int cA = 0; cA < iColsA; ++cA )
            {
                float* t_pMatU_A_full_cur = t_pMatU_A_full + (cA*p_iLeadFieldRows);

                int cB = threadIdx.z;
                while(cB < p_iColsB)
                {
                    float* U_B_cur = U_B + (cB*p_iLeadFieldRows);
                    float* t_pCor_cur = Cor + (cA*p_iColsB + cB);
                    *t_pCor_cur = 0;

                    cuScalarProductY_shared( U_B_cur, t_pMatU_A_full_cur, p_iLeadFieldRows, t_pCacheY, t_pCor_cur );
                    __syncthreads();
//                     for( int rAB = 0; rAB < p_iLeadFieldRows; ++rAB)
//                         *t_pCor_cur += U_B_cur[rAB]*t_pMatU_A_full_cur[rAB];

                    cB += blockDim.z;
                }
                __syncthreads();
            }
        }
        else//ToDo to debug
        {
            //C = U_A^T * U_B
            for( int cB = 0; cB < p_iColsB; ++cB )
            {
                float* U_B_cur = U_B + (cB*p_iLeadFieldRows);

                int cA = threadIdx.z;
                while(cA < iColsA)
                {
                    float* t_pMatU_A_full_cur = t_pMatU_A_full + (cA*p_iLeadFieldRows);
                    float* t_pCor_cur = Cor + (cB*iColsA+cA);
                    *t_pCor_cur = 0;

                    cuScalarProductY_shared( t_pMatU_A_full_cur, U_B_cur, p_iLeadFieldRows, t_pCacheY, t_pCor_cur );
                    __syncthreads();
//                     for( int rAB = 0; rAB < p_iLeadFieldRows; ++rAB)
//                         *t_pCor_cur += t_pMatU_A_full_cur[rAB]*U_B_cur[rAB];

                    cA += blockDim.z;
                }
                __syncthreads();
            }
        }

    //Part 3
        int rows = p_iColsB;
        int cols = iColsA;

        if (p_iColsB < iColsA)
        {
            rows = iColsA;
            cols = p_iColsB;
        }

        //cols are maximal iColsA = 6.
        //That's why we can use w and shared cache again. they are 6 width
        cuSVD_W_shared( Cor, rows, cols, w, t_pSVDCache, t_pCacheYZ);
//        cuSVD_W( Cor, rows, cols, w);
        __syncthreads();
        
        p_pRoh[t_iCombIdx] = 0;
        if(threadIdx.y == 0 && threadIdx.z == 0)
            for(int i = 0; i < cols; ++i)
                if (p_pRoh[t_iCombIdx] < w[i])
                    p_pRoh[t_iCombIdx] = w[i];

        __syncthreads();

         t_iCombIdx += gridDim.x*t_iNumPairMatsPerBlock;
    }

    __syncthreads();
}






__device__ void cuPowellOffset( float* p_pMatProjLeadField, //Input
                                    int p_iLeadFieldRows,
                                    float* p_pMatProj_G,
                                    int p_iIdx1, int p_iIdx2)
{

}


//*************************************************************************************************************

__device__ int cuPowellOffset(int p_iRow, int p_iNumPoints)
{

    return p_iRow*p_iNumPoints - (( (p_iRow-1)*p_iRow) / 2); //triangular series 1 3 6 10 ... = (num_pairs*(num_pairs+1))/2

}


//*************************************************************************************************************

__global__ void cuPowellIdxVec(int p_iRow, int p_iNumPoints, int* p_pVecElements)
{

    //     if(p_pVecElements != NULL)
    //         delete[] p_pVecElements;
    // 
    //     p_pVecElements = new int(p_iNumPoints);

    //col combination index
    int t_iIdx = threadIdx.x + blockIdx.x * blockDim.x;


    while (t_iIdx < p_iRow) {
        p_pVecElements[t_iIdx] = cuPowellOffset(t_iIdx+1,p_iNumPoints)-(p_iNumPoints-p_iRow);
        t_iIdx += gridDim.x*blockDim.x;
    }

    //row combination index
    int off = cuPowellOffset(p_iRow,p_iNumPoints);
    int length = p_iNumPoints - p_iRow;

    t_iIdx = threadIdx.x + blockIdx.x * blockDim.x;

    while (t_iIdx+p_iRow < p_iRow+length) {
        p_pVecElements[t_iIdx+p_iRow] = off+t_iIdx;
        t_iIdx += gridDim.x*blockDim.x;
    }
}


//*************************************************************************************************************

__global__ void PowellRapMusicSubcorr(  float* p_pMatProjLeadField, //Input
                                        int p_iLeadFieldRows,
                                        int p_iLeadFieldCols,
                                        int* p_pPairIdxCombinations, //Combination
                                        int* p_pRowIndezes,
                                        int p_iNumOfDipoles,

                                        float* p_pMatU_B, //[rowsA x colsB] //from kernel part2
                                        int p_iColsB,
                                        float* p_pRoh )
{
    const int t_iSizePairMat = p_iLeadFieldRows * 6;

    const int t_iPairCols = 6;
    const int iColsA = 6;
    const int t_iSizeSVDCache = t_iPairCols+1+1;

    const int t_iSizeCorMat = iColsA * p_iColsB;

    //Create all Pair Mats in shared mem
    extern __shared__ float t_pSharedMem[];
    
    //Split shared Memory
    float* t_pMatProj_G_all = t_pSharedMem;//size = t_iSizePairMat*blockDim.x
    float* t_pW_all = t_pMatProj_G_all + t_iSizePairMat*blockDim.x;//size = t_iPairCols*blockDim.x
    float* t_pCor_all = t_pW_all + t_iPairCols*blockDim.x;//size = t_iSizeCorMat*blockDim.x
    float* t_pCacheYZ_all = t_pCor_all + t_iSizeCorMat*blockDim.x;//size = blockDim.y*blockDim.z*blockDim.x
    float* t_pSVDCache_all = t_pCacheYZ_all + blockDim.y*blockDim.z*blockDim.x;//size = t_iSizeSVDCache*blockDim.x

    //Split YZ Cache
    float* t_pCacheYZ = t_pCacheYZ_all +(threadIdx.x*blockDim.y*blockDim.z);
    float* t_pCacheY = t_pCacheYZ + threadIdx.z*blockDim.y;

    int t_iNumPairMatsPerBlock = blockDim.x;//Number of threads in x
    int t_iCurrentBlock = blockIdx.x;

    int t_iCombIdx = threadIdx.x + t_iCurrentBlock * t_iNumPairMatsPerBlock;


    while (t_iCombIdx < p_iNumOfDipoles) {
        int t_iCurrentIdx = p_pRowIndezes[t_iCombIdx];

        int idx1 = p_pPairIdxCombinations[2*t_iCurrentIdx];//*3, 3x more cols -> x y z
        int idx2 = p_pPairIdxCombinations[2*t_iCurrentIdx+1];//*3, 3x more cols -> x y z

        float* t_pMatProj_G = t_pMatProj_G_all + t_iSizePairMat*threadIdx.x;

        cuGetLeadFieldPair( p_pMatProjLeadField, //Input
                            p_iLeadFieldRows,
                            t_pMatProj_G,
                            idx1, idx2);
        __syncthreads();

    //Part 1
        //führe svd auf paarmat aus
        float* w = t_pW_all + t_iPairCols*threadIdx.x;
        float* t_pSVDCache = t_pSVDCache_all + t_iSizeSVDCache*threadIdx.x;

        cuSVD_UW_shared( t_pMatProj_G, p_iLeadFieldRows, t_iPairCols, w, t_pSVDCache, t_pCacheYZ);
//        cuSVD_UW( t_pMatProj_G, p_iLeadFieldRows, t_iPairCols, w);
        __syncthreads();

        float* t_pMatU_A_full = t_pMatProj_G;

        //if once a singularvalue is smaller than epsilon = 10^-5 the following values are also smaller
        // -> because Singular values are ordered
        int t_iRank = threadIdx.z;
        while(t_iRank < t_iPairCols)
        {
            if (w[t_iRank] < 0.00001f)//set Eigenvectors with corresponding 0 eigenvalue to 0
            {
                int c = threadIdx.y;
                while(c < p_iLeadFieldRows)
                {
                    t_pMatU_A_full[t_iRank*p_iLeadFieldRows + c] = 0.0f; 
                    c += blockDim.y;
                }
            }
            t_iRank += blockDim.z;
        }
        __syncthreads();


    //Part 2
        float* U_B = p_pMatU_B;

        float* Cor = t_pCor_all + t_iSizeCorMat*threadIdx.x;//new float[t_iSizeCorMat];//p_pMatCor+(t_iCombIdx * t_iSizeCorMat);

        //lt. Mosher 1998: C = U_A^T * U_B

        //Cor.cols() >= Cor.rows() == U_B.cols > U_A.cols
        if(p_iColsB >= iColsA) //Bug ID 1 - fixed! changed from > to >=
        {
            //C = U_B^T * U_A
            for( int cA = 0; cA < iColsA; ++cA )
            {
                float* t_pMatU_A_full_cur = t_pMatU_A_full + (cA*p_iLeadFieldRows);

                int cB = threadIdx.z;
                while(cB < p_iColsB)
                {
                    float* U_B_cur = U_B + (cB*p_iLeadFieldRows);
                    float* t_pCor_cur = Cor + (cA*p_iColsB + cB);
                    *t_pCor_cur = 0;

                    cuScalarProductY_shared( U_B_cur, t_pMatU_A_full_cur, p_iLeadFieldRows, t_pCacheY, t_pCor_cur );
                    __syncthreads();
//                     for( int rAB = 0; rAB < p_iLeadFieldRows; ++rAB)
//                         *t_pCor_cur += U_B_cur[rAB]*t_pMatU_A_full_cur[rAB];

                    cB += blockDim.z;
                }
                __syncthreads();
            }
        }
        else//ToDo to debug
        {
            //C = U_A^T * U_B
            for( int cB = 0; cB < p_iColsB; ++cB )
            {
                float* U_B_cur = U_B + (cB*p_iLeadFieldRows);

                int cA = threadIdx.z;
                while(cA < iColsA)
                {
                    float* t_pMatU_A_full_cur = t_pMatU_A_full + (cA*p_iLeadFieldRows);
                    float* t_pCor_cur = Cor + (cB*iColsA+cA);
                    *t_pCor_cur = 0;

                    cuScalarProductY_shared( t_pMatU_A_full_cur, U_B_cur, p_iLeadFieldRows, t_pCacheY, t_pCor_cur );
                    __syncthreads();
//                     for( int rAB = 0; rAB < p_iLeadFieldRows; ++rAB)
//                         *t_pCor_cur += t_pMatU_A_full_cur[rAB]*U_B_cur[rAB];

                    cA += blockDim.z;
                }
                __syncthreads();
            }
        }

    //Part 3
        int rows = p_iColsB;
        int cols = iColsA;

        if (p_iColsB < iColsA)
        {
            rows = iColsA;
            cols = p_iColsB;
        }

        //cols are maximal iColsA = 6.
        //That's why we can use w and shared cache again. they are 6 width
        cuSVD_W_shared( Cor, rows, cols, w, t_pSVDCache, t_pCacheYZ);
//        cuSVD_W( Cor, rows, cols, w);
        __syncthreads();
        
        p_pRoh[t_iCurrentIdx] = 0;
        if(threadIdx.y == 0 && threadIdx.z == 0)
            for(int i = 0; i < cols; ++i)
                if (p_pRoh[t_iCurrentIdx] < w[i])
                    p_pRoh[t_iCurrentIdx] = w[i];

        __syncthreads();

         t_iCombIdx += gridDim.x*t_iNumPairMatsPerBlock;
    }

    __syncthreads();
}

}//Namespace