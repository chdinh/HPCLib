//=============================================================================================================
/**
* @file     RapMusic_kernel.cuh
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

#ifndef RAPMUSIC_KERNEL_CUH
#define RAPMUSIC_KERNEL_CUH


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
// Kernels and Device Functions
//=============================================================================================================

//=========================================================================================================
/**
* Calculates the combination indices Idx1 and Idx2 of n points.\n
*   [  (0,0)   (0,1)   (0,2)  ...  (0,n-1)   (0,n)\n
*              (1,1)   (1,2)  ...  (1,n-1)   (1,n)\n
*                      (2,2)  ...  (2,n-1)   (2,n)\n
*\n
*                                  (n-1,n-1) (n-1,n)\n
*                                            (n,n)]
*
*
* @param[in] p_iPoints  The number of points n which are combined with each other.
* @param[in] p_iCurIdx  The current combination index (between 0 and nchoosek(n+1,2))
* @param[out] p_iIdx1   The resulting index 1.
* @param[out] p_iIdx2   The resulting index 2.
*/
extern __device__ void cuGetPointPair(  const int p_iPoints,
                                        const int p_iCurIdx,
                                        int *p_pIdx1, int *p_pIdx2);


//=============================================================================================================
/**
* Pre-Calculates the Lead Field index combinations to search for a two dipole independent topography
* (IT = source).
*
* @param[in] p_iNumPoints   The number of Lead Field points -> for dimension check
* @param[in] p_iNumCombinations The number of pair index combinations.
* @param[out] p_pPairIdxCombinations    The destination which contains pointer of index combinations of
*                                       Lead Field indices -> Number of pointers = Combination (number of grid
*                                       points over 2 = Num + 1 C 2)
*/
extern __global__ void cuCalcPairCombinations(  int p_iNumPoints,
                                                int p_iNumCombinations,
                                                int* p_pPairIdxCombinations);



extern __global__ void cuCalcU_B(   float* p_pMatProj_Phi_s,
                                    int p_iProj_Phi_sRows,
                                    int p_iProj_Phi_sCols,
                                    float* p_pMatU_B,
                                    int* p_pRank);



//=============================================================================================================
/**
* ToDo
*
* @param[in] p_pMatProjLeadField    [p_iLeadFieldRows x p_iLeadFieldCols] ToDo
* @param[in] p_iLeadFieldRows   ToDo
* @param[out] p_pMatProj_G  [p_iLeadFieldRows x 6]
* @param[in] p_iIdx1    ToDo
* @param[in] p_iIdx2    ToDo
*/
extern __device__ void cuGetLeadFieldPair(  float* p_pMatProjLeadField, //Input
                                            int p_iLeadFieldRows,
                                            float* p_pMatProj_G,
                                            int p_iIdx1, int p_iIdx2);


//=============================================================================================================
/**
* ToDo
*
* @param[in] p_pMatProjLeadField    [p_iLeadFieldRows x p_iLeadFieldCols] ToDo
* @param[in] p_iLeadFieldRows   ToDo
* @param[in] p_iLeadFieldCols   ToDo
* @param[in] p_pPairIdxCombinations    ToDo
* @param[in] p_iNumOfCombinations    ToDo
* @param[in] p_pMatU_B  ToDo
* @param[in] p_iColsB   ToDo
* @param[out] p_pRoh    ToDo
*/
extern __global__ void RapMusicSubcorr( float* p_pMatProjLeadField, //Input
                                        int p_iLeadFieldRows,
                                        int p_iLeadFieldCols, 
                                        int* p_pPairIdxCombinations, //Combination
                                        int p_iNumOfCombinations,

                                        float* p_pMatU_B, //[rowsA x colsB] //from kernel part2
                                        int p_iColsB,
                                        float* p_pRoh);





extern __device__ int cuPowellOffset(int p_iRow, int p_iNumPoints);

extern __global__ void cuPowellIdxVec(int p_iRow, int p_iNumPoints, int* p_pVecElements);

extern __global__ void PowellRapMusicSubcorr(   float* p_pMatProjLeadField, //Input
                                                int p_iLeadFieldRows,
                                                int p_iLeadFieldCols,
                                                int* p_pPairIdxCombinations, //Combination
                                                int* p_pRowIndezes,
                                                int p_iNumOfDipoles,

                                                float* p_pMatU_B, //[rowsA x colsB] //from kernel part2
                                                int p_iColsB,
                                                float* p_pRoh );




}//Namespace

#endif
