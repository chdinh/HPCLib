//=============================================================================================================
/**
* @file     rapmusic_cuda.cuh
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


#ifndef RAPMUSIC_CUDA_CUH
#define RAPMUSIC_CUDA_CUH


//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "../../cpp/include/hpcmatrix.h"


//*************************************************************************************************************
//=============================================================================================================
// CUDA INCLUDES
//=============================================================================================================

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include "../include/cuhpcmatrix.cuh"


//*************************************************************************************************************
//=============================================================================================================
// STL INCLUDES
//=============================================================================================================

#include <iostream>


//*************************************************************************************************************
//=============================================================================================================
// DEFINE NAMESPACE HPCLib
//=============================================================================================================

namespace HPCLib
{

//*************************************************************************************************************
//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================

template<class T>
class Model;

// template<class T>
// class cuHPCMatrix;

template<class T>
class RapDipoles;

class CudaDevice;


//*************************************************************************************************************
//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================


//=============================================================================================================
/**
* DECLARE CLASS RapMusic_Cuda
*
* @brief ToDo
*/
    //template<class T>
class RapMusic_Cuda
{
public:
    RapMusic_Cuda();
    ~RapMusic_Cuda();

    //####Chris####
    bool initRAPMusic(  HPCLib::CudaDevice* p_pDeviceInfo,
                        HPCLib::Model<float>* p_pModel,
                        bool p_bSparsed = false,
                        int p_iN = 2,
                        double p_dThr = 0.5);

    bool initRAPMusic(  HPCLib::CudaDevice* p_pDeviceInfo,
                        HPCMatrix<float>* p_pMatLeadField,
                        HPCMatrix<float>* p_pMatGrid,
                        int p_iN, double p_dThr);

    bool calcRapMusic(  HPCMatrix<float>* p_pMatMeasurement, RapDipoles<float>*& p_pRapDipoles);


    bool calcPowellRAPMusic(  HPCMatrix<float>* p_pMatMeasurement, RapDipoles<float>*& p_pRapDipoles);



    //=========================================================================================================
    /**
    * Calculates Combination of (n over 2 -> nchoosek(n,2) with)
    *
    * @param n number of things to combinate with (n over 2)
    * @return the number of combinations
    */
    int nchoose2(int n);


protected:

// 	//template <T>
// 	static __global__ void getMaxCorPair( )
// 	{
// 	}//float/*T*/ *p_pMatProjLeadField, /*Pair*/int **p_ppPairIdxCombinations, float/*T*/ *p_pMatU_B );

// 	__device__ void subcorr(float/*T*/ *p_pMatProjLeadField,);


    //=========================================================================================================
    /**
    * Computes the signal subspace Phi_s out of the measurement F.
    *
    * @param[in] p_pMatMeasurement  The current measured data to process (for best performance it should have
                                    the dimension channels x samples with samples = number of channels)
    * @param[out] p_pMatPhi_s   The calculated signal subspace.
    * @return   The rank of the measurement F (named r lt. Mosher 1998, 1999)
    */
    int calcPhi_s(const HPCMatrix<float>& p_pMatMeasurement, cuHPCMatrix<float>* &p_dev_pMatPhi_s);


    //=========================================================================================================
    /**
    * Computes the subspace correlation between the projected G_rho and the projected signal subspace Phi_s, as
    * well as the resulting direction.
    * For speed-up: we calculate the decomposition of the projected Phi_s before this function. So the argument
    * for this function is U_B instead of Phi_s.
    *
    * @param[in] p_matProj_G    The projected Lead Field combination. This is a m x 6 matrix composed of the
    *                           Lead Field combination of two Points for all m channels and 3 orthogonal
    *                           components (x y z).
    * @param[in] p_matU_B    The matrix U is the subspace projection of the orthogonal projected Phi_s
    * @param[out] p_vec_phi_k_1 Returns the orientation for a correlated dipole pair.
                                (phi_x1, phi_y1, phi_z1, phi_x2, phi_y2, phi_z2)
    * @return   The maximal correlation c_1 of the subspace correlation of the current projected Lead Field
    *           combination and the projected measurement.
    */
    float subcorr(HPCMatrix<float>& p_matProj_G, HPCMatrix<float>& p_matU_B, HPCMatrix<float>& p_vec_phi_k_1);





    //=========================================================================================================
    /**
    * Calculates the accumulated manifold vectors A_{k1}
    *
    * @param[in] p_matG_k_1 The Lead Field combination for the currently best correlated pair.
    * @param[in] p_matPhi_k_1   Is equal to u_k_1 in the paper and it is the direction of the currently best
    *                           correlated pair.
    * @param[in] p_iIdxk_1  The current position in the manifold vector array A_k_1
    * @param[out] p_matA_k_1    The array of the manifold vectors.
    */
    void calcA_k_1( const HPCMatrix<float>& p_matG_k_1,
                    const HPCMatrix<float>& p_matPhi_k_1,
                    const int p_iIdxk_1,
                    HPCMatrix<float>& p_matA_k_1);


    //=========================================================================================================
    /**
    * Calculates the orthogonal projector Phi_A_k_1 like in the paper Mosher 1999 (13)
    *
    * @param[in] p_matA_k_1 The array of the manifold vectors.
    * @param[out] p_matOrthProj The orthogonal projector.
    */
    void calcOrthProj(const HPCMatrix<float>& p_matA_k_1, HPCMatrix<float>& p_matOrthProj);




    void getLeadFieldPair(  HPCMatrix<float>& p_matLeadField,
                            HPCMatrix<float>& p_matLeadField_Pair,
                            int p_iIdx1, int p_iIdx2);

private:
    //CUDA Device properties
    int m_iMultiProcessorCount;                 /**< Number of multiprocessors on device */
    const int m_iMaxBlocksPerMultiProcessor;    /**< Maximum number of resident blocks per multiprocessor */
    int m_iWarpSize;                            /**< Warp size in threads */
    int m_iMaxThreadsPerMultiProcessor;         /**< Maximum resident threads per multiprocessor */
    int m_iSharedMemoryPerMultiProcessor;       /**< Shared memory available per block in bytes */

    //Sizes
    const int m_iPairCols;


    //device memory ptr



    cuHPCMatrix<float>* m_dev_pLeadFieldMat;

    //host memory ptr
    float *m_host_pLeadFieldMat;//Attention also used by ####Chris####
    HPCMatrix<float>* m_pMatLeadField;  /**< Holds the given Lead Field.*/ //ToDo replace m_host_pLeadFieldMat with this data pointer

    float *m_host_pResultMat;


    //RAP MUSIC stuff

    HPCMatrix<float>* m_pMatGrid;   /**< Holds the Grid which corresponds to the Lead Field. This variable can be
                                         NULL, then no coordinates are assigned to the located dipoles.*/


    int m_iN;               /**< Number of Sources to find*/
    double m_dThreshold;    /**< Threshold which defines the minimal correlation. Is the correlation of
                                 the found dipole pair smaller as this threshold than the RAP MUSIC
                                 calculation is stopped. */

    //Measurement Set-Up
    int m_iNumGridPoints;               /**< Number of Grid points. */
    int m_iNumChannels;                 /**< Number of channels */
    int m_iNumLeadFieldCombinations;    /**< Number of Lead Filed combinations (grid points + 1 over 2)*/

    //Preallocated combinations accessing by idx1 = m_pPairIdxCombinations[2*i]; idx2 = m_pPairIdxCombinations[(2*i)+1];
    thrust::device_vector<int> *m_dev_pVecPairIdxCombinations; //Idx11, Idx12, Idx21, Idx22, Idx31, ..., IdxN1, IdxN2
    int *m_dev_pPairIdxCombinations;

    bool m_bIsInit; /**< Wether the algorithm is initialized. */


    //=========================================================================================================
    /**
    * Returns the rank r of a singular value matrix based on non-zero singular values
    * (singular value > epsilon = 10^-5) 
    *
    * @param[in] p_matSigma diagonal matrix which contains the Singular values (Dimension n x n)
    * @return The rank r.
    */
    inline int getRank(const float* p_matSigma, int p_iSize)
    {
        int t_iRank;
        //if once a singularvalue is smaller than epsilon = 10^-5 the following values are also smaller
        // -> because Singular values are ordered
        for(t_iRank = p_iSize-1; t_iRank > 0; t_iRank--)
            if (p_matSigma[t_iRank] > 0.00001)
                break;

        t_iRank++;//rank corresponding to epsilon

        return t_iRank;
    }


    inline int getRank(HPCMatrix<float>* p_matSigma)
    {
        int t_iRank;
        //if once a singularvalue is smaller than epsilon = 10^-5 the following values are also smaller
        // -> because Singular values are ordered
        for(t_iRank = p_matSigma->rows()-1; t_iRank > 0; t_iRank--)
            if (p_matSigma->data()[t_iRank] > 0.00001)
                break;

        t_iRank++;//rank corresponding to epsilon

        return t_iRank;
    }

        
    //=========================================================================================================
    /**
    * Performs F * F^Transposed, is used when n > m 
    *
    * @param[in] p_matF The matrix which should be transformed.
    * @return F * F^Transposed (we call it FFT ;))
    */
    inline HPCMatrix<float> makeSquareMat(const HPCMatrix<float>& p_matF)
    {
        cuHPCMatrix<float> t_dev_matF(p_matF);

        cuHPCMatrix<float> t_dev_matResult(p_matF.rows(), p_matF.rows());

        t_dev_matResult.cuHPCMatMult('N','T',t_dev_matF,t_dev_matF);//FF^T


        //Make rectangular - p_matF*p_matF^T
        //MatrixXT FFT = p_matF*p_matF.transpose();

        return t_dev_matResult.toHPCMatrix();//p_matF*p_matF.transpose();
    }

};


} // NAMESPACE

#endif // RAPMUSIC_CUDA_CUH