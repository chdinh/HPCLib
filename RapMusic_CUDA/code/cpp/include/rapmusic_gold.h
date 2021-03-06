//=============================================================================================================
/**
* @file     rapmusic_gold.h
* @author   Christoph Dinh <christoph.dinh@live.de>
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
* @brief    Contains the CPU implementation (gold standard) of the RAP MUSIC algorithm.
*
*/

#ifndef RAPMUSIC_GOLD_H
#define RAPMUSIC_GOLD_H


//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#ifdef _OPENMP
#include <omp.h>
#endif


//*************************************************************************************************************
//=============================================================================================================
// STL INCLUDES
//=============================================================================================================

#include <iostream>
#include <time.h>


//*************************************************************************************************************
//=============================================================================================================
// EIGEN 3.0 INCLUDES
//=============================================================================================================

#include "../Eigen/Core"
#include "../Eigen/SVD"
#include "../Eigen/LU"


//*************************************************************************************************************
//=============================================================================================================
// DEFINE NAMESPACE HPCLib
//=============================================================================================================

namespace HPCLib
{


//*************************************************************************************************************
//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================


//*************************************************************************************************************
//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================

template<class T>
class Model;

template<class T>
class HPCMatrix;

template<class T>
class RapDipoles;


//*************************************************************************************************************
//=============================================================================================================
// SOME DEFINES
//=============================================================================================================

#define NOT_TRANSPOSED   0  /**< Defines NOT_TRANSPOSED */
#define IS_TRANSPOSED   1   /**< Defines IS_TRANSPOSED */


//=============================================================================================================
/**
* Declares a pair structure for index combinations used in RAP MUSIC algorithm.
*/
typedef struct Pair
{
    int x1; /**< Index one of the pair. */
    int x2; /**< Index two of the pair. */
} Pair;


//=============================================================================================================
/**
* @brief    The RapMusic class provides the RAP MUSIC Algorithm CPU implementation.
*
* ToDo Detailed description
*/
template<class T>
class RapMusic
{

//*************************************************************************************************************
//=============================================================================================================
// TYPEDEFS
//=============================================================================================================

typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> MatrixXT;  /**< Defines Eigen::Matrix<T, Eigen::Dynamic,
                                                                         Eigen::Dynamic> as MatrixXT type. */
typedef Eigen::Matrix<T, Eigen::Dynamic, 6> MatrixX6T;              /**< Defines Eigen::Matrix<T, Eigen::Dynamic,
                                                                         6> as MatrixX6T type. */
typedef Eigen::Matrix<T, 6, Eigen::Dynamic> Matrix6XT;              /**< Defines Eigen::Matrix<T, 6,
                                                                         Eigen::Dynamic> as Matrix6XT type. */
typedef Eigen::Matrix<T, 6, 6> Matrix6T;                            /**< Defines Eigen::Matrix<T, 6, 6>
                                                                         as Matrix6T type. */
typedef Eigen::Matrix<T, Eigen::Dynamic, 1> VectorXT;               /**< Defines Eigen::Matrix<T, Eigen::Dynamic,
                                                                         1> as VectorXT type. */
typedef Eigen::Matrix<T, 6, 1> Vector6T;                            /**< Defines Eigen::Matrix<T, 6, 1>
                                                                         as Vector6T type. */

public:

    //=========================================================================================================
    /**
    * Default constructor creates an empty RapMusic algorithm which still needs to be initialized.
    */
    RapMusic();


    //=========================================================================================================
    /**
    * Constructor which initializes the RapMusic algorithm with the given model.
    *
    * @param[in] p_pModel  The model which contains the Lead Field matrix and its corresponding Grid matrix.
    * @param[in] p_bSparsed    True when sparse matrices should be used.
    * @param[in] p_iN   The number (default 2) of uncorrelated sources, which should be found. Starting with
    *                   the strongest.
    * @param[in] p_dThr    The correlation threshold (default 0.5) at which the search for sources stops.
    */
    RapMusic(Model<T>* p_pModel, bool p_bSparsed, int p_iN = 2, double p_dThr = 0.5);


    //=========================================================================================================
    /**
    * Constructor which initializes the RapMusic algorithm with the given Lead Field and optional Grid matrix.
    *
    * @param[in] p_pMatLeadField   The Lead Field Matrix (m x n*3) with rows = channels m;
    *                               and cols = grid n*3 (x1 y1 z1; x2 y2 z2..).
    * @param[in] p_pMatGrid    The Grid Matrix (n x 3)(default NULL) with rows = grid points n;
    *                           and cols = coordinates (x y z).
    * @param[in] p_iN   The number (default 2) of uncorrelated sources, which should be found. Starting with
    *                   the strongest.
    * @param[in] p_dThr    The correlation threshold (default 0.5) at which the search for sources stops.
    */
    RapMusic(   HPCMatrix<T>* p_pMatLeadField,
                HPCMatrix<T>* p_pMatGrid = NULL,
                int p_iN = 2, double p_dThr = 0.5);


    //=========================================================================================================
    /**
    * Th destructor destroys the RapMusic object properly by performing a cleanup.
    */
    ~RapMusic();


    //=========================================================================================================
    /**
    * Initializes the RAP MUSIC algorithm with the given model.
    *
    * @param[in] p_pModel  The model which contains the Lead Field matrix and its corresponding Grid matrix.
    * @param[in] p_bSparsed    True when sparse matrices should be used.
    * @param[in] p_iN   The number (default 2) of uncorrelated sources, which should be found. Starting with
    *                   the strongest.
    * @param[in] p_dThr    The correlation threshold (default 0.5) at which the search for sources stops.
    * @return   true if successful initialized, false otherwise.
    */
    bool initRAPMusic(Model<T>* p_pModel, bool p_bSparsed = false, int p_iN = 2, double p_dThr = 0.5);


    //=========================================================================================================
    /**
    * Initializes the RAP MUSIC algorithm with the given Lead Field and optional Grid matrix.
    *
    * @param[in] p_pMatLeadField    The Lead Field Matrix (m x n*3) with rows = channels m;
    *                               and cols = grid n*3 (x1 y1 z1; x2 y2 z2..).
    * @param[in] p_pMatGrid     The Grid Matrix (n x 3)(default NULL) with rows = grid points n;
    *                           and cols = coordinates (x y z).
    * @param[in] p_iN   The number (default 2) of uncorrelated sources, which should be found. Starting with
    *                   the strongest.
    * @param[in] p_dThr     The correlation threshold (default 0.5) at which the search for sources stops.
    * @return   True if successful initialized, false otherwise.
    */
    bool initRAPMusic(  HPCMatrix<T>* p_pMatLeadField,
                        HPCMatrix<T>* p_pMatGrid = NULL,
                        int p_iN = 2, double p_dThr = 0.5);


    //=========================================================================================================
    /**
    * Runs the RAP Music algorithm with given Lead Field and the current Measurement matrix.
    *
    * @param[in] p_pMatMeasurement The Mesaurement matrix (m x s) with rows = channels m;
    *                               and cols = sampling time points s
    * @param[out] p_pRapDipoles     Reference to pointer which contains the located and correlated Dipoles
    *                               after the calculation has finished. When a Grid was specified it contains
    *                               the coordinates to the locations.
    * @return   True if calculation was successful, false otherwise.
    */
    bool calcRAPMusic(HPCMatrix<T>* p_pMatMeasurement, RapDipoles<T>* &p_pRapDipoles);

    //=========================================================================================================
    /**
    * Runs the Powell search accelearetd RAP Music algorithm with given Lead Field and the current Measurement
    * matrix.
    *
    * @param[in] p_pMatMeasurement The Mesaurement matrix (m x s) with rows = channels m;
    *                               and cols = sampling time points s
    * @param[out] p_pRapDipoles     Reference to pointer which contains the located and correlated Dipoles
    *                               after the calculation has finished. When a Grid was specified it contains
    *                               the coordinates to the locations.
    * @return   True if calculation was successful, false otherwise.
    */
    bool calcPowellRAPMusic(HPCMatrix<T>* p_pMatMeasurement, RapDipoles<T>* &p_pRapDipoles);


    int PowellOffset(int p_iRow, int p_iNumPoints);


    void PowellIdxVec(int p_iRow, int p_iNumPoints, Eigen::VectorXi& p_pVecElements);


    //=========================================================================================================
    /**
    * Calculates the combination of n over 2 (nchoosek(n,2))
    *
    * @param[in] n  The number of elements which should be combined with each other (n over 2)
    * @return   The number of combinations
    */
    int nchoose2(int n);


protected:
    
    //=========================================================================================================
    /**
    * Computes the signal subspace Phi_s out of the measurement F.
    *
    * @param[in] p_pMatMeasurement  The current measured data to process (for best performance it should have
                                    the dimension channels x samples with samples = number of channels)
    * @param[out] p_pMatPhi_s   The calculated signal subspace.
    * @return   The rank of the measurement F (named r lt. Mosher 1998, 1999)
    */
    int calcPhi_s(const MatrixXT& p_pMatMeasurement, MatrixXT* &p_pMatPhi_s);


    //=========================================================================================================
    /**
    * Computes the subspace correlation between the projected G_rho and the projected signal subspace Phi_s.
    * For speed-up: we calculate the decomposition of the projected Phi_s before this function. So the argument
    * for this function is U_B instead of Phi_s.
    *
    * @param[in] p_matProj_G    The projected Lead Field combination. This is a m x 6 matrix composed of the
    *                           Lead Field combination of two Points for all m channels and 3 orthogonal
    *                           components (x y z).
    * @param[in] p_matU_B    The matrix U is the subspace projection of the orthogonal projected Phi_s
    * @return   The maximal correlation c_1 of the subspace correlation of the current projected Lead Field
    *           combination and the projected measurement.
    */
    T subcorr(MatrixX6T& p_matProj_G, const MatrixXT& p_matU_B);


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
    T subcorr(MatrixX6T& p_matProj_G, const MatrixXT& p_matU_B, Vector6T& p_vec_phi_k_1);


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
    void calcA_k_1( const MatrixX6T& p_matG_k_1,
                    const Vector6T& p_matPhi_k_1,
                    const int p_iIdxk_1,
                    MatrixXT& p_matA_k_1);


    //=========================================================================================================
    /**
    * Calculates the orthogonal projector Phi_A_k_1 like in the paper Mosher 1999 (13)
    *
    * @param[in] p_matA_k_1 The array of the manifold vectors.
    * @param[out] p_matOrthProj The orthogonal projector.
    */
    void calcOrthProj(const MatrixXT& p_matA_k_1, MatrixXT& p_matOrthProj);


    //=========================================================================================================
    /**
    * Pre-Calculates the Lead Field index combinations to search for a two dipole independent topography
    * (IT = source).
    *
    * @param[in] p_iNumPoints   The number of Lead Field points -> for dimension check
    * @param[in] p_iNumCombinations The number of pair index combinations.
    * @param[out] p_ppPairIdxCombinations   The destination which contains pointer to pointer of index
    *                                       combinations of Lead Field indices -> Number of pointers =
    *                                       Combination (number of grid points over 2 = Num + 1 C 2)
    */
    void calcPairCombinations(  const int p_iNumPoints,
                                const int p_iNumCombinations,
                                Pair** p_ppPairIdxCombinations);


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
    void getPointPair(const int p_iPoints, const int p_iCurIdx, int &p_iIdx1, int &p_iIdx2);


    //=========================================================================================================
    /**
    * Returns a Lead Field pair for the given indices 
    *
    * @param[in] p_matLeadField The Lead Field matrix.
    * @param[out] p_matLeadField_Pair   Lead Field combination (dimension: m x 6)
    * @param[in] p_iIdx1 first Lead Field index point
    * @param[in] p_iIdx2 second Lead Field index point
    */
    void getLeadFieldPair(const MatrixXT& p_matLeadField,
                                MatrixX6T& p_matLeadField_Pair,
                                int p_iIdx1, int p_iIdx2);


private:

    const Eigen::Map<MatrixXT>* m_pMappedMatLeadField;  /**< Holds the given Lead Field which is mapped
                                                             to an Eigen matrix structure.*/

    HPCMatrix<T>* m_pMatGrid;   /**< Holds the Grid which corresponds to the Lead Field. This variable can be
                                     NULL, then no coordinates are assigned to the located dipoles.*/

    int m_iN;               /**< Number of Sources to find*/
    double m_dThreshold;    /**< Threshold which defines the minimal correlation. Is the correlation of
                                 the found dipole pair smaller as this threshold than the RAP MUSIC
                                 calculation is stopped. */

    int m_iNumGridPoints;               /**< Number of Grid points. */
    int m_iNumChannels;                 /**< Number of channels */
    int m_iNumLeadFieldCombinations;    /**< Number of Lead Filed combinations (grid points + 1 over 2)*/

    Pair** m_ppPairIdxCombinations; /**< Index combination vector with grid pair indices. */

    int m_iMaxNumThreads;   /**< Number of available CPU threads. */

    bool m_bIsInit; /**< Wether the algorithm is initialized. */


    //=========================================================================================================
    /**
    * Returns the rank r of a singular value matrix based on non-zero singular values
    * (singular value > epsilon = 10^-5) 
    *
    * @param[in] p_matSigma diagonal matrix which contains the Singular values (Dimension n x n)
    * @return The rank r.
    */
    inline int getRank(const MatrixXT& p_matSigma)
    {
        int t_iRank;
        //if once a singularvalue is smaller than epsilon = 10^-5 the following values are also smaller
        // -> because Singular values are ordered
        for(t_iRank = p_matSigma.rows()-1; t_iRank > 0; t_iRank--)
            if (p_matSigma(t_iRank, t_iRank) > 0.00001)
                break;

        t_iRank++;//rank corresponding to epsilon

        return t_iRank;
    }


    //=========================================================================================================
    /**
    * lt. Mosher 1998 -> Only Retain those Components of U_A and U_B that correspond to nonzero singular values
    * for U_A and U_B the number of columns corresponds to their ranks
    *
    * @param[in] p_Mat  The Matrix which should be reduced to its rank.
    * @param[in] p_matSigma_src The singular values of the matrix
    * @param[out] p_matFull_Rank    The corresponding full rank matrix.
    * @param[in] type   Whether p_Mat is transposed, than rows and columns are changed.
    */
    inline int useFullRank( const MatrixXT& p_Mat,
                            const MatrixXT& p_matSigma_src, 
                            MatrixXT& p_matFull_Rank,
                            int type = NOT_TRANSPOSED)
    {
        int rank = getRank(p_matSigma_src);

        if (type == NOT_TRANSPOSED)
            p_matFull_Rank = p_Mat.block(0,0,p_Mat.rows(),rank);
        else
            p_matFull_Rank = p_Mat.block(0,0,rank,p_Mat.cols());

        return rank;
    }
    
    //=========================================================================================================
    /**
    * Performs F * F^Transposed, is used when n > m 
    *
    * @param[in] p_matF The matrix which should be transformed.
    * @return F * F^Transposed (we call it FFT ;))
    */
    inline MatrixXT makeSquareMat(const MatrixXT& p_matF)
    {
        //Make rectangular - p_matF*p_matF^T
        //MatrixXT FFT = p_matF*p_matF.transpose();

        return p_matF*p_matF.transpose();
    }
};

} //NAMESPACE

//Make the template definition visible to compiler in the first point of instantiation
#include "../src/rapmusic_gold.cpp"

#endif // RAPMUSIC_GOLD_H