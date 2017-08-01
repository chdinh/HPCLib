//=============================================================================================================
/**
* @file     rapmusic_gold.cpp
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
* @brief    Continuation of the RapMusic CPU template declaration.
*
*/

#ifndef RAPMUSIC_GOLD_SOURCES //Because this cpp is part of the header -> template
#define RAPMUSIC_GOLD_SOURCES


//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "../include/rapmusic_gold.h"

#include "../include/model.h"
#include "../include/hpcmatrix.h"
#include "../include/rapdipoles.h"


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
// USED NAMESPACES
//=============================================================================================================


//*************************************************************************************************************
//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

template <class T>
RapMusic<T>::RapMusic()
: m_pMappedMatLeadField(NULL)
, m_pMatGrid(NULL)
, m_iN(0)
, m_dThreshold(0)
, m_iNumGridPoints(0)
, m_iNumChannels(0)
, m_iNumLeadFieldCombinations(0)
, m_ppPairIdxCombinations(NULL)
, m_iMaxNumThreads(1)
, m_bIsInit(false)
{

}


//*************************************************************************************************************

template <class T>
RapMusic<T>::RapMusic(Model<T>* p_pModel, bool p_bSparsed, int p_iN, double p_dThr)
: m_pMappedMatLeadField(NULL)
, m_pMatGrid(NULL)
, m_iN(0)
, m_dThreshold(0)
, m_iNumGridPoints(0)
, m_iNumChannels(0)
, m_iNumLeadFieldCombinations(0)
, m_ppPairIdxCombinations(NULL)
, m_iMaxNumThreads(1)
, m_bIsInit(false)
{
    //Init
    initRAPMusic(p_pModel, p_bSparsed, p_iN, p_dThr);
}


//*************************************************************************************************************

template <class T>
RapMusic<T>::RapMusic(HPCMatrix<T>* p_pMatLeadField, HPCMatrix<T>* p_pMatGrid, int p_iN, double p_dThr)
: m_pMappedMatLeadField(NULL)
, m_pMatGrid(NULL)
, m_iN(0)
, m_dThreshold(0)
, m_iNumGridPoints(0)
, m_iNumChannels(0)
, m_iNumLeadFieldCombinations(0)
, m_ppPairIdxCombinations(NULL)
, m_iMaxNumThreads(1)
, m_bIsInit(false)
{
    //Init
    initRAPMusic(p_pMatLeadField, p_pMatGrid, p_iN, p_dThr);
}



//*************************************************************************************************************

template <class T>
RapMusic<T>::~RapMusic()
{
    if(m_pMappedMatLeadField != NULL)
        delete m_pMappedMatLeadField;//ToDo debug this -> make sure that Lead Field is still available

    if(m_ppPairIdxCombinations != NULL)
        free(m_ppPairIdxCombinations);

    //delete []m_ppPairIdxCombinations;
}


//*************************************************************************************************************

template <class T>
bool RapMusic<T>::initRAPMusic(Model<T>* p_pModel, bool p_bSparsed, int p_iN, double p_dThr)
{
    return initRAPMusic(p_bSparsed ? p_pModel->getSparsedLeadFieldMat() : p_pModel->getLeadFieldMat(),
                        p_bSparsed ? p_pModel->getSparsedGridMat() : p_pModel->getGridMat(),
                        p_iN, p_dThr);
}


//*************************************************************************************************************

template <class T>
bool RapMusic<T>::initRAPMusic(	HPCMatrix<T>* p_pMatLeadField,
                                HPCMatrix<T>* p_pMatGrid,
                                int p_iN, double p_dThr)
{
    //Get available thread number
    #ifdef _OPENMP
        std::cout << "OpenMP enabled" << std::endl;
        m_iMaxNumThreads = omp_get_max_threads();
    #else
        std::cout << "OpenMP disabled (to enable it: VS2010->Project Properties->C/C++->Language, then modify OpenMP Support)" << std::endl;
        m_iMaxNumThreads = 1;
    #endif
        std::cout << "Available Threats: " << m_iMaxNumThreads << std::endl << std::endl;

    //Initialize RAP MUSIC
    std::cout << "##### Initialization RAP MUSIC started ######\n\n";

    m_iN = p_iN;
    m_dThreshold = p_dThr;

    //Grid check
    if(p_pMatGrid != NULL)
    {
        if ( p_pMatGrid->rows() != p_pMatLeadField->cols() / 3 )
        {
            std::cout << "Grid does not fit to given Lead Field!\n";
            return false;
        }
    }

    m_pMatGrid = p_pMatGrid;


    //Lead Field check
    if ( p_pMatLeadField->cols() % 3 != 0 )
    {
        std::cout << "Lead Field is not associated with a 3D grid!\n";
        return false;
    }

    m_iNumGridPoints = p_pMatLeadField->cols()/3;

    m_iNumChannels = p_pMatLeadField->rows();

    m_pMappedMatLeadField = new Eigen::Map<MatrixXT>
        (   p_pMatLeadField->data(),
            p_pMatLeadField->rows(),
            p_pMatLeadField->cols() );

    //##### Calc lead field combination #####

    std::cout << "Calculate lead field combinations. \n";

    m_iNumLeadFieldCombinations = nchoose2(m_iNumGridPoints+1);

    m_ppPairIdxCombinations = (Pair **)malloc(m_iNumLeadFieldCombinations * sizeof(Pair *));

    calcPairCombinations(m_iNumGridPoints, m_iNumLeadFieldCombinations, m_ppPairIdxCombinations);

    std::cout << "Lead Field combinations calculated. \n\n";

    //##### Calc lead field combination end #####

    std::cout << "Number of grid points: " << m_iNumGridPoints << "\n\n";

    std::cout << "Number of combinated points: " << m_iNumLeadFieldCombinations << "\n\n";

    std::cout << "Number of sources to find: " << m_iN << "\n\n";

    std::cout << "Threshold: " << m_dThreshold << "\n\n";

    //Init end

    std::cout << "##### Initialization RAP MUSIC completed ######\n\n\n";

    m_bIsInit = true;

    return m_bIsInit;
}


//*************************************************************************************************************

template <class T>
bool RapMusic<T>::calcRAPMusic(HPCMatrix<T>* p_pMatMeasurement, RapDipoles<T>* &p_pRapDipoles)
{
    //if not initialized -> break
    if(!m_bIsInit)
    {
        std::cout << "RAP MUSIC wasn't initialized!"; //ToDo: catch this earlier
        return false;
    }

    //Test if data are correct
    if(p_pMatMeasurement->rows() != m_iNumChannels)
    {
        std::cout << "Lead Field channels do not fit to number of measurement channels!"; //ToDo: catch this earlier
        return false;
    }

    //Inits
    //Stop the time for benchmark purpose
    clock_t start, end;
    start = clock();

    //Map HPCMatrix to Eigen Matrix
    Eigen::Map<MatrixXT>
        t_MappedMatMeasurement(	p_pMatMeasurement->data(),
        p_pMatMeasurement->rows(),
        p_pMatMeasurement->cols() );

    //Calculate the signal subspace (t_pMatPhi_s)
    MatrixXT* t_pMatPhi_s = NULL;//(m_iNumChannels, m_iN < t_r ? m_iN : t_r);
    int t_r = calcPhi_s(/*(MatrixXT)*/t_MappedMatMeasurement, t_pMatPhi_s);

    int t_iMaxSearch = m_iN < t_r ? m_iN : t_r; //The smallest of Rank and Iterations

    if (t_r < m_iN)
    {
        std::cout << "Warning: Rank " << t_r << " of the measurement data is smaller than the " << m_iN;
        std::cout << " sources to find." << std::endl;
        std::cout << "         Searching now for " << t_iMaxSearch << " correlated sources.";
        std::cout << std::endl << std::endl;
    }

    //Create Orthogonal Projector
    //OrthProj
    MatrixXT t_matOrthProj(m_iNumChannels,m_iNumChannels);
    t_matOrthProj.setIdentity();

    //A_k_1
    MatrixXT t_matA_k_1(m_iNumChannels, t_iMaxSearch);
    t_matA_k_1.setZero();

    if (m_pMatGrid != NULL)
    {
        if(p_pRapDipoles != NULL)
            p_pRapDipoles->initRapDipoles(m_pMatGrid);
        else
            p_pRapDipoles = new RapDipoles<T>(m_pMatGrid);
    }
    else
    {
        if(p_pRapDipoles != NULL)
            delete p_pRapDipoles;

        p_pRapDipoles = new RapDipoles<T>();
    }

    std::cout << "##### Calculation of RAP MUSIC started ######\n\n";

    MatrixXT t_matProj_Phi_s(t_matOrthProj.rows(), t_pMatPhi_s->cols());
    //new Version: Calculate projection before
    MatrixXT t_matProj_LeadField(m_pMappedMatLeadField->rows(), m_pMappedMatLeadField->cols());

    for(int r = 0; r < t_iMaxSearch ; ++r)
    {
        t_matProj_Phi_s = t_matOrthProj*(*t_pMatPhi_s);

        //new Version: Calculating Projection before
        t_matProj_LeadField = t_matOrthProj * (*m_pMappedMatLeadField);//Subtract the found sources from the current found source

        //###First Option###
        //Step 1: lt. Mosher 1998 -> Maybe tmp_Proj_Phi_S is already orthogonal -> so no SVD needed -> U_B = tmp_Proj_Phi_S;
        Eigen::JacobiSVD< MatrixXT > t_svdProj_Phi_S(t_matProj_Phi_s, Eigen::ComputeThinU);
        MatrixXT t_matU_B;
        useFullRank(t_svdProj_Phi_S.matrixU(), t_svdProj_Phi_S.singularValues().asDiagonal(), t_matU_B);

        //Inits
        VectorXT t_vecRoh(m_iNumLeadFieldCombinations,1);
        t_vecRoh.setZero();

        //subcorr benchmark
        //Stop the time
        clock_t start_subcorr, end_subcorr;
        start_subcorr = clock();

        //Multithreading correlation calculation
        #ifdef _OPENMP
        #pragma omp parallel num_threads(m_iMaxNumThreads)
        #endif
        {
        #ifdef _OPENMP
        #pragma omp for
        #endif
            for(int i = 0; i < m_iNumLeadFieldCombinations; i++)
            {
                //new Version: calculate matrix multiplication before
                //Create Lead Field combinations -> It would be better to use a pointer construction, to increase performance
                MatrixX6T t_matProj_G(t_matProj_LeadField.rows(),6);

                int idx1 = m_ppPairIdxCombinations[i]->x1;
                int idx2 = m_ppPairIdxCombinations[i]->x2;

                getLeadFieldPair(t_matProj_LeadField, t_matProj_G, idx1, idx2);

                t_vecRoh(i) = (T)subcorr(t_matProj_G, t_matU_B);//t_vecRoh holds the correlations roh_k
            }
        }


//         if(r==0)
//         {
//             std::fstream filestr;
//             std::stringstream filename;
//             filename << "Roh_gold.txt";
// 
//             filestr.open ( filename.str().c_str(), std::fstream::out);
//             for(int i = 0; i < m_iNumLeadFieldCombinations; ++i)
//             {
//               filestr << t_vecRoh(i) << "\n";
//             }
//             filestr.close();
// 
//             //exit(0);
//         }


        //subcorr benchmark
        end_subcorr = clock();

        float t_fSubcorrElapsedTime = ( (float)(end_subcorr-start_subcorr) / (float)CLOCKS_PER_SEC ) * 1000.0f;
        std::cout << "Time Elapsed: " << t_fSubcorrElapsedTime << " ms" << std::endl;

        //Find the maximum of correlation - can't put this in the for loop because it's running in different threads.
        T t_val_roh_k;

        VectorXT::Index t_iMaxIdx;

        t_val_roh_k = t_vecRoh.maxCoeff(&t_iMaxIdx);//p_vecCor = ^roh_k

        //get positions in sparsed leadfield from index combinations;
        int t_iIdx1 = m_ppPairIdxCombinations[t_iMaxIdx]->x1;
        int t_iIdx2 = m_ppPairIdxCombinations[t_iMaxIdx]->x2;

        // (Idx+1) because of MATLAB positions -> starting with 1 not with 0
        std::cout << "Iteration: " << r+1 << " of " << t_iMaxSearch
            << "; Correlation: " << t_val_roh_k<< "; Position (Idx+1): " << t_iIdx1+1 << " - " << t_iIdx2+1 <<"\n\n";

        //Calculations with the max correlated dipole pair G_k_1 -> ToDo Obsolet when taking direkt Projected Lead Field
        MatrixX6T t_matG_k_1(m_pMappedMatLeadField->rows(),6);
        getLeadFieldPair(*m_pMappedMatLeadField, t_matG_k_1, t_iIdx1, t_iIdx2);

        MatrixX6T t_matProj_G_k_1(t_matOrthProj.rows(), t_matG_k_1.cols());
        t_matProj_G_k_1 = t_matOrthProj * t_matG_k_1;//Subtract the found sources from the current found source
//         MatrixX6T t_matProj_G_k_1(t_matProj_LeadField.rows(), 6);
//         getLeadFieldPair(t_matProj_LeadField, t_matProj_G_k_1, t_iIdx1, t_iIdx2);

        //Calculate source direction
        //source direction (p_pMatPhi) for current source r (phi_k_1)
        Vector6T t_vec_phi_k_1(6, 1);
        subcorr(t_matProj_G_k_1, t_matU_B, t_vec_phi_k_1);//Correlate the current source to calculate the direction

        //Set return values
        p_pRapDipoles->insertSource(t_iIdx1, t_iIdx2, t_vec_phi_k_1.data(), t_val_roh_k);

        //Stop Searching when Correlation is smaller then the Threshold
        if (t_val_roh_k < m_dThreshold)
        {
            std::cout << "Searching stopped, last correlation " << t_val_roh_k;
            std::cout << " is smaller then the given threshold " << m_dThreshold << std::endl << std::endl;
            break;
        }

        //Calculate A_k_1 = [a_theta_1..a_theta_k_1] matrix for subtraction of found source
        calcA_k_1(t_matG_k_1, t_vec_phi_k_1, r, t_matA_k_1);

        //Calculate new orthogonal Projector (Pi_k_1)
        calcOrthProj(t_matA_k_1, t_matOrthProj);

        //garbage collecting
        //ToDo
    }

    std::cout << "##### Calculation of RAP MUSIC completed ######"<< std::endl << std::endl << std::endl;

    end = clock();

    float t_fElapsedTime = ( (float)(end-start) / (float)CLOCKS_PER_SEC ) * 1000.0f;
    std::cout << "Total Time Elapsed: " << t_fElapsedTime << " ms" << std::endl << std::endl;

    //garbage collecting
    delete t_pMatPhi_s;

    return true;
}


//*************************************************************************************************************

template <class T>
bool RapMusic<T>::calcPowellRAPMusic(HPCMatrix<T>* p_pMatMeasurement, RapDipoles<T>* &p_pRapDipoles)
{
    //if not initialized -> break
    if(!m_bIsInit)
    {
        std::cout << "RAP MUSIC wasn't initialized!"; //ToDo: catch this earlier
        return false;
    }

    //Test if data are correct
    if(p_pMatMeasurement->rows() != m_iNumChannels)
    {
        std::cout << "Lead Field channels do not fit to number of measurement channels!"; //ToDo: catch this earlier
        return false;
    }

    //Inits
    //Stop the time for benchmark purpose
    clock_t start, end;
    start = clock();

    //Map HPCMatrix to Eigen Matrix
    Eigen::Map<MatrixXT>
        t_MappedMatMeasurement(	p_pMatMeasurement->data(),
        p_pMatMeasurement->rows(),
        p_pMatMeasurement->cols() );

    //Calculate the signal subspace (t_pMatPhi_s)
    MatrixXT* t_pMatPhi_s = NULL;//(m_iNumChannels, m_iN < t_r ? m_iN : t_r);
    int t_r = calcPhi_s(/*(MatrixXT)*/t_MappedMatMeasurement, t_pMatPhi_s);

    int t_iMaxSearch = m_iN < t_r ? m_iN : t_r; //The smallest of Rank and Iterations

    if (t_r < m_iN)
    {
        std::cout << "Warning: Rank " << t_r << " of the measurement data is smaller than the " << m_iN;
        std::cout << " sources to find." << std::endl;
        std::cout << "         Searching now for " << t_iMaxSearch << " correlated sources.";
        std::cout << std::endl << std::endl;
    }

    //Create Orthogonal Projector
    //OrthProj
    MatrixXT t_matOrthProj(m_iNumChannels,m_iNumChannels);
    t_matOrthProj.setIdentity();

    //A_k_1
    MatrixXT t_matA_k_1(m_iNumChannels, t_iMaxSearch);
    t_matA_k_1.setZero();

    if (m_pMatGrid != NULL)
    {
        if(p_pRapDipoles != NULL)
            p_pRapDipoles->initRapDipoles(m_pMatGrid);
        else
            p_pRapDipoles = new RapDipoles<T>(m_pMatGrid);
    }
    else
    {
        if(p_pRapDipoles != NULL)
            delete p_pRapDipoles;

        p_pRapDipoles = new RapDipoles<T>();
    }

    std::cout << "##### Calculation of RAP MUSIC started ######\n\n";

    MatrixXT t_matProj_Phi_s(t_matOrthProj.rows(), t_pMatPhi_s->cols());
    //new Version: Calculate projection before
    MatrixXT t_matProj_LeadField(m_pMappedMatLeadField->rows(), m_pMappedMatLeadField->cols());

    for(int r = 0; r < t_iMaxSearch ; ++r)
    {
        t_matProj_Phi_s = t_matOrthProj*(*t_pMatPhi_s);

        //new Version: Calculating Projection before
        t_matProj_LeadField = t_matOrthProj * (*m_pMappedMatLeadField);//Subtract the found sources from the current found source

        //###First Option###
        //Step 1: lt. Mosher 1998 -> Maybe tmp_Proj_Phi_S is already orthogonal -> so no SVD needed -> U_B = tmp_Proj_Phi_S;
        Eigen::JacobiSVD< MatrixXT > t_svdProj_Phi_S(t_matProj_Phi_s, Eigen::ComputeThinU);
        MatrixXT t_matU_B;
        useFullRank(t_svdProj_Phi_S.matrixU(), t_svdProj_Phi_S.singularValues().asDiagonal(), t_matU_B);

        //Inits
        VectorXT t_vecRoh(m_iNumLeadFieldCombinations,1);
        t_vecRoh.setZero();

        //subcorr benchmark
        //Stop the time
        clock_t start_subcorr, end_subcorr;
        start_subcorr = clock();

        T t_val_roh_k;

        //Powell
        int t_iCurrentRow = 2;

        int t_iIdx1 = -1;
        int t_iIdx2 = -1;

        int t_iMaxIdx_old = -1;

        int t_iMaxFound = 0;

        Eigen::VectorXi t_pVecIdxElements(m_iNumGridPoints);

        PowellIdxVec(t_iCurrentRow, m_iNumGridPoints, t_pVecIdxElements);

        int t_iNumVecElements = m_iNumGridPoints;

        while(t_iMaxFound == 0)
        {

            //Multithreading correlation calculation
            #ifdef _OPENMP
            #pragma omp parallel num_threads(m_iMaxNumThreads)
            #endif
            {
            #ifdef _OPENMP
            #pragma omp for
            #endif
                for(int i = 0; i < t_iNumVecElements; i++)
                {
                    int k = t_pVecIdxElements(i);
                    //new Version: calculate matrix multiplication before
                    //Create Lead Field combinations -> It would be better to use a pointer construction, to increase performance
                    MatrixX6T t_matProj_G(t_matProj_LeadField.rows(),6);

                    int idx1 = m_ppPairIdxCombinations[k]->x1;
                    int idx2 = m_ppPairIdxCombinations[k]->x2;

                    getLeadFieldPair(t_matProj_LeadField, t_matProj_G, idx1, idx2);

                    t_vecRoh(k) = (T)subcorr(t_matProj_G, t_matU_B);//t_vecRoh holds the correlations roh_k
                }
            }

    //         if(r==0)
    //         {
    //             std::fstream filestr;
    //             std::stringstream filename;
    //             filename << "Roh_gold.txt";
    // 
    //             filestr.open ( filename.str().c_str(), std::fstream::out);
    //             for(int i = 0; i < m_iNumLeadFieldCombinations; ++i)
    //             {
    //               filestr << t_vecRoh(i) << "\n";
    //             }
    //             filestr.close();
    // 
    //             //exit(0);
    //         }

            //Find the maximum of correlation - can't put this in the for loop because it's running in different threads.

            VectorXT::Index t_iMaxIdx;

            t_val_roh_k = t_vecRoh.maxCoeff(&t_iMaxIdx);//p_vecCor = ^roh_k

            if((int)t_iMaxIdx == t_iMaxIdx_old)
            {
                t_iMaxFound = 1;
                break;
            }
            else
            {
                t_iMaxIdx_old = t_iMaxIdx;
                //get positions in sparsed leadfield from index combinations;
                t_iIdx1 = m_ppPairIdxCombinations[t_iMaxIdx]->x1;
                t_iIdx2 = m_ppPairIdxCombinations[t_iMaxIdx]->x2;
            }


            //set new index
            if(t_iIdx1 == t_iCurrentRow)
                t_iCurrentRow = t_iIdx2;
            else
                t_iCurrentRow = t_iIdx1;

            PowellIdxVec(t_iCurrentRow, m_iNumGridPoints, t_pVecIdxElements);
        }

        //subcorr benchmark
        end_subcorr = clock();

        float t_fSubcorrElapsedTime = ( (float)(end_subcorr-start_subcorr) / (float)CLOCKS_PER_SEC ) * 1000.0f;
        std::cout << "Time Elapsed: " << t_fSubcorrElapsedTime << " ms" << std::endl;


        // (Idx+1) because of MATLAB positions -> starting with 1 not with 0
        std::cout << "Iteration: " << r+1 << " of " << t_iMaxSearch
            << "; Correlation: " << t_val_roh_k<< "; Position (Idx+1): " << t_iIdx1+1 << " - " << t_iIdx2+1 <<"\n\n";

        //Calculations with the max correlated dipole pair G_k_1
        MatrixX6T t_matG_k_1(m_pMappedMatLeadField->rows(),6);
        getLeadFieldPair(*m_pMappedMatLeadField, t_matG_k_1, t_iIdx1, t_iIdx2);

        MatrixX6T t_matProj_G_k_1(t_matOrthProj.rows(), t_matG_k_1.cols());
        t_matProj_G_k_1 = t_matOrthProj * t_matG_k_1;//Subtract the found sources from the current found source
//         MatrixX6T t_matProj_G_k_1(t_matProj_LeadField.rows(), 6);
//         getLeadFieldPair(t_matProj_LeadField, t_matProj_G_k_1, t_iIdx1, t_iIdx2);

        //Calculate source direction
        //source direction (p_pMatPhi) for current source r (phi_k_1)
        Vector6T t_vec_phi_k_1(6, 1);
        subcorr(t_matProj_G_k_1, t_matU_B, t_vec_phi_k_1);//Correlate the current source to calculate the direction

        //Set return values
        p_pRapDipoles->insertSource(t_iIdx1, t_iIdx2, t_vec_phi_k_1.data(), t_val_roh_k);

        //Stop Searching when Correlation is smaller then the Threshold
        if (t_val_roh_k < m_dThreshold)
        {
            std::cout << "Searching stopped, last correlation " << t_val_roh_k;
            std::cout << " is smaller then the given threshold " << m_dThreshold << std::endl << std::endl;
            break;
        }

        //Calculate A_k_1 = [a_theta_1..a_theta_k_1] matrix for subtraction of found source
        calcA_k_1(t_matG_k_1, t_vec_phi_k_1, r, t_matA_k_1);

        //Calculate new orthogonal Projector (Pi_k_1)
        calcOrthProj(t_matA_k_1, t_matOrthProj);

        //garbage collecting
        //ToDo
    }

    std::cout << "##### Calculation of RAP MUSIC completed ######"<< std::endl << std::endl << std::endl;

    end = clock();

    float t_fElapsedTime = ( (float)(end-start) / (float)CLOCKS_PER_SEC ) * 1000.0f;
    std::cout << "Total Time Elapsed: " << t_fElapsedTime << " ms" << std::endl << std::endl;

    //garbage collecting
    delete t_pMatPhi_s;

    return true;
}


//*************************************************************************************************************

template <class T>
int RapMusic<T>::PowellOffset(int p_iRow, int p_iNumPoints)
{

    return p_iRow*p_iNumPoints - (( (p_iRow-1)*p_iRow) / 2); //triangular series 1 3 6 10 ... = (num_pairs*(num_pairs+1))/2

}


//*************************************************************************************************************

template <class T>
void RapMusic<T>::PowellIdxVec(int p_iRow, int p_iNumPoints, Eigen::VectorXi& p_pVecElements)
{

    //     if(p_pVecElements != NULL)
    //         delete[] p_pVecElements;
    // 
    //     p_pVecElements = new int(p_iNumPoints);


    //col combination index
    for(int i = 0; i <= p_iRow; ++i)//=p_iNumPoints-1
        p_pVecElements(i) = PowellOffset(i+1,p_iNumPoints)-(p_iNumPoints-p_iRow);


    //row combination index
    int off = PowellOffset(p_iRow,p_iNumPoints);
    int length = p_iNumPoints - p_iRow;
    int k=0;
    for(int i = p_iRow; i < p_iRow+length; ++i)//=p_iNumPoints-1
    {
        p_pVecElements(i) = off+k;
        k = k + 1;
    }
}


//*************************************************************************************************************

template <class T>
int RapMusic<T>::nchoose2(int n)
{

    //nchoosek(n, k) with k = 2, equals n*(n-1)*0.5

    int t_iNumOfCombination = (int)(n*(n-1)*0.5);

    return t_iNumOfCombination;
}


//*************************************************************************************************************

template <class T>
int RapMusic<T>::calcPhi_s(const MatrixXT& p_pMatMeasurement, MatrixXT* &p_pMatPhi_s)
{
    //Calculate p_pMatPhi_s
    MatrixXT t_matF;//t_matF = makeSquareMat(p_pMatMeasurement); //FF^T -> ToDo Check this
    if (p_pMatMeasurement.cols() > p_pMatMeasurement.rows())
        t_matF = makeSquareMat(p_pMatMeasurement); //FF^T
    else
        t_matF = MatrixXT(p_pMatMeasurement);

    Eigen::JacobiSVD<MatrixXT> t_svdF(t_matF, Eigen::ComputeThinU);

    int t_r = getRank(t_svdF.singularValues().asDiagonal());

    int t_iCols = t_r;//t_r < m_iN ? m_iN : t_r;

    if (p_pMatPhi_s != NULL)
        delete p_pMatPhi_s;

    //m_iNumChannels has to be equal to t_svdF.matrixU().rows()
    p_pMatPhi_s = new MatrixXT(m_iNumChannels, t_iCols);

    //assign the signal subspace
    memcpy(p_pMatPhi_s->data(), t_svdF.matrixU().data(), sizeof(T) * m_iNumChannels * t_iCols);

    return t_r;
}


//*************************************************************************************************************

template <class T>
T RapMusic<T>::subcorr(	MatrixX6T& p_matProj_G, const MatrixXT& p_matU_B)
{
    //Orthogonalisierungstest wegen performance weggelassen -> ohne is es viel schneller

    Matrix6T t_matSigma_A(6, 6);
    Matrix6XT t_matU_A_T(6, p_matProj_G.rows()); //rows and cols are changed, because of CV_SVD_U_T
    
    Eigen::JacobiSVD<MatrixXT> t_svdProj_G(p_matProj_G, Eigen::ComputeThinU);

    t_matSigma_A = t_svdProj_G.singularValues().asDiagonal();
    t_matU_A_T = t_svdProj_G.matrixU().transpose();
    
    //lt. Mosher 1998 ToDo: Only Retain those Components of U_A and U_B that correspond to nonzero singular values
    //for U_A and U_B the number of columns corresponds to their ranks
    MatrixXT t_matU_A_T_full;
    //reduce to rank only when directions aren't calculated, otherwise use the full t_matU_A_T

    useFullRank(t_matU_A_T, t_matSigma_A, t_matU_A_T_full, IS_TRANSPOSED);//lt. Mosher 1998: Only Retain those Components of U_A that correspond to nonzero singular values -> for U_A the number of columns corresponds to their ranks

    MatrixXT t_matCor(t_matU_A_T_full.rows(), p_matU_B.cols());
    
    //Step 2: compute the subspace correlation
    t_matCor = t_matU_A_T_full*p_matU_B;//lt. Mosher 1998: C = U_A^T * U_B

    VectorXT t_vecSigma_C;
    
    if (t_matCor.cols() > t_matCor.rows())
    {
        MatrixXT t_matCor_H(t_matCor.cols(), t_matCor.rows());
        t_matCor_H = t_matCor.adjoint(); //for complex it has to be adjunct
        //ToDo -> use instead adjointInPlace

        Eigen::JacobiSVD<MatrixXT> t_svdCor_H(t_matCor_H);

        t_vecSigma_C = t_svdCor_H.singularValues();
    }
    else
    {
        Eigen::JacobiSVD<MatrixXT> t_svdCor(t_matCor);

        t_vecSigma_C = t_svdCor.singularValues();
    }

    //Step 3
    T t_dRetSigma_C;
    t_dRetSigma_C = t_vecSigma_C(0); //Take only the correlation of the first principal components

    //garbage collecting
    //ToDo

    return t_dRetSigma_C;
}


//*************************************************************************************************************
template <class T>
T RapMusic<T>::subcorr(MatrixX6T& p_matProj_G, const MatrixXT& p_matU_B, Vector6T& p_vec_phi_k_1)
{
    //Orthogonalisierungstest wegen performance weggelassen -> ohne is es viel schneller

    Matrix6T sigma_A(6, 6);
    Matrix6XT U_A_T(6, p_matProj_G.rows()); //rows and cols are changed, because of CV_SVD_U_T
    Matrix6T V_A(6, 6);

    Eigen::JacobiSVD<MatrixXT> svdOfProj_G(p_matProj_G, Eigen::ComputeThinU | Eigen::ComputeThinV);

    sigma_A = svdOfProj_G.singularValues().asDiagonal();
    U_A_T = svdOfProj_G.matrixU().transpose();
    V_A = svdOfProj_G.matrixV();

    //lt. Mosher 1998 ToDo: Only Retain those Components of U_A and U_B that correspond to nonzero singular values
    //for U_A and U_B the number of columns corresponds to their ranks
    //-> reduce to rank only when directions aren't calculated, otherwise use the full U_A_T

    Matrix6XT t_matCor(6, p_matU_B.cols());

    //Step 2: compute the subspace correlation
    t_matCor = U_A_T*p_matU_B;//lt. Mosher 1998: C = U_A^T * U_B


    VectorXT sigma_C;

    //Step 4
    Matrix6XT U_C;

    if (t_matCor.cols() > t_matCor.rows())
    {
        MatrixX6T Cor_H(t_matCor.cols(), 6);
        Cor_H = t_matCor.adjoint(); //for complex it has to be adjunct

        Eigen::JacobiSVD<MatrixXT> svdOfCor_H(Cor_H, Eigen::ComputeThinV);

        U_C = svdOfCor_H.matrixV(); //because t_matCor Hermitesch U and V are exchanged
        sigma_C = svdOfCor_H.singularValues();
    }
    else
    {
        Eigen::JacobiSVD<MatrixXT> svdOfCor(t_matCor, Eigen::ComputeThinU);

        U_C = svdOfCor.matrixU();
        sigma_C = svdOfCor.singularValues();
    }

    Matrix6T sigma_a_inv;
    sigma_a_inv = sigma_A.inverse();

    Matrix6XT X;
    X = (V_A*sigma_a_inv)*U_C;//X = V_A*Sigma_A^-1*U_C

    Vector6T X_max;//only for the maximum c - so instead of X->cols use 1
    X_max = X.col(0);

    T norm_X = 1/(X_max.norm());

    //Multiply a scalar with an Array -> linear transform
    p_vec_phi_k_1 = X_max*norm_X;//u1 = x1/||x1|| this is the orientation

    //garbage collecting
    //ToDo

    //Step 3
    T ret_sigma_C;
    ret_sigma_C = sigma_C(0); //Take only the correlation of the first principal components

    //garbage collecting
    //ToDo

    return ret_sigma_C;
}


//*************************************************************************************************************

template <class T>
void RapMusic<T>::calcA_k_1(	const MatrixX6T& p_matG_k_1,
                                const Vector6T& p_matPhi_k_1,
                                const int p_iIdxk_1,
                                MatrixXT& p_matA_k_1)
{
    //Calculate A_k_1 = [a_theta_1..a_theta_k_1] matrix for subtraction of found source
    VectorXT t_vec_a_theta_k_1(p_matG_k_1.rows(),1);

    t_vec_a_theta_k_1 = p_matG_k_1*p_matPhi_k_1; // a_theta_k_1 = G_k_1*phi_k_1   this corresponds to the normalized signal component in subspace r

    p_matA_k_1.block(0,p_iIdxk_1,p_matA_k_1.rows(),1) = t_vec_a_theta_k_1;
}


//*************************************************************************************************************

template <class T>
void RapMusic<T>::calcOrthProj(const MatrixXT& p_matA_k_1, MatrixXT& p_matOrthProj)
{
    //Calculate OrthProj=I-A_k_1*(A_k_1'*A_k_1)^-1*A_k_1' //Wetterling -> A_k_1 = Gain

    MatrixXT t_matA_k_1_tmp(p_matA_k_1.cols(), p_matA_k_1.cols());
    t_matA_k_1_tmp = p_matA_k_1.adjoint()*p_matA_k_1;//A_k_1'*A_k_1 = A_k_1_tmp -> A_k_1' has to be adjoint for complex


    int t_size = t_matA_k_1_tmp.cols();

    while (!t_matA_k_1_tmp.block(0,0,t_size,t_size).fullPivLu().isInvertible())
    {
        --t_size;
    }

    MatrixXT t_matA_k_1_tmp_inv(t_matA_k_1_tmp.rows(), t_matA_k_1_tmp.cols());
    t_matA_k_1_tmp_inv.setZero();

    t_matA_k_1_tmp_inv.block(0,0,t_size,t_size) = t_matA_k_1_tmp.block(0,0,t_size,t_size).inverse();//(A_k_1_tmp)^-1 = A_k_1_tmp_inv


    t_matA_k_1_tmp = MatrixXT::Zero(p_matA_k_1.rows(), p_matA_k_1.cols());

    t_matA_k_1_tmp = p_matA_k_1*t_matA_k_1_tmp_inv;//(A_k_1*A_k_1_tmp_inv) = A_k_1_tmp


    MatrixXT t_matA_k_1_tmp2(p_matA_k_1.rows(), p_matA_k_1.rows());
    
    t_matA_k_1_tmp2 = t_matA_k_1_tmp*p_matA_k_1.adjoint();//(A_k_1_tmp)*A_k_1' -> here A_k_1' is only transposed - it has to be adjoint


    MatrixXT I(m_iNumChannels,m_iNumChannels);
    I.setIdentity();

    p_matOrthProj = I-t_matA_k_1_tmp2; //OrthProj=I-A_k_1*(A_k_1'*A_k_1)^-1*A_k_1';

    //garbage collecting
    //ToDo
}


//*************************************************************************************************************

template <class T>
void RapMusic<T>::calcPairCombinations( const int p_iNumPoints,
                                            const int p_iNumCombinations,
                                            Pair** p_ppPairIdxCombinations)
{
    int idx1 = 0;
    int idx2 = 0;

    //Process Code in {m_max_num_threads} threads -> When compile with Intel Compiler -> probably obsolete
    #ifdef _OPENMP
    #pragma omp parallel num_threads(m_iMaxNumThreads) private(idx1, idx2)
    #endif
    {
    #ifdef _OPENMP
    #pragma omp for
    #endif
        for (int i = 0; i < p_iNumCombinations; ++i)
        {
            getPointPair(p_iNumPoints, i, idx1, idx2);

            Pair* t_pairCombination = new Pair();
            t_pairCombination->x1 = idx1;
            t_pairCombination->x2 = idx2;

            p_ppPairIdxCombinations[i] = t_pairCombination;
        }
    }
}


//*************************************************************************************************************

template <class T>
void RapMusic<T>::getPointPair(const int p_iPoints, const int p_iCurIdx, int &p_iIdx1, int &p_iIdx2)
{
    int ii = p_iPoints*(p_iPoints+1)/2-1-p_iCurIdx;
    int K = (int)floor((sqrt((double)(8*ii+1))-1)/2);
    
    p_iIdx1 = p_iPoints-1-K;

    p_iIdx2 = (p_iCurIdx-p_iPoints*(p_iPoints+1)/2 + (K+1)*(K+2)/2)+p_iIdx1;
}


//*************************************************************************************************************
//ToDo don't make a real copy
template <class T>
void RapMusic<T>::getLeadFieldPair(	const MatrixXT& p_matLeadField,
                                    MatrixX6T& p_matLeadField_Pair,
                                    int p_iIdx1, int p_iIdx2)
{
    p_matLeadField_Pair.block(0,0,p_matLeadField.rows(),3) = p_matLeadField.block(0, p_iIdx1*3, p_matLeadField.rows(), 3);

    p_matLeadField_Pair.block(0,3,p_matLeadField.rows(),3) = p_matLeadField.block(0, p_iIdx2*3, p_matLeadField.rows(), 3);
}


}//Namespace

#endif //RAPMUSIC_GOLD_SOURCES