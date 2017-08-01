//=============================================================================================================
/**
* @file     rapdipoles.cpp
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
* @brief    Continuation of the RapDipoles template declaration.
*
*/

#ifndef RAPDIPOLES_SOURCES //Because this cpp is part of the header -> template
#define RAPDIPOLES_SOURCES


//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "../include/rapdipoles.h"

#include "../include/model.h"

#include "../include/hpcmatrix.h"


//*************************************************************************************************************
//=============================================================================================================
// STL INCLUDES
//=============================================================================================================


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
RapDipoles<T>::RapDipoles()
: m_pMatGrid(NULL)
, m_bGridInitialized(false)
{

}


//*************************************************************************************************************

template <class T>
RapDipoles<T>::RapDipoles(HPCMatrix<T>* p_pMatGrid)
: m_pMatGrid(NULL)
, m_bGridInitialized(false)
{
    initRapDipoles(p_pMatGrid);
}


//*************************************************************************************************************

template <class T>
RapDipoles<T>::~RapDipoles()
{
    clean();
}


//*************************************************************************************************************

template <class T>
bool RapDipoles<T>::initRapDipoles(HPCMatrix<T>* p_pMatGrid)
{
    clean();

    m_pMatGrid = p_pMatGrid;

    m_bGridInitialized = true;

    return m_bGridInitialized;
}


//*************************************************************************************************************

template <class T>
void RapDipoles<T>::insertSource(	int p_iDipoleIdx1, int p_iDipoleIdx2,
                                const T *p_vec_phi_k_1, 
                                T p_valCor)
{
    RapDipolePair<T>* t_pRapDipolePair = new RapDipolePair<T>();

    t_pRapDipolePair->m_iIdx1 = p_iDipoleIdx1; //p_iDipoleIdx1+1 because of MATLAB index
    t_pRapDipolePair->m_iIdx2 = p_iDipoleIdx2;


    t_pRapDipolePair->m_Dipole1.x() = m_bGridInitialized ? (*m_pMatGrid)(p_iDipoleIdx1, 0) : 0;
    t_pRapDipolePair->m_Dipole1.y() = m_bGridInitialized ? (*m_pMatGrid)(p_iDipoleIdx1, 1) : 0;
    t_pRapDipolePair->m_Dipole1.z() = m_bGridInitialized ? (*m_pMatGrid)(p_iDipoleIdx1, 2) : 0;

    t_pRapDipolePair->m_Dipole2.x() = m_bGridInitialized ? (*m_pMatGrid)(p_iDipoleIdx2, 0) : 0;
    t_pRapDipolePair->m_Dipole2.y() = m_bGridInitialized ? (*m_pMatGrid)(p_iDipoleIdx2, 1) : 0;
    t_pRapDipolePair->m_Dipole2.z() = m_bGridInitialized ? (*m_pMatGrid)(p_iDipoleIdx2, 2) : 0;


    t_pRapDipolePair->m_Dipole1.phi_x() = p_vec_phi_k_1[0];
    t_pRapDipolePair->m_Dipole1.phi_y() = p_vec_phi_k_1[1];
    t_pRapDipolePair->m_Dipole1.phi_z() = p_vec_phi_k_1[2];

    t_pRapDipolePair->m_Dipole2.phi_x() = p_vec_phi_k_1[3];
    t_pRapDipolePair->m_Dipole2.phi_y() = p_vec_phi_k_1[4];
    t_pRapDipolePair->m_Dipole2.phi_z() = p_vec_phi_k_1[5];


    t_pRapDipolePair->m_vCorrelation = p_valCor;

    m_vecRapDipolePairs.push_back(t_pRapDipolePair);
}



//*************************************************************************************************************

template <class T>
void RapDipoles<T>::print()
{
    //Output
    std::cout << "Results of RAPMusic: \n\n";

    printf( "Source Pairs (Sparsed Grid Number {Coordinates}):\n" );
    for(unsigned int i = 0; i < m_vecRapDipolePairs.size(); i++ )
    {
        printf( "\n Source %i:\n\n", i+1 );

        printf( "     - Correlation: %f \n",m_vecRapDipolePairs[i]->m_vCorrelation);

        printf( "\n   Dipole Point:\n\n");
        if(m_bGridInitialized)
        {
            printf("     - Position:    %i {x = %4.2f, y = %4.2f, z = %4.2f}\n\n", m_vecRapDipolePairs[i]->m_iIdx1+1, m_vecRapDipolePairs[i]->m_Dipole1.x(), m_vecRapDipolePairs[i]->m_Dipole1.y(), m_vecRapDipolePairs[i]->m_Dipole1.z());
        }
        else
        {
            printf( "     - Position:    %i {Coordinates can't be displayed - Grid is not initialized}\n\n", m_vecRapDipolePairs[i]->m_iIdx1+1);
        }
        printf( "     - Orientation: %f %f %f \n", m_vecRapDipolePairs[i]->m_Dipole1.phi_x(), m_vecRapDipolePairs[i]->m_Dipole1.phi_y(), m_vecRapDipolePairs[i]->m_Dipole1.phi_z());


        printf( "\n   Correlated Dipole Point:\n\n");
        if(m_bGridInitialized)
        {
            printf("     - Position:    %i {x = %4.2f, y = %4.2f, z = %4.2f}\n\n", m_vecRapDipolePairs[i]->m_iIdx2+1, m_vecRapDipolePairs[i]->m_Dipole2.x(), m_vecRapDipolePairs[i]->m_Dipole2.y(), m_vecRapDipolePairs[i]->m_Dipole2.z());
        }
        else
        {
            printf("     - Position:    %i {Coordinates can't be displayed - Grid is not initialized}\n\n", m_vecRapDipolePairs[i]->m_iIdx2+1);
        }
        printf("     - Orientation: %f %f %f \n", m_vecRapDipolePairs[i]->m_Dipole2.phi_x(), m_vecRapDipolePairs[i]->m_Dipole2.phi_y(), m_vecRapDipolePairs[i]->m_Dipole2.phi_z());

        printf( "\n" );
    }
    printf( "\n" );


//     printf( "Dipole Orientations (phi):\n" );
//     for(unsigned int i = 0; i < m_vecRapDipolePairs.size(); i++ )
//     {
//         printf( "\n Source %i:\n", i+1 );
// 
//         printf( "\n   Dipole Point:\n     %f %f %f ", m_vecRapDipolePairs[i]->m_Dipole1.phi_x(), m_vecRapDipolePairs[i]->m_Dipole1.phi_y(), m_vecRapDipolePairs[i]->m_Dipole1.phi_z());
// 
//         printf("\n\n   Correlated Dipole Point:\n     %f %f %f ", m_vecRapDipolePairs[i]->m_Dipole2.phi_x(), m_vecRapDipolePairs[i]->m_Dipole2.phi_y(), m_vecRapDipolePairs[i]->m_Dipole2.phi_z());
// 
//         printf( "\n" );
//     }
//     printf( "\n\n" );


//     printf( "Source Correlations:\n" );
//     for(unsigned int i = 0; i < m_vecRapDipolePairs.size(); i++ )
//     {
//         printf( "\n Source %i:\n", i+1 );
// 
//         printf( "     %f ",m_vecRapDipolePairs[i]->m_vCorrelation);
// 
//         printf( "\n" );
//     }
//     printf( "\n" );
}


//*************************************************************************************************************

template <class T>
void RapDipoles<T>::clean()
{
    m_vecRapDipolePairs.clear();

    m_pMatGrid = NULL;

    m_bGridInitialized = false;
}


}//Namespace


#endif //RAPDIPOLES_SOURCES