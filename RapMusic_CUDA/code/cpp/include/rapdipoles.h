//=============================================================================================================
/**
* @file     rapdipoles.h
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
* @brief    Contains the declaration of the RapDipoles class.
*
*/


#ifndef RAPDIPOLES_H
#define RAPDIPOLES_H


//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "dipole.h"


//*************************************************************************************************************
//=============================================================================================================
// STL INCLUDES
//=============================================================================================================

#include <iostream>
#include <vector>


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
class HPCMatrix;


//*************************************************************************************************************
//=============================================================================================================
// SOME DEFINES
//=============================================================================================================

//=============================================================================================================
/**
* Declares a RapDipolePair structure consisting of two correlated dipoles which are the result of the RAP MUSIC
* searching algorithm.
*/
template<typename T>
struct RapDipolePair
{
    int m_iIdx1;            /**< Index of dipole one. */
    Dipole<T> m_Dipole1;    /**< Dipole one. */

    int m_iIdx2;            /**< Index of dipole two. */
    Dipole<T> m_Dipole2;    /**< Dipole two. */

    T m_vCorrelation;     /**< Correlation of the dipole pair. */
};


//=============================================================================================================
/**
* @brief    The RapDipoles class provides a list of correlated dipoles which are the result of the RAP MUSIC
*           algorithm.
*
* ToDo Detailed description
*/
template<class T>
class RapDipoles
{
public:

    //=========================================================================================================
    /**
    * Default constructor which creates a RapDipoles object without an associated Grid matrix. Though no
    * coordinates are assigned to the located sources.
    */
    RapDipoles();

    //=========================================================================================================
    /**
    * Constructor which creates a RapDipoles object with a Grid matrix. To located dipoles coordinates are
    * assigned.
    *
    * @param[in] p_pMatGrid The Grid matrix (n x 3) with coordinates (x y z). This grid must fit to the LeadField
    *                       which is used in the RAP MUSIC algorithm (n = number of grid points).
    */
    RapDipoles(HPCMatrix<T>* p_pMatGrid);

 
    //=========================================================================================================
    /**
    * The destructor destroys the RapDipoles object properly by performing a cleanup. The member function
    * clean() is called.
    */
    virtual ~RapDipoles();


    //=========================================================================================================
    /**
    * Initializes the RapDipoles with a Grid matrix.
    *
    * @param[in] p_pMatGrid The Grid matrix (n x 3) with coordinates (x y z). This grid must fit to the LeadField
    *                       which is used in the RAP MUSIC algorithm (n = number of grid points).
    * @return   True if successful initialized, false otherwise.
    */
    bool initRapDipoles(HPCMatrix<T>* p_pMatGrid);


    //=========================================================================================================
    /**
    * Returns whether the grid is initialized.
    *
    * @return   True when grid is available.
    */
    bool gridInitialized()
    {
        return m_bGridInitialized;
    }


    //=========================================================================================================
    /**
    * Adds a new correlated dipole pair to th RapDipoles. This function is called by the RAP MUSIC Algorithm.
    *
    * @param[in] p_iDipoleIdx1  Index (Lead Field grid index) of the first dipole.
    * @param[in] p_iDipoleIdx2  Index (Lead Field grid index) of the second dipole.
    * @param[in] p_vec_phi_k_1  Array of the dipole directories (phi_x1, phi_y1, phi_z1, phi_x2, phi_y2, phi_z2).
    * @param[in] p_valCor   Correlation value of the dipole pair.
    */
    void insertSource(  int p_iDipoleIdx1, int p_iDipoleIdx2,
                        const T *p_vec_phi_k_1,
                        T p_valCor);

    
    //=========================================================================================================
    /**
    * Displays the RapDipoles
    */
    void print();
    

    //=========================================================================================================
    /**
    * Performs a cleanup of the RapDipoles object.
    */
    void clean();

private:
    HPCMatrix<T>* m_pMatGrid;   /**< Holds the Grid of coordinates (n x 3); n = number of grid points. */

    bool m_bGridInitialized;    /**< Whether the grid is initialized. */

    std::vector<RapDipolePair<T>*> m_vecRapDipolePairs; /**< A vector of the located correlated dipoles. */

}; // class

} // NAMESPACE

//Make the template definition visible to compiler in the first point of instantiation
#include "../src/rapdipoles.cpp"

#endif // RAPDIPOLES_H