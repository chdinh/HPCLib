//=============================================================================================================
/**
* @file     cuhpcmatrix.cu
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
* @brief    ToDo Documentation...
*
*/

#ifndef CUHPCMATRIX_SOURCES //Because this cpp is part of the header -> template
#define CUHPCMATRIX_SOURCES

//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "../include/cuhpcmatrix.cuh"


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
cuHPCMatrix<T>::cuHPCMatrix(int p_iRows, int p_iCols, bool p_bColumnMajor)
: m_bColumnMajor(p_bColumnMajor)
, m_dev_pData(NULL)
, m_bCreated(false)
{
    init(p_iRows, p_iCols); 
}


//*************************************************************************************************************

template <class T>
cuHPCMatrix<T>::cuHPCMatrix(const HPCMatrix<T> &p_matHPCMatrix)
: m_bColumnMajor(p_matHPCMatrix.m_bColumnMajor)
, m_dev_pData(NULL)
, m_bCreated(false)
{
    init(p_matHPCMatrix.rows(),p_matHPCMatrix.cols());

    int sz = p_matHPCMatrix.rows()*p_matHPCMatrix.cols() ;

    //Set device to identity
    HANDLE_ERROR( cudaMemcpy(   m_dev_pData,
                                p_matHPCMatrix.m_pData,
                                sizeof(T) * sz,
                                cudaMemcpyHostToDevice ) );
}


//*************************************************************************************************************

template <class T>
cuHPCMatrix<T>::~cuHPCMatrix()
{
    clean();
}


//*************************************************************************************************************

template <class T>
void cuHPCMatrix<T>::clean()
{
    m_iRows = 0;
    m_iCols = 0;
    m_iSize = 0;

    if(m_dev_pData != NULL && m_bCreated)
    {
        HANDLE_ERROR( cudaFree( m_dev_pData ) );
        m_dev_pData = NULL;
    }

    m_bCreated = false;
}


//*************************************************************************************************************

template <class T>
std::ostream& cuHPCMatrix<T>::print(std::ostream& os)
{
    HPCMatrix<T> t_Mat = HPCLib::toHPCMatrix(*this);
    int i, j;
    const int r = t_Mat.rows();
    const int c = t_Mat.cols();


    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++){
            os <<  std::setw(2) << t_Mat.elem(i,j) << ' ' ;
        }
        os << '\n';
    }

//     if ( m_bColumnMajor )
//         for (j = 0; j < c; j++)
//         {
//             for (i = 0; i < r; i++) {
//                 os <<  std::setw(2) << t_Mat.elem(i,j) << ' ';
//             }
//             os << '\n';
//         }
//     else
//         for (i = 0; i < r; i++)
//         {
//             for (j = 0; j < c; j++){
//                 os <<  std::setw(2) << t_Mat.elem(i,j) << ' ' ;
//             }
//             os << '\n';
//         }

        return os;
}


}//Namespace

#endif //CUHPCMATRIX_SOURCES