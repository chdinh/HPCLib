//=============================================================================================================
/**
* @file     eigeninterface.h
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

#ifndef EIGENINTERFACE_H
#define EIGENINTERFACE_H


//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "../include/hpcmatrix.h"

//*************************************************************************************************************
//=============================================================================================================
// DEFINE NAMESPACE HPCLib
//=============================================================================================================

namespace HPCLib
{

class SVD
{
public:
    SVD(HPCMatrix<float>& matrix, unsigned int computationOptions = 0);

    ~SVD();

    HPCMatrix<float>* matrixU()
    {
        return m_pMatrixU;
    }

    HPCMatrix<float>* matrixV()
    {
        return m_pMatrixV;
    }

    HPCMatrix<float>* singularValues()
    {
        return m_pSingularValues;
    }



private:
    HPCMatrix<float>* m_pMatrixU;
    HPCMatrix<float>* m_pMatrixV;
    HPCMatrix<float>* m_pSingularValues;
};


class LU
{
public:
    LU(HPCMatrix<float>* matrix);

    ~LU();

    HPCMatrix<float> invert();

private:
    HPCMatrix<float>* m_pMatrix;
};




}//Namespace



#endif // EIGENINTERFACE_H