//=============================================================================================================
/**
* @file     eigeninterface.cpp
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

//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "../include/eigeninterface.h"


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

SVD::SVD(HPCMatrix<float>& matrix, unsigned int computationOptions)
: m_pMatrixU(NULL)
, m_pMatrixV(NULL)
, m_pSingularValues(NULL)
{
    Eigen::Map< Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> >
        MappedMat( matrix.data(), matrix.rows(), matrix.cols() );

    if(computationOptions == 0)
    {
        Eigen::JacobiSVD< Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> > 
            eigenSVD(MappedMat);

        m_pSingularValues = new HPCMatrix<float>(eigenSVD.singularValues().rows(), eigenSVD.singularValues().cols());

        memcpy( m_pSingularValues->data(),
            eigenSVD.singularValues().data(),
            eigenSVD.singularValues().rows()*eigenSVD.singularValues().cols()*sizeof(float));
    }
    else if(computationOptions == 1)
    {
        Eigen::JacobiSVD< Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> > 
            eigenSVD(MappedMat, Eigen::ComputeThinU);

        m_pMatrixU = new HPCMatrix<float>(eigenSVD.matrixU().rows(), eigenSVD.matrixU().cols());

        memcpy( m_pMatrixU->data(),
            eigenSVD.matrixU().data(),
            eigenSVD.matrixU().rows()*eigenSVD.matrixU().cols()*sizeof(float));

        m_pSingularValues = new HPCMatrix<float>(eigenSVD.singularValues().rows(), eigenSVD.singularValues().cols());

        memcpy( m_pSingularValues->data(),
            eigenSVD.singularValues().data(),
            eigenSVD.singularValues().rows()*eigenSVD.singularValues().cols()*sizeof(float));
    }
    else if(computationOptions == 2)
    {
        Eigen::JacobiSVD< Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> > 
            eigenSVD(MappedMat, Eigen::ComputeThinV);

        m_pMatrixV = new HPCMatrix<float>(eigenSVD.matrixV().rows(), eigenSVD.matrixV().cols());

        memcpy( m_pMatrixV->data(),
            eigenSVD.matrixV().data(),
            eigenSVD.matrixV().rows()*eigenSVD.matrixV().cols()*sizeof(float));


        m_pSingularValues = new HPCMatrix<float>(eigenSVD.singularValues().rows(), eigenSVD.singularValues().cols());

        memcpy( m_pSingularValues->data(),
            eigenSVD.singularValues().data(),
            eigenSVD.singularValues().rows()*eigenSVD.singularValues().cols()*sizeof(float));
    }
    else if(computationOptions == 3)
    {
        Eigen::JacobiSVD< Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> > 
            eigenSVD(MappedMat, Eigen::ComputeThinU | Eigen::ComputeThinV);

        m_pMatrixU = new HPCMatrix<float>(eigenSVD.matrixU().rows(), eigenSVD.matrixU().cols());

        memcpy( m_pMatrixU->data(),
            eigenSVD.matrixU().data(),
            eigenSVD.matrixU().rows()*eigenSVD.matrixU().cols()*sizeof(float));


        m_pMatrixV = new HPCMatrix<float>(eigenSVD.matrixV().rows(), eigenSVD.matrixV().cols());

        memcpy( m_pMatrixV->data(),
            eigenSVD.matrixV().data(),
            eigenSVD.matrixV().rows()*eigenSVD.matrixV().cols()*sizeof(float));


        m_pSingularValues = new HPCMatrix<float>(eigenSVD.singularValues().rows(), eigenSVD.singularValues().cols());

        memcpy( m_pSingularValues->data(),
            eigenSVD.singularValues().data(),
            eigenSVD.singularValues().rows()*eigenSVD.singularValues().cols()*sizeof(float));
    }
}

SVD::~SVD()
{
    if(m_pMatrixU != NULL)
        delete m_pMatrixU;

    if(m_pMatrixV != NULL)
        delete m_pMatrixV;

    if(m_pSingularValues != NULL)
        delete m_pSingularValues;
}












LU::LU(HPCMatrix<float>* matrix)
: m_pMatrix(matrix)
{

}

LU::~LU()
{

}

HPCMatrix<float> LU::invert()
{
    Eigen::Map< Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> >
        MappedMat( m_pMatrix->data(), m_pMatrix->rows(), m_pMatrix->cols() );

    HPCMatrix<float> t_matInv(m_pMatrix->rows(), m_pMatrix->cols());

    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> t_invMat = MappedMat.inverse();

    memcpy( t_matInv.data(), t_invMat.data(),
        (t_matInv.size()*sizeof(float)));

    return t_matInv;

}


}//Namespace