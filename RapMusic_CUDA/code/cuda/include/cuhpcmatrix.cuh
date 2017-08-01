//=============================================================================================================
/**
* @file     cuhpcmatrix.cuh
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


#ifndef CUHPCMATRIX_CUH
#define CUHPCMATRIX_CUH


//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include <iostream>


#include "../../cpp/include/hpcmatrix.h"


//*************************************************************************************************************
//=============================================================================================================
// CUDA INCLUDES
//=============================================================================================================

#include "handle_error.cuh"
#include <cublas.h>



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
// FORWARD DECLERATIONS
//=============================================================================================================

template <class T> class cuHPCMatrix;


// Predefining every friend functions 
// This is required by latest ISO C++ draft
//template <class T> std::istream& operator>>(std::istream& is, cuHPCMatrix<T>& ary);
template <class T> std::ostream& operator<<(std::ostream& os, cuHPCMatrix<T>& ary);


template<class T>
void cuHPCMatMult(char TA, char TB, cuHPCMatrix<T> &A, cuHPCMatrix<T> &B, cuHPCMatrix<T> &C);


template<class T>
HPCMatrix<T> toHPCMatrix(cuHPCMatrix<T> &mat);


template<class T>
void initCUHPCMatrix(cuHPCMatrix<T>& mat, const int m, const int n);



//=============================================================================================================
/**
* DECLARE CLASS CudaDevice
*
* @brief ToDo
*/
template<class T>
class cuHPCMatrix
{
public:

    //=========================================================================================================
    /**
    * Constructor which initializes the device cuHPCMatrix with given size and allocates needed memory.
    *
    * @param[in] p_iRows    The number of rows.
    * @param[in] p_iCols    The number of columns.
    * @param[in] p_bColumnMajor Whether Matrix should be stored in column major (true by default).
    */
    cuHPCMatrix(int p_iRows, int p_iCols, bool p_bColumnMajor = true);


//     //=========================================================================================================
//     /**
//     * Copy Constructor
//     *
//     * @param[in] p_matHPCMatrix The matrix which should be copied.
//     */
//     cuHPCMatrix(const cuHPCMatrix<T>& p_matCUHPCMatrix);


    //=========================================================================================================
    /**
    * Copy Constructor which copies a HPCMatrix.
    *
    * @param[in] p_matHPCMatrix The matrix which should be copied.
    */
    cuHPCMatrix(const HPCMatrix<T>& p_matHPCMatrix);
 

    //=========================================================================================================
    /**
    * dtor
    * Do garbage collecting
    */
    virtual ~cuHPCMatrix();
    

    //=========================================================================================================
    /**
    * Performs a cleanup of the HPCMatrix object.
    */
    void clean();


    //=========================================================================================================
    /**
    * Whether the HPCMatrix is stored in column major.
    *
    * @return True when data are stored in column major.
    */
    bool isColumnMajor()
    {
        return m_bColumnMajor;
    }


    //=========================================================================================================
    /**
    * The data of the matrix are stored in an array.
    *
    * @return The data array of the matrix.
    */
    T* data()
    {
        return m_dev_pData;
    }
    

    //=========================================================================================================
    /**
    * The size of the data array of the matrix is calculated by rows*cols.
    *
    * @return The size of the data array.
    */
    inline int size() const
    { 
        return m_iSize;
    }


    //=========================================================================================================
    /**
    * The number of rows of the matrix.
    *
    * @return Matrix rows.
    */
    inline int rows() const
    { 
        return m_iRows;
    }


    //=========================================================================================================
    /**
    * The number of columns of the matrix.
    *
    * @return Matrix columns.
    */
    inline int cols() const
    { 
        return m_iCols;
    }



    //=========================================================================================================

    friend HPCMatrix<T> toHPCMatrix<>(cuHPCMatrix<T> &mat);

    HPCMatrix<T> toHPCMatrix() { return HPCLib::toHPCMatrix(*this); }

    //=========================================================================================================
    /**
    * Attaches the matrix to an output stream and returns the stream.
    *
    * @param[in] os  The out stream to which the matrix is attached.
    * @return The output stream with the attached matrix.
    */
    std::ostream& print(std::ostream& os);

//     //=========================================================================================================
//     /**
//     * Implements the operator (>>) which performs an input operation on a stream. A matrix can be stored
//     * through a input stream to a matrix. Returns the input stream with the hold matrix to be attached.
//     * Example: srcMat >> dstMat;
//     *
//     * @param[in] is The input stream which holds the matrix to be stored.
//     * @param[in] ary    The matrix to which the stream should be assigned.
//     * @return The input stream which hold the matrix to be attached.
//     */
//     friend std::istream& operator>>  <>(std::istream& is, cuHPCMatrix<T>& ary);
    //=========================================================================================================
    /**
    * Implements the operator (<<) which performs an output operation on a matrix by calling the member
    * function print(). The matrix is attached to the output stream. Example: std::cout << matrix;
    *
    * @param[in] os  The out stream to which the matrix is attached.
    * @param[in] ary  The matrix which should be attached to the output stream.
    * @return The output stream with the attached matrix.
    */
    friend std::ostream& operator<<  <>(std::ostream& os, cuHPCMatrix<T>& ary);



    friend void cuHPCMatMult<>(char TA, char TB, cuHPCMatrix<T> &A, cuHPCMatrix<T> &B, cuHPCMatrix<T> &C);

    void cuHPCMatMult(char TA, char TB, cuHPCMatrix<T> &A, cuHPCMatrix<T> &B) { HPCLib::cuHPCMatMult(TA, TB, A, B, *this); }


    //=========================================================================================================
    /**
    * Initializes a HPCMatrix with the given number of rows (m) and columns (n). It also allocates the needed
    * memory.
    *
    * @param[out] mat  Matrix which should be initialized.
    * @param[in] m  Number of rows.
    * @param[in] n  Number of cols.
    */
    friend void initCUHPCMatrix<>(cuHPCMatrix<T> &mat, const int m,const int n);


protected:

    //=========================================================================================================
    /**
    * Initializes the matrix with the given number of rows (m) and columns (n) by calling HPCMatrix().
    *
    * @param[in] m  Number of rows.
    * @param[in] n  Number of cols.
    */
    void init(const int m = 1, const int n = 1) 
    {
        initCUHPCMatrix(*this,m,n); 
    }

private:
    bool m_bColumnMajor;    /**< Holds whether matrix is stored in column major. */

    int m_iRows;    /**< Number of rows. */
    int m_iCols;    /**< Number of columns. */

    int m_iSize;    /**< Size of data array, which is calculated by rows*columns. */

    T* m_dev_pData;     /**< The device data array. */

    bool m_bCreated;    /**< True when data array was created by this objected false otherwise, when false
                            data array is not deleted when cleanup is performed */

}; //Class


//*************************************************************************************************************

template <class T>
HPCMatrix<T> toHPCMatrix(cuHPCMatrix<T> &p_devMat)
{
    HPCMatrix<T> t_matHPCMatrix(p_devMat.rows(),p_devMat.cols());

    HANDLE_ERROR( cudaMemcpy( t_matHPCMatrix.data(), p_devMat.data(),
                            p_devMat.size()*sizeof(T),
                            cudaMemcpyDeviceToHost ) );

    return t_matHPCMatrix;
}


//*************************************************************************************************************

template <class T>
void cuHPCMatMult(char TA, char TB, cuHPCMatrix<T> &A, cuHPCMatrix<T> &B, cuHPCMatrix<T> &C)
{
    if( TA == 'N' && TB == 'N') //C = A * B;
    {
        if(A.rows() != C.rows() || B.cols() != C.cols() || A.cols() != B.rows())
        {
            //ToDo throw error
            std::cout << "cuHPCMatMult 'N', 'N': Matrix Dimensions does not fit!" << std::endl;
            return;
        }

        cublasSgemm(TA, TB,
                A.rows(), /* leading dimension */
                B.cols(), /* trailing   "      */
                A.cols(), /* inner      "      */
                1, /* alpha */
                A.data(), A.rows(),//Column_Major
                B.data(), B.rows(),//Column_Major
                0, /* beta */
                C.data(), C.rows());//Column_Major
        HANDLE_ERROR( cudaThreadSynchronize() );
    }
    else if( TA == 'N' && TB == 'T')//C = A * B^T;
    {
        if(A.rows() != C.rows() || B.rows() != C.cols() || A.cols() != B.cols())
        {
            //ToDo throw error
            std::cout << "cuHPCMatMult 'N', 'T': Matrix Dimensions does not fit!" << std::endl;
            return;
        }

        cublasSgemm(TA, TB,
                A.rows(), /* leading dimension */
                B.rows(), /* trailing   "      */
                A.cols(), /* inner      "      */
                1, /* alpha */
                A.data(), A.rows(),//Column_Major
                B.data(), B.rows(),//Column_Major
                0, /* beta */
                C.data(), C.rows());//Column_Major
        HANDLE_ERROR( cudaThreadSynchronize() );
    }
}


//*************************************************************************************************************

template <class T>
void initCUHPCMatrix(cuHPCMatrix<T> &mat, const int m,const int n)
{
    mat.clean();

    if ( m <= 0  || n <= 0 )
    {
        //todo throw error
        return ;
    }

    mat.m_iRows = m;     mat.m_iCols = n;

    mat.m_iSize = mat.m_iRows*mat.m_iCols ;


    //allocate the memory on the GPU
    HANDLE_ERROR( cudaMalloc( (void**)&mat.m_dev_pData, sizeof(T) * mat.m_iSize ) );
    mat.m_bCreated = true ;

//     int i;
// 
//     //Map array to 2d array
//     if (mat.m_bColumnMajor)
//     {
//         mat.m_pVMat = new T*[mat.m_iCols];
//         for(i = mat.m_iCols-1; i >= 0; --i)
//             mat.m_pVMat[i] = &mat.m_pData[i*mat.m_iRows] ;
//     }
//     else
//     {
//         mat.m_pVMat = new T*[mat.m_iRows];
//         for(i = mat.m_iRows-1; i >= 0; --i)
//             mat.m_pVMat[i] = &mat.m_pData[i*mat.m_iCols] ;
//     }
}


//*************************************************************************************************************

// template <class T>
// std::istream& operator>> (std::istream& is, cuHPCMatrix<T>& a)
// {
//     int i, j;
//     int r = a.rows(), c = a.cols();
// 
//     if ( a.m_bColumnMajor )
//         for (j = 0; j < c; j++)
//             for (i = 0; i < r; i++)
//                 is >>  a.elem(i,j) ;
//     else
//         for (i = 0; i < r; i++)
//             for (j = 0; j < c; j++)
//                 is >>  a.elem(i,j) ;
// 
//     return is;    
// }


//*************************************************************************************************************

template <class T>
std::ostream& operator<<(std::ostream& os, cuHPCMatrix<T>& mat)
{
    return mat.print(os); 
}



}// NAMESPACE

//Make the template definition visible to compiler in the first point of instantiation
#include "../src/cuhpcmatrix.cu"

#endif // CUHPCMATRIX_CUH