//=============================================================================================================
/**
* @file     cuhpcvalue.cuh
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


#ifndef CUHPCVALUE_CUH
#define CUHPCVALUE_CUH


//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include <iostream>
#include <iomanip>

#include "handle_error.cuh"




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

template <class T> class cuHPCValue;


// Predefining every friend functions 
// This is required by latest ISO C++ draft
//template <class T> std::istream& operator>>(std::istream& is, cuHPCValue<T>& ary);
template <class T> std::ostream& operator<<(std::ostream& os, cuHPCValue<T>& ary);



template<class T>
T toHostValue(cuHPCValue<T> &p_HPCVal);

template<class T>
void initCUHPCValue(cuHPCValue<T>& p_HPCVal);



//=============================================================================================================
/**
* DECLARE CLASS CudaDevice
*
* @brief ToDo
*/
template<class T>
class cuHPCValue
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
    cuHPCValue(T p_Val);
 

    //=========================================================================================================
    /**
    * dtor
    * Do garbage collecting
    */
    virtual ~cuHPCValue();
    

    //=========================================================================================================
    /**
    * Performs a cleanup of the HPCMatrix object.
    */
    void clean();



    //=========================================================================================================
    /**
    * The data of the matrix are stored in an array.
    *
    * @return The data array of the matrix.
    */
    T* data()
    {
        return m_dev_pValue;
    }


    //=========================================================================================================

    friend T toHostValue<>(cuHPCValue<T> &p_HPCVal);

    T toHostValue() { return HPCLib::toHostValue(*this); }

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
    friend std::ostream& operator<<  <>(std::ostream& os, cuHPCValue<T>& val);


    //=========================================================================================================
    /**
    * Initializes a HPCMatrix with the given number of rows (m) and columns (n). It also allocates the needed
    * memory.
    *
    * @param[out] mat  Matrix which should be initialized.
    * @param[in] m  Number of rows.
    * @param[in] n  Number of cols.
    */
    friend void initCUHPCValue<>(cuHPCValue<T> &p_HPCVal);


protected:

    //=========================================================================================================
    /**
    * Initializes the matrix with the given number of rows (m) and columns (n) by calling HPCMatrix().
    *
    * @param[in] m  Number of rows.
    * @param[in] n  Number of cols.
    */
    void init() 
    {
        initCUHPCValue(*this); 
    }

private:
    T* m_dev_pValue;     /**< The device data array. */

    bool m_bCreated;    /**< True when data array was created by this objected false otherwise, when false
                            data array is not deleted when cleanup is performed */

}; //Class


//*************************************************************************************************************

template <class T>
T toHostValue(cuHPCValue<T> &p_HPCVal)
{
    T t_Val;

    HANDLE_ERROR( cudaMemcpy( &t_Val, p_HPCVal.data(),
                            sizeof(T),
                            cudaMemcpyDeviceToHost ) );

    return t_Val;
}

//*************************************************************************************************************

template <class T>
void initCUHPCValue(cuHPCValue<T> &p_HPCVal)
{
    p_HPCVal.clean();

    //allocate the memory on the GPU
    HANDLE_ERROR( cudaMalloc( (void**)&p_HPCVal.m_dev_pValue, sizeof(T) ) );
    p_HPCVal.m_bCreated = true ;
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
std::ostream& operator<<(std::ostream& os, cuHPCValue<T>& val)
{
    return val.print(os); 
}



}// NAMESPACE

//Make the template definition visible to compiler in the first point of instantiation
#include "../src/cuhpcvalue.cu"

#endif // CUHPCVALUE_CUH