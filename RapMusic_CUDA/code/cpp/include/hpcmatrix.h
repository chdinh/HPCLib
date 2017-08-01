//=============================================================================================================
/**
* @file     hpcmatrix.h
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
* @brief    Contains the declaration of the HPCMatrix class.
*
*/


#ifndef HPCMATRIX_H
#define HPCMATRIX_H


//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "nomenclature.h"
#include "error.h"


//*************************************************************************************************************
//=============================================================================================================
// STL INCLUDES
//=============================================================================================================

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>


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

template <class T> class HPCMatrix ;

// Predefining every friend functions 
// This is required by latest ISO C++ draft
template <class T> std::istream& operator>>(std::istream& is, HPCMatrix<T>& ary);
template <class T> std::ostream& operator<<(std::ostream& os, const HPCMatrix<T>& ary);


template<class T>
void initHPCMatrix(HPCMatrix<T>& mat, const int m, const int n);



template<class T>
void resizeKeepHostMatrix(HPCMatrix<T>& ary, const int nr, const int nc);

template <class T> HPCLib::HPCMatrix<T> operator+(const HPCLib::HPCMatrix<T>&,const HPCLib::HPCMatrix<T>&);
template <class T> HPCLib::HPCMatrix<T> operator-(const HPCLib::HPCMatrix<T>&,const HPCLib::HPCMatrix<T>&);
template <class T> HPCLib::HPCMatrix<T> operator*(const HPCLib::HPCMatrix<T>&,const HPCLib::HPCMatrix<T>&);
template <class T> HPCLib::HPCMatrix<T> operator*(const double,const HPCLib::HPCMatrix<T>&);
// template <class T> HPCLib::HPCMatrix<T> operator*(const std::complex<T>&,const HPCLib::HPCMatrix<T>&);
// template <class T> HPCLib::Vector<T> operator*(const HPCLib::HPCMatrix<T>&,const HPCLib::Vector<T>&);
template <class T> int operator==(const HPCLib::HPCMatrix<T>&,const HPCLib::HPCMatrix<T>&);
template <class T> int operator!=(const HPCLib::HPCMatrix<T>& a,const HPCLib::HPCMatrix<T>& b) ;

// template <> HPCLib::HPCMatrix<Complex>  operator*(const double d, const HPCLib::HPCMatrix<Complex> &a);
// template <> HPCLib::HPCMatrix<Complex>  operator*(const Complex &d, const HPCLib::HPCMatrix<Complex> &a);





//=============================================================================================================
/**
* @brief    The HPCMatrix class provides a matrix implementation for general use in the HPC lib.
*
* ToDo Detailed description
*/
template<class T>
class HPCMatrix
{
public:

    template<class T> friend class cuHPCMatrix;

    //=========================================================================================================
    /**
    * Default constructor
    */
    HPCMatrix();

    //=========================================================================================================
    /**
    * Constructor which initializes the HPCMatrix with given size and allocates needed memory.
    *
    * @param[in] p_iRows    The number of rows.
    * @param[in] p_iCols    The number of columns.
    * @param[in] p_bColumnMajor Whether Matrix should be stored in column major (true by default).
    */
    HPCMatrix(int p_iRows, int p_iCols, bool p_bColumnMajor = true);

    //=========================================================================================================
    /**
    * Copy Constructor
    *
    * @param[in] p_matHPCMatrix The matrix which should be copied.
    */
    HPCMatrix(const HPCMatrix<T>& p_matHPCMatrix);

    //=========================================================================================================
    /**
    * Constructor
    *
    * @param[in] p_pData    A pointer to the data of the matrix (they should be allocated with new).
    *                       !!Pay attention when HPCMatrix is destroyed this data are only deleted,when cleanup
    *                       is set true!!
    * @param[in] p_iRows    The number of rows.
    * @param[in] p_iCols    The number of columns.
    * @param[in] p_bCleanup Whether the data pointer should be deleted when HPCMatrix is destroyed.
    * @param[in] p_bColumnMajor Whether the data array of p_pData is stored in column major (true by default).
    */
    HPCMatrix(T *p_pData, int p_iRows, int p_iCols, bool p_bCleanup = false, bool p_bColumnMajor = true);

 
    //=========================================================================================================
    /**
    * The destructor destroys the HPCMatrix object properly by performing a cleanup. The member function
    * clean() is called.
    */
    virtual ~HPCMatrix();


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
        return m_pData;
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



    void resize(const int m, const int n); 
    void resize(const HPCMatrix<T>& A) { resize(A.rows(),A.cols()) ; }//attention all data are erased - use resizeKeep instead
    void resizeKeep(const int nr, const int nc) { resizeKeepHostMatrix(*this,nr,nc) ; }


    //=========================================================================================================
    /**
    * Assigns a value to all elements of a matrix
    *
    * @param val value which should be assigned.
    */
    void reset(const T val = 0.0);


    HPCMatrix<T>& operator=(const HPCMatrix<T>& mat);

    T operator=(const T v) { reset((T)0); diag(v); return v; }//assigns the value only to the diagonal
    //T operator=(const T val) { reset(val); return val; }//assins the value to the whole matrix


    T* operator[](const int i) { return m_pVMat[i]; } 
    T* operator[](const int i) const { return m_pVMat[i];}

    //=========================================================================================================
    /**
    * Implements the operator Matrix(i,j) which accesses a specific element by reference.
    *
    * @param[in] i  The element row.
    * @param[in] j  The element column.
    * @return A reference to the element at position (row, column) in this matrix.
    */
    T& operator()(const int i,const int j)  { return elem(i,j); }
    //=========================================================================================================
    /**
    * Implements the operator Matrix(i,j) which accesses a specific element. Keeps the object unchanged.
    *
    * @param[in] i  The element row.
    * @param[in] j  The element column.
    * @return An element at position (row, column) in this matrix.
    */
    T  operator()(const int i,const int j) const { return elem(i,j); }


    //=========================================================================================================
    /**
    * Accesses a specific element by reference.
    *
    * @param[in] i  The element row.
    * @param[in] j  The element column.
    * @return A reference to the element at position (row, column) in this matrix.
    */
    T& elem(const int i,const int j) 
    { return m_bColumnMajor ? m_pVMat[j][i] : m_pVMat[i][j]; }  // no error message are generated if the index
                                                                // are out of range
    //=========================================================================================================
    /**
    * Accesses a specific element. Keeps the object unchanged.
    *
    * @param[in] i  The element row.
    * @param[in] j  The element column.
    * @return An element at position (row, column) in this matrix.
    */
    T  elem(const int i,const int j) const
    { return m_bColumnMajor ? m_pVMat[j][i] : m_pVMat[i][j]; }  // no error message are generated if the index
                                                                // are out of range


    void submatrix(int i, int j, HPCMatrix<T>&); 
    void as(int rw, int cl, HPCMatrix<T>&) ;
    HPCMatrix<T> get(int rw, int cl, int nr, int nc) const ;


    //=========================================================================================================
    /**
    * Attaches the matrix to an output stream and returns the stream.
    *
    * @param[in] os  The out stream to which the matrix is attached.
    * @return The output stream with the attached matrix.
    */
    std::ostream& print(std::ostream& os) const;

    //=========================================================================================================
    /**
    * Implements the operator (>>) which performs an input operation on a stream. A matrix can be stored
    * through a input stream to a matrix. Returns the input stream with the hold matrix to be attached.
    * Example: srcMat >> dstMat;
    *
    * @param[in] is The input stream which holds the matrix to be stored.
    * @param[in] ary    The matrix to which the stream should be assigned.
    * @return The input stream which hold the matrix to be attached.
    */
    friend std::istream& operator>>  <>(std::istream& is, HPCMatrix<T>& ary);
    //=========================================================================================================
    /**
    * Implements the operator (<<) which performs an output operation on a matrix by calling the member
    * function print(). The matrix is attached to the output stream. Example: std::cout << matrix;
    *
    * @param[in] os  The out stream to which the matrix is attached.
    * @param[in] ary  The matrix which should be attached to the output stream.
    * @return The output stream with the attached matrix.
    */
    friend std::ostream& operator<<  <>(std::ostream& os, const HPCMatrix<T>& ary);



    // Mathematical functions
    HPCMatrix<T>& operator+=(double d);
    HPCMatrix<T>& operator-=(double d);
    HPCMatrix<T>& operator*=(double d);
    HPCMatrix<T>& operator/=(double d);

    HPCMatrix<T>& operator+=(const HPCMatrix<T>&);
    HPCMatrix<T>& operator-=(const HPCMatrix<T>&);



    friend HPCMatrix<T> operator+  <>(const HPCMatrix<T>&, const HPCMatrix<T>&);
    friend HPCMatrix<T> operator-  <>(const HPCMatrix<T>&, const HPCMatrix<T>&);
    friend HPCMatrix<T> operator*  <>(const HPCMatrix<T>&, const HPCMatrix<T>&);
    friend HPCMatrix<T> operator*  <>(const double, const HPCMatrix<T>&);
//    friend HostMatrix<T> operator*  <>(const std::complex<T>&, const HostMatrix<T>&);
//    friend Vector<T> operator*  <>(const HostMatrix<T>&, const Vector<T>&);
    friend int operator==  <>(const HPCMatrix<T>&, const HPCMatrix<T>&);
    friend int operator!=  <>(const HPCMatrix<T>& a, const HPCMatrix<T>& b) ;



    HPCMatrix<T> herm() const ;
    HPCMatrix<T> transpose() const ;
    HPCMatrix<T> flop() const ;
    T trace() const ;

    double norm(void) ;
    void diag(const T fv);
    //  Vector<T> getDiag(); 

    void qSort() ;


    void setIdentity()
    { this->diag(1); }

    void setZero()
    { this->reset(0); }


    //=========================================================================================================
    /**
    * Initializes a HPCMatrix with the given number of rows (m) and columns (n). It also allocates the needed
    * memory.
    *
    * @param[out] mat  Matrix which should be initialized.
    * @param[in] m  Number of rows.
    * @param[in] n  Number of cols.
    */
    friend void initHPCMatrix<>(HPCMatrix<T> &mat, const int m,const int n);

    friend void resizeKeepHostMatrix<>(HPCMatrix<T>&, const int, const int);

protected:

    //=========================================================================================================
    /**
    * Initializes the matrix with the given number of rows (m) and columns (n) by calling HPCMatrix().
    *
    * @param[in] m  Number of rows.
    * @param[in] n  Number of cols.
    */
    void init(const int m = 1, const int n = 1) 
    { initHPCMatrix(*this,m,n); }

private:
    bool m_bColumnMajor;    /**< Holds whether matrix is stored in column major. */

    int m_iRows;    /**< Number of rows. */
    int m_iCols;    /**< Number of columns. */

    int m_iSize;    /**< Size of data array, which is calculated by rows*columns. */

    T *m_pData;     /**< The data array. */
    T **m_pVMat;    /**< Holds the mapping from 1D data array to a 2D matrix. */

    bool m_bCreated;    /**< True when data array was created by this objected false otherwise, when false
                             data array is not deleted when cleanup is performed */

};//Class


//*************************************************************************************************************

//=============================================================================================================
/**
* Initializes the data of the given pointer.
*
* @param[out] pointer   Pointer to data array which should be initialized.
* @param[in] size   Size of the data pointer.
*/
template <class T>
inline void initMemory(T* pointer, int size){
    T* p1;
    p1 = pointer-1;
    for(int i = size; i > 0; --i)
        *(++p1) = T() ;  
}


//*************************************************************************************************************

template <class T>
void initHPCMatrix(HPCMatrix<T> &mat, const int m,const int n)
{
    mat.clean();

    if ( m <= 0  || n <= 0 )
    {
        //todo throw error
        return ;
    }

    mat.m_iRows = m;     mat.m_iCols = n;

    mat.m_iSize = mat.m_iRows*mat.m_iCols ;


    mat.m_pData = new T [mat.m_iSize];
    mat.m_bCreated = true ;

    initMemory(mat.m_pData,mat.m_iSize);

    int i;

    //Map array to 2d array
    if (mat.m_bColumnMajor)
    {
        mat.m_pVMat = new T*[mat.m_iCols];
        for(i = mat.m_iCols-1; i >= 0; --i)
            mat.m_pVMat[i] = &mat.m_pData[i*mat.m_iRows] ;
    }
    else
    {
        mat.m_pVMat = new T*[mat.m_iRows];
        for(i = mat.m_iRows-1; i >= 0; --i)
            mat.m_pVMat[i] = &mat.m_pData[i*mat.m_iCols] ;
    }
}



//*************************************************************************************************************

template <class T>
void resizeKeepHostMatrix(HPCMatrix<T> &mat,const int nr,const int nc)
{

    if(nr == mat.m_iRows && nc == mat.m_iCols){
        return ;
    }

    T *mn ;
    T *p,*pn ;

    mn = new T[nr*nc] ;

    int i,j ;

    //this has to be tested
    if(mat.m_bColumnMajor)
    {
        for(j=0 ;j < std::min(nc,mat.m_iCols); j++){
            p = &mat.m_pData[j*mat.m_iRows] -1;
            pn = &mn[j*nr] -1 ;

            for(i=0; i < std::min(nr,mat.m_iRows); i++)
                *(++pn) = *(++p) ;

            for(i=mat.m_iRows;i<nr;i++)
                *(++pn) = T() ;
        }

        for(j=mat.m_iCols;j<nc;j++){
            pn = &mn[j*nr]-1 ;
            for(i=0;i<nr;i++)
                *(++pn) = T() ;
        }
    }
    else
    {
        for(i=0;i<minimum(nr,mat.m_iRows);i++){
            p = &mat.m_pData[i*mat.m_iCols] -1;
            pn = &mn[i*nc] -1 ;
            for(j=0;j<minimum(nc,mat.m_iCols);j++){
                *(++pn) = *(++p) ;
            }
            for(j=mat.m_iCols;j<nc;j++)
                *(++pn) = T() ;
        }

        for(i=mat.m_iRows;i<nr;i++){
            pn = &mn[i*nc]-1 ;
            for(j=0;j<nc;j++)
                *(++pn) = T() ;
        }
    }

    mat.m_iRows = nr ;
    mat.m_iCols = nc ;

    if(mat.m_pData && mat.m_bCreated)
        delete []mat.m_pData ;
    mat.m_pData = mn ;
    if(mat.m_pVMat)
        delete []mat.m_pVMat ;

    if(mat.m_bColumnMajor)
    {
        mat.m_pVMat = new T* [mat.m_iCols] ;
        for(i=0;i<mat.cz;++i)
            mat.m_pVMat[i] = &mat.m[i*mat.m_iRows] ;
    }
    else
    {
        mat.m_pVMat = new T* [mat.m_iRows] ;
        for(i=0;i<mat.m_iRows;++i)
            mat.m_pVMat[i] = &mat.m[i*mat.m_iCols] ;
    }
}



//*************************************************************************************************************

template <class T>
std::istream& operator>> (std::istream& is, HPCMatrix<T>& a)
{
    int i, j;
    int r = a.rows(), c = a.cols();

    if ( a.m_bColumnMajor )
        for (j = 0; j < c; j++)
            for (i = 0; i < r; i++)
                is >>  a.elem(i,j) ;
    else
        for (i = 0; i < r; i++)
            for (j = 0; j < c; j++)
                is >>  a.elem(i,j) ;

    return is;    
}


//*************************************************************************************************************

template <class T>
std::ostream& operator<<(std::ostream& os, const HPCMatrix<T>& mat)
{
    return mat.print(os) ; 
}












//Matrix Mathematical Operations
//*************************************************************************************************************

template <class T>
HPCMatrix<T> operator+(const HPCMatrix<T> &a,const HPCMatrix<T> &b)
{
    HPCMatrix<T> sum(a) ;
    sum += b ;
    return sum;
}


//*************************************************************************************************************

template <class T>
HPCMatrix<T> operator-(const HPCMatrix<T> &a,const HPCMatrix<T> &b)
{
    HPCMatrix<T> diff(a) ;
    diff -= b ;
    return diff;
}


//*************************************************************************************************************

template <class T>
HPCMatrix<T> operator*(const HPCMatrix<T> &a,const HPCMatrix<T> &b)
{
    if ( a.cols() != b.rows() )
    {
        #ifdef USE_EXCEPTION
            throw WrongSize2D(a.rows(),a.cols(),b.rows(),b.cols()) ;
        #else
            Error error("HPCMatrix<T> operator*(HPCMatrix<T>&,HPCMatrix<T>&)");
            error << "HPCMatrix<T> a * HPCMatrix<T> b incommensurate, a.cols() = " << a.cols() << ", b.rows() = " << b.rows() << std::endl ;
            error.fatal() ;
        #endif
    }

    if ( a.m_bColumnMajor != b.m_bColumnMajor )
    {
        Error error("HPCMatrix<T> operator*(HPCMatrix<T>&,HPCMatrix<T>&)");
        error << "HPCMatrix<T> a and HPCMatrix<T> b do not have the same column-major" << std::endl ;
        error.fatal() ;
    }

    int i, j, k, row=a.rows(), col=b.cols(),size = a.cols();
    HPCMatrix<T> prod(row,col);
    T zero = (T)0;

    if ( a.m_bColumnMajor )
    {
        for(i=row-1;i>=0;--i)
            for(j=size-1;j>=0;--j)
                if(a(i,j) != zero){
                    for(k=col-1;k>=0; --k)
                        prod(i,k) += a(i,j)* b(j,k) ;
                }
    }
    else
    {
        T *pptr,*aptr,*bptr ;
        aptr = a.m_pData ;
        for (i = 0; i < row; ++i)
            for (j = 0; j < size; ++j){
                if ( *aptr != zero )
                {
                    pptr = prod[i]-1;
                    bptr = b[j]-1;
                    for (k = col; k > 0; --k){
                        *(++pptr) += *aptr * *(++bptr);
                    }
                }
                ++aptr;
            }
    }
    return prod;
}


//*************************************************************************************************************

template <class T>
HPCMatrix<T>  operator*(const double d, const HPCMatrix<T> &a)
{
    int i, size=a.rows()*a.cols() ;
    HPCMatrix<T> b(a.rows(),a.cols());

    T *bptr,*aptr ;
    bptr = b.m_pData - 1 ;
    aptr = a.m_pData - 1 ;
    for (i = size; i > 0; --i)
        *(++bptr) = (T)(d * (*(++aptr))) ;
    return b;
}

// 
// //*************************************************************************************************************
// 
// template <class T>
// HPCMatrix<T>  operator*(const std::complex<T> &d, const HPCMatrix<T> &a)
// {
//     int i, size=a.rows()*a.cols() ;
//     HPCMatrix<T> b(a.rows(),a.cols());
// 
//     T *bptr,*aptr ;
//     bptr = b.m_pData - 1 ;
//     aptr = a.m_pData - 1 ;
//     for (i = size; i > 0; --i)
//         *(++bptr) = d.real() * *(++aptr) ;
//     return b;
// }


//*************************************************************************************************************

// template <class T>
// Vector<T> operator*(const HPCMatrix<T> &a, const Vector<T> &x)
// {
// 	if ( a.cols() != x.size() )
// 	{
// 		#ifdef USE_EXCEPTION
// 			throw WrongSize2D(a.rows(),a.cols(),x.size(),1);
// 		#else
// 			Error error("HPCMatrix<T> operator*(HPCMatrix<T>& a,Vector<T>& b)");
// 			error << "a * b incommensurate, a.cols() = " << a.cols() << ", b.rows() = " << x.size() << endl ;
// 			error.fatal() ;
// 		#endif
// 	}
// 
// 	int i, k, row = a.rows(), size = a.cols();
// 	Vector<T> prod(row);
// 	T zero = (T)0;
// 
// 	T *pptr,*aptr,*xptr ;
// 	aptr = a.m_pData - 1 ;
// 	pptr = &prod[0] ;
// 	for (i = row; i > 0; --i){
// 		xptr = x.memory()-1 ;
// 		for (k = size, *pptr = zero; k > 0 ; --k)
// 			*pptr += *(++aptr) * *(++xptr) ;
// 		++pptr ;
// 	}
// 
// 	return prod;
// }














//*************************************************************************************************************

template <class T>
int operator==(const HPCMatrix<T> &a,const HPCMatrix<T> &b)
{
    if ( a.rows() != b.rows() || a.cols() != b.cols() )
    {
        #ifdef USE_EXCEPTION
            throw WrongSize2D(a.rows(),a.cols(),b.rows(),b.cols()) ;
        #else
            Error error("operator==(const HostMatrix<T>&,const HostMatrix<T>&)");
            if ( b.rows() != a.rows() )
                error << "Matrices are of diferent size, a.rows() = " << a.rows() << " and b.rows() = " << b.rows() << endl ;
            if ( b.cols() != a.cols())
                error << "Matrices are of diferent size, a.cols() = " << a.cols() << " and b.cols() = " << b.cols() << endl ;
            error.fatal() ;
        #endif
    }

    int i, j, row = a.rows(), col = a.cols();
    int l = 1;

    for (i = 0; i < row; ++i)
        for (j = 0; j < col; ++j)
            l = l && ( a.elem(i,j) == b.elem(i,j) );

    return l;
}


} // NAMESPACE


//Make the template definition visible to compiler in the first point of instantiation
#include "../src/hpcmatrix.cpp"

#endif // HPCMATRIX_H