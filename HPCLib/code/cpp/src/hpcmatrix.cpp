//=============================================================================================================
/**
* @file     hpcmatrix.cpp
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
* @brief    Continuation of the HPCMatrix template declaration.
*
*/

#ifndef HPCMATRIX_SOURCES //Because this cpp is part of the header -> template
#define HPCMATRIX_SOURCES

//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

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
HPCMatrix<T>::HPCMatrix()
: m_bColumnMajor(false)
, m_pData(NULL)
, m_pVMat(NULL)
, m_bCreated(false)
{
    init(0,0);
}


//*************************************************************************************************************

template <class T>
HPCMatrix<T>::HPCMatrix(int p_iRows, int p_iCols, bool p_bColumnMajor)
: m_bColumnMajor(p_bColumnMajor)
, m_pData(NULL)
, m_pVMat(NULL)
, m_bCreated(false)
{
    init(p_iRows, p_iCols); 
}


//*************************************************************************************************************

template <class T>
HPCMatrix<T>::HPCMatrix(const HPCMatrix<T> &p_matHPCMatrix)
: m_bColumnMajor(p_matHPCMatrix.m_bColumnMajor)
, m_pData(NULL)
, m_pVMat(NULL)
, m_bCreated(false)
{
    init(p_matHPCMatrix.rows(),p_matHPCMatrix.cols());

    T *p1,*pa ;
    int sz = p_matHPCMatrix.rows()*p_matHPCMatrix.cols() ;
    p1 = m_pData-1 ;
    pa = p_matHPCMatrix.m_pData-1 ;

    for (int i = sz; i > 0 ; --i)
        *(++p1) = *(++pa) ;
}


//*************************************************************************************************************

template <class T>
HPCMatrix<T>::HPCMatrix(T *p_tData, int p_iRows, int p_iCols, bool p_bCleanup, bool p_bColumnMajor)
: m_bColumnMajor(p_bColumnMajor)
, m_pData(NULL)
, m_pVMat(NULL)
, m_bCreated(p_bCleanup)
{
    m_iRows = p_iRows;
    m_iCols = p_iCols;

    m_iSize = p_iRows * p_iCols;

    m_pData = p_tData;

    //Map array to 2d array
    if(m_bColumnMajor)
    {
        m_pVMat = new T* [m_iCols] ;
        for(int i=m_iCols-1;i>=0;--i)
            m_pVMat[i] = &m_pData[i*m_iRows] ;
    }
    else
    {
        m_pVMat = new T* [m_iRows] ;
        for(int i=m_iRows-1;i>=0;--i)
            m_pVMat[i] = &m_pData[i*m_iCols] ;
    }
}


//*************************************************************************************************************

template <class T>
HPCMatrix<T>::~HPCMatrix()
{
    clean();
}



//*************************************************************************************************************

template <class T>
void HPCMatrix<T>::resize(const int r, const int c)
{
/*    if(rows()>1 || cols()>1){
        if(m_pData && m_bCreated)
            delete []m_pData ;
        if(m_pVMat)
            delete []m_pVMat ;
    }
    else{
        if(m_pData && m_bCreated)
            delete []m_pData ;
        if(m_pVMat)
            delete []m_pVMat ;
    }
	*/
    init(r, c);
}


//*************************************************************************************************************

template <class T>
void HPCMatrix<T>::reset(const T v)
{
    T *p1 ;
    p1 = m_pData-1 ;
    const int size = rows()*cols() ;
    for(int i = size; i > 0; --i)
        *(++p1) = v ;
}


//*************************************************************************************************************

template <class T>
HPCMatrix<T>& HPCMatrix<T>::operator=(const HPCMatrix<T> &mat)
{
    int i;

    if ( this == &mat )
        return *this;

    m_bColumnMajor = mat.m_bColumnMajor;

    if(rows() != mat.rows() || cols() != mat.cols()){
        resize(mat.rows(),mat.cols()) ;
    }

    int sze = mat.rows()*mat.cols() ;

    T *p1,*pa ;
    p1 = m_pData-1 ;
    pa = mat.m_pData-1 ;

    for (i = sze; i > 0; --i)
        *(++p1) = *(++pa) ;

    return *this;
}


//*************************************************************************************************************

template <class T>
void HPCMatrix<T>::clean()
{
    m_iRows = 0;
    m_iCols = 0;
    m_iSize = 0;

    if(m_pData != NULL && m_bCreated)
	{
        delete []m_pData;
		m_pData = NULL;
	}

    if(m_pVMat != NULL)
	{
        delete []m_pVMat;
		m_pVMat = NULL;
	}

    m_bCreated = false;
}


//*************************************************************************************************************

template <class T>
std::ostream& HPCMatrix<T>::print(std::ostream& os) const
{
    int i, j;
    const int r = this->rows();
    const int c = this->cols();


    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++){
            os <<  std::setw(2) << elem(i,j) << ' ' ;
        }
        os << '\n';
    }

//     if ( m_bColumnMajor )
//         for (j = 0; j < c; j++)
//         {
//             for (i = 0; i < r; i++) {
//                 os <<  std::setw(2) << elem(i,j) << ' ';
//             }
//             os << '\n';
//         }
//     else
//         for (i = 0; i < r; i++)
//         {
//             for (j = 0; j < c; j++){
//                 os <<  std::setw(2) << elem(i,j) << ' ' ;
//             }
//             os << '\n';
//         }

        return os;
}




//MATRIX Stuff
//*************************************************************************************************************

template <class T>
void HPCMatrix<T>::submatrix(int sr, int sc, HPCMatrix<T> &a)
{
    int rwz,coz,i,j;

    if ( rows() % a.rows() != 0 || cols() % a.cols() != 0 || rows() < a.rows() || cols() < a.cols() )
    {
        #ifdef USE_EXCEPTION
            throw WrongSize2D(rows(),cols(),a.rows(),a.cols()) ;
        #else
            Error error("Matrix<T>::submatrix");
            error << "Matrix and submatrix incommensurate" ;
            error.fatal() ;
        #endif
    }

    if ( sr >= rows()/a.rows() || sr < 0 || sc >= cols()/a.cols() || sc < 0 )
    {
        #ifdef USE_EXCEPTION
            throw OutOfBound2D(sr,sc,0,rows()/a.rows()-1,0,cols()/a.cols()-1) ;
        #else
            Error error("Matrix<T>::submatrix");
            error << "Submatrix location out of bounds.\nrowblock " << sr << ", " << rows()/a.rows() << " colblock " << sc << ", " << a.cols() << std::endl ;
            error.fatal() ;
        #endif
    }
    rwz = sr*a.rows();
    coz = sc*a.cols();

    if ( m_bColumnMajor )
    {
        for ( i = a.rows()-1; i >= 0; --i )
            for(j=a.cols()-1;j>=0;--j)
                elem(i+rwz,j+coz) = a(i,j) ;
    }
    else
    {
        T *ptr, *aptr ;
        aptr = a.m_pData - 1;
        for ( i = a.rows()-1; i >= 0; --i )
        {
            ptr = &m_pData[(i+rwz)*cols()+coz]-1 ;
            for ( j = a.cols(); j > 0; --j)
                *(++ptr) = *(++aptr) ;
        }
    }
}


//*************************************************************************************************************

template <class T> 
void HPCMatrix<T>::as(int rw, int cl, HPCMatrix<T>& a)
{
    // Assign matrix a to this matrix at (i,j)
    int i, j;

    if ( (rw + a.rows()) > rows() || ( cl + a.cols()) > cols()) {
        #ifdef USE_EXCEPTION
            throw MatrixErr();
        #else
            Error error("HPCMatrix<T>::as") ;
            error << "Matrix A will not fit in this Matrix at " << rw << ", " << cl << std::endl ;
            error.fatal() ;
        #endif
    }

    if ( m_bColumnMajor )
    {
    for(i=0;i<a.rows();++i)
        for(j=0;j<a.cols();++j)
            elem(i+rw,j+cl) = a(i,j) ;
    }
    else
    {
        T *pptr,*aptr ;
        aptr = a.m_pData-1 ;
        for ( i = 0; i<a.rows(); ++i) {
            pptr = &m_pData[(i+rw)*cols()+cl]-1 ;
            for ( j = 0; j < a.cols(); ++j)
                *(++pptr) = *(++aptr);
        }
    }
}


//*************************************************************************************************************

template <class T> 
HPCMatrix<T> HPCMatrix<T>::get(int rw, int cl, int nr, int nc) const
{
    HPCMatrix<T> getmat(nr,nc) ;
    if ( (rw+nr) > rows() || (cl+nc) > cols()) {
        #ifdef USE_EXCEPTION
            throw MatrixErr();
        #else
            Error error("Matrix<T>::get") ;
            error << "Matrix of size "<< nr << ", " << nc << " not available at " << rw << ", " << cl << std::endl ;
            error.fatal() ;
        #endif
    }

    int i, j;


    if ( m_bColumnMajor )
    {
        for(i=0;i<nr;++i)
            for(j=0;j<nc;++j)
                getmat(i,j) = elem(i+rw,j+cl) ;
    }
    else
    {
        T *pptr,*aptr ;
        aptr = getmat.m_pData-1;
        for (i = 0; i < nr; ++i) {
            pptr = &m_pData[(i+rw)*cols()+cl]-1 ;
            for ( j = 0; j < nc; ++j)
                *(++aptr) = *(++pptr) ;
        }
    }

    return getmat ;
}

//Matrix Mathematical functions
//*************************************************************************************************************

template <class T> 
HPCMatrix<T>& HPCMatrix<T>::operator+=(double d)
{
    T *p1 ;
    p1 = m_pData-1 ;
    const int size = rows()*cols() ;
    for(int i = size; i > 0; --i)
        *(++p1) += d ;  
    return *this ;
}


//*************************************************************************************************************

template <class T> 
HPCMatrix<T>& HPCMatrix<T>::operator-=(double d)
{
    T *p1 ;
    p1 = m_pData-1 ;
    const int size = rows()*cols() ;
    for(int i=size; i>0; --i)
        *(++p1) -= d ;  
    return *this ;
}


//*************************************************************************************************************

template <class T> 
HPCMatrix<T>& HPCMatrix<T>::operator*=(double d)
{
    T *p1 ;
    p1 = m_pData-1 ;
    const int size = rows()*cols() ;
    for(int i=size; i>0; --i)
        *(++p1) *= d ;  
    return *this ;
}


//*************************************************************************************************************

template <class T> 
HPCMatrix<T>& HPCMatrix<T>::operator/=(double d)
{
    T *p1 ;
    p1 = m_pData-1 ;
    const int size = rows()*cols() ;
    for(int i=size; i>0; --i)
        *(++p1) /= d ;  
    return *this ;
}


//*************************************************************************************************************

template <class T> 
HPCMatrix<T>& HPCMatrix<T>::operator+=(const HPCMatrix<T> &mat)
{
    if ( mat.rows() != rows() || mat.cols() != cols() )
    {
        #ifdef USE_EXCEPTION
            throw WrongSize2D(rows(),cols(),mat.rows(),mat.cols());
        #else
            Error error("HPCMatrix<T>::operator+=") ;
            if ( rows() != mat.rows() )
                error << "Matrices are of diferent size, a.rows() = " << rows() << " and b.rows() = " << mat.rows() << std::endl ;
            if ( cols() != mat.cols())
                error << "Matrices are of diferent size, a.cols() = " << cols() << " and b.cols() = " << mat.cols() << std::endl ;
            error.fatal() ;
        #endif
    }

    int i, sze ;
    T *aptr,*sptr ;
    aptr = mat.m_pData - 1 ;
    sptr = m_pData - 1 ;
    sze = rows()*cols() ;
    for (i = sze; i > 0; --i){
        *(++sptr) += *(++aptr) ;
    }
    return *this ;
}


//*************************************************************************************************************

template <class T> 
HPCMatrix<T>& HPCMatrix<T>::operator-=(const HPCMatrix<T> &mat)
{
    if ( mat.rows() != rows() || mat.cols() != cols() )
    {
        #ifdef USE_EXCEPTION
            throw WrongSize2D(rows(),cols(),mat.rows(),mat.cols());
        #else
            Error error("HPCMatrix<T>::operator-=") ;
            if ( rows() != mat.rows() )
                error << "Matrices are of diferent size, a.rows() = " << rows() << " and b.rows() = " << mat.rows() << std::endl ;
            if ( cols() != mat.cols())
                error << "Matrices are of diferent size, a.cols() = " << cols() << " and b.cols() = " << mat.cols() << std::endl ;
            error.fatal() ;
        #endif
    }

    int i, size;
    T *aptr,*sptr ;
    aptr = mat.m_pData - 1 ;
    sptr = m_pData - 1 ;
    size = rows()*cols() ;
    for (i = size; i > 0; --i){
        *(++sptr) -= *(++aptr) ;
    }
    return *this ;
}


//Matrix Math
//*************************************************************************************************************

template <class T>
HPCMatrix<T> HPCMatrix<T>::herm() const
{
    int i, j, r = cols(), c = rows();
    HPCMatrix<T> adj(r,c);

    for (i = 0; i < r; ++i)
        for (j = 0; j < c; ++j)
            adj.elem(i,j) = elem(j,i) ;

    return adj;
}


//*************************************************************************************************************

template <class T>
HPCMatrix<T> HPCMatrix<T>::transpose() const
{                                       
    // same as hermitian for real Matrix<T>
    int i, j;
    const int& r = cols();
    const int& c = rows();
    HPCMatrix<T> adj(r,c);

    for (i = r-1; i >=0; --i)
        for (j = c-1; j >=0; --j)
            adj.elem(i,j) = elem(j,i) ;

    return adj; 
}


//*************************************************************************************************************

template <class T>
HPCMatrix<T> HPCMatrix<T>::flop() const
{                                       
    HPCMatrix<T> f(rows(),cols()) ;

    for(int i=rows()-1;i>=0;--i)
        for(int j=cols()-1;j>=0;--j)
            f(i,j) = elem(i,cols()-j-1);

    return f; 
}


//*************************************************************************************************************

template <class T>
T HPCMatrix<T>::trace() const
{
    int size = rows();
    T sum = (T)0;

    if ( size > cols() )
        size = cols();

    for (int d = 0; d < size; ++d)
        sum += elem(d,d) ;

    return sum;
}


//*************************************************************************************************************

template <class T> 
double HPCMatrix<T>::norm(void){
    int i,j ;
    double sum, maxsum;
    int init=0 ;
    T *pptr ;
    pptr = m_pData-1 ;
    maxsum = 0 ; // Silence the warning message
    for(i=0;i<rows();++i)
    {
        sum = 0 ;
        for ( j = 0; j < cols(); ++j) 
            sum += *(++pptr) ;
        if(init)
            maxsum = (maxsum>sum) ? maxsum : sum;
        else{
            maxsum = sum ;
            init = 1;
        }
    }
    return maxsum;
}


//*************************************************************************************************************

template <class T>
void HPCMatrix<T>::diag(const T a)
{
    int i, iend;

    iend = std::min(rows(),cols());

    for (i = iend-1; i >=0; --i)
        elem(i,i) = a;

}


//*************************************************************************************************************

//template <class T>
//Vector<T> Matrix<T>::getDiag(){
// 	int i, iend;
// 	Vector<T> vec(std::min(rows(),cols())) ;
// 	iend = std::min(rows(),cols());
// 	for (i = iend-1; i >=0; --i)
// 		vec[i] = elem(i,i);
// 	return vec ;
//}


}//Namespace

#endif //HPCMATRIX_SOURCES