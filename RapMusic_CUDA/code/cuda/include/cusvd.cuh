//=============================================================================================================
/**
* @file     cusvd.cuh
* @author   Christoph Dinh <christoph.dinh@live.de>; Johannes Ruehle <johannes.ruehle@tu-ilmenau.de>
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
* @brief    Contains the svd gpu implementations.
*
*/

#ifndef CUSVD_CUH
#define CUSVD_CUH


//*************************************************************************************************************
//=============================================================================================================
// CUDA INCLUDES
//=============================================================================================================

#include "cuhelpers.cuh"


//*************************************************************************************************************
//=============================================================================================================
// DEFINE NAMESPACE HPCLib
//=============================================================================================================

namespace HPCLib
{

//*************************************************************************************************************
//=============================================================================================================
// Kernels and Device Functions
//=============================================================================================================


    //=============================================================================================================
/**
* Calculates the SVD on the GPU.
*
* CODE of 3rd EDITION OF Numerical Recipes in C (2007)
* Given a matrix u[0..m-1][0..n-1] (u == a), this routine computes its singular value decomposition, A =
* U*W*V^T. Thematrix U replaces a on output. The diagonal matrix of singular values W is output
* as a vector w[0..n-1]. The matrix V (not the transpose V^T ) is output as v[0..n-1][0..n-1].
*
* @param[in,out] u  The column-major array of the input matrix [m x n]. After calculation has finished
*                   this variable holds the return matrix u which contains the Eigenvectors.
* @param[in] m  The number of rows of u.
* @param[in] n  The number of columns of u.
* @param[out] w  The singular values of the input matrix in a row vector [n x 1].
* @param[out] v  The matrix v [n x n] of the singular value decomposition of the input matrix.
*/
__device__ void cuSVD(  float *u,   /* [m x n ]*/
                        int m,      /* rows */
                        int n,      /* columns */
                        float *w,   /* [nx1]*/
                        float *v    /* [nxn]*/      )
{
    //double
    //double eps = 1.19e-16; //see page 1165, Numerical Recipes in C 3.Edition (2007); // eps = numeric_limits<Doub>::epsilon();
    //float
    float eps = 1.19e-7;

    bool flag;
    int i,its,j,jj,k,l,nm;
    float anorm,c,f,g,h,s,scale,x,y,z;

    float* rv1 = new float[n];//This is working since CUDA 4.0, on devices with compute capatibility 2.x

    //Part 1
    g = scale = anorm = 0.0; //Householder reduction to bidiagonal form.
    for (i=0;i<n;i++) {
        l=i+2;
        rv1[i]=scale*g;
        g=s=scale=0.0;
        if (i < m) {
            for (k=i;k<m;k++) scale += abs(u[k+i*m]);//ToDo parallel
            if (scale != 0.0) {
                for (k=i;k<m;k++) {//ToDo parallel
                    u[k+i*m] /= scale;
                    s += u[k+i*m]*u[k+i*m];
                }
                f=u[i+i*m];
                g = -SIGN(sqrt(s),f);
                h=f*g-s;
                u[i+i*m]=f-g;
                for (j=l-1;j<n;j++) {//ToDo parallel - nonsens cause n is too small
                    for (s=0.0,k=i;k<m;k++) s += u[k+i*m]*u[k+j*m];
                    f=s/h;
                    for (k=i;k<m;k++) u[k+j*m] += f*u[k+i*m];
                }
                for (k=i;k<m;k++) u[k+i*m] *= scale;//ToDo parallel
            }
        }
        w[i]=scale *g;
        g=s=scale=0.0;
        if (i+1 <= m && i+1 != n) {
            for (k=l-1;k<n;k++) scale += abs(u[i+k*m]);
            if (scale != 0.0) {
                for (k=l-1;k<n;k++) {
                    u[i+k*m] /= scale;
                    s += u[i+k*m]*u[i+k*m];

                }
                f=u[i + m*(l-1) ];
                g = -SIGN(sqrt(s),f);
                h=f*g-s;
                u[i + m*(l-1)]=f-g;
                for (k=l-1;k<n;k++) rv1[k]=u[i + m*k]/h;
                for (j=l-1;j<m;j++) {
                    for (s=0.0,k=l-1;k<n;k++) s += u[j + m*k]*u[i + m*k];
                    for (k=l-1;k<n;k++) u[j + m*k] += s*rv1[k];
                }
                for (k=l-1;k<n;k++) u[i + m*k] *= scale;
            }
        }
        anorm = DMAX(anorm,(abs(w[i])+abs(rv1[i])));
    }
    //Part 2
    for (i=n-1;i>=0;i--) {// Accumulation of right-hand transformations.
        if (i < n-1) {
            if (g != 0.0) {
                for (j=l;j<n;j++) //Double division to avoid possible underow.
                    v[j + n*i]=(u[i + m*j]/u[i + m*l])/g;
                for (j=l;j<n;j++) {
                    for (s=0.0,k=l;k<n;k++) s += u[i + m*k]*v[k + n*j];
                    for (k=l;k<n;k++) v[k + n*j] += s*v[k + n*i];
                }
            }
            for (j=l;j<n;j++) v[i + n*j]=v[j + n*i]=0.0;
        }
        v[i + n*i]=1.0;
        g=rv1[i];
        l=i;
    }
    //Part 3
    for (i=IMIN(m,n)-1;i>=0;i--) { //Accumulation of left-hand transformations.
        l=i+1;
        g=w[i];
        for (j=l;j<n;j++) u[i + m*j]=0.0;
        if (g != 0.0) {
            g=1.0/g;
            for (j=l;j<n;j++) {
                for (s=0.0,k=l;k<m;k++) s += u[k + m*i]*u[k + m*j];
                f=(s/u[i + m*i])*g;
                for (k=i;k<m;k++) u[k + m*j] += f*u[k + m*i];
            }
            for (j=i;j<m;j++) u[j + m*i] *= g;
        } else for (j=i;j<m;j++) u[j + m*i]=0.0;
        ++u[i + m*i];
    }
    //Part 4
    for (k=n-1;k>=0;k--) { //Diagonalization of the bidiagonal form: Loop over
        for (its=0;its<30;its++) { //singular values, and over allowed iterations.
            flag=true;
            for (l=k;l>=0;l--) {// Test for splitting.
                nm=l-1;
                if (l == 0 || abs(rv1[l]) <= eps*anorm) {
                    flag=false;
                    break;
                }
                if (abs(w[nm]) <= eps*anorm) break;
            }
            if (flag) {
                c=0.0; //Cancellation of rv1[l], if l > 0.
                    s=1.0;
                for (i=l;i<k+1;i++) {
                    f=s*rv1[i];
                    rv1[i]=c*rv1[i];
                    if (abs(f) <= eps*anorm) break;
                    g=w[i];
                    h=pythag(f,g);
                    w[i]=h;
                    h=1.0/h;
                    c=g*h;
                    s = -f*h;
                    for (j=0;j<m;j++) {
                        y=u[j + m*nm];
                        z=u[j + m*i];
                        u[j + m*nm]=y*c+z*s;
                        u[j + m*i]=z*c-y*s;
                    }
                }
            }
            z=w[k];
            if (l == k) { //Convergence.
                if (z < 0.0) { //Singular value is made nonnegative.
                    w[k] = -z;
                    for (j=0;j<n;j++) v[j + n*k] = -v[j + n*k];
                }
                break;
            }
//			if (its == 29) throw("no convergence in 30 svdcmp iterations");
            x=w[l]; //Shift from bottom 2-by-2 minor.
            nm=k-1;
            y=w[nm];
            g=rv1[nm];
            h=rv1[k];
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g=pythag(f,1.0);
            f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
            c=s=1.0;// Next QR transformation:
            for (j=l;j<=nm;j++) {
                i=j+1;
                g=rv1[i];
                y=w[i];
                h=s*g;
                g=c*g;
                z=pythag(f,h);
                rv1[j]=z;
                c=f/z;
                s=h/z;
                f=x*c+g*s;
                g=g*c-x*s;
                h=y*s;
                y *= c;
                for (jj=0;jj<n;jj++) {
                    x=v[jj + n*j];
                    z=v[jj + n*i];
                    v[jj + n*j]=x*c+z*s;
                    v[jj + n*i]=z*c-x*s;
                }
                z=pythag(f,h);
                w[j]=z;// Rotation can be arbitrary if z D 0.
                if (z) {
                    z=1.0/z;
                    c=f*z;
                    s=h*z;
                }
                f=c*g+s*y;
                x=c*y-s*g;
                for (jj=0;jj<m;jj++) {
                    y=u[jj + m*j];
                    z=u[jj + m*i];
                    u[jj + m*j]=y*c+z*s;
                    u[jj + m*i]=z*c-y*s;
                }
            }
            rv1[l]=0.0;
            rv1[k]=f;
            w[k]=x;
        }
    }

    delete [] rv1;
}


//=============================================================================================================
/**
* Calculates only the U an W Matrix of the SVD on the GPU.
*
* CODE of 3rd EDITION OF Numerical Recipes in C (2007)
* Given a matrix u[0..m-1][0..n-1] (u == a), this routine computes its singular value decomposition, A =
* U*W*V^T. Thematrix U replaces a on output. The diagonal matrix of singular values W is output
* as a vector w[0..n-1]. Thematrix V (not the transpose V^T ) is output as v[0..n-1][0..n-1].
*
* @param[in,out] u  The column-major array of the input matrix [m x n]. After calculation has finished
*                   this variable holds the return matrix u which contains the Eigenvectors.
* @param[in] m  The number of rows of u.
* @param[in] n  The number of columns of u.
* @param[out] w  The singular values of the input matrix in a row vector [n x 1].
*/
__device__ void cuSVD_UW(   float *u,   /* [m x n ]*/
                            int m,      /* rows */
                            int n,      /* columns */
                            float *w    /* [nx1]*/)
{
    //double
    //double eps = 1.19e-16; //see page 1165, Numerical Recipes in C 3.Edition (2007); // eps = numeric_limits<Doub>::epsilon();
    //float
    float eps = 1.19e-7;

    bool flag;
    int i,its,j,jj,k,l,nm;
    float anorm,c,f,g,h,s,scale,x,y,z;
    
    float* rv1 = new float[n];//This is working since CUDA 4.0, on devices with compute capatibility 2.x

    //Part 1
    g = scale = anorm = 0.0; //Householder reduction to bidiagonal form.
    for (i=0;i<n;i++) {
        l=i+2;
        rv1[i]=scale*g;
        g=s=scale=0.0;
        if (i < m) {
            for (k=i;k<m;k++) scale += abs(u[k+i*m]);
            if (scale != 0.0) {
                for (k=i;k<m;k++) {
                    u[k+i*m] /= scale;
                    s += u[k+i*m]*u[k+i*m];
                }
                f=u[i+i*m];
                g = -SIGN(sqrt(s),f);
                h=f*g-s;
                u[i+i*m]=f-g;
                for (j=l-1;j<n;j++) {
                    for (s=0.0,k=i;k<m;k++) s += u[k+i*m]*u[k+j*m];
                    f=s/h;
                    for (k=i;k<m;k++) u[k+j*m] += f*u[k+i*m];
                }
                for (k=i;k<m;k++) u[k+i*m] *= scale;
            }
        }
        w[i]=scale *g;
        g=s=scale=0.0;
        if (i+1 <= m && i+1 != n) {
            for (k=l-1;k<n;k++) scale += abs(u[i+k*m]);
            if (scale != 0.0) {
                for (k=l-1;k<n;k++) {
                    u[i+k*m] /= scale;
                    s += u[i+k*m]*u[i+k*m];

                }
                f=u[i + m*(l-1) ];
                g = -SIGN(sqrt(s),f);
                h=f*g-s;
                u[i + m*(l-1)]=f-g;
                for (k=l-1;k<n;k++) rv1[k]=u[i + m*k]/h;
                for (j=l-1;j<m;j++) {
                    for (s=0.0,k=l-1;k<n;k++) s += u[j + m*k]*u[i + m*k];
                    for (k=l-1;k<n;k++) u[j + m*k] += s*rv1[k];
                }
                for (k=l-1;k<n;k++) u[i + m*k] *= scale;
            }
        }
        anorm = DMAX(anorm,(abs(w[i])+abs(rv1[i])));
    }
    //Part 3
    for (i=IMIN(m,n)-1;i>=0;i--) { //Accumulation of left-hand transformations.
        l=i+1;
        g=w[i];
        for (j=l;j<n;j++) u[i + m*j]=0.0;
        if (g != 0.0) {
            g=1.0/g;
            for (j=l;j<n;j++) {
                for (s=0.0,k=l;k<m;k++) s += u[k + m*i]*u[k + m*j];
                f=(s/u[i + m*i])*g;
                for (k=i;k<m;k++) u[k + m*j] += f*u[k + m*i];
            }
            for (j=i;j<m;j++) u[j + m*i] *= g;
        } else for (j=i;j<m;j++) u[j + m*i]=0.0;
        ++u[i + m*i];
    }
    //Part 4
    for (k=n-1;k>=0;k--) { //Diagonalization of the bidiagonal form: Loop over
        for (its=0;its<30;its++) { //singular values, and over allowed iterations.
            flag=true;
            for (l=k;l>=0;l--) {// Test for splitting.
                nm=l-1;
                if (l == 0 || abs(rv1[l]) <= eps*anorm) {
                    flag=false;
                    break;
                }
                if (abs(w[nm]) <= eps*anorm) break;
            }
            if (flag) {
                c=0.0; //Cancellation of rv1[l], if l > 0.
                    s=1.0;
                for (i=l;i<k+1;i++) {
                    f=s*rv1[i];
                    rv1[i]=c*rv1[i];
                    if (abs(f) <= eps*anorm) break;
                    g=w[i];
                    h=pythag(f,g);
                    w[i]=h;
                    h=1.0/h;
                    c=g*h;
                    s = -f*h;
                    for (j=0;j<m;j++) {
                        y=u[j + m*nm];
                        z=u[j + m*i];
                        u[j + m*nm]=y*c+z*s;
                        u[j + m*i]=z*c-y*s;
                    }
                }
            }
            z=w[k];
            if (l == k) { //Convergence.
                if (z < 0.0) { //Singular value is made nonnegative.
                    w[k] = -z;
//                    for (j=0;j<n;j++) v[j + n*k] = -v[j + n*k];
                }
                break;
            }
//          if (its == 29) throw("no convergence in 30 svdcmp iterations");
            x=w[l]; //Shift from bottom 2-by-2 minor.
            nm=k-1;
            y=w[nm];
            g=rv1[nm];
            h=rv1[k];
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g=pythag(f,1.0);
            f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
            c=s=1.0;// Next QR transformation:
            for (j=l;j<=nm;j++) {
                i=j+1;
                g=rv1[i];
                y=w[i];
                h=s*g;
                g=c*g;
                z=pythag(f,h);
                rv1[j]=z;
                c=f/z;
                s=h/z;
                f=x*c+g*s;
                g=g*c-x*s;
                h=y*s;
                y *= c;
                z=pythag(f,h);
                w[j]=z;// Rotation can be arbitrary if z D 0.
                if (z) {
                    z=1.0/z;
                    c=f*z;
                    s=h*z;
                }
                f=c*g+s*y;
                x=c*y-s*g;
                for (jj=0;jj<m;jj++) {
                    y=u[jj + m*j];
                    z=u[jj + m*i];
                    u[jj + m*j]=y*c+z*s;
                    u[jj + m*i]=z*c-y*s;
                }
            }
            rv1[l]=0.0;
            rv1[k]=f;
            w[k]=x;
        }
    }

    delete [] rv1;
}


//=============================================================================================================
/**
* Calculates only the W Matrix of the SVD on the GPU. Attention the input Matrix is still changed.
*
* CODE of 3rd EDITION OF Numerical Recipes in C (2007)
* Given a matrix u[0..m-1][0..n-1] (u == a), this routine computes its singular value decomposition, A =
* U*W*V^T. Thematrix U replaces a on output. The diagonal matrix of singular values W is output
* as a vector w[0..n-1]. Thematrix V (not the transpose V^T ) is output as v[0..n-1][0..n-1].
*
* @param[in,out] u  The column-major array of the input matrix [m x n]. After calculation has finished
*                   this variable holds the return matrix u which contains the Eigenvectors.
* @param[in] m  The number of rows of u.
* @param[in] n  The number of columns of u.
* @param[out] w The singular values of the input matrix in a row vector [n x 1].
*/
__device__ void cuSVD_W(    float *u,   /* [m x n ]*/
                            int m,      /* rows */
                            int n,      /* columns */
                            float *w    /* [nx1]*/)
{
    //double
    //double eps = 1.19e-16; //see page 1165, Numerical Recipes in C 3.Edition (2007); // eps = numeric_limits<Doub>::epsilon();
    //float
    float eps = 1.19e-7;

    bool flag;
    int i,its,j,k,l,nm;
    float anorm,c,f,g,h,s,scale,x,y,z;
    
    float* rv1 = new float[n];//This is working since CUDA 4.0, on devices with compute capatibility 2.x

    //Part 1
    g = scale = anorm = 0.0; //Householder reduction to bidiagonal form.
    for (i=0;i<n;i++) {
        l=i+2;
        rv1[i]=scale*g;
        g=s=scale=0.0;
        if (i < m) {
            for (k=i;k<m;k++) scale += abs(u[k+i*m]);
            if (scale != 0.0) {
                for (k=i;k<m;k++) {
                    u[k+i*m] /= scale;
                    s += u[k+i*m]*u[k+i*m];
                }
                f=u[i+i*m];
                g = -SIGN(sqrt(s),f);
                h=f*g-s;
                u[i+i*m]=f-g;
                for (j=l-1;j<n;j++) {
                    for (s=0.0,k=i;k<m;k++) s += u[k+i*m]*u[k+j*m];
                    f=s/h;
                    for (k=i;k<m;k++) u[k+j*m] += f*u[k+i*m];
                }
                for (k=i;k<m;k++) u[k+i*m] *= scale;
            }
        }
        w[i]=scale *g;
        g=s=scale=0.0;
        if (i+1 <= m && i+1 != n) {
            for (k=l-1;k<n;k++) scale += abs(u[i+k*m]);
            if (scale != 0.0) {
                for (k=l-1;k<n;k++) {
                    u[i+k*m] /= scale;
                    s += u[i+k*m]*u[i+k*m];

                }
                f=u[i + m*(l-1) ];
                g = -SIGN(sqrt(s),f);
                h=f*g-s;
                u[i + m*(l-1)]=f-g;
                for (k=l-1;k<n;k++) rv1[k]=u[i + m*k]/h;
                for (j=l-1;j<m;j++) {
                    for (s=0.0,k=l-1;k<n;k++) s += u[j + m*k]*u[i + m*k];
                    for (k=l-1;k<n;k++) u[j + m*k] += s*rv1[k];
                }
                for (k=l-1;k<n;k++) u[i + m*k] *= scale;
            }
        }
        anorm = DMAX(anorm,(abs(w[i])+abs(rv1[i])));
    }
    //Part 4
    for (k=n-1;k>=0;k--) { //Diagonalization of the bidiagonal form: Loop over
        for (its=0;its<30;its++) { //singular values, and over allowed iterations.
            flag=true;
            for (l=k;l>=0;l--) {// Test for splitting.
                nm=l-1;
                if (l == 0 || abs(rv1[l]) <= eps*anorm) {
                    flag=false;
                    break;
                }
                if (abs(w[nm]) <= eps*anorm) break;
            }
            if (flag) {
                c=0.0; //Cancellation of rv1[l], if l > 0.
                    s=1.0;
                for (i=l;i<k+1;i++) {
                    f=s*rv1[i];
                    rv1[i]=c*rv1[i];
                    if (abs(f) <= eps*anorm) break;
                    g=w[i];
                    h=pythag(f,g);
                    w[i]=h;
                    h=1.0/h;
                    c=g*h;
                    s = -f*h;
                }
            }
            z=w[k];
            if (l == k) { //Convergence.
                if (z < 0.0) { //Singular value is made nonnegative.
                    w[k] = -z;
                }
                break;
            }
//          if (its == 29) throw("no convergence in 30 svdcmp iterations");
            x=w[l]; //Shift from bottom 2-by-2 minor.
            nm=k-1;
            y=w[nm];
            g=rv1[nm];
            h=rv1[k];
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g=pythag(f,1.0);
            f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
            c=s=1.0;// Next QR transformation:
            for (j=l;j<=nm;j++) {
                i=j+1;
                g=rv1[i];
                y=w[i];
                h=s*g;
                g=c*g;
                z=pythag(f,h);
                rv1[j]=z;
                c=f/z;
                s=h/z;
                f=x*c+g*s;
                g=g*c-x*s;
                h=y*s;
                y *= c;
                z=pythag(f,h);
                w[j]=z;// Rotation can be arbitrary if z D 0.
                if (z) {
                    z=1.0/z;
                    c=f*z;
                    s=h*z;
                }
                f=c*g+s*y;
                x=c*y-s*g;
            }
            rv1[l]=0.0;
            rv1[k]=f;
            w[k]=x;
        }
    }

    delete [] rv1;
}


//=============================================================================================================
/**


* This function provides a higly parallelized device SVD.
* The parallelization is done with the help of y and z threads. Because reductions are performed y & z must be
* a power of 2, respectively.
* This function calculates the whole SVD on the GPU.
*
* Given a matrix u[0..m-1][0..n-1] (u == a), this routine computes its singular value decomposition, A =
* U*W*V^T. The matrix U replaces a on output. The diagonal matrix of singular values W is output
* as a vector w[0..n-1]. The matrix V (not the transpose V^T ) is output as v[0..n-1][0..n-1].
*
* @param[in,out] u  The column-major array of the input matrix [m x n]. After calculation has finished
*                   this variable holds the return matrix u which contains the Eigenvectors.
* @param[in] m  The number of rows of u.
* @param[in] n  The number of columns of u.
* @param[out] w  The singular values of the input matrix in a row vector [n x 1].
* @param[out] v  The matrix v [n x n] of the singular value decomposition of the input matrix.
* @param shared_cache   Intern working array of size [nx1] in the shared memory for intern storage.
* @param t_pCacheYZ     Intern working array in the shared memory of size [blockDim.yxblockDim.z].
*/
__device__ void cuSVD_shared(   float *u,   /* [m x n ]*/
                                int m,      /* rows */
                                int n,      /* columns */
                                float *w,   /* [nx1]*/
                                float *v,   /* [nxn]*/
                                float *shared_cache, /* [nx1] */
                                float *t_pCacheYZ)
{
    //double
    //double eps = 1.19e-16; //see page 1165, Numerical Recipes in C 3.Edition (2007); // eps = numeric_limits<Doub>::epsilon();
    //float
    float eps = 1.19e-7;

    bool flag;
    int i,its,j,jj,k,l,nm;
    float anorm,c,f,g,h,s,x,y,z;

    float* rv1 = shared_cache;//size = n
    float* t_pScale = rv1+n;
    float* t_pS = t_pScale + 1;

//    float* t_pCacheY = t_pCacheYZ + threadIdx.z*blockDim.y;
    float* t_pCacheZ = t_pCacheYZ + threadIdx.y*blockDim.z;

    //Part 1
    g = *t_pScale = anorm = 0.0; //Householder reduction to bidiagonal form.
    for (i=0;i<n;i++) {
        l=i+2;
        rv1[i]=*t_pScale*g;
        g=*t_pS=*t_pScale=0.0;
        if (i < m) {

       /*9*/cuReduceAbsYZ_shared( u+(i+i*m), m-i-1, t_pCacheYZ, t_pScale);
            //for (k=i;k<m;k++) *t_pScale += abs(u[k+i*m]);
            __syncthreads();

            if (*t_pScale != 0.0) {

           /*8*/cuVecScaleSquareReduceYZ_shared( t_pScale, u+i+i*m, m-i, t_pCacheYZ, t_pS );
                __syncthreads();
                //for (k=i;k<m;k++) {//ToDo parallel
                //    u[k+i*m] /= *t_pScale;
                //    s += u[k+i*m]*u[k+i*m];
                //}

                f=u[i+i*m];
                g = -SIGN(sqrt(*t_pS),f);
                h=f*g-*t_pS;
                u[i+i*m]=f-g;

                for (j=l-1;j<n;j++) {//ToDo parallel - nonsens cause n is too small
              /*10*/cuScalarProductYZ_shared( u+(i+i*m), u+(i+j*m), m-i, t_pCacheYZ, t_pS );
                    __syncthreads();
                    //for (s=0.0,k=i;k<m;k++) s += u[k+i*m]*u[k+j*m];

                    f=*t_pS/h;

               /*1*/k = i + threadIdx.y+blockDim.y*threadIdx.z;
                    while(k < m)
                    {
                        u[k+j*m] += f*u[k+i*m];
                        k += blockDim.y*blockDim.z;
                    }//for (k=i;k<m;k++) u[k+j*m] += f*u[k+i*m];//same like /*7*/
                    __syncthreads();
                }

           /*2*/k = i + threadIdx.y+blockDim.y*threadIdx.z;
                while(k < m)
                {
                    u[k+i*m] *= *t_pScale;
                    k += blockDim.y*blockDim.z;
                }//for (k=i;k<m;k++) u[k+i*m] *= *t_pScale;//ToDo parallel
                __syncthreads();
            }
        }
        w[i]=*t_pScale *g;
        g=*t_pS=*t_pScale=0.0;
        if (i+1 <= m && i+1 != n) {
            
      /*11*/cuReduceAbsYZ_shared( u+(i+(l-1)*m), n-(l-1), t_pCacheYZ, t_pScale);
            //for (k=l-1;k<n;k++) *t_pScale += abs(u[i+k*m]);
            __syncthreads();

            if (*t_pScale != 0.0) {

          /*12*/cuVecScaleSquareReduceYZ_incr_shared( t_pScale, u+i+(l-1)*m, n-(l-1), t_pCacheYZ, t_pS, m);
                __syncthreads();
                //for (k=l-1;k<n;k++) {
                //    u[i+k*m] /= *t_pScale;
                //    s += u[i+k*m]*u[i+k*m];
                //}

                f=u[i + m*(l-1) ];
                g = -SIGN(sqrt(*t_pS),f);
                h=f*g-*t_pS;
                u[i + m*(l-1)]=f-g;

          /*13*/k = l-1 + threadIdx.y+blockDim.y*threadIdx.z;
                while(k < n)
                {
                    rv1[k]=u[i + m*k]/h;
                    k += blockDim.y*blockDim.z;
                }//for (k=l-1;k<n;k++) rv1[k]=u[i + m*k]/h;
                __syncthreads();

           /*3*/j = l-1 + threadIdx.y;
                while(j < m)
                {
              /*15*/cuScalarProductAddZ_incr_shared(u+j+m*(l-1), u+i+m*(l-1), n-(l-1), t_pCacheZ, t_pS, m, rv1+(l-1));
                    //for (s=0.0,k=l-1;k<n;k++) s += u[j + m*k]*u[i + m*k];
                    //for (k=l-1;k<n;k++) u[j + m*k] += s*rv1[k];
                    __syncthreads();

                    j += blockDim.y;
                }//for (j=l-1;j<m;j++) {
                 //    for (s=0.0,k=l-1;k<n;k++) s += u[j + m*k]*u[i + m*k];
                 //    for (k=l-1;k<n;k++) u[j + m*k] += s*rv1[k];
                 //}
                __syncthreads();

          /*14*/k = l-1 + threadIdx.y+blockDim.y*threadIdx.z;
                while(k < n)
                {
                    u[i + m*k] *= *t_pScale;
                    k += blockDim.y*blockDim.z;
                }//for (k=l-1;k<n;k++) u[i + m*k] *= *t_pScale;
            }
        }
        anorm = DMAX(anorm,(abs(w[i])+abs(rv1[i])));
    }
    //Part 2
    for (i=n-1;i>=0;i--) {// Accumulation of right-hand transformations.
        if (i < n-1) {
            if (g != 0.0) {
                for (j=l;j<n;j++) //Double division to avoid possible underow.
                    v[j + n*i]=(u[i + m*j]/u[i + m*l])/g;
                for (j=l;j<n;j++) {
                    for (s=0.0,k=l;k<n;k++) s += u[i + m*k]*v[k + n*j];
                    for (k=l;k<n;k++) v[k + n*j] += s*v[k + n*i];
                }
            }
            for (j=l;j<n;j++) v[i + n*j]=v[j + n*i]=0.0;
        }
        v[i + n*i]=1.0;
        g=rv1[i];
        l=i;
    }
    //Part 3
    for (i=IMIN(m,n)-1;i>=0;i--) { //Accumulation of left-hand transformations.
        l=i+1;
        g=w[i];
        for (j=l;j<n;j++) u[i + m*j]=0.0;
        if (g != 0.0) {
            g=1.0/g;
            for (j=l;j<n;j++) {

          /*17*/cuScalarProductYZ_shared( u+(l+m*i), u+(l+m*j), m-l, t_pCacheYZ, t_pS );
                //for (s=0.0,k=l;k<m;k++) s += u[k + m*i]*u[k + m*j];//scalarprod y
                __syncthreads();

                f=(*t_pS/u[i + m*i])*g;

           /*7*/k = i + threadIdx.y+blockDim.y*threadIdx.z;
                while(k < m)
                {
                    u[k + m*j] += f*u[k + m*i];
                    k += blockDim.y*blockDim.z;
                }//for (k=i;k<m;k++) u[k + m*j] += f*u[k + m*i];//same like /*1*/
                __syncthreads();
            }
       /*4*/j = i + threadIdx.y+blockDim.y*threadIdx.z;
            while(j < m)
            {
                u[j + m*i] *= g;
                j += blockDim.y*blockDim.z;
            }//for (j=i;j<m;j++) u[j + m*i] *= g;
            __syncthreads();
        }
        else {
       /*4*/j = i + threadIdx.y+blockDim.y*threadIdx.z;
            while(j < m)
            {
                u[j + m*i]=0.0;
                j += blockDim.y*blockDim.z;
            }//for (j=i;j<m;j++) u[j + m*i]=0.0;
            __syncthreads();
        }
        ++u[i + m*i];
    }
    //Part 4
    for (k=n-1;k>=0;k--) { //Diagonalization of the bidiagonal form: Loop over
        for (its=0;its<30;its++) { //singular values, and over allowed iterations.
            flag=true;
            for (l=k;l>=0;l--) {// Test for splitting.
                nm=l-1;
                if (l == 0 || abs(rv1[l]) <= eps*anorm) {
                    flag=false;
                    break;
                }
                if (abs(w[nm]) <= eps*anorm) break;
            }
            if (flag) {
                c=0.0; //Cancellation of rv1[l], if l > 0.
                    s=1.0;

          /*16*/i = l + threadIdx.z;
                while(i < k+1)
                {
//                for (i=l;i<k+1;i++) {
                    f=s*rv1[i];
                    rv1[i]=c*rv1[i];
                    if (abs(f) <= eps*anorm) break;
                    g=w[i];
                    h=pythag(f,g);
                    w[i]=h;
                    h=1.0/h;
                    c=g*h;
                    s = -f*h;

               /*5*/j = threadIdx.y;
                    while(j < m)
                    {
                        y=u[j + m*nm];
                        z=u[j + m*i];
                        u[j + m*nm]=y*c+z*s;
                        u[j + m*i]=z*c-y*s;
                        j += blockDim.y;
                    }//for (j=0;j<m;j++) {
                     //    y=u[j + m*nm];
                     //    z=u[j + m*i];
                     //    u[j + m*nm]=y*c+z*s;
                     //    u[j + m*i]=z*c-y*s;
                     //}
                    __syncthreads();

                    i += blockDim.z;
                }
            }
            z=w[k];
            if (l == k) { //Convergence.
                if (z < 0.0) { //Singular value is made nonnegative.
                    w[k] = -z;
                    for (j=0;j<n;j++) v[j + n*k] = -v[j + n*k];
                }
                break;
            }
//          if (its == 29) throw("no convergence in 30 svdcmp iterations");
            x=w[l]; //Shift from bottom 2-by-2 minor.
            nm=k-1;
            y=w[nm];
            g=rv1[nm];
            h=rv1[k];
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g=pythag(f,1.0);
            f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
            c=s=1.0;// Next QR transformation:
            for (j=l;j<=nm;j++) {
                i=j+1;
                g=rv1[i];
                y=w[i];
                h=s*g;
                g=c*g;
                z=pythag(f,h);
                rv1[j]=z;
                c=f/z;
                s=h/z;
                f=x*c+g*s;
                g=g*c-x*s;
                h=y*s;
                y *= c;
          /*6b*/jj = threadIdx.y+blockDim.y*threadIdx.z;
                while(jj < m)
                {
                    x=v[jj + n*j];
                    z=v[jj + n*i];
                    v[jj + n*j]=x*c+z*s;
                    v[jj + n*i]=z*c-x*s;
                    jj += blockDim.y*blockDim.z;
                }//for (jj=0;jj<n;jj++) {
                 //    x=v[jj + n*j];
                 //    z=v[jj + n*i];
                 //    v[jj + n*j]=x*c+z*s;
                 //    v[jj + n*i]=z*c-x*s;
                 //}
                __syncthreads();
                z=pythag(f,h);
                w[j]=z;// Rotation can be arbitrary if z D 0.
                if (z) {
                    z=1.0/z;
                    c=f*z;
                    s=h*z;
                }
                f=c*g+s*y;
                x=c*y-s*g;
           /*6*/jj = threadIdx.y+blockDim.y*threadIdx.z;
                while(jj < m)
                {
                    y=u[jj + m*j];
                    z=u[jj + m*i];
                    u[jj + m*j]=y*c+z*s;
                    u[jj + m*i]=z*c-y*s;
                    jj += blockDim.y*blockDim.z;
                }//for (jj=0;jj<m;jj++) {
                 //    y=u[jj + m*j];
                 //    z=u[jj + m*i];
                 //    u[jj + m*j]=y*c+z*s;
                 //    u[jj + m*i]=z*c-y*s;
                 //}
                __syncthreads();
            }
            rv1[l]=0.0;
            rv1[k]=f;
            w[k]=x;
        }
    }
}


//=============================================================================================================
/**
* This function provides a higly parallelized device SVD.
* The parallelization is done with the help of y and z threads. Because reductions are performed y & z must be
* a power of 2, respectively.
* This function calculates only the U an W Matrix of the SVD on the GPU.
*
* Given a matrix u[0..m-1][0..n-1] (u == a), this routine computes its singular value decomposition, A =
* U*W*V^T. The matrix U replaces a on output. The diagonal matrix of singular values W is output
* as a vector w[0..n-1].
*
* @param[in,out] u  The column-major array of the input matrix [m x n]. After calculation has finished
*                   this variable holds the return matrix u which contains the Eigenvectors.
* @param[in] m  The number of rows of u.
* @param[in] n  The number of columns of u.
* @param[out] w  The singular values of the input matrix in a row vector [n x 1].
* @param shared_cache   Intern working array of size [nx1] in the shared memory for intern storage.
* @param t_pCacheYZ     Intern working array in the shared memory of size [blockDim.yxblockDim.z].
*/
__device__ void cuSVD_UW_shared(    float *u,   /* [m x n ]*/
                                    int m,      /* rows */
                                    int n,      /* columns */
                                    float *w,   /* [nx1]*/
                                    float *shared_cache, /* [nx1] */
                                    float *t_pCacheYZ)
{
    //double
    //double eps = 1.19e-16; //see page 1165, Numerical Recipes in C 3.Edition (2007); // eps = numeric_limits<Doub>::epsilon();
    //float
    float eps = 1.19e-7;

    bool flag;
    int i,its,j,jj,k,l,nm;
    float anorm,c,f,g,h,s,x,y,z;
    
    float* rv1 = shared_cache;//size = n
    float* t_pScale = rv1 + n;
    float* t_pS = t_pScale + 1;

//    float* t_pCacheY = t_pCacheYZ + threadIdx.z*blockDim.y;
    float* t_pCacheZ = t_pCacheYZ + threadIdx.y*blockDim.z;

    //Part 1
    g = *t_pScale = anorm = 0.0; //Householder reduction to bidiagonal form.
    for (i=0;i<n;i++) {
        l=i+2;
        rv1[i]=*t_pScale*g;
        g=*t_pS=*t_pScale=0.0;
        if (i < m) {

       /*9*/cuReduceAbsYZ_shared( u+(i+i*m), m-i-1, t_pCacheYZ, t_pScale);
            //for (k=i;k<m;k++) *t_pScale += abs(u[k+i*m]);
            __syncthreads();

            if (*t_pScale != 0.0) {

           /*8*/cuVecScaleSquareReduceYZ_shared( t_pScale, u+i+i*m, m-i, t_pCacheYZ, t_pS );
                __syncthreads();
                //for (k=i;k<m;k++) {
                //    u[k+i*m] /= *t_pScale;
                //    s += u[k+i*m]*u[k+i*m];
                //}

                f=u[i+i*m];
                g = -SIGN(sqrt(*t_pS),f);
                h=f*g-*t_pS;
                u[i+i*m]=f-g;

                for (j=l-1;j<n;j++) {//ToDo parallel with z
              /*10*/cuScalarProductYZ_shared( u+(i+i*m), u+(i+j*m), m-i, t_pCacheYZ, t_pS );
                    __syncthreads();
                    //for (s=0.0,k=i;k<m;k++) s += u[k+i*m]*u[k+j*m];

                    f=*t_pS/h;

               /*1*/k = i + threadIdx.y+blockDim.y*threadIdx.z;
                    while(k < m)
                    {
                        u[k+j*m] += f*u[k+i*m];
                        k += blockDim.y*blockDim.z;
                    }//for (k=i;k<m;k++) u[k+j*m] += f*u[k+i*m];
                    __syncthreads();
                }

           /*2*/k = i + threadIdx.y+blockDim.y*threadIdx.z;
                while(k < m)
                {
                    u[k+i*m] *= *t_pScale;
                    k += blockDim.y*blockDim.z;
                }//for (k=i;k<m;k++) u[k+i*m] *= *t_pScale;
                __syncthreads();
            }
        }
        w[i]=*t_pScale *g;
        g=*t_pS=*t_pScale=0.0;
        if (i+1 <= m && i+1 != n) {

      /*11*/cuReduceAbsYZ_shared( u+(i+(l-1)*m), n-(l-1), t_pCacheYZ, t_pScale);
            //for (k=l-1;k<n;k++) *t_pScale += abs(u[i+k*m]);
            __syncthreads();

            if (*t_pScale != 0.0) {

          /*12*/cuVecScaleSquareReduceYZ_incr_shared( t_pScale, u+i+(l-1)*m, n-(l-1), t_pCacheYZ, t_pS, m);
                __syncthreads();
                //for (k=l-1;k<n;k++) {
                //    u[i+k*m] /= *t_pScale;
                //    s += u[i+k*m]*u[i+k*m];
                //}

                f=u[i + m*(l-1) ];
                g = -SIGN(sqrt(*t_pS),f);
                h=f*g-*t_pS;
                u[i + m*(l-1)]=f-g;

          /*13*/k = l-1 + threadIdx.y+blockDim.y*threadIdx.z;
                while(k < n)
                {
                    rv1[k]=u[i + m*k]/h;
                    k += blockDim.y*blockDim.z;
                }//for (k=l-1;k<n;k++) rv1[k]=u[i + m*k]/h;
                __syncthreads();

           /*3*/j = l-1 + threadIdx.y;
                while(j < m)
                {
              /*15*/cuScalarProductAddZ_incr_shared(u+j+m*(l-1), u+i+m*(l-1), n-(l-1), t_pCacheZ, t_pS, m, rv1+(l-1));
                    //for (s=0.0,k=l-1;k<n;k++) s += u[j + m*k]*u[i + m*k];
                    //for (k=l-1;k<n;k++) u[j + m*k] += s*rv1[k];
                    __syncthreads();

                    j += blockDim.y;
                }//for (j=l-1;j<m;j++) {
                 //    for (s=0.0,k=l-1;k<n;k++) s += u[j + m*k]*u[i + m*k];
                 //    for (k=l-1;k<n;k++) u[j + m*k] += s*rv1[k];
                 //}
                __syncthreads();

          /*14*/k = l-1 + threadIdx.y+blockDim.y*threadIdx.z;
                while(k < n)
                {
                    u[i + m*k] *= *t_pScale;
                    k += blockDim.y*blockDim.z;
                }//for (k=l-1;k<n;k++) u[i + m*k] *= *t_pScale;
                __syncthreads();
            }
        }
        anorm = DMAX(anorm,(abs(w[i])+abs(rv1[i])));
    }
    //Part 3
    for (i=IMIN(m,n)-1;i>=0;i--) { //Accumulation of left-hand transformations.
        l=i+1;
        g=w[i];
        for (j=l;j<n;j++) u[i + m*j]=0.0;
        if (g != 0.0) {
            g=1.0/g;
            for (j=l;j<n;j++) {

          /*17*/cuScalarProductYZ_shared( u+(l+m*i), u+(l+m*j), m-l, t_pCacheYZ, t_pS );
                //for (s=0.0,k=l;k<m;k++) s += u[k + m*i]*u[k + m*j];//scalarprod y
                __syncthreads();

                f=(*t_pS/u[i + m*i])*g;

           /*7*/k = i + threadIdx.y+blockDim.y*threadIdx.z;
                while(k < m)
                {
                    u[k + m*j] += f*u[k + m*i];
                    k += blockDim.y*blockDim.z;
                }//for (k=i;k<m;k++) u[k + m*j] += f*u[k + m*i];//same like /*1*/
                __syncthreads();
            }
       /*4*/j = i + threadIdx.y+blockDim.y*threadIdx.z;
            while(j < m)
            {
                u[j + m*i] *= g;
                j += blockDim.y*blockDim.z;
            }//for (j=i;j<m;j++) u[j + m*i] *= g;
            __syncthreads();
        }
        else {
       /*4*/j = i + threadIdx.y+blockDim.y*threadIdx.z;
            while(j < m)
            {
                u[j + m*i]=0.0;
                j += blockDim.y*blockDim.z;
            }//for (j=i;j<m;j++) u[j + m*i]=0.0;
            __syncthreads();
        }
        ++u[i + m*i];
    }
    //Part 4
    for (k=n-1;k>=0;k--) { //Diagonalization of the bidiagonal form: Loop over
        for (its=0;its<30;its++) { //singular values, and over allowed iterations.
            flag=true;
            for (l=k;l>=0;l--) {// Test for splitting.
                nm=l-1;
                if (l == 0 || abs(rv1[l]) <= eps*anorm) {
                    flag=false;
                    break;
                }
                if (abs(w[nm]) <= eps*anorm) break;
            }
            if (flag) {
                c=0.0; //Cancellation of rv1[l], if l > 0.
                    s=1.0;

          /*16*/i = l + threadIdx.z;
                while(i < k+1)
                {
//                for (i=l;i<k+1;i++) {
                    f=s*rv1[i];
                    rv1[i]=c*rv1[i];
                    if (abs(f) <= eps*anorm) break;
                    g=w[i];
                    h=pythag(f,g);
                    w[i]=h;
                    h=1.0/h;
                    c=g*h;
                    s = -f*h;

               /*5*/j = threadIdx.y;
                    while(j < m)
                    {
                        y=u[j + m*nm];
                        z=u[j + m*i];
                        u[j + m*nm]=y*c+z*s;
                        u[j + m*i]=z*c-y*s;
                        j += blockDim.y;
                    }//for (j=0;j<m;j++) {
                     //    y=u[j + m*nm];
                     //    z=u[j + m*i];
                     //    u[j + m*nm]=y*c+z*s;
                     //    u[j + m*i]=z*c-y*s;
                     //}
                    __syncthreads();

                    i += blockDim.z;
                }
            }
            z=w[k];
            if (l == k) { //Convergence.
                if (z < 0.0) { //Singular value is made nonnegative.
                    w[k] = -z;
//                    for (j=0;j<n;j++) v[j + n*k] = -v[j + n*k];
                }
                break;
            }
//          if (its == 29) throw("no convergence in 30 svdcmp iterations");
            x=w[l]; //Shift from bottom 2-by-2 minor.
            nm=k-1;
            y=w[nm];
            g=rv1[nm];
            h=rv1[k];
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g=pythag(f,1.0);
            f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
            c=s=1.0;// Next QR transformation:
            for (j=l;j<=nm;j++) {
                i=j+1;
                g=rv1[i];
                y=w[i];
                h=s*g;
                g=c*g;
                z=pythag(f,h);
                rv1[j]=z;
                c=f/z;
                s=h/z;
                f=x*c+g*s;
                g=g*c-x*s;
                h=y*s;
                y *= c;
                z=pythag(f,h);
                w[j]=z;// Rotation can be arbitrary if z D 0.
                if (z) {
                    z=1.0/z;
                    c=f*z;
                    s=h*z;
                }
                f=c*g+s*y;
                x=c*y-s*g;
           /*6*/jj = threadIdx.y+blockDim.y*threadIdx.z;
                while(jj < m)
                {
                    y=u[jj + m*j];
                    z=u[jj + m*i];
                    u[jj + m*j]=y*c+z*s;
                    u[jj + m*i]=z*c-y*s;
                    jj += blockDim.y*blockDim.z;
                }//for (jj=0;jj<m;jj++) {
                 //    y=u[jj + m*j];
                 //    z=u[jj + m*i];
                 //    u[jj + m*j]=y*c+z*s;
                 //    u[jj + m*i]=z*c-y*s;
                 //}
                __syncthreads();
            }
            rv1[l]=0.0;
            rv1[k]=f;
            w[k]=x;
        }
    }
}


//=============================================================================================================
/**
* This function provides a higly parallelized device SVD.
* The parallelization is done with the help of y and z threads. Because reductions are performed y & z must be
* a power of 2, respectively.
* This function calculates only the W Matrix of the SVD on the GPU. !Attention the input Matrix is still changed!
*
* Given a matrix u[0..m-1][0..n-1] (u == a), this routine computes its singular value decomposition, A =
* U*W*V^T. The diagonal matrix of singular values W is output as a vector w[0..n-1].
*
* @param[in,out] u  The column-major array of the input matrix [m x n]. After calculation has finished
*                   this variable holds the return matrix u which contains the Eigenvectors.
* @param[in] m  The number of rows of u.
* @param[in] n  The number of columns of u.
* @param[out] w The singular values of the input matrix in a row vector [n x 1].
* @param shared_cache   Intern working array in the shared memory of size [nx1] for intern storage.
* @param t_pCacheYZ     Intern working array in the shared memory of size [blockDim.yxblockDim.z].
*/
__device__ void cuSVD_W_shared( float *u,   /* [m x n ]*/
                                int m,      /* rows */
                                int n,      /* columns */
                                float *w,   /* [nx1]*/
                                float *shared_cache, /* [nx1] */
                                float *t_pCacheYZ)
{
    //double
    //double eps = 1.19e-16; //see page 1165, Numerical Recipes in C 3.Edition (2007); // eps = numeric_limits<Doub>::epsilon();
    //float
    float eps = 1.19e-7;

    bool flag;
    int i,its,j,k,l,nm;
    float anorm,c,f,g,h,s,x,y,z;
    
    float* rv1 = shared_cache;//size = n
    float* t_pScale = rv1+n;
    float* t_pS = t_pScale + 1;


//    float* t_pCacheY = t_pCacheYZ + threadIdx.z*blockDim.y;
    float* t_pCacheZ = t_pCacheYZ + threadIdx.y*blockDim.z;

    //Part 1
    g = *t_pScale = anorm = 0.0; //Householder reduction to bidiagonal form.
    for (i=0;i<n;i++) {
        l=i+2;
        rv1[i]=*t_pScale*g;
        g=*t_pS=*t_pScale=0.0;
        if (i < m) {

       /*9*/cuReduceAbsYZ_shared( u+(i+i*m), m-i-1, t_pCacheYZ, t_pScale);
            //for (k=i;k<m;k++) *t_pScale += abs(u[k+i*m]);
            __syncthreads();

            if (*t_pScale != 0.0) {

           /*8*/cuVecScaleSquareReduceYZ_shared( t_pScale, u+i+i*m, m-i, t_pCacheYZ, t_pS );
                __syncthreads();
                //for (k=i;k<m;k++) {
                //    u[k+i*m] /= *t_pScale;
                //    s += u[k+i*m]*u[k+i*m];
                //}

                f=u[i+i*m];
                g = -SIGN(sqrt(*t_pS),f);
                h=f*g-*t_pS;
                u[i+i*m]=f-g;

                for (j=l-1;j<n;j++) {
              /*10*/cuScalarProductYZ_shared( u+(i+i*m), u+(i+j*m), m-i, t_pCacheYZ, t_pS );
                    __syncthreads();
                    //for (s=0.0,k=i;k<m;k++) s += u[k+i*m]*u[k+j*m];

                    f=*t_pS/h;

               /*1*/k = i + threadIdx.y+blockDim.y*threadIdx.z;
                    while(k < m)
                    {
                        u[k+j*m] += f*u[k+i*m];
                        k += blockDim.y*blockDim.z;
                    }//for (k=i;k<m;k++) u[k+j*m] += f*u[k+i*m];
                    __syncthreads();
                }

           /*2*/k = i + threadIdx.y+blockDim.y*threadIdx.z;
                while(k < m)
                {
                    u[k+i*m] *= *t_pScale;
                    k += blockDim.y*blockDim.z;
                }//for (k=i;k<m;k++) u[k+i*m] *= *t_pScale;
                __syncthreads();
            }
        }
        w[i]=*t_pScale *g;
        g=*t_pS=*t_pScale=0.0;
        if (i+1 <= m && i+1 != n) {

      /*11*/cuReduceAbsYZ_shared( u+(i+(l-1)*m), n-(l-1), t_pCacheYZ, t_pScale);
            //for (k=l-1;k<n;k++) *t_pScale += abs(u[i+k*m]);
            __syncthreads();

            if (*t_pScale != 0.0) {
          /*12*/cuVecScaleSquareReduceYZ_incr_shared( t_pScale, u+i+(l-1)*m, n-(l-1), t_pCacheYZ, t_pS, m);
                __syncthreads();
                //for (k=l-1;k<n;k++) {
                //    u[i+k*m] /= *t_pScale;
                //    s += u[i+k*m]*u[i+k*m];
                //}

                f=u[i + m*(l-1) ];
                g = -SIGN(sqrt(*t_pS),f);
                h=f*g-*t_pS;
                u[i + m*(l-1)]=f-g;

          /*13*/k = l-1 + threadIdx.y+blockDim.y*threadIdx.z;
                while(k < n)
                {
                    rv1[k]=u[i + m*k]/h;
                    k += blockDim.y*blockDim.z;
                }//for (k=l-1;k<n;k++) rv1[k]=u[i + m*k]/h;

           /*3*/j = l-1 + threadIdx.y;
                while(j < m)//ToDo put this into a device function
                {
              /*15*/cuScalarProductAddZ_incr_shared(u+j+m*(l-1), u+i+m*(l-1), n-(l-1), t_pCacheZ, t_pS, m, rv1+(l-1));
                    //for (s=0.0,k=l-1;k<n;k++) s += u[j + m*k]*u[i + m*k];
                    //for (k=l-1;k<n;k++) u[j + m*k] += s*rv1[k];
                    __syncthreads();

                    j += blockDim.y;
                }//for (j=l-1;j<m;j++) {
                 //    for (s=0.0,k=l-1;k<n;k++) s += u[j + m*k]*u[i + m*k];
                 //    for (k=l-1;k<n;k++) u[j + m*k] += s*rv1[k];
                 //}
                __syncthreads();

          /*14*/k = l-1 + threadIdx.y+blockDim.y*threadIdx.z;
                while(k < n)
                {
                    u[i + m*k] *= *t_pScale;
                    k += blockDim.y*blockDim.z;
                }//for (k=l-1;k<n;k++) u[i + m*k] *= *t_pScale;
                __syncthreads();
            }
        }
        anorm = DMAX(anorm,(abs(w[i])+abs(rv1[i])));
    }
    //Part 4
    for (k=n-1;k>=0;k--) { //Diagonalization of the bidiagonal form: Loop over
        for (its=0;its<30;its++) { //singular values, and over allowed iterations.
            flag=true;
            for (l=k;l>=0;l--) {// Test for splitting.
                nm=l-1;
                if (l == 0 || abs(rv1[l]) <= eps*anorm) {
                    flag=false;
                    break;
                }
                if (abs(w[nm]) <= eps*anorm) break;
            }
            if (flag) {
                c=0.0; //Cancellation of rv1[l], if l > 0.
                    s=1.0;

          /*16*/i = l + threadIdx.y+blockDim.y*threadIdx.z; //Different than other SVDs
                while(i < k+1)
                {
//                for (i=l;i<k+1;i++) {
                    f=s*rv1[i];
                    rv1[i]=c*rv1[i];
                    if (abs(f) <= eps*anorm) break;
                    g=w[i];
                    h=pythag(f,g);
                    w[i]=h;
                    h=1.0/h;
                    c=g*h;
                    s = -f*h;

                    i += blockDim.y*blockDim.z;
                }
            }
            z=w[k];
            if (l == k) { //Convergence.
                if (z < 0.0) { //Singular value is made nonnegative.
                    w[k] = -z;
                }
                break;
            }
//          if (its == 29) throw("no convergence in 30 svdcmp iterations");
            x=w[l]; //Shift from bottom 2-by-2 minor.
            nm=k-1;
            y=w[nm];
            g=rv1[nm];
            h=rv1[k];
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g=pythag(f,1.0);
            f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
            c=s=1.0;// Next QR transformation:
            for (j=l;j<=nm;j++) {
                i=j+1;
                g=rv1[i];
                y=w[i];
                h=s*g;
                g=c*g;
                z=pythag(f,h);
                rv1[j]=z;
                c=f/z;
                s=h/z;
                f=x*c+g*s;
                g=g*c-x*s;
                h=y*s;
                y *= c;
                z=pythag(f,h);
                w[j]=z;// Rotation can be arbitrary if z D 0.
                if (z) {
                    z=1.0/z;
                    c=f*z;
                    s=h*z;
                }
                f=c*g+s*y;
                x=c*y-s*g;
            }
            rv1[l]=0.0;
            rv1[k]=f;
            w[k]=x;
        }
    }
}

}//Namespace

#endif
