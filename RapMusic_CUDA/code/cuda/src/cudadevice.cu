//=============================================================================================================
/**
* @file     cudadevice.cu
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

#include "../include/cudadevice.cuh"


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

CudaDevice::CudaDevice()
: m_iDeviceCount(0)
, m_pCudaDeviceProp(NULL)
{

    cudaError_t err = cudaGetDeviceCount( &m_iDeviceCount );

    if(err)
        std::cout << std::endl << "No CUDA Device detected" << std::endl << std::endl;


    if(m_iDeviceCount > 0)
    {
        m_pCudaDeviceProp = new cudaDeviceProp[m_iDeviceCount];

        for (int i = 0; i < m_iDeviceCount; i++)
            HANDLE_ERROR( cudaGetDeviceProperties( &m_pCudaDeviceProp[i], i ) );

        m_iSelectedDevice = 0;
    }
}


//*************************************************************************************************************

CudaDevice::~CudaDevice()
{
    if(m_pCudaDeviceProp != NULL) 
        delete m_pCudaDeviceProp;
}


//*************************************************************************************************************

bool CudaDevice::cudaDeviceAvailable()
{
    if(m_iDeviceCount > 0)
        return true;

    return false;
}


//*************************************************************************************************************

bool CudaDevice::cuda2xDeviceAvailable()
{
    for (int i = 0; i< m_iDeviceCount; i++)
        if(m_pCudaDeviceProp[i].major >= 2)
            return true;

    return false;
}


//*************************************************************************************************************

void CudaDevice::printDeviceProperties()
{
    if(m_iDeviceCount > 0)
    {
        for (int i = 0; i< m_iDeviceCount; i++) {
            printf( "   --- General Information for device %d ---\n", i );
            printf( "Name:  %s\n", m_pCudaDeviceProp[i].name );
            printf( "Compute capability:  %d.%d\n", m_pCudaDeviceProp[i].major, m_pCudaDeviceProp[i].minor );
            printf( "Clock rate:  %d\n", m_pCudaDeviceProp[i].clockRate );
            printf( "Device copy overlap:  " );
            if (m_pCudaDeviceProp[i].deviceOverlap)
                printf( "Enabled\n" );
            else
                printf( "Disabled\n");
            printf( "Device can map host memory:  " );
            if (m_pCudaDeviceProp[i].canMapHostMemory)
                printf( "Enabled\n" );
            else
                printf( "Disabled\n");
            printf( "Kernel execution timeout :  " );
            if (m_pCudaDeviceProp[i].kernelExecTimeoutEnabled)
                printf( "Enabled\n" );
            else
                printf( "Disabled\n" );

            printf( "   --- Memory Information for device %d ---\n", i );
            printf( "Total global mem:  %ld\n", m_pCudaDeviceProp[i].totalGlobalMem );
            printf( "Total constant Mem:  %ld\n", m_pCudaDeviceProp[i].totalConstMem );
            printf( "Max mem pitch:  %ld\n", m_pCudaDeviceProp[i].memPitch );
            printf( "Texture Alignment:  %ld\n", m_pCudaDeviceProp[i].textureAlignment );

            printf( "   --- MP Information for device %d ---\n", i );
            printf( "Multiprocessor count:  %d\n",
                        m_pCudaDeviceProp[i].multiProcessorCount );
            printf( "Shared mem per mp:  %ld\n", m_pCudaDeviceProp[i].sharedMemPerBlock );
            printf( "Registers per mp:  %d\n", m_pCudaDeviceProp[i].regsPerBlock );
            printf( "Threads in warp:  %d\n", m_pCudaDeviceProp[i].warpSize );
            printf( "Max threads per block:  %d\n",
                        m_pCudaDeviceProp[i].maxThreadsPerBlock );
            printf( "Max thread dimensions:  (%d, %d, %d)\n",
                        m_pCudaDeviceProp[i].maxThreadsDim[0], m_pCudaDeviceProp[i].maxThreadsDim[1],
                        m_pCudaDeviceProp[i].maxThreadsDim[2] );
            printf( "Max grid dimensions:  (%d, %d, %d)\n",
                        m_pCudaDeviceProp[i].maxGridSize[0], m_pCudaDeviceProp[i].maxGridSize[1],
                        m_pCudaDeviceProp[i].maxGridSize[2] );
            printf( "\n" );
        }

    }
    else
    {
        printf( "No Cuda Device available.\n");
    }
}


}//Namespace