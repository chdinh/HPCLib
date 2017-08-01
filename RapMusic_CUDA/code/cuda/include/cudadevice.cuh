//=============================================================================================================
/**
* @file		cudadevice.cuh
* @author	Christoph Dinh <christoph.dinh@live.de>;
* @version	1.0
* @date		March, 2011
*
* @section	LICENSE
*
* Copyright (C) 2011 Christoph Dinh. All rights reserved.
*
* No part of this program may be photocopied, reproduced,
* or translated to another program language without the
* prior written consent of the author.
*
*
* @brief	ToDo Documentation...
*
*/


#ifndef CUDADEVICE_CUH
#define CUDADEVICE_CUH


//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include <iostream>

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


//=============================================================================================================
/**
* DECLARE CLASS CudaDevice
*
* @brief ToDo
*/
class CudaDevice
{
public:

	//=========================================================================================================
	/**
	* Default constructor
	*/
	CudaDevice();

 
 	//=========================================================================================================
 	/**
 	* dtor
 	* Do garbage collecting
 	*/
 	virtual ~CudaDevice();


	//=========================================================================================================
	/**
	* Initializes the HostMatrix
	*
	* @return True if a CUDA device is available
	*/
	bool cudaDeviceAvailable();

	//=========================================================================================================
	/**
	* Initializes the HostMatrix
	*
	* @return True if a CUDA 2.0 device is available
	*/
	bool cuda2xDeviceAvailable();

    //=========================================================================================================
    /**
    * Initializes the HostMatrix
    *
    * @return True if a CUDA 2.0 device is available
    */
    cudaDeviceProp getSelectedDeviceProperties()
    {
        return m_pCudaDeviceProp[m_iSelectedDevice];
    }

	//=========================================================================================================
	/**
	* Initializes the HostMatrix
	*
	* @return True if a CUDA device is available
	*/
	void printDeviceProperties();

protected:

private:
    int m_iDeviceCount;

    int m_iSelectedDevice;

    cudaDeviceProp*  m_pCudaDeviceProp;


}; //Class



}// NAMESPACE

#endif // CUDADEVICE_CUH