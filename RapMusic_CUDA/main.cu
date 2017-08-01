//=============================================================================================================
/**
* @file     main.cu
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
* @brief    Implements the main() application function.
*
*/

//*************************************************************************************************************
//=============================================================================================================
// CPP INCLUDES
//=============================================================================================================

#include "code/cpp/include/model.h"
#include "code/cpp/include/rapinterface.h"

#include "code/cpp/include/rapdipoles.h"
#include "code/cpp/include/measurement.h"


//*************************************************************************************************************
//=============================================================================================================
// CUDA INCLUDES
//=============================================================================================================

#include "code/cuda/include/cudadevice.cuh"

#include "code/cuda/include/rapmusic_cuda.cuh"


//*************************************************************************************************************
//=============================================================================================================
// STL INCLUDES
//=============================================================================================================

#include <iostream>


//*************************************************************************************************************
//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================


//*************************************************************************************************************
//=============================================================================================================
// MAIN
//=============================================================================================================

//=============================================================================================================
/**
* The function main marks the entry point of the program.
* By default, main has the storage class extern.
*
* @param [in] argc (argument count) is an integer that indicates how many arguments were entered on the command line when the program was started.
* @param [in] argv (argument vector) is an array of pointers to arrays of character objects. The array objects are null-terminated strings, representing the arguments that were entered on the command line when the program was started.
* @return the value that was set to exit() (which is 0 if exit() is called via quit()).
*/
int main(int argc, char *argv[])
{
    //create model
    HPCLib::Model<float> *t_pModel = new HPCLib::Model<float>(".");

    //read Measurement
    HPCLib::Measurement<float> *t_pMeasurement = new HPCLib::Measurement<float>(".");

    //Get GPU Device Parameters
    HPCLib::CudaDevice* t_pCudaDevice = new HPCLib::CudaDevice();

    //When CUDA Device is available do CUDA stuff
    if(t_pCudaDevice->cuda2xDeviceAvailable())
    {
        std::cout << "================================================================" << std::endl;
        std::cout << "======================== RAP MUSIC CUDA ========================" << std::endl;
        std::cout << "================================================================" << std::endl << std::endl;

        std::cout << std::endl << "CUDA 2.x device available" << std::endl << std::endl;
        t_pCudaDevice->printDeviceProperties();
        std::cout << std::endl << std::endl;

        //----------------------------------------------------------------------------------
        HPCLib::RapMusic_Cuda *t_pRapCuda = new HPCLib::RapMusic_Cuda();
    
        //init memorys etc for cuda run
        t_pRapCuda->initRAPMusic(t_pCudaDevice, t_pModel, NULL, 7); //give handle of Model_interface to the init-function

        //Destination for Return Dipoles
        HPCLib::RapDipoles<float> *t_pDipoles = NULL;

        //run the rap music algorithms using cuda
        
        //t_pRapCuda->calcRapMusic(t_pMeasurement->getMeasurementMat(), t_pDipoles);

        t_pRapCuda->calcPowellRAPMusic(t_pMeasurement->getMeasurementMat(), t_pDipoles);

        //get correlation results
        t_pDipoles->print();


        //clean up dipoles
        if(t_pDipoles != NULL)
            delete t_pDipoles;

        //clean up
        delete t_pRapCuda;
    }
    else
    {
        std::cout << "Skipping RAP MUSIC CUDA, because no CUDA 2.x compatible device available" << std::endl << std::endl;
    }
    //clean up
    delete t_pCudaDevice;

    //RAP MUSIC gold -> CPU Implementation
    std::cout << std::endl << std::endl << std::endl;
    std::cout << "================================================================" << std::endl;
    std::cout << "======================== RAP MUSIC GOLD ========================" << std::endl;
    std::cout << "================================================================" << std::endl << std::endl;
    HPCLib::RapInterface* t_pRapInterface = new HPCLib::RapInterface( t_pModel );

    t_pRapInterface->calcRapMusic_gold( t_pMeasurement );

    std::cout << std::endl << std::endl << std::endl;
    std::cout << "================================================================" << std::endl;
    std::cout << "===================== Powell RAP MUSIC GOLD ====================" << std::endl;
    std::cout << "================================================================" << std::endl << std::endl;
    t_pRapInterface->calcPowellRapMusic_gold( t_pMeasurement );

    //clean up rap interface
    delete t_pRapInterface;

    //----------------------------------------------------------

    //clean up measurement
    delete t_pMeasurement;
    
    //clean up model
    delete t_pModel;


    return 0;
}
