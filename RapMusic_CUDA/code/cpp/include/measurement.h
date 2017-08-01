//=============================================================================================================
/**
* @file		measurement.h
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


#ifndef MEASUREMENT_H
#define MEASUREMENT_H


//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "hpcmatrix.h"
#include "nomenclature.h"


//*************************************************************************************************************
//=============================================================================================================
// STL INCLUDES
//=============================================================================================================

#include <vector>
#include <string>
#include <iostream>
#include <deque>
	

//*************************************************************************************************************
//=============================================================================================================
// QT INCLUDES
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


//=============================================================================================================
/**
* DECLARE CLASS Measurement
*
* @brief The Measurement class loads and holds a static measurement + ToDo Online measurement
*/
template<class T>
class Measurement
{
public:

	//=========================================================================================================
	/**
	* Default constructor
	*/
	Measurement();

	//=========================================================================================================
	/**
	* ctor
	*
	* @param p_sDataPath the path to the directory which contains the data folder.
	*/
	Measurement(std::string p_sDataPath);

 
 	//=========================================================================================================
 	/**
 	* dtor
 	* Do garbage collecting
 	*/
 	virtual ~Measurement();


	//=========================================================================================================
	/**
	* Initializes the Model
	*
	* @param p_sDataPath the path to the directory which contains the resource folder ("data").
	* @return True if loaded successfully, false otherwise
	*/
	bool initMeasurement(std::string p_sDataPath);


	//=========================================================================================================
	/**
	* Reads the configuration file
	*
	* @param p_sFilePath the path to the file which contains the configuration file.
	* @return True if loaded successfully, false otherwise
	*/
	bool readConfigFile(std::string p_sFilePath);

	
	//=========================================================================================================
	/**
	* Reads a txt matrix file
	*
	* @param p_sFilePath the path to the matrix file 
	* @param p_pMat reference to the destination matrix pointer which is new instantiated inside this function.
	* @param p_enumMatType the Matrix type which should be loaded.
	* @return true if loaded successfully, false otherwise
	*/
	bool readFileToMatrix(const std::string p_sFilePath, 
						HPCMatrix<T>* &p_pMat, 
						MATRIX_TYPE p_enumMatType = _default);


	//=========================================================================================================
	/**
	* Is Model initialized.
	*
	* @return true when model is initialized
	*/
	bool isInitialized() const
	{
		return m_bIsInitialized;
	}

	HPCMatrix<T>* getMeasurementMat()
	{
		return m_pHPCMatMeasurement;
	}

	T* getMeasurementData()
	{
		return m_pHPCMatMeasurement->data();
	}

protected:

private:
	bool m_bIsInitialized;		/**< ToDo Documentation */

	unsigned int m_uiNumChannels;			/**< ToDo Documentation */

	std::string m_sConfigDir;		/**< ToDo Documentation */
	std::string m_sConfigFile;		/**< ToDo Documentation */

	std::string m_sResDir;					/**< ToDo Documentation */
	std::string m_sResFileMeasurement;		/**< ToDo Documentation */

	HPCMatrix<T>* m_pHPCMatMeasurement;				/**< ToDo Documentation */
};

} // NAMESPACE

//Make the template definition visible to compiler in the first point of instantiation
#include "../src/measurement.cpp"

#endif // MODEL_H