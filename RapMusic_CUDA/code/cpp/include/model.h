//=============================================================================================================
/**
* @file		model.h
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


#ifndef MODEL_H
#define MODEL_H


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
* DECLARE CLASS Model
*
* @brief The Model class loads and holds static data (e.g. Lead field, Grid, SparsedGrid)
*/
template<class T>
class Model
{
public:

	//=========================================================================================================
	/**
	* Default constructor
	*/
	Model();

	//=========================================================================================================
	/**
	* ctor
	*
	* @param p_sDataPath the path to the directory which contains the data folder.
	*/
	Model(std::string p_sDataPath);

 
 	//=========================================================================================================
 	/**
 	* dtor
 	* Do garbage collecting
 	*/
 	virtual ~Model();


	//=========================================================================================================
	/**
	* Initializes the Model
	*
	* @param p_sDataPath the path to the directory which contains the resource folder ("data").
	* @return True if loaded successfully, false otherwise
	*/
	bool initModel(std::string p_sDataPath);


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
	* returns the number of grid points
	*
	* @return the number of grid points
	*/
	inline unsigned int getNumGridPoints() const
	{return m_uiNumGridPoints;};

	//=========================================================================================================
	/**
	* returns the number of sparsed grid points
	*
	* @return the number of sparsed grid points
	*/
	inline unsigned int getNumSparsedGridPoints() const
	{return m_uiNumSparsedGridPoints;};


	//=========================================================================================================
	/**
	* returns the number of channels
	*
	* @return the number of channels
	*/
	inline unsigned int getNumChannels() const
	{return m_uiNumChannels;};


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

	//LeadField
	HPCMatrix<T>* getLeadFieldMat()
	{
		return m_pMatLeadField;
	}

	T* getLeadFieldData()
	{
		return m_pMatLeadField->data();
	}

	//SparsedLeadField
	HPCMatrix<T>* getSparsedLeadFieldMat()
	{
		return m_pMatLeadField;
	}

	T* getSparsedLeadFieldData()
	{
		return m_pMatLeadField->data();
	}

	//Grid
	HPCMatrix<T>* getGridMat()
	{
		return m_pMatGrid;
	}

	T* getGridData()
	{
		return m_pMatGrid->data();
	}

	//SparsedGrid
	HPCMatrix<T>* getSparsedGridMat()
	{
		return m_pMatSparsedGrid;
	}

	T* getSparsedGridData()
	{
		return m_pMatSparsedGrid->data();
	}

protected:

private:
	bool m_bIsInitialized;		/**< ToDo Documentation */

	unsigned int m_uiNumChannels;			/**< ToDo Documentation */
	unsigned int m_uiNumGridPoints;			/**< ToDo Documentation */
	unsigned int m_uiNumSparsedGridPoints;	/**< ToDo Documentation */

	std::string m_sConfigDir;		/**< ToDo Documentation */
	std::string m_sConfigFile;		/**< ToDo Documentation */

	std::string m_sResDir;					/**< ToDo Documentation */
	std::string m_sResFileLeadField;		/**< ToDo Documentation */
	std::string m_sResFileGrid;				/**< ToDo Documentation */
	std::string m_sResFileSparsedGrid;		/**< ToDo Documentation */
	std::string m_sResFileSparsedChannels;	/**< ToDo Documentation */

	//T
	HPCMatrix<T>* m_pMatGrid;				/**< ToDo Documentation */ //Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> instead of MatrixXf cause of template stuff
	HPCMatrix<T>* m_pMatLeadField;		/**< ToDo Documentation */

	HPCMatrix<int>* m_pVecSparsedGrid;		/**< ToDo Documentation */
	HPCMatrix<int>* m_pVecSparsedChannels;	/**< ToDo Documentation */

	HPCMatrix<T>* m_pMatSparsedLeadField;		/**< ToDo Documentation */
	HPCMatrix<T>* m_pMatSparsedGrid;		/**< ToDo Documentation */
};

} // NAMESPACE

//Make the template definition visible to compiler in the first point of instantiation
#include "../src/model.cpp"

#endif // MODEL_H