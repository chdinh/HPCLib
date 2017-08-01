//=============================================================================================================
/**
* @file		model.cpp
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

#ifndef MODEL_SOURCES //Because this cpp is part of the header -> template
#define MODEL_SOURCES


//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "../include/model.h"


//*************************************************************************************************************
//=============================================================================================================
// STL INCLUDES
//=============================================================================================================

#include <iostream>
#include <fstream>

#include <vector>


//*************************************************************************************************************
//=============================================================================================================
// DEFINE NAMESPACE RAPMusicOpenCV
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
Model<T>::Model()
: m_bIsInitialized(false)
, m_uiNumChannels(0)
, m_uiNumGridPoints(0)
, m_uiNumSparsedGridPoints(0)
, m_sConfigDir("")
, m_sConfigFile("")
, m_sResDir("")
, m_sResFileLeadField ("")
, m_sResFileGrid ("")
, m_sResFileSparsedGrid ("")
, m_sResFileSparsedChannels ("")
, m_pMatGrid(NULL)
, m_pMatLeadField(NULL)
, m_pVecSparsedGrid(NULL)
, m_pVecSparsedChannels(NULL)
, m_pMatSparsedLeadField(NULL)
{

}


//*************************************************************************************************************

template <class T>
Model<T>::Model(std::string p_sDataPath)
: m_bIsInitialized(false)
, m_uiNumChannels(0)
, m_uiNumGridPoints(0)
, m_uiNumSparsedGridPoints(0)
, m_sConfigDir("")
, m_sConfigFile("")
, m_sResDir("")
, m_sResFileLeadField ("")
, m_sResFileGrid ("")
, m_sResFileSparsedGrid ("")
, m_sResFileSparsedChannels ("")
, m_pMatGrid(NULL)
, m_pMatLeadField(NULL)
, m_pVecSparsedGrid(NULL)
, m_pVecSparsedChannels(NULL)
, m_pMatSparsedLeadField(NULL)
{
    initModel(p_sDataPath);
}


//*************************************************************************************************************

template <class T>
Model<T>::~Model()
{
	if(m_pMatGrid != NULL)
		delete m_pMatGrid;
	if(m_pMatLeadField != NULL)
		delete m_pMatLeadField;

	if(m_pVecSparsedGrid != NULL)
		delete m_pVecSparsedGrid;
	if(m_pVecSparsedChannels != NULL)
		delete m_pVecSparsedChannels;

	if(m_pMatSparsedLeadField != NULL)
		delete m_pMatSparsedLeadField;

}


//*************************************************************************************************************

template <class T>
bool Model<T>::initModel(std::string p_sDataPath)
{
	//Configuration Directory
	m_sConfigDir = p_sDataPath+"/config/";
	//Path to Configuration File
	m_sConfigFile = m_sConfigDir+"rap_music.cfg";
	//Resource Directory
	m_sResDir = p_sDataPath+"/data/";

	//Set by default true, when an error occures that false.
	m_bIsInitialized = true;

	//read the configuration file
	if(!readConfigFile(m_sConfigFile))
	{
		std::cout << "#####ERROR##### Configuration file could not be read!\n\n";
		m_bIsInitialized = false;
		return m_bIsInitialized;
	}

	//Loading LeadField
	if(!readFileToMatrix( m_sResFileLeadField, m_pMatLeadField, LEADFIELD))
	{
		std::cout << "#####ERROR##### Lead Field not read correctly!\n\n";
		m_bIsInitialized = false;
		return m_bIsInitialized;
	}

	//Loading Grid
	if(!readFileToMatrix( m_sResFileGrid, m_pMatGrid, GRID))
	{
		std::cout << "#####ERROR##### Grid not read correctly!\n\n";
		m_bIsInitialized = false;
		return m_bIsInitialized;
	}
	else
	{
		m_uiNumGridPoints = m_pMatGrid->rows();
	}

	//Loading SparsedGrid
	HPCMatrix<T>* t_pMatSparsedGrid = NULL;
	if(!readFileToMatrix( m_sResFileSparsedGrid, t_pMatSparsedGrid, SPARSEDGRID))
	{
		std::cout << "#####ERROR##### Sparsed grid not read correctly!\n\n";

		if (t_pMatSparsedGrid != NULL)
			delete t_pMatSparsedGrid;

		m_bIsInitialized = false;
		return m_bIsInitialized;
	}
	else
	{
		m_pVecSparsedGrid = new HPCMatrix<int>(t_pMatSparsedGrid->rows(),1);
		for(int i = 0; i < t_pMatSparsedGrid->rows(); ++i)
			m_pVecSparsedGrid->data()[i] = (int)(t_pMatSparsedGrid->data()[i]) - 1;// -1 because of MATLAB positions -> starting with 1 not with 0
		if (t_pMatSparsedGrid != NULL)
			delete t_pMatSparsedGrid;

		m_uiNumSparsedGridPoints = m_pVecSparsedGrid->rows();
	}


	//######## Check data Integrity ########

	//Number of channel is correct?
	if(m_pMatLeadField->rows() != m_uiNumChannels)
	{
		std::cout << "#####ERROR##### Number of lead field channels don't fit to the given channel number of the configuration file!\n\n";
		m_bIsInitialized = false;
		return m_bIsInitialized;
	}

	//Numbers of grid points is correct?

	//lead field is 3D?
	if(m_pMatLeadField->cols()%3)
	{
		std::cout << "#####ERROR##### Lead field is not 3D!\n\n";
		m_bIsInitialized = false;
		return m_bIsInitialized;
	}
	//Grid is 3D?
	if(m_pMatGrid->cols() != 3)
	{
		std::cout << "#####ERROR##### Grid is not 3D!\n\n";
		m_bIsInitialized = false;
		return m_bIsInitialized;
	}
	//number of lead field points fit to number of grid points?
	if(m_pMatLeadField->cols()/3 != m_pMatGrid->rows())
	{
		std::cout << "#####ERROR##### Lead field points does not fit to Grid!\n\n";
		m_bIsInitialized = false;
		return m_bIsInitialized;
	}
	else
	{
		m_uiNumGridPoints = m_pMatGrid->rows();
	}
	//SparsedGrid is a selection of points -> 1D?
	if(m_pVecSparsedGrid->cols() != 1)
	{
		std::cout << "#####ERROR##### Sparsed Grid is not 1D!\n\n";
		m_bIsInitialized = false;
		return m_bIsInitialized;
	}
	//Number of SparsedGrid points is <= Number of Grid points?
	if(m_pVecSparsedGrid->rows() > m_pMatGrid->rows())
	{
		std::cout << "#####ERROR##### SparsedGrid contains more points than the Grid!\n\n";
		m_bIsInitialized = false;
		return m_bIsInitialized;
	}

	return m_bIsInitialized;
}


//*************************************************************************************************************

template <class T>
bool Model<T>::readConfigFile(std::string p_sConfigFilePath)
{
	std::ifstream t_ifs( p_sConfigFilePath.data(), std::ifstream::in );

	if(!t_ifs){
		std::cout << "#####ERROR##### Cannot open ConfigFile: " << p_sConfigFilePath << std::endl;
		return false;
	}
	else
		std::cout << " <<<Start loading Configuration File " << p_sConfigFilePath << ">>> " << std::endl;


	std::string t_sFileCont((std::istreambuf_iterator<char>(t_ifs)), std::istreambuf_iterator<char>());

	t_ifs.close();


	std::string t_sLine;
	const std::string t_sLineDelimiter = "\n";

	const std::string t_sLineSubDelimiter = " : ";

	std::vector<std::string> t_vecTokens;

	//Read first line of the config file
	std::string::size_type size_tLastPos = t_sFileCont.find_first_not_of(t_sLineDelimiter, 0);	// Find first "non-delimiter".
	std::string::size_type size_tPos     = t_sFileCont.find_first_of(t_sLineDelimiter, size_tLastPos);	// Find first "delimiter".

	while (std::string::npos != size_tPos || std::string::npos != size_tLastPos)
	{
		t_sLine = t_sFileCont.substr(size_tLastPos, size_tPos - size_tLastPos);

		//split the line
		std::string::size_type size_tSubPos	= t_sLine.find_first_of(t_sLineSubDelimiter, 0);

		t_vecTokens.push_back(t_sLine.substr(0, size_tSubPos));
		t_vecTokens.push_back(t_sLine.substr(size_tSubPos+t_sLineSubDelimiter.size(), t_sLine.size()));


		if("channels" == t_vecTokens.at(0))
		{
			m_uiNumChannels = std::atoi(t_vecTokens.at(1).data());
			std::cout << "  -Number of Channels: " << m_uiNumChannels << std::endl;
		}
		else if("leadfield_file" == t_vecTokens.at(0))
		{
			m_sResFileLeadField = m_sResDir+t_vecTokens.at(1);
			std::cout << "  -LeadField file: " << m_sResFileLeadField << std::endl;
		}
		else if("grid_file" == t_vecTokens.at(0))
		{
			m_sResFileGrid = m_sResDir+t_vecTokens.at(1);
			std::cout << "  -Grid file: " << m_sResFileGrid << std::endl;
		}
		else if("sparsed_grid_file" == t_vecTokens.at(0))
		{
			m_sResFileSparsedGrid = m_sResDir+t_vecTokens.at(1);
			std::cout << "  -SparsedGrid file:  " << m_sResFileSparsedGrid << std::endl;
		}

		t_vecTokens.clear();

		// Next line
		size_tLastPos = t_sFileCont.find_first_not_of(t_sLineDelimiter, size_tPos);
		size_tPos = t_sFileCont.find_first_of(t_sLineDelimiter, size_tLastPos);
	}

	//Parameter Check
	if(m_uiNumChannels == 0)
	{
		std::cout << "#####ERROR##### Channels not read correctly!\n\n";
		return false;
	}
	if(m_sResFileLeadField == "")
	{
		std::cout << "#####ERROR##### Leadfield not read correctly!\n\n";
		return false;
	}
	if(m_sResFileGrid == "")
	{
		std::cout << "#####ERROR##### Grid not read correctly!\n\n";
		return false;
	}
	if(m_sResFileSparsedGrid == "")
	{
		std::cout << "#####ERROR##### SparsedGrid not read correctly!\n\n";
		return false;
	}
	std::cout << " <<<ConfigFile loaded>>> " << std::endl << std::endl;

	return true;
}



//*************************************************************************************************************

template <class T>
bool Model<T>::readFileToMatrix(const std::string p_sFilePath,
							HPCMatrix<T>* &p_pMat,
							const MATRIX_TYPE p_enumMatType)
{

	//read file to string
	std::ifstream t_ifs( p_sFilePath.data(), std::ifstream::in );

	if(!t_ifs){
		std::cout << "#####ERROR##### Cannot open Matrix File: " << p_sFilePath << std::endl;
		return false;
	}
	else
		std::cout << " <<<Start loading Matrix File " << p_sFilePath << ">>> " << std::endl;


	std::string t_sFileCont((std::istreambuf_iterator<char>(t_ifs)), std::istreambuf_iterator<char>());
	t_ifs.close();


	//parse the string to a matrix
	std::string t_sLine;
	const std::string t_sLineDelimiter = "\n";

	const std::string t_sLineSubDelimiter = " : ";
	const std::string t_sDimDelimiter = "x";
	const std::string t_sNumDelimiter = "#";

	std::vector<std::string> t_vecTokens;

	//Read first line of matrix file
	std::string::size_type size_tLastPos = t_sFileCont.find_first_not_of(t_sLineDelimiter, 0);	// Find first "non-delimiter".
	std::string::size_type size_tPos     = t_sFileCont.find_first_of(t_sLineDelimiter, size_tLastPos);	// Find first "delimiter".

	int m = 0;
	int n = 0;

	while (std::string::npos != size_tPos || std::string::npos != size_tLastPos)
	{
		t_sLine = t_sFileCont.substr(size_tLastPos, size_tPos - size_tLastPos);

		//split the line
		std::string::size_type size_tSubPos	= t_sLine.find_first_of(t_sLineSubDelimiter, 0);

		t_vecTokens.push_back(t_sLine.substr(0, size_tSubPos));
		t_vecTokens.push_back(t_sLine.substr(size_tSubPos+t_sLineSubDelimiter.size(), t_sLine.size()));


		if("Type" == t_vecTokens.at(0))
		{
			if (p_enumMatType != _default)
			{
				if ("LeadField" == t_vecTokens.at(1))
				{
					if (p_enumMatType != LEADFIELD)
					{
						std::cout << "#####ERROR#####  File contains a " << t_vecTokens.at(1) << ", which was not requested" << std::endl;
						return false;
					}
				}
				else if ("Grid" == t_vecTokens.at(1))
				{
					if (p_enumMatType != GRID)
					{
						std::cout << "#####ERROR##### File contains a " << t_vecTokens.at(1) << ", which was not requested" << std::endl;
						return false;
					}
				}
				else if ("SparsedGrid" == t_vecTokens.at(1))
				{
					if (p_enumMatType != SPARSEDGRID)
					{
						std::cout << "#####ERROR##### File contains a " << t_vecTokens.at(1) << ", which was not requested" << std::endl;
						return false;
					}
				}
				else if ("Measurement" == t_vecTokens.at(1))
				{
					if (p_enumMatType != MEASUREMENT)
					{
						std::cout << "#####ERROR##### File contains a " << t_vecTokens.at(1) << ", which was not requested" << std::endl;
						return false;
					}
				}
				else
				{
					std::cout << "#####ERROR##### Matrix type " << t_vecTokens.at(1) << " unknown" << std::endl;
					return false;
				}
			}

			std::cout << "  -Type: " << t_vecTokens.at(1) << std::endl;
		}
		else if("Dimension" == t_vecTokens.at(0))
		{
			//get position where dimension is split by "x"
			std::string::size_type dimPos = t_vecTokens.at(1).find_first_of(t_sDimDelimiter, 0);

			//Extract Dimension
			m = atoi(t_vecTokens.at(1).substr(0, dimPos).data());
			n = atoi(t_vecTokens.at(1).substr(dimPos+t_sDimDelimiter.size(), t_vecTokens.at(1).size()).data());
			
			std::cout << "  -Dimension: " << m << " x " << n << std::endl;
		}
		else if("Matrix" == t_vecTokens.at(0))
		{
			if ("begin" == t_vecTokens.at(1))
			{
				std::cout << "  -Matrix begin " << std::endl;
				if (!m && !n)
				{
					std::cout << "#####ERROR##### Dimension not set" << std::endl;
					return false;
				}

				if(p_pMat != NULL)
					delete p_pMat;

				p_pMat = new HPCMatrix<T>(m, n);

				for (int i = 0; i < m; i++)
				{
					//new line range
					size_tLastPos = t_sFileCont.find_first_not_of(t_sLineDelimiter, size_tPos);
					size_tPos = t_sFileCont.find_first_of(t_sLineDelimiter, size_tLastPos);

					t_sLine = t_sFileCont.substr(size_tLastPos, size_tPos - size_tLastPos);

					// initial Number range
					std::string::size_type size_tLastNumPos = t_sLine.find_first_not_of(t_sNumDelimiter, 0);
					std::string::size_type size_tNumPos     = t_sLine.find_first_of(t_sNumDelimiter, size_tLastNumPos);

					int j = 0;//column count
					while (std::string::npos != size_tNumPos || std::string::npos != size_tLastNumPos)
					{
// 						//read in row major
// 						p_pMat->data()[i*n+j] = (T)atof(t_sLine.substr(size_tLastNumPos, size_tNumPos - size_tLastNumPos).data());

						//read in column major
						p_pMat->data()[i+j*m] = (T)atof(t_sLine.substr(size_tLastNumPos, size_tNumPos - size_tLastNumPos).data());

						// New number range
						size_tLastNumPos = t_sLine.find_first_not_of(t_sNumDelimiter, size_tNumPos);
						size_tNumPos     = t_sLine.find_first_of(t_sNumDelimiter, size_tLastNumPos);
						++j;
					}
				}
			} 
			else
			{
				std::cout << "  -Matrix end " << std::endl;
			}
		}


		t_vecTokens.clear();

		// Next line
		size_tLastPos = t_sFileCont.find_first_not_of(t_sLineDelimiter, size_tPos);
		size_tPos = t_sFileCont.find_first_of(t_sLineDelimiter, size_tLastPos);
	}

	std::cout << " <<<Matrix loaded>>> " << std::endl << std::endl;

	return true;
}


}//Namespace


#endif //MODEL_SOURCES
