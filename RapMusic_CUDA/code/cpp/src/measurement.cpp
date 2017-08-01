//=============================================================================================================
/**
* @file		measurement.cpp
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

#ifndef MEASUREMENT_SOURCES //Because this cpp is part of the header -> template
#define MEASUREMENT_SOURCES


//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "../include/measurement.h"


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
Measurement<T>::Measurement()
: m_bIsInitialized(false)
, m_uiNumChannels (0)
, m_sConfigDir("")
, m_sConfigFile("")
, m_sResDir("")
, m_sResFileMeasurement("")
, m_pHPCMatMeasurement(NULL)
{

}


//*************************************************************************************************************

template <class T>
Measurement<T>::Measurement(std::string p_sDataPath)
: m_bIsInitialized(false)
, m_uiNumChannels (0)
, m_sConfigDir("")
, m_sConfigFile("")
, m_sResDir("")
, m_sResFileMeasurement("")
, m_pHPCMatMeasurement(NULL)
{
	initMeasurement(p_sDataPath);
}


//*************************************************************************************************************

template <class T>
Measurement<T>::~Measurement()
{
	if(m_pHPCMatMeasurement != NULL)
		delete m_pHPCMatMeasurement;
}


//*************************************************************************************************************

template <class T>
bool Measurement<T>::initMeasurement(std::string p_sDataPath)
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

	//Loading Measurement to HPCMatrix
	if(!readFileToMatrix( m_sResFileMeasurement, m_pHPCMatMeasurement, MEASUREMENT))
	{
		std::cout << "#####ERROR##### Lead Field not read correctly!\n\n";
		m_bIsInitialized = false;
		return m_bIsInitialized;
	}

	//######## Check data Integrity ########

	//Number of channel is correct?
	if(m_pHPCMatMeasurement->rows() != m_uiNumChannels)
	{
		std::cout << "#####ERROR##### Number of Measurement channels don't fit the given channel number of the configuration file!\n\n";
		m_bIsInitialized = false;
		return m_bIsInitialized;
	}

	return m_bIsInitialized;
}


//*************************************************************************************************************

template <class T>
bool Measurement<T>::readConfigFile(std::string p_sConfigFilePath)
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
		else if("measurement_file" == t_vecTokens.at(0))
		{
			m_sResFileMeasurement = m_sResDir+t_vecTokens.at(1);
			std::cout << "  -Measurement file: " << m_sResFileMeasurement << std::endl;
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
	if(m_sResFileMeasurement == "")
	{
		std::cout << "#####ERROR##### Measurement not read correctly!\n\n";
		return false;
	}
	std::cout << " <<<ConfigFile loaded>>> " << std::endl << std::endl;

	return true;
}


//*************************************************************************************************************

template <class T>
bool Measurement<T>::readFileToMatrix(const std::string p_sFilePath,
							HPCMatrix<T>* &p_pMat,
							MATRIX_TYPE p_enumMatType)
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


#endif //MEASUREMENT_SOURCES
