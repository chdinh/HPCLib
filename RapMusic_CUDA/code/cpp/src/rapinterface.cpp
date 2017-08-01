//=============================================================================================================
/**
* @file		rapinterface.cpp
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

//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "../include/rapinterface.h"

#include "../include/model.h"
#include "../include/measurement.h"
#include "../include/rapmusic_gold.h"
#include "../include/rapdipoles.h"



#include "../include/dipole.h"



//*************************************************************************************************************
//=============================================================================================================
// STL INCLUDES
//=============================================================================================================

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


//*************************************************************************************************************
//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

RapInterface::RapInterface(Model<float>* p_pModel)
:// m_pModel (new Model<float>(p_sDataPath))
//, m_pMeasurement(new Measurement<float>(p_sDataPath))
 m_pRapMusic (new RapMusic<float>(p_pModel, false, 7))
//, m_pRapMusic (new RapMusic<float>(m_pModel->getLeadFieldMat(), m_pModel->getGridMat(), 4))
//, m_pRapMusic (new RapMusic<float>(m_pModel->getLeadFieldMat(), NULL, 4))
{

}


//*************************************************************************************************************

RapInterface::~RapInterface()
{
// 	if(m_pModel != NULL)
// 		delete m_pModel;
// 	if(m_pMeasurement != NULL)
// 		delete m_pMeasurement;
	if(m_pRapMusic != NULL)
		delete m_pRapMusic;
}


//*************************************************************************************************************

bool RapInterface::calcRapMusic_gold(Measurement<float>* p_pMeasurement)
{ 
	RapDipoles<float> *t_pDipoles = NULL;

	m_pRapMusic->calcRAPMusic( p_pMeasurement->getMeasurementMat(), t_pDipoles);

	t_pDipoles->print();

	if(t_pDipoles != NULL)
		delete t_pDipoles;

	return true;
}


//*************************************************************************************************************

bool RapInterface::calcPowellRapMusic_gold(Measurement<float>* p_pMeasurement)
{ 
	RapDipoles<float> *t_pDipoles = NULL;

	m_pRapMusic->calcPowellRAPMusic( p_pMeasurement->getMeasurementMat(), t_pDipoles);

	t_pDipoles->print();

	if(t_pDipoles != NULL)
		delete t_pDipoles;

	return true;
}


}//Namespace