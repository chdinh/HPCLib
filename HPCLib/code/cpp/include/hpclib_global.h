//=============================================================================================================
/**
* @file     hpclib_global.h
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
* @brief    ToDo
*
*/


#ifndef HPCLIB_GLOBAL_H
#define HPCLIB_GLOBAL_H


//*************************************************************************************************************
//=============================================================================================================
// PREPROCESSOR DEFINES
//=============================================================================================================

//#define HPCLIB_EXPORT   __declspec(dllimport) //-> for static linking
#define HPCLIB_EXPORT   __declspec(dllexport)

// #ifdef HPCLib_LIB
// # define HPCLIB_EXPORT __declspec( dllexport )
// #else
// # define HPCLIB_EXPORT __declspec( dllimport ) //-> for static linking
// #endif

#endif // HPCLIB_GLOBAL_H
