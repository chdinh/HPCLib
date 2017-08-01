//=============================================================================================================
/**
* @file     nomenclature.h
* @author   Christoph Dinh <christoph.dinh@live.de>
* @version  1.0
* @date     March, 2011
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
* @brief    Contains the nomenclature: enumerations and definitions of new types.
*
*/

#ifndef NOMENCLATURE_H
#define NOMENCLATURE_H


//*************************************************************************************************************
//=============================================================================================================
// TYPEDEF
//=============================================================================================================

typedef unsigned char  U8;  /**< Defines unsigned char as U8 type. */
typedef          char  S8;  /**< Defines char as S8 type. */
typedef unsigned short U16; /**< Defines unsigned short as U16 type. */
typedef          short S16; /**< Defines short as S16 type. */
typedef unsigned int   U32; /**< Defines unsigned int as U32 type. */
typedef          int   S32; /**< Defines int as S32 type. */


//*************************************************************************************************************
//=============================================================================================================
// DEFINE NAMESPACE HPCLib
//=============================================================================================================

namespace HPCLib
{
    //=========================================================================================================
    /**
    * Matrix type enumeration.
    */
    enum MATRIX_TYPE{
        LEADFIELD,      /**< Lead Field matrix type */
        GRID,           /**< Grid matrix type */
        SPARSEDGRID,    /**< Sparsed grid matrix type */
        MEASUREMENT,    /**< Measurement matrix type */
        _default        /**< default type */
    };

}//NAMESPACE

#endif // NOMENCLATURE_H