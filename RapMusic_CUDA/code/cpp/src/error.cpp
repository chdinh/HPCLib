//=============================================================================================================
/**
* @file		error.cpp
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

#include "../include/error.h"

#include <string.h>

//*************************************************************************************************************
//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace HPCLib;



//*************************************************************************************************************
//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================


Error::Error(const char *title)
: prog(0)
{
	int l = strlen( title ) + 1;

	prog = new char[l];
	strcpy_s(prog,l, title);

	clear() ;
}

/*!
  \brief reports an error message
  The usual call to this function is made through the 
               $\ll$ operator, {\em i.e.} 
	       \begin{verbatim}
	         Error error("routine name") ;
		 error << "An error occured " ;
		 error.report() ;
	       \end{verbatim}
  \param msg  the message to report
  \return 
  \warning
*/
void Error::report(const char *msg)
{
	if ( msg == 0 )		// message must have been sent via <<
		std::cerr << str() ;
	else
		std::cerr << msg;

	std::cerr << '\n';  
}

/*!
  \brief reports a warning error message
  The usual call to this function is made through the 
               $\ll$ operator, {\em i.e.} 
	       \begin{verbatim}
	         Error error("routine name") ;
		 error << "An error occured " ;
		 error.warning() ;
	       \end{verbatim}
  \param msg  the message to report
  \return 
  \warning
*/
void Error::warning(const char *msg)
{
	std::cerr << "\nRoutine: " << prog << "\nWarning: ";

	report( msg );
}

/*!
  \brief reports a fatal error
  Reports a fatal error. If the \verb.DEBUG_MATRIX. flag has been
               set then the routine starts an infinte loop. An exit(1) is
	       executed otherwise.

	       The usual call to this function is made through the 
               $\ll$ operator, {\em i.e.} 
	       \begin{verbatim}
	         Error error("routine name") ;
		 error << "A fatal error occured " ;
		 error.fatal() ;
	       \end{verbatim}

  \param msg  a message to report
  \return 
  \warning
*/
void Error::fatal(const char *msg)
{
  std::cerr << "\nRoutine: " << prog << "\nFatal error: ";
  
  report( msg );
  exit(1);
}

/*!
  \brief reports a memory allocation error
  Reports a memory allocation error. If the DEBUG_PLIB
  flag has been set then the routine starts an infinte loop. 
  An exit(1) is executed otherwise.
	       
  The usual call to this function is made through the 
   operator>> , {\em i.e.} 
  \begin{verbatim}
	         Error error("variable name") ;
		 float* var = new float[bigSize] ;
		 error.memory(var) ;
		\end{verbatim}

  \param p  the pointer to test
  \return 
  \warning
*/
void Error::memory(const void *p)
{
  if ( p == 0)
    {
      std::cerr << "\nRoutine: " << prog << " Memory allocation error\n";
      exit(1);
    }  
}


