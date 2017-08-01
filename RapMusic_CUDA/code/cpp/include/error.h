//=============================================================================================================
/**
* @file		error.h
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

#ifndef ERROR_H
#define ERROR_H


//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================



//*************************************************************************************************************
//=============================================================================================================
// STL INCLUDES
//=============================================================================================================

#include <stdlib.h>
#include <iostream>
#include <sstream>


//*************************************************************************************************************
//=============================================================================================================
// TYPEDEFS
//=============================================================================================================

typedef std::ostringstream ErrorStream;


//*************************************************************************************************************
//=============================================================================================================
// DEFINE NAMESPACE HPCLib
//=============================================================================================================

namespace HPCLib {

	struct MatrixErr { 
		MatrixErr() { print_debug(); }
		void print_debug() { 
			#ifdef VERBOSE_EXCEPTION
				print();
			#else
				;
			#endif
		}
		virtual void print() { std::cerr << "Matrix error.\n" ; }
	};


//*************************************************************************************************************

	struct MatrixInputError : public MatrixErr {
		MatrixInputError() { print_debug();}
		virtual void print(){
			std::cerr << "One of the input value is not in appropriate.\n";
		}
	};


//*************************************************************************************************************

	struct OutOfBound : public MatrixInputError {
		int i ;
		int s,e ;
		OutOfBound(int index, int from, int to): i(index), s(from), e(to) { print_debug(); }
		virtual void print() { 
			std::cerr << "Out of bound error, trying to access " << i << 
				" but the valid range is [ " << s << "," << e << "]\n" ; 
		}
	};


//*************************************************************************************************************

	struct OutOfBound2D : public MatrixInputError {
		int i,j ;
		int s_i,e_i ;
		int s_j,e_j ;
		OutOfBound2D(int I, int J, int fI, int tI, int fJ, int tJ): i(I), j(J), s_i(fI), e_i(tI), s_j(fJ), e_j(tJ) { print_debug(); } 
		virtual void print() { 
			std::cerr << "Out of bound error, trying to access (" << i << ',' << j <<
		") but the valid range is ([ " << s_i << "," << e_i << "], [" <<
		s_j << ',' << e_j << "])\n" ;
		}
	};


//*************************************************************************************************************

	struct WrongSize : public MatrixInputError {
		int s1,s2 ;
		WrongSize(int a, int b) : s1(a), s2(b) { print_debug();}
		virtual void print(){
			std::cerr << "The vector sizes  " << s1 << " and " << s2 << " are incompatible.\n" ;
		}
	}; 


//*************************************************************************************************************

	struct WrongSize2D : public MatrixInputError {
		int rows,cols ;
		int bad_rows, bad_cols ;
		WrongSize2D(int r, int c, int br, int bc) : rows(r), cols(c), bad_rows(br), bad_cols(bc) { print_debug();}
		virtual void print(){
			std::cerr << "The matrix sizes  (" << rows << " x " << cols << ") and (" << bad_rows << " x " << bad_cols << ") are incompatible.\n" ;
		}
	}; 
  
//=============================================================================================================
/**
*  A class to print and handle error messages
*
* @brief The HostMatrix class...
*/
class Error : public ErrorStream
{
public:
	Error(): ErrorStream(), prog(0) {}
	Error(const char *s);
	~Error(){ if (prog) delete []prog ; }

	void warning(const char *msg = 0);
	void nonfatal(const char *msg = 0) { warning(msg); }
	void fatal(const char * = 0 );
	void memory(const void * = 0 );

private:
	char* prog;
	void report(const char *msg = NULL);
};

} // end namespace

#endif

