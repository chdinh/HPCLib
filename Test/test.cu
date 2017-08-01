//=============================================================================================================
/**
* @file     test.cu
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

#include "../HPCLib/code/cpp/include/MathFuncsDLL.h"

//#include "../HPCLib/code/cpp/include/hpcmatrix.h"


//*************************************************************************************************************
//=============================================================================================================
// CPP INCLUDES
//=============================================================================================================


//*************************************************************************************************************
//=============================================================================================================
// CUDA INCLUDES
//=============================================================================================================


//*************************************************************************************************************
//=============================================================================================================
// STL INCLUDES
//=============================================================================================================

#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <fstream>


//*************************************************************************************************************
//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

int offset(int row, int numPoints)
{

return row*numPoints - (( (row-1)*row) / 2); //triangular series 1 3 6 10 ... = (num_pairs*(num_pairs+1))/2

}




void idxVec(int row, int numPoints, int*& elementsVec)
{

//     if(elementsVec != NULL)
//         delete[] elementsVec;
// 
//     elementsVec = new int(numPoints);


    //col combination index
    for(int i = 0; i <= row; ++i)//=numPoints-1
        elementsVec[i] = offset(i+1,numPoints)-(numPoints-row);


    //row combination index
    int off = offset(row,numPoints);
    int length = numPoints - row;
    int k=0;
    for(int i = row; i < row+length; ++i)//=numPoints-1
    {
        elementsVec[i] = off+k;
        k = k + 1;
    }
}



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
    double a;
    a = MathFuncs::MyMathFuncs::Add(10,12);

    std::cout << "Ergebnis: " << a << std::endl;


    //create combinations
    const int numPoints = 121;
    const int numCombinations = ((numPoints+1)*numPoints)/2;
    int combis[numCombinations*2];

    int k = 0;
    for(int i = 0; i < numPoints; ++i)
    {
        for(int j = i; j < numPoints; ++j)
        {
            combis[k*2] = i;
            combis[k*2+1] = j;
            k++;
        }
    }

    //load subcorr simulation
    std::fstream filestr;
    filestr.open ("test.txt", std::fstream::in | std::fstream::out | std::fstream::app);

    double val;




    double subcorr[numCombinations];
    double roh[numCombinations];

    for(int i = 0; i < numCombinations; ++i)
    {
//        std::cout << combis[i*2] << ", " << combis[i*2+1] << std::endl;

        /* generate secret number: */
            
        filestr >> val;
        subcorr[i] = val;
        roh[i] = 0;

//        std::cout << subcorr[i] << std::endl << std::endl;
    }

   filestr.close();






    int current_row = 2;
    int current_combination[2] = {-1, -1};

    double max_val = 0;
    int new_max_idx = -1;
    int max_idx = -1;

    int maxfound = 0;

    int* elVec = new int[numPoints];

    idxVec(current_row, numPoints, elVec);

    int elVec_length = numPoints;

//     std::cout << "IndexVector" << std::endl;
//     for(int i = 0; i < numPoints; ++i)
//     {
//         std::cout << elVec[i] << std::endl;
// 
//     }


    while(maxfound == 0)
    {
        std::cout << "Curr Combi: " << current_combination[0] << ", " << current_combination[1] << std::endl;

        //Subcorr simulation
        for(int i = 0; i < elVec_length; ++i)//length(sel_ele)-1
        {
            k = elVec[i];

            //actual calculate subcorr with combis
            roh[k] = subcorr[k];
        }


        //find maximum
        max_val = 0;
        for(int i = 0; i < numCombinations; ++i)
        {
            if(roh[i] > max_val)
            {
                max_val = roh[i];
                new_max_idx = i;
            }
        }

        if(new_max_idx == max_idx)
        {
            std::cout << "Result: " << current_combination[0] << ", " << current_combination[1] << std::endl;
            maxfound = 1;
            break;
        }
        else
        {
            max_idx = new_max_idx;
            current_combination[0] = combis[2*max_idx];
            current_combination[1] = combis[2*max_idx+1];
        }



        //set new index
        if(current_combination[0] == current_row)
            current_row = current_combination[1];
        else
            current_row = current_combination[0];


        idxVec(current_row, numPoints, elVec);

    }







//     HPCLib::HPCMatrix<float> myMatrix(2,1);
// 
//     myMatrix(0,0) = 4;
// 
//     myMatrix(1,0) = 4;
// 
// 
//     std::cout << "Matrix: " << myMatrix << std::endl;

    return 0;
}
