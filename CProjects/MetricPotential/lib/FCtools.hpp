#ifndef H_FCTOOLS_H
#define H_FCTOOLS_H

#include "Vec.h"
#include <string>
#include <cstdio>
#include <iostream>
#include <fstream>
#include "FCmatrix.h"

using namespace std;

class FCtools 
{
  public:
    static vector<string> ParseStringToVector(string);    
    // expects a file in text format
    static vector<string> ReadFileToStrings(string);

    // expects a file in text format
    static vector<Vec> ReadFileTToVecs(string);

    //expects a file in binary format
    //static vector<Vec> ReadFileBToVecs(string);
    
    // prints a matrix to a text file
    static void PrintMatrixToFileT(const FCmatrix&, string, string format = "%f");

    // prints a matrix in binary file
    //static void PrintMatrixToFileB(const FCmatrix&, string);
    
};

#endif
    
    
