#ifndef HFCMATRIXFULLh
#define HFCMATRIXFULLh
#include "Vec.h"
#include "FCmatrix.h"
#include <cstdio>
#include <vector>
#include <string>
using namespace std;

class FCmatrixFull : public FCmatrix {
    
    public:
        FCmatrixFull(double**, int, int); //matrix fm x fn
        FCmatrixFull(double* fM = nullptr, int fm = 0, int fn = 0);
        FCmatrixFull(vector<Vec>);
        FCmatrixFull(const FCmatrix&);
        int GetFM() const;
        int GetFN() const;
        Vec GetRow(int) const;
        Vec GetCol(int) const;
        vector<Vec> GetM() const;
        double GetRowMax(int);
        double GetColMax(int);
        double Determinant();
        Vec& operator[] (int);
        FCmatrixFull operator=(const FCmatrix&);
        FCmatrixFull operator+(const FCmatrix&); // add 2 matrices of any kind
        FCmatrixFull operator-(const FCmatrix&); // sub 2 matrices of any kind
        FCmatrixFull operator*(const FCmatrix&); // mul 2 matrices of any kind
        FCmatrixFull operator*(double lambda); // mul matrix of any kind by scalar
        FCmatrixFull operator*(const Vec&); // mul matrix by Vec
	Vec dot(const Vec&);
	
	string PrintToString(string format = "%f") const;
        void Print(string format = "%f") const;
        void PrintIndices();
        int Getrowindices(int);
        int* Getcolindices();
        void swapRows(int,int);
        void swapCols(int,int);
        FCmatrixFull Transpose();
        FCmatrixFull EliminateRow(int);
        FCmatrixFull EliminateCol(int);

        ~FCmatrixFull();
    
    protected:
        int* rowindices = NULL; // row indices (0,1,2,...)
        int* colindices = NULL; // col indices (0,1,2,...)
};

#endif
