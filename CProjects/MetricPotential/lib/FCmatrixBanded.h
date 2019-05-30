#ifndef HFCMATRIXBANDEDh
#define HFCMATRIXBANDEDh
#include "Vec.h"
#include "FCmatrix.h"
#include "FCmatrixFull.h"
#include <vector>
#include <string>
using namespace std;

class FCmatrixBanded : public FCmatrix {
    
    public:
        FCmatrixBanded() : FCmatrix() {;}
        //BandedMatrix(Vec, Vec, Vec);
        //BandedMatrix(double*a0 = nullptr, double*b0 = nullptr, double*c0 = nullptr);
        FCmatrixBanded(vector<Vec>);
        //BandedMatrix(double **);
        FCmatrixBanded(const FCmatrix&);
        int GetFM() const;
        int GetFN() const;
        Vec GetVec(int) const;
        Vec GetRow(int) const;
        Vec GetCol(int) const;
        void Seti(int, int, double);
        vector<Vec> GetM() const;
        double GetRowMax(int);
        double GetColMax(int);
        Vec& operator[] (int);
        double Determinant();
        FCmatrixBanded operator=(const FCmatrix&);
        FCmatrixBanded operator+(const FCmatrixBanded&); // add 2 matrices of any kind
        FCmatrixBanded operator-(const FCmatrixBanded&); // sub 2 matrices of any kind
        FCmatrixFull operator*(const FCmatrix&); // mul 2 matrices of any kind
        FCmatrixBanded operator*(double lambda); // mul matrix of any kind by scalar
        FCmatrixFull operator*(const Vec&); // mul matrix by Vec
        FCmatrixBanded Transpose();
        FCmatrixFull BandedToFCmatrixFull();
        void PrintVecs();
        void PrintFull();
	string PrintToString(string format = "%f") const;
        void Print(string format = "%f") const;
        FCmatrixBanded GetL();
        FCmatrixBanded GetU();

        ~FCmatrixBanded();

    private:
        FCmatrixBanded LUdecomposition(FCmatrixBanded&);
};

#endif
