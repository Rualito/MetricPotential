#ifndef HFCMATRIXh
#define HFCMATRIXh
#include "Vec.h"
#include <vector>
#include <string>
using namespace std;

class FCmatrix {
    
    public:
        FCmatrix(double**, int, int); //matrix fm x fn
        FCmatrix(double* fM = nullptr, int fm = 0, int fn = 0);
        FCmatrix(const vector<Vec>);
        virtual string Getclassname() const;

	virtual int GetFM() const;
        virtual int GetFN() const;
        virtual Vec& operator[] (int) = 0;
        virtual Vec GetRow(int) const = 0 ; // retrieve row i
        virtual Vec GetCol(int) const = 0; // retrieve column i
        virtual vector<Vec> GetM() const = 0;
        virtual double Determinant() = 0;
	virtual string PrintToString(string format = "%f") const = 0;
        virtual void Print(string format = "%f") const = 0;
	virtual double GetRowMax(int)  = 0;
        virtual double GetColMax(int) = 0;
        virtual ~FCmatrix();

	// its diferent from GetM, this isnt supposed to be overloaded
	virtual vector<Vec> GetRawM() const;
	
    protected:
        vector<Vec> M;
        string classname;
};



#endif
