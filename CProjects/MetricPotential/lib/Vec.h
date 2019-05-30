#ifndef HVECH
#define HVECH
#include <iostream>
#include <string>
using namespace std;

class Vec {

    public:
        Vec(int n = 0, double *a = nullptr);
        Vec(int, double);
        Vec(const Vec&);
        ~Vec();
        void SetEntries (int, double*);
        Vec operator = (const Vec&);
        Vec operator + (const Vec&); 
        Vec operator += (const Vec&);
        Vec operator - (const Vec&);
        Vec operator -= (const Vec&);
        Vec operator - ();
        Vec operator + ();
        Vec operator * (const Vec&);
        Vec operator * (double);
	double& operator[](int) const;
	int size() const;
        double dot(const Vec&);
        void Print(string format="%f") const;
	string PrintToString(string format="%f") const;
	void swap(int, int);
        void Seti(int, double); //coloca o double na posição 
	
    private:
        int N;
        double *entries;
};

#endif
