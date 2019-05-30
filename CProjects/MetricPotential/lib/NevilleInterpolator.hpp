#ifndef H_NEVILLEINTER_H
#define H_NEVILLEINTER_H

#include "DataPoints.hpp"
#include "TF1.h"
#include <string>

class NevilleInterpolator : public DataPoints 
{
  public:
    NevilleInterpolator(int N=0, double *xs=NULL, double *ys=NULL);
    ~NevilleInterpolator();

    double Interpolate(double);

    void SetFunction(TF1* f);
    TF1* GetInterpolationFunction();

    void Draw(double step = 0);
    //string PrintToString(string format);
    void Print(string filename="", string format = "%f");

  private:
    double fInterpolator(double *fx, double *par);
    TF1* FInterpolator;
    TF1* F0;
};

	
#endif
