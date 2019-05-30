#ifndef H_EMD_H
#define H_EMD_H

#include "DataPoints.hpp"
#include "IMF.hpp"
#include "Vec.h"
#include "Spline3Interpolator.hpp"

class EMD : public DataPoints
{
  public:
    EMD();
    EMD(int, double*, double*);
    ~EMD();

    void SetVarianceParameters(double maxMed, double maxVAR, double medVRate);
    
    bool operator+=(int);
    IMF& operator[](int); // returns IMF

    // Drawing functions
    void DrawIMFs(int perPad = 0);

    // imf's size
    int size() const;
    
    // Calculates all of the IMF extras
    void UpdateAll();
    
  private:
    vector<IMF> IMFs;
    double maxMedVar;
    double maxVar;
    double medVarRate;

    bool SiftingProcess();
    
    bool CheckVariance(const IMF&);
};    

#endif
