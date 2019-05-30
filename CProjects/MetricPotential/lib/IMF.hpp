#ifndef H_IMF_H
#define H_IMF_H

#include "Vec.h"
#include "DataPoints.hpp"
#include "Spline3Interpolator.hpp"
#include <cstdint>
#include <cmath>

// one-hot
#define DRAW_SPLINEMAX 1
#define DRAW_SPLINEMIN 2
#define DRAW_SPLINEMED 4
#define DRAW_SPLINEAMP 8
#define DRAW_FREQUENCY 16
#define DRAW_INTRINSIC 32
#define DRAW_VARIANCE 64
#define DRAW_OLDINTRINSIC 128
#define DRAW_DATA 256
#define DRAW_MEAN 512

class IMF : public DataPoints
{
  public:
    IMF();
    IMF(int n, double* xp=NULL, double* yp=NULL);
    IMF(const IMF&);
    ~IMF();
    
    double* GetVariance() const;
    
    // double[N] with current 'IMF'
    double* GetIntrinsic() const;
    
    double GetCorrelation(double*);
    
    // draws the residual, last refined
    void DrawIntrinsic();
    // draws the signal, the maxSpline, the minSpline, the medSpline
    void DrawSplines(double step = 0.01);
    
    bool operator++();
    const IMF& operator=(const IMF&);
    double operator[](int) const; // residual on point i

    double GetAmplitude(double);
    double* GetFrequency();
    double* GetMean();
    
    int GetMaxCount();
    int GetMinCount();
    int GetZeroCount();

    // calculates the extra parameters if queued
    void CalculateExtras();

    TGraph* GetDrawable(int what, int color = -1, double step = 0) const;
  private:    
    void CopyParameters(const IMF& imf);

    // Calcultes essencial parameters and queues extra
    void Update();

    // these are required for sifting
    void CalculateIntrinsic();
    void CalculateMean();
    void CalculateVariance();
    
    // these are extra parameters that dont need to be updated in every loop 
    void CalculateFrequency();
    void CountExtrema();
    void CalculateCorrelation(double*);
    
    double SplineMed(double) const;
    double SplineAmp(double) const;
    
    // interpolates the array to Spline<CODE>
    // CODE_MAX -> SplineMax  
    bool CalculateSpline(double*, uint8_t);
    
    vector<int> GetIndexes(double*, uint8_t) const; 
    int GetCount(double*, uint8_t) const;
    
    // vector with y points, x is always constant
    double* intrinsic, *old_intrinsic;
    double* variance;
    double* mean;
    double *frequency;
    
    bool queueCalculateExtras;
    
    Spline3Interpolator *SplineMax, *SplineMin;
    double  correlation;

    int maxCount, minCount, zeroCount;
    
    int index;
    static int indexCount;
};


#endif
