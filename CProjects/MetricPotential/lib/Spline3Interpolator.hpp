#ifndef H_SPLINE3INTERPOLATOR_H
#define H_SPLINE3INTERPOLATOR_H

#include "DataPoints.hpp"
#include "FCmatrixBanded.h"
#include "FCmatrixFull.h"
#include "EqSolver.h"
#include "TF1.h"
#include <string>

using namespace std;

class Spline3Interpolator : public DataPoints {
  public:
    Spline3Interpolator(int fN=0,
			double *fx=NULL,
			double *fy=NULL);
    Spline3Interpolator(const Spline3Interpolator&);

    ~Spline3Interpolator();

    //void Draw(double step = 0);
    double Interpolate(double) const;

    double operator()(double);

    TF1* GetFInterpolator();
    
  private:
    void SetCurvatureLines();
    double fInterpolator(double *fx, double *par);
    
    TF1* FInterpolator;
    double *K;

    int index;
    static int indexCount;
};

#endif
