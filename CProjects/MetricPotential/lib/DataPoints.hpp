#ifndef H_DATAPOINTS_H
#define H_DATAPOINTS_H

#include "cFCgraphics.h"
#include "TGraph.h"
#include "FCtools.hpp"
#include "FCmatrixFull.h"

class DataPoints
{
  public:
    DataPoints();
    DataPoints(int n, double* xp = NULL, double* yp = NULL);
    virtual ~DataPoints();
    virtual void Draw(double step=0.01);
    virtual double Interpolate(double)const;
    virtual void Print(string filename = "", string format = "%f");
    
//    virtual TGraph& GetPointGraph();
//    virtual TGraph& GetInterpolatedGraph(double step = 0.01);

    virtual double* GetArrayX();
    virtual double* GetArrayY();
    virtual int size() const;
    
    void StartGraphics();
    
  protected:
    bool graphics;
    int N;
    double *x, *y;
    cFCgraphics G;
};

#endif
