#ifndef H_FCRAND_H
#define H_FCRAND_H

#include <cstdlib>
#include "Vec.h"
#include <ctime>
#include "TF1.h"

//#include "benchmark.hpp"

class FCrand
{
  public:
    FCrand(int seed=-1);
    ~FCrand();
    
    double GetRandom(double min=0, double max=1);
    Vec GetRandomArray(int N=0, double* min=NULL, double* max=NULL);

}; 

#endif
