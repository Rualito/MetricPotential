#ifndef H_ODESOLVER_H
#define H_ODESOLVER_H

#include <vector>
#include "TFormula.h"
#include "FCmatrixFull.h"
#include "ODEpoint.hpp"

using namespace std;

class ODEsolver
{
  public:
    // functions and how many
    ODEsolver();
    ODEsolver(TFormula*, int N = 0);
    ODEsolver(vector<TFormula>);
    ~ODEsolver();

    void SetODEfunctions(vector<TFormula>);

    // x range and number of subintervals 
    vector<ODEpoint> Euler(const ODEpoint& p0, double xmin, double xmax, int Nstep);
    vector<ODEpoint> RK2(const ODEpoint& p0, double step, double time);
    vector<ODEpoint> RK4(const ODEpoint& p0, double step, double time);
    

  private:
    vector<TFormula> F;
    ODEpoint EulerIterator(const ODEpoint&, double step);
    ODEpoint RK2Iterator(const ODEpoint&, double step);
    ODEpoint RK4Iterator(const ODEpoint&, double step);
};

#endif
