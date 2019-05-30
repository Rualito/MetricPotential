#ifndef H_INTEGRATORMC_H
#define H_INTEGRATORMC_H

#include "TF1.h"
#include "FCrand.hpp"
#include "Integrator.h"

class IntegratorMC : public Integrator 
{
  public:
    
    IntegratorMC(double fx0 = 0, double fx1= 0, TF1* fp=NULL);
    ~IntegratorMC();
    
    // simple integration
    void IntegralMC(int N, double& result, double& error);
    
    void IntegralMC(TFormula *F,int N, double& result, double& error, double *xmin, double *xmax);
	
    // importance sampling
    void IntegralMCIS(int N, double& result, double& error, TF1* pdf, TF1* iPdf);
    
    //void IntegralMCIS(double* xmin, double* xmax, int N, double& result, double& error, TF1* pdf, TF1* iPdf);
    // acception rejection
    
    void IntegralMCAR(int N, double& result, double& error, TF1* q = NULL, double qIval = 0); // q is the polyfunction, qIval is its integral in [xmin, xmax]
    //void IntegralMCAR(double xmin, double xmax, double& result, double& error);
    
    //void IntegralMCAR(double* xmin, double* xmax, double& result, double& error, TF1* q=NULL);
};

#endif
