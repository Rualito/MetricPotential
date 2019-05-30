#ifndef H_INTEGRATOR_H
#define H_INTEGRATOR_H

#include "TF1.h"

class Integrator 
{
    public:
        Integrator(double fx0=0, double fx1=0, TF1* fp=NULL);
        ~Integrator();

	void SetIntegrandFunction(TF1*);
	
	void TrapezoidalRule(int, double&, double&, double x0p = 0, double x1p = 0);
        void SimpsonRule(int, double&, double&, double x0p = 0, double x1p = 0);
    
    protected:
        double x0;
        double x1;
        TF1* F;
};

#endif
