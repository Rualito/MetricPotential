#include "Integrator.h"

//#define DEBUG
#include "Debug.hpp"

//#define BENCHMARK
#include "BenchMark.hpp"

Integrator::Integrator(double fx0, double fx1, TF1* fp) : x0(fx0), x1(fx1){

    if(fp){
	F = new TF1(*fp);
    } else {
	F = new TF1("nullfunc", "0");
    }
    
}

void Integrator::SetIntegrandFunction(TF1* fp)
{
    ECHOS;
    if(fp){
	if(F){
	    delete F;
	}

	F = new TF1(*fp);
    }
    ECHOF;
}

Integrator::~Integrator() {
    ECHOS;
    delete F;
    ECHOF;
}

void Integrator::TrapezoidalRule(int n, double& Integral, double& Error, double x0p, double x1p) {
    ECHOS;

    if(x0p == 0 && x1p == 0){
	x0p = x0;
	x1p = x1;
    }
    
    Integral = F->Eval(x0p);
    Error = 0;
    double h = (x1p - x0p)/n;
    
    for(int i = 1; i < n; i++)
    {
        Integral += 2*F->Eval(x0p + i*h);
        Error += -(h*h*h/12) * F->Derivative2(x0p + (2*i-1)*h/2);
    }
    
    Integral += F->Eval(x1p);
    Integral *= h/2;

    ECHOF;
}

void Integrator::SimpsonRule(int n, double& Integral, double& Error, double x0p, double x1p) {
    ECHOS;
    
    if(x0p == 0 && x1p == 0){
	x0p = x0;
	x1p = x1;
    }
    
    Integral = F->Eval(x0p);
    Error = 0;
    double h = (x1p - x0p)/n;
    int teste = 0;
    
    if(n%2 == 1)
    {
        n = n-1;
        teste = 1;
    }

    for(int i = 1; i < n; i++)
    {
        if(i%2 == 1)
        {
            Integral += 4*F->Eval(x0p + i*h);
        }
        else
        {
            Integral += 2*F->Eval(x0p + i*h);
        }
	//printf(">> int %f\n", Integral);
    }

    double xi, xf;
    for(int i = 0; i < n/2; i++)
    {
        xi = x0p + 2*i*h;
        xf = x0p + 2*(i+1)*h;
        Error += -(h*h*h*h*h/90) * (F->Derivative3(xf) + F->Derivative3(xi))/(xf - xi);
    }

    if(teste == 1)
    {
        Integral += F->Eval(x0 + (n-1)*h) + F->Eval(x0 + n*h);
        xi = x0 + 2*n*h;
        xf = x0 + 2*(n+1)*h;
        Error += -(h*h*h*h*h/90) * (F->Derivative3(xf) + F->Derivative3(xi))/(xf - xi);
    }
    else
    {
        Integral += F->Eval(x1p);
    }
    
    Integral *= h/3;

    ECHOF;
}
