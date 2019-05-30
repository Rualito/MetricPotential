#include "IntegratorMC.hpp"

#include <cmath>

#define ABS(x) (x>0?x:-x)

//x#define DEBUG
//#define COUNT
#include "Debug.hpp"

//#define BENCHMARK
#include "BenchMark.hpp"
#include <cstdlib>

#define RAND rand()/((double) RAND_MAX)
#define RAND_(min, max) (max-min)*(rand()/((double) RAND_MAX)) + min

IntegratorMC::IntegratorMC(double fx0, double fx1, TF1*fp) : Integrator(fx0, fx1, fp){;}

IntegratorMC::~IntegratorMC(){;}


void IntegratorMC::IntegralMC(int N, double& result, double& error)
{
    ECHOS;
    result = 0;
    error = 0;

    BENCHMARK_START(0);
    for(int i = 0; i<N; ++i){
	double randfx = F->Eval(RAND_(x0, x1)); 
	result += randfx;
	error += randfx*randfx;
    }
    BENCHMARK_END(0);

    result *= (x1-x0)/N;
    // result*N = (xmax-xmin); (xmax-xmin)/sqrt(N) = result*sqrt(N)
    error = (x1-x0)*sqrt(error - N*result*result)/N;
    ECHOF;
}

void IntegratorMC::IntegralMC(TFormula *F,int N, double& result, double& error, double *xmin, double *xmax){

    int npar = F->GetNdim();

    result = 0;
    error = 0;
    
    double randx[npar];
    double val;

    for(int i = 0; i<N; ++i){
	
	for(int j = 0; j<npar; ++j){
	    randx[j] = RAND_(xmin[j], xmax[j]);
	    //printf(">< %f\n", randx[j]);
	}
	
	val = F->EvalPar(randx, NULL);
	//printf("> %f \n", val);
	error += val*val; 
	result += val;
    }
    
    double volume = 1; 

    for(int i = 0; i<npar; i++){
	volume *= xmax[i]-xmin[i];
    }

    result *= volume/N;
    error = volume*sqrt(error - N*result*result)/N;

}
void IntegratorMC::IntegralMCIS(int N, double& result, double& error, TF1* pdf, TF1* iPdf)
{
    ECHOS;
    result = 0;
    error = 0;

    double x, y, randfx;

    BENCHMARK_START(0);
    for(int i = 0; i<N; ++i){
        y = RAND;
	x = iPdf->Eval(y);
	//printf(">>x: %f, F: %f, ipdf %f\n", x, F->Eval(x), iPdf->Eval(x));
	randfx = F->Eval(x)/pdf->Eval(x);
	
	result += randfx;
	error += randfx*randfx;
    }
    BENCHMARK_END(0);
    
    error = sqrt(N*error - result*result)/N;
    
    result /=N;
    
    ECHOF;
}

void IntegratorMC::IntegralMCAR(int N, double& result, double& error, TF1* q, double qIval)
{
    ECHOS;

    result = 0;
    error = 0;
    
    int count = 0;
    double x,u;
    
    double fmax = F->GetMaximum(x0, x1);
    if(!q){
	
	BENCHMARK_START(0);
	for(int i = 0; i<N; ++i){
	    x = RAND_(x0, x1);
	    u = RAND_(0, fmax);
	    
	    if(ABS(u) <= ABS(F->Eval(x))){
		++count;
	    }
	}
	BENCHMARK_END(0);
	
	result = (x1-x0)*fmax/N;
	error = result*sqrt(count*(1-((double)count)/N));
	result *= count;
    } else {

	BENCHMARK_START(1);
	for(int i = 0; i<N; ++i){
	    x = RAND_(x0, x1);
	    u = RAND;
	    
	    if( u*q->Eval(x) <= F->Eval(x) ){
		++count;
	    }
	}
	BENCHMARK_END(1);

	result = (x1-x0)*fmax/N;
	error = result*sqrt(count*(1-((double) count)/N));
	result *= count;
    }
    
    ECHOF;
}
    
    
    
    
    
    
