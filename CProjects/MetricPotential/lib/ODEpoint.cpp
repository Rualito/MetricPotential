#include "ODEpoint.hpp"

//#define BENCHMARK
#include "BenchMark.hpp"

//#define DEBUG
//#define COUNT 
#include "Debug.hpp"

ODEpoint::ODEpoint(double* a, int nDIM) : x(a), nDim(nDIM)
{
    ECHOS;
    if(nDim != 0){
	x = new double[nDIM+1];
	
	// +1 because of the independent variable
	for(int i = 0; i<nDIM+1; i++){
	    x[i] = a[i];
	}
    }
  
    ECHOF;
}

ODEpoint::ODEpoint(const ODEpoint& points) 
{
    ECHOS;
    
    if(this != &points){

	if(!points.x){
	    x = NULL;
	    nDim = 0;
	    ECHOF;
	    return;
	}
	
	x = new double[points.nDim+1];
	nDim = points.nDim;
	
	for(int i=0; i<=nDim; i++){
	    x[i] = points.x[i];
	}	
    }

    ECHOF;
}   

ODEpoint::~ODEpoint()
{
    ECHOS;
    
    if(x){
	delete[] x;
    }

    ECHOF;
}

double ODEpoint::operator[](int i) const
{
    ECHOS;
    return x[i];
}

double& ODEpoint::operator[](int i)
{
    ECHOS;
    return x[i];
}

double* ODEpoint::GetX() const
{
    return x;
}

int ODEpoint::dim() const{
    return nDim;
}

const ODEpoint& ODEpoint::operator=(const ODEpoint& p)
{
    ECHOS;
    if(&p != this){
	if(p.nDim != nDim){
	    delete[] x;
	    
	    nDim = p.nDim;
	    x = new double[nDim+1];
	}

	for(int i = 0; i<nDim+1; i++) {
	    x[i] = p[i];
	}
    }

    ECHOF;
    return *this;
}

	

	
    
