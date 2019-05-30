#include "NevilleInterpolator.hpp"
#define DEBUG

NevilleInterpolator::NevilleInterpolator(int N, double *xs, double *ys): DataPoints(N, xs, ys)
{
    #ifdef DEBUG
    printf("\n[NevilleInterpolator::NevilleInterpolator] START\n");
    #endif
    
    FInterpolator = new TF1("fInterpolator", this, &NevilleInterpolator::fInterpolator, -0.1, 3.1, 0, "NevilleInterpolator", "fInterpolator");
    
    F0 = NULL;
    
    #ifdef DEBUG
    //DataPoints::Print();
    printf("\n[NevilleInterpolator::NevilleInterpolator] END\n");
    #endif
}

void NevilleInterpolator::SetFunction(TF1* f)
{
    (*F0) = (*f);
}

double NevilleInterpolator::Interpolate(double xv)
{
    double *yp = new double[N];
    
    for(int i = 0; i<N; i++)
    {
	    yp[i] = y[i];
    }

    for(int k = 1; k<N; k++)
    {
        for(int i  = 0; i<N-k; i++)
        {
            yp[i] = ((xv-x[i+k])*yp[i] - (xv-x[i])*yp[i+1])/(x[i]-x[i+k]);
        }
    }

    double A = yp[0];

    delete[] yp;

    return A;
}
	
void NevilleInterpolator::Draw(double step)
{
    DataPoints::Draw(step);
}

TF1* NevilleInterpolator::GetInterpolationFunction()
{
    return FInterpolator;
}

// string PrintToString(string format)
// {
//     return "";
// }

void NevilleInterpolator::Print(string filename, string format)
{
    DataPoints::Print(filename, format);
}

double NevilleInterpolator::fInterpolator(double *fx, double *par)
{
    return Interpolate(fx[0]);
}

NevilleInterpolator::~NevilleInterpolator()
{
    delete FInterpolator;
    delete F0;
    
}


