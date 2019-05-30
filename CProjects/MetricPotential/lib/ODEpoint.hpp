#ifndef H_ODEPOINT_H
#define H_ODEPOINT_H

#include <cstdio>
#include <cstdlib>

using namespace std;

class ODEpoint
{
  public:
    // nDim = number of dependent vars, a[nDim+1]
    ODEpoint(double*a=NULL, int nDIM=0);
    ODEpoint(const ODEpoint&);
    ~ODEpoint();

    double operator[](int) const;
    double& operator[](int);

    const ODEpoint& operator=(const ODEpoint&);
    double* GetX() const;

    int dim() const;
    
  private:
    double *x;
    int nDim;
    /* x[0] indep var
     * x[1] first dep var
     * ...
     * x[nDim] last dep var
     */
    
};

#endif
    
