#include "ODEsolver.hpp"

//#define DEBUG
//#define COUNT
#include "Debug.hpp"

//#define BENCHMARK
#include "BenchMark.hpp"

ODEsolver::ODEsolver(){;}
      
ODEsolver::ODEsolver(TFormula *f, int N) {    
    ECHOS;
    
    for(int i = 0; i<N; i++){
	F.push_back(f[i]);
    }
    
    ECHOF; 
}

ODEsolver::ODEsolver(vector<TFormula> VF) {
    F = VF;
}

ODEsolver::~ODEsolver() {;}

void ODEsolver::SetODEfunctions(vector<TFormula> VF) {
    F.clear();
    F = VF;
}

vector<ODEpoint> ODEsolver::Euler(const ODEpoint& p0, double xmin, double xmax, int Nstep) {
    ECHOS;
    double step = (xmax-xmin)/Nstep;
    vector<ODEpoint> v;
    COUNT_P(1);
    v.push_back(p0);
    COUNT_P(1);
    for(int i = 1; i<Nstep; i++){
	v.push_back(EulerIterator(v.back(), step));
    }
    ECHOF;
    return v;
}

ODEpoint ODEsolver::EulerIterator(const ODEpoint& a, double step) {
    ECHOS;
    double x[F.size()+1];
    COUNT_P(1);
    x[0] = a[0] + step;
    COUNT_P(1);
    for(int i = 0; i<(int) F.size(); i++){
	COUNT_P(2);
	x[i+1] = a[i+1] + step*F[i].EvalPar(a.GetX(), NULL);
    }
    
    ECHOF;
    return ODEpoint(x, (int) F.size());
}

vector<ODEpoint> ODEsolver::RK2(const ODEpoint& p0, double step, double time) {
    ECHOS;
    
    int Nstep = (int)(time/step);
    vector<ODEpoint> v;
    
    v.push_back(p0);
    
    for(int i = 1; i<Nstep; i++){
	v.push_back(RK2Iterator(v.back(), step));
    }
    
    ECHOF;
    return v;
}

ODEpoint ODEsolver::RK2Iterator(const ODEpoint& a, double step) {
    ECHOS;
    double x[F.size()+1];
    //double K1[F.size()], K2[F.size()];
    double k1[F.size()];
    double k2;
    
    x[0] = a[0] + step;


    double v[F.size()+1];
    v[0] = a[0] + step/2;
    //printf(">> ");
    for(int i = 0; i<(int) F.size(); i++){
	k1[i] = step * F[i].EvalPar(a.GetX(), NULL);
	v[i+1] = a[i+1] + k1[i]/2;
	//printf(" (%f %f)", k1[i], v[i]);
    }

    for(int i=0; i<(int)F.size(); i++){
	k2 = step * F[i].EvalPar(v, NULL);
	x[i+1] = a[i+1] + k2;
	//printf(" [%f %f]", a[i+1], k2);
    }
    //printf("\n");
    
    ECHOF;
    return ODEpoint(x, (int) F.size());
}

vector<ODEpoint> ODEsolver::RK4(const ODEpoint& p0, double step, double time)
{
    ECHOS;
    
    int Nstep = (int)(time/step);
    vector<ODEpoint> v(Nstep+1);
    
    v[0] = p0;
    
    for(int i = 1; i<Nstep; i++){
	v[i] = RK4Iterator(v[i-1], step);
    }
    
    ECHOF;
    return v;
}

ODEpoint ODEsolver::RK4Iterator(const ODEpoint& a, double step)
{
    ECHOS;
    int dim = F.size();
    double K[4][dim];
    
    for(int i= 0; i<F.size(); i++){
	K[0][i] = step * F[i].EvalPar(a.GetX(), NULL);
    }

    for(int j = 1; j<4;j++){
	double vec[dim+1];

	vec[0] = a[0] +  step / (3-(j+1)/2);

	for(int i = 1; i<F.size()+1; i++){
	    vec[i] = a[i] + K[j-1][i] / (3-(j+1)/2);
	}

	for(int i = 0; i < F.size(); i++){
	    K[j][i] = step * F[i].EvalPar(vec, NULL);
	}
    }

    
    double v[dim+1];
    v[0] = a[0] + step;

    //printf("<>--<>--<>\n");
    for(int i = 0 ; i<(int) F.size()+1; i++){
	//printf(">> %12.3e %12.3e %12.3e %12.3e\n", K[0][i], K[1][i], K[2][i], K[3][i]);
	v[i+1] = a[i+1] + (K[0][i] + 2*K[1][i] + 2*K[2][i] + K[3][i])/6;
    }
    //if(COUNT_B(4) > 100) exit(-1);

    ECHOF;
    return ODEpoint(v, dim);
}
