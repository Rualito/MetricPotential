#include "Spline3Interpolator.hpp"

#define POW3(x) ((x)*(x)*(x))
#define ABS(x) (x>0? (x) :-(x))

//#define BENCHMARK
#include "BenchMark.hpp"

//#define DEBUG
#include "Debug.hpp"

//#define DEBUG_VAL

int Spline3Interpolator::indexCount = 0;

Spline3Interpolator::Spline3Interpolator(int fN, double *fx, double *fy) : DataPoints(fN,fx,fy) {
    ECHOS;
    //DataPoints::Print();
    //F0=NULL;
    
    index = indexCount++;
   
    string name = "SplineInterpoltator ";  
    name+=index;

    if(fN && fx && fy){
	FInterpolator = new TF1(name.c_str(), this,
				&Spline3Interpolator::fInterpolator,
				x[0]-0.1 ,x[N-1]+0.1, 0);
	K = new double[N];
    
	SetCurvatureLines(); //define segment interpolators
    }
    ECHOF;
}

Spline3Interpolator::Spline3Interpolator(const Spline3Interpolator& spl) : DataPoints(spl.N, spl.x, spl.y){
    ECHOS;
    FInterpolator = new TF1(*spl.FInterpolator);

    index = indexCount++;
    
    K = new double[N];
    
    for(int i = 0; i<N; i++){
	K[i] = spl.K[i];
    }

    ECHOF;
}

void Spline3Interpolator::SetCurvatureLines() {
    // define tri-diagonal matrix and axrray of constants
    ECHOS;
    // diagonal size
    int n = N-2;

    //printf("mat size >> %d \n", n);
    
    // diagonal
    double diag[n];
    
    // diagonal sides
    double up[n-1]; // up

    double low[n-1]; // down
    
    // b vector
    double *b = new double[n];
 
    BENCHMARK_START(0);
    
    // the number of elements in x is diferent from K
    for(int i = 0; i<n; i++){ // loop on rows
	
	double up_i = x[i+1]-x[i+2];  // upper diagonal element
	double diag_i = 2.*(x[i]-x[i+2]); //diag
	double low_i = x[i]-x[i+1]; // lower diagonal element
	
	if(i < n-1){
	    // upper diagonal
	    up[i] = up_i;
	}
	if(i > 0){
	    //lower diagonal
	    low[i-1] = low_i;
	} 

	diag[i] = diag_i;
	
	b[i] = 6 * ((y[i]-y[i+1])/(x[i]-x[i+1])
		   - (y[i+1]-y[i+2])/(x[i+1]-x[i+2]));
    }
    BENCHMARK_END(0);
	
    vector<Vec> vv;
    
    Vec B(n, b);
    
    vv.push_back(Vec(n-1, up));
    vv.push_back(Vec(n, diag));
    vv.push_back(Vec(n-1, low));
    
    FCmatrixBanded mat(vv);
    BENCHMARK_END(0);
    //FCmatrixBanded bM(vv);
    
    EqSolver eq(mat, B);
    BENCHMARK_END(0);
    Vec s = eq.TridiagonalSolver();
    BENCHMARK_END(0);
    #ifdef DEBUG_VAL
    printf("vec\n");
    
    B.Print();
    
    printf("Banded: \n");
    
    mat.PrintFull();

    printf("Solution: \n");
    s.Print();
    #endif
    
    for(int i = 0; i<n; i++){
	K[i+1] = s[i];
    }
    
    K[0] = 0;
    K[n+1] = 0;
    
    BENCHMARK_END(0);
    
    ECHOF;
}
double Spline3Interpolator::Interpolate(double fx) const {
// detect in wich segment is x
    //return 0;
    int i;
    for (i=0; i<=N; ++i) {
	if ((fx-x[i])<0.) break;
    } //upper bound returned
    if (i==0 || i==N+1){ // out of range
	return 0.;
    }

    i--;

    double val = (K[i]/6.) * ( POW3(fx-x[i+1])/(x[i]-x[i+1]) - (fx-x[i+1])*(x[i] - x[i+1]))
      - (K[i+1]/6.) * ( POW3(fx-x[i])/(x[i]-x[i+1]) - (fx-x[i])*(x[i]-x[i+1])) +
	(y[i]*(fx-x[i+1]) - y[i+1]*(fx-x[i]))/(x[i]-x[i+1]); 

    if(isnan(val)){
        return 0.;
    }
    return val; 
}

double Spline3Interpolator::operator()(double fx){
    return Interpolate(fx);
}

TF1* Spline3Interpolator::GetFInterpolator()
{
    return FInterpolator;
}

Spline3Interpolator::~Spline3Interpolator()
{
    ECHOS;
    if(FInterpolator)
	delete FInterpolator;
    if(K)
	delete[] K;
    ECHOF;
}

double Spline3Interpolator::fInterpolator(double *fx, double *par)
{
    return Interpolate(fx[0]);
}
/*
void Spline3Interpolator::Draw(double step)
{
    ECHOS;
    DataPoints::Draw(step);
    ECHOF;
}

*/    
    
