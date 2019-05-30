#include "FCrand.hpp"

//#define BENCHMARK

FCrand::FCrand(int seed)
{
    if(seed < 0){
	srand(time(NULL));
    } else{
	srand(seed);
    }
    
}

FCrand::~FCrand()
{}   

double FCrand::GetRandom(double min, double max){
    return (max-min)*(rand()/((double) RAND_MAX)) + min;
}

Vec FCrand::GetRandomArray(int N, double* min, double* max)
{
    double v[N];
    if(min && max){
	for(int i = 0; i<N; i++){
	    v[i] = GetRandom(min[i], max[i]);
	}
    } else{
	for(int i = 0; i<N; i++){
	    v[i] = GetRandom(0, 1);
	}
    }

    return Vec(N, v);
}



