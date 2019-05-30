#include "Vec.h"

//#define DEBUG
//#define COUNT
// channel 10 counts how many times Vec(Vec&) is called
#include "Debug.hpp"
#include "BenchMark.hpp"

using namespace std;

Vec::Vec(int n, double *a) : N(n) {
    ECHOS;
    entries = new double [N];
    
    if(a){
	for (int i = 0; i<N; i++){
	    entries[i] = a[i];
	}
    } else {
	for (int i = 0; i<N; i++){
	    entries[i] = 0;
	}
    }
    
    ECHOF;
    /*cout << endl;
    for(int i = 0; i < N; i++)
    {
        cout << entries[i] << "  ";
    }
    cout << endl;*/
    
    //cout << endl << "Constructor 1" << endl;
}

Vec::Vec(int n, double c) : N(n) {
    ECHOS;
    //cout << endl << "Entrei" << endl;
    if(N != 0)
    {
        entries = new double [N];
        for(int i = 0; i < n; i++)
        {
            entries[i] = c;
        }
    }
    //cout << endl << "Constructor 3" << endl;
    //cout << "Saí!" << endl;
    ECHOF;
}

Vec::Vec(const Vec& p) {
    ECHOS;
    COUNT_B(10);
    if(&p != this)
    {
	//p.Print();
        N = p.N;
        entries = new double [N];
        
        for (int i = 0; i < p.N; i++) {
            entries[i] = p.entries[i];
        }
    }
    //cout << endl << "Copy Constructor" << endl;
    ECHOF;
}

void Vec::SetEntries (int n, double* a) {
    ECHOS;
    N = n;

    if(entries != nullptr)
    {
        delete[] entries;
        entries = nullptr;
    }

    entries = new double [N];
    
    for (int i = 0; i < N; i++) {
        entries[i] = a[i];
    }
    ECHOF;
    //cout << endl << "SetEntries " << endl;
}

string Vec::PrintToString(string format) const
{
    //ECHOS;
    string temp = "";
    try{
	for(int i = 0; i<N; i++){
	    char c[32];
	    sprintf(c, format.c_str(), entries[i]);

	    temp += " " + string(c);
	}
    } catch (const exception& e){
	printf("String size too big, trying again with scientific notation - %%5.2e\n");
	temp = "";
	
	format = "%5.2e";
	
        for(int i = 0; i<N; i++){
	    char c[32];
	    sprintf(c, format.c_str(), entries[i]);

	    temp += " " + string(c);
	}
    }
    

    //ECHOF;
    return temp;
}

void Vec::Print(string format) const{

    // for(int i = 0; i<N; i++){
    // 	char c[24];
    // 	sprintf(c, format.c_str(), entries[i]);

    // 	printf()
    // }
    cout << PrintToString(format) << '\n';
}
    
Vec Vec::operator = (const Vec& V0) {

    if(this != &V0)
    {
        N = V0.N;

        if(entries!= nullptr)
        {
            delete[] entries;
            entries = nullptr;
        }
	N = V0.N;
	    
	entries = new double [N];
	
        for (int i = 0; i < N; i++) 
        {
            entries[i] = V0.entries[i];
        }
    }

    return *this;
}

double& Vec :: operator [] (int i) const {

    return entries[i];
}

Vec Vec::operator + (const Vec& V0) {

    Vec V(*this);

    if(V0.size() == (*this).size())
    {
        for(int i = 0; i < N; i++)
        {
            V.entries[i] += V0.entries[i];
        }
    }
    else
    {
        cout << endl << "Os vetores não têm dimensões apropriadas para a soma!" << endl;
    }

    return V;
}

Vec Vec::operator += (const Vec& V0) {

    if(V0.size() == (*this).size())
    {
        for(int i = 0; i < N; i++)
        {
            entries[i] += V0.entries[i];
        }
    }
    else
    {
        cout << endl << "Os vetores não têm dimensões apropriadas para a soma!" << endl;
    }

    return *this;
}

Vec Vec::operator - (const Vec& V0) {

    Vec V(*this);

    if(V0.size() == (*this).size())
    {
        for(int i = 0; i < N; i++)
        {
            V.entries[i] -= V0.entries[i];
        }
    }
    else
    {
        cout << endl << "Os vetores não têm dimensões apropriadas para a subtração!" << endl;
    }

    return V;
}

Vec Vec::operator -= (const Vec& V0) {

    if(V0.size() == (*this).size())
    {
        for(int i = 0; i < N; i++)
        {
            entries[i] -= V0.entries[i];
        }
    }
    else
    {
        cout << endl << "Os vetores não têm dimensões apropriadas para a subtração!" << endl;
    }

    return *this;
}

Vec Vec::operator - () {

    Vec V0 = (*this);
    for(int i = 0; i < N; i++)
    {
        V0.entries[i] = -entries[i];
    }

    return V0;
}

Vec Vec::operator + () {

    return *this;
}

Vec Vec::operator * (const Vec& V0) {

    Vec V(*this);

    if(V0.size() == (*this).size())
    {
        for(int i = 0; i < N; i++)
        {
            V.entries[i] *= V0.entries[i];
        }
    }
    else
    {
        cout << endl << "Os vetores não têm dimensões apropriadas para a multiplicação escalar!" << endl;
    }

    return V;
}

Vec Vec::operator * (double lambda) {

    Vec V(*this);

    for(int i = 0; i < N; i++)
    {
        V.entries[i] *= lambda;
    }

    return V;
}

int Vec::size() const {

    return N;
}

double Vec::dot(const Vec& V0) {

    if( N != V0.N)
    {
        cout << endl << "As dimensões dos vetores são inapropriadas para o produto interno!" << endl;
    }

    double mul = 0;

    for(int i = 0; i < N; i++)
    {
        mul += entries[i]*V0.entries[i];
    }
    //cout << mul << endl;

    return mul;
}

void Vec::swap(int i, int j) {

    double x = entries[i];
    entries[i] = entries[j];
    entries[j] = x;
}

void Vec::Seti(int i, double d) {

    entries[i] = d;
}

Vec::~Vec() {
    
    if (entries != nullptr)
    {
        delete[] entries;
    }  
}
