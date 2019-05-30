#include "FCmatrix.h"
#include "FCmatrixBanded.h"
#include "FCmatrixFull.h"
#include <vector>
#include "Vec.h"

//#define DEBUG
#include "Debug.hpp"

using namespace std;

#define ABS(x) ((x)>0?(x) :(-x))

FCmatrixBanded::FCmatrixBanded(vector<Vec> M) : FCmatrix(M) {
    
    if(M.size() != 3)
    {
        cout << "A matriz não foi bem construída!" << endl 
             << "A matriz Banded é cosntruída com recurso a 3 vetores que correspondem âs entradas das três diagonais centrais da matriz!" << endl;
    }

    if(M[0].size() != M[2].size() || M[0].size() != M[1].size() - 1)
    {
        cout << "A matriz não foi bem construída!" << endl << "A matriz Banded é quadrada e cosntruída com recurso a 3 vetores:" << endl
             << "Os Vecs a e c têm de ter a mesma dimensão, igual a n-1 elementos, e o Vec b tem de ter n elementos!" << endl << endl;
    }

    classname = "FCmatrixBanded";
}

FCmatrixBanded::FCmatrixBanded(const FCmatrix& Mx) {
    ECHOS;
    // this casts any matrix into a BandedMatrix, no matter if it was originally banded
    classname = "FCmatrixBanded";
    
    int N = Mx.GetFM();
    Vec v0 (N-1, 0.);
    Vec v1 (N, 0.);
    Vec v2 (N-1, 0.);
    vector<Vec> mm = Mx.GetRawM();    
    
    if(Mx.Getclassname() == "FCmatrixBanded"){
        M.push_back(mm[0]);
        M.push_back(mm[1]);
        M.push_back(mm[2]);
    }
    
    if(Mx.GetFM() != Mx.GetFN()){
	throw runtime_error("Wrong matrix dimensions");
    }
    
    if(Mx.Getclassname() == "FCmatrixFull")
    {
	//printf("Checking for FULL\n");
        //ver se é possível converter a Full para Banded
        //Não pode ter elementos não nulo fora das três diagonais centrais

        for(int i=0; i < N-1; i++)
        {
            for(int j=i+2; j < N-1; j++)
            {
                if(Mx.GetM()[i][j] != 0)
                {
                    cout << endl << "A FCmatrixFull fornecida tem valores não nulos foras das três diagonais centrais!"
                                 << endl << "A matriz FCmatrixBanded não pode, por isso, er construída!" << endl;
                    exit(1);
                }
            }

            for(int j=0; j < i-1; j++)
            {
                if(mm[i][j] != 0)
                {
                    cout << endl << "A FCmatrixFull fornecida tem valores não nulos foras das três diagonais centrais!"
                                 << endl << "A matriz FCmatrixBanded não pode, por isso, er construída!" << endl;
                    exit(1);
                }
            }            
        }
	    
        for(int i=0; i < N-1; i++)
        {
            v0.Seti(i, mm[i][i+1]);
            v2.Seti(i, mm[i+1][i]);

            v1.Seti(i, mm[i][i]);
        }

        v1.Seti(N-1, mm[N-1][N-1]);
        M.push_back(v0);
        M.push_back(v1);
        M.push_back(v2);
    }
    ECHOF;
}

int FCmatrixBanded::GetFM() const {

    return M[1].size();
}

int FCmatrixBanded::GetFN() const{

    return M[1].size();
}

Vec FCmatrixBanded::GetVec(int i) const {

    if(0 <= i && i <= 2)
    {
        return M[i];
    }
    else
    {
        cout << endl << "Só existem 3 Vecs!" << endl;
        Vec V;
        return V;
    }
}

vector<Vec> FCmatrixBanded::GetM() const {

    vector<Vec> M0;

    for(int i = 0; i < GetFM(); i++)
    {
      M0.push_back(GetRow(i));
    }

    return M0;
}

Vec FCmatrixBanded::GetRow(int i) const {

    if(i < 0 || i >= M[1].size())
    {
        cout << endl << "Essa linha não existe na matriz!" << endl;
        Vec V;
        return V;
    }
    else
    {
        Vec row( (int) M[1].size(), 0.);

        for(int j = 0; j < row.size(); j++)
        {
	  //row.Seti(j, (ABS(j-i)<=1?M[i-j+1][(j==i+1?j-1:j]:0));
	  // Does the same as this (in theory)
            if(j == i)
            {
                row.Seti(j, M[1][j]);
            }
            else if(j == i-1)
            {
                row.Seti(j, M[2][j]);
            }
            else if(j == i+1)
            {
                row.Seti(j, M[0][j-1]);
            }
            else
            {
                row.Seti(j, 0);
	    }   
        }    

        return row;
    }   
}

Vec FCmatrixBanded::GetCol(int i) const {

    if(i < 0 || i >= M[1].size())
    {
        cout << endl << "Essa coluna não existe na matriz!" << endl;
        Vec V;
        return V;
    }
    else
    {
        Vec col( (int) M[1].size(), 0.);

        for(int j = 0; j < col.size(); j++)
        {
	  //col.Seti(j, ( ABS(j-i)<=1 ? M[j-i+1][ (j==i+1?j-1:j) ] : 0) );

	  //Does the same as this
            if(j == i)
            {
                col.Seti(j, M[1][j]);
            }
            else if(j == i-1)
            {
                col.Seti(j, M[0][j]);
            }
            else if(j == i+1)
            {
                col.Seti(j, M[2][j-1]);
            }  
            else
            {
                col.Seti(j, 0);
	    }       
	}    

        return col;
    }   
}

void FCmatrixBanded::Seti(int i, int j, double ent) {

    if(i < 0 || i >= 3 || j < 0 || j >= (*this).GetFN() )
    {
        cout << endl << "Os valores introduzidos são inadequados!" << endl;
    }
    else
    {
        M[i].Seti(j, ent);
    }
}

double FCmatrixBanded::GetRowMax(int i) {

    if(i < 0 || i >= M[1].size())
    {
        cout << endl << "Essa linha não existe na matriz!" << endl;
        return 0;
    }
    
    double MAX = M[1][i];

    if(i - 1 >= 0 && M[2][i-1] > MAX)
    {
        MAX = M[2][i-1];
    }
    
    if(M[0][i] > MAX)
    {
        MAX = M[0][i];
    }

    return MAX;
}

double FCmatrixBanded::GetColMax(int i) {

    if(i < 0 || i >= M[1].size())
    {
        cout << endl << "Essa coluna não existe na matriz!" << endl;
        return 0;
    }
    
    double MAX = M[1][i];

    if(i - 1 >= 0 && M[0][i-1] > MAX)
    {
        MAX = M[0][i-1];
    }
    
    if(M[2][i] > MAX)
    {
        MAX = M[2][i];
    }

    return MAX;
}

void FCmatrixBanded::PrintVecs() {

    cout << endl;
    cout << "Vec 0: ";
    for(int i = 0; i < M[0].size(); i++)
    {
        cout << M[0][i] << "  ";
    }
    cout << endl;
    cout << "Vec 1: ";
    for(int i = 0; i < M[1].size(); i++)
    {
        cout << M[1][i] << "  ";
    }
    cout << endl;
    cout << "Vec 2: ";
    for(int i = 0; i < M[2].size(); i++)
    {
        cout << M[2][i] << "  ";
    }
    cout << endl << endl;
}

void FCmatrixBanded::PrintFull() {

    cout<<endl;
    Print();
    cout<<endl;
}

Vec& FCmatrixBanded::operator[] (int n) {

    if(n < 0 || n >= (*this).GetFM())
    {
        cout << endl << "Essa linha não existe!" <<  endl;
        exit(1);
    }
    else
    {
        Vec *V = new Vec((*this).GetRow(n));
        return *V;
    }
}

double FCmatrixBanded::Determinant() {
    
    FCmatrixFull M0(*this);
    return M0.Determinant();;
}

FCmatrixBanded FCmatrixBanded::operator=(const FCmatrix& Mx) {

    if(&Mx != this)
    {
        M.clear();
        // this casts any matrix into a BandedMatrix, no matter if it was originally banded
        classname = "FCmatrixBanded";
        
        int N = Mx.GetFM();
        Vec v0 (N-1, 0.);
        Vec v1 (N, 0.);
        Vec v2 (N-1, 0.);

        if(Mx.Getclassname() == "FCmatrixBanded")
        {
            for(int i=0; i < N-1; i++)
            {
                v0.Seti(i, Mx.GetM()[i][i+1]);
                v2.Seti(i, Mx.GetM()[i+1][i]);

                v1.Seti(i, Mx.GetM()[i][i]);
            }

            v1.Seti(N-1, Mx.GetM()[N-1][N-1]);
            M.push_back(v0);
            M.push_back(v1);
            M.push_back(v2);
        }

        if(Mx.GetFM() != Mx.GetFN()){
        throw runtime_error("Wrong matrix dimensions");
        }
        
        if(Mx.Getclassname() == "FCmatrixFull")
        {
            //ver se é possível converter a Full para Banded
            //Não pode ter elementos não nulo fora das três diagonais centrais

            for(int i=0; i < N-1; i++)
            {
                for(int j=i+2; j < N-1; j++)
                {
                    if(Mx.GetM()[i][j] != 0)
                    {
                        cout << endl << "A FCmatrixFull fornecida tem valores não nulos foras das três diagonais centrais!"
                                    << endl << "A matriz FCmatrixBanded não pode, por isso, er construída!" << endl;
                        exit(1);
                    }
                }

                for(int j=0; j < i-1; j++)
                {
                    if(Mx.GetM()[i][j] != 0)
                    {
                        cout << endl << "A FCmatrixFull fornecida tem valores não nulos foras das três diagonais centrais!"
                                    << endl << "A matriz FCmatrixBanded não pode, por isso, er construída!" << endl;
                        exit(1);
                    }
                }            
            }
            
            for(int i=0; i < N-1; i++)
            {
                v0.Seti(i, Mx.GetM()[i][i+1]);
                v2.Seti(i, Mx.GetM()[i+1][i]);

                v1.Seti(i, Mx.GetM()[i][i]);
            }

            v1.Seti(N-1, Mx.GetM()[N-1][N-1]);
            M.push_back(v0);
            M.push_back(v1);
            M.push_back(v2);
        }
    }


    return (*this);
}

FCmatrixBanded FCmatrixBanded::operator+ (const FCmatrixBanded& Mx) {

    if( M[1].size() != Mx.GetFM() )
    {
        cout << "As matrizes não têm dimensões apropriadas para a soma!" << endl;
        return *this;
    }

    FCmatrixBanded M0 = *this;

    for(int i = 0; i < M[1].size(); i++)
    {
        M0[1].Seti(i, M[1][i] + Mx.GetVec(1)[i]);
    }

    for(int i = 0; i < M[1].size()-1; i++)
    {
        M0[0].Seti(i, M[0][i] + Mx.GetVec(0)[i]);
        M0[2].Seti(i, M[2][i] + Mx.GetVec(2)[i]);
    }

    return M0;  
}

FCmatrixBanded FCmatrixBanded::operator-(const FCmatrixBanded& Mx) {
    
    if( M[1].size() != Mx.GetFM() )
    {
        cout << "As matrizes não têm dimensões apropriadas para a subtração!" << endl;
        return (*this);
    }

    FCmatrixBanded M0 = *this;

    for(int i = 0; i < M[1].size(); i++)
    {
        M0[1].Seti(i, M[1][i] - Mx.GetVec(1)[i]);
    }

    for(int i = 0; i < M[1].size()-1; i++)
    {
        M0[0].Seti(i, M[0][i] - Mx.GetVec(0)[i]);
        M0[2].Seti(i, M[2][i] - Mx.GetVec(2)[i]);
    }

    return M0;  
}

FCmatrixFull FCmatrixBanded::operator*(const FCmatrix& Mx) {

    vector<Vec> M0;

    if( M[1].size() != Mx.GetFM() || M[1].size() != Mx.GetFN() )
    {
        cout << "As matrizes não têm dimensões apropriadas para a multiplicação!" << endl;
        M0 = (*this).GetM();
        FCmatrixFull NEW(M0);
        return NEW;
    }

    double mul;
    Vec V0((*this).GetFM(), 0.);
    for(int i = 0; i < (*this).GetFM(); i++)
    {
        for(int j = 0; j < (*this).GetFM(); j++)
        {
            mul = (*this).GetRow(i).dot(Mx.GetCol(j));
            V0.Seti(j, mul);
        }
        M0.push_back(V0);
    }

    //cout << endl << "Hello" << endl;
    FCmatrixFull NEW(M0);

    #ifdef DEBUG
    printf("\n[BandedMatrix::operator* FCmatrix] END\n\n");
    #endif
    
    return NEW;
}

FCmatrixFull FCmatrixBanded::operator*(const Vec& V) {

    vector<Vec> M0;

    if (V.size() != (*this).GetFM())
    {
        M0 = (*this).GetM();
        FCmatrixFull NEW(M0);
        
        cout << "A matriz e o vetor não têm dimensões que permitam a multiplicação!" << endl;
        
        return NEW;
    }

    Vec V0(1, 0.);
    double mul;
    for(int i = 0; i < (*this).GetFM(); i++)
    {
        mul = (*this).GetRow(i).dot(V);
        V0.Seti(0, mul);
        M0.push_back(V0);        
    }

    FCmatrixFull NEW(M0);

    #ifdef DEBUG
    printf("\n[BandedMatrix::operator* Vec] END\n\n");
    #endif
    
    return NEW;
}

FCmatrixBanded FCmatrixBanded::operator*(double lambda) {

    vector<Vec> M0;

    for(int i = 0; i < 3; i++)
    {
        M0.push_back((*this).GetVec(i) * lambda);
    }

    FCmatrixBanded NEW(M0);
    
    return NEW;
}

FCmatrixBanded FCmatrixBanded::Transpose() {

    double c;
    for(int i = 0; i < (*this).GetFM() - 1; i++)
    {
        c = M[0][i];
        M[0].Seti(i, M[2][i]);
        M[2].Seti(i, c);
    }

    return *this;
}

FCmatrixFull FCmatrixBanded::BandedToFCmatrixFull() {

    vector<Vec> M0;

    M0 = (*this).GetM();
    FCmatrixFull NEW(M0);
    return NEW;
}

string FCmatrixBanded::PrintToString(string format) const {
    string temp = "";
    for(int i = 0; i<GetFM(); i++){
        temp += GetRow(i).PrintToString(format) + "\n";
    }
    
    return temp;
}

void FCmatrixBanded::Print(string format) const {
    for(int i = 0; i<GetFM(); i++){
	GetRow(i).Print(format);
    }
}

FCmatrixBanded::~FCmatrixBanded() {

    M.clear();
    //cout << endl << "Destructor de FCmatrixFull" << endl;
}
