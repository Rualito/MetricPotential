#include "FCmatrixFull.h"
using namespace std;

//#define DEBUG
#include "Debug.hpp"
#include "BenchMark.hpp"

FCmatrixFull::FCmatrixFull(double** fM, int fm, int fn) : FCmatrix(fM, fm, fn) {

    ECHOS;
    
    classname = "FCmatrixFull";
    rowindices = new int [fm];
    colindices = new int [fn];

    for(int i = 0; i < fm; i++)
    {
        rowindices[i] = i;
    }

    for(int i = 0; i < fn; i++)
    {
        colindices[i] = i;
    }

    ECHOF;
} 

FCmatrixFull::FCmatrixFull(double* fM, int fm, int fn) : FCmatrix(fM, fm, fn) {
    ECHOS;
    
    classname = "FCmatrixFull";

    rowindices = new int [fm];
    colindices = new int [fn];

    for(int i = 0; i < fm; i++)
    {
        rowindices[i] = i;
    }

    for(int i = 0; i < fn; i++)
    {
        colindices[i] = i;
    }

    ECHOF;
}

FCmatrixFull::FCmatrixFull(vector<Vec> V) : FCmatrix(V) {
    ECHOS;
    classname = "FCmatrixFull";

    int fm = V.size();
    int fn = (V[0]).size();
    
    rowindices = new int [fm];
    colindices = new int [fn];

    for(int i = 0; i < fm; i++)
    {
        rowindices[i] = i;
    }

    for(int i = 0; i < fn; i++)
    {
        colindices[i] = i;
    }

    ECHOF;
}

FCmatrixFull::FCmatrixFull(const FCmatrix& Mx) : FCmatrix(Mx.GetM()) {
    ECHOS;
    classname = "FCmatrixFull";

    int fm = Mx.GetFM();
    int fn = Mx.GetFN();

    rowindices = new int [fm];
    colindices = new int [fn];

    for(int i = 0; i < fm; i++)
    {
        rowindices[i] = i;
    }

    for(int i = 0; i < fn; i++)
    {
        colindices[i] = i;
    }

    ECHOF;
}    

int FCmatrixFull::GetFM() const {

    return M.size();
}

int FCmatrixFull::GetFN() const {

    return M[0].size();
}

Vec FCmatrixFull::GetRow(int i) const {
    ECHOS;
    
    if(i >= M.size() || i < 0)
    {
        cout << endl << "Essa linha não existe na matriz!" << endl;
        Vec V;
	ECHOF;
	return V;
    }
    else
    {
	ECHOF;
        return M[rowindices[i]];
    }
}

Vec FCmatrixFull::GetCol(int i) const {
    ECHOS;
    Vec V;
    double a[M.size()];

    if(i >= (*this).GetFN() || i < 0)
    {
        cout << endl << "Essa coluna não existe na matriz!" << endl;
	ECHOF;
        return V;
    }
    
    for(int j = 0; j < M.size(); j++){
	a[j] = (M[j])[colindices[i]];
    }
    
    V.SetEntries(M.size(), a);
    ECHOF;
    return V;
}

vector<Vec> FCmatrixFull::GetM() const {
    
    return M;
}

double FCmatrixFull::GetRowMax(int i) {

    if(i >= M.size() || i < 0)
    {
        cout << endl << "Essa linha não existe na matriz!" << endl;
        return 0;
    }
    
    
    i = rowindices[i];
    int SIZE = (*M.begin()).size();
    double MAX = (M[i])[0];

    for(int k = 0; k < SIZE; k++)
    {
	if(MAX < (M[i])[k])
	{
	    MAX = (M[i])[k];
     	}
    }
    
    return MAX;
}

double FCmatrixFull::GetColMax(int j) {

    if(j >= (*this).GetFN() || j < 0)
    {
        cout << endl << "Essa coluna não existe na matriz!" << endl;
        return 0;
    }
    else
    {
        j = colindices[j];
        double MAX = (M[0])[j];

        for(int i = 0; i < M.size(); i++)
        {
            if(MAX < (M[i])[j])
            {
                MAX = (M[i])[j];
            }
        }

        return MAX;
    }
}

double FCmatrixFull::Determinant() {

    ECHOS;
    BENCHMARK_START(0);
    //Determinantes específicos para n = 1, 2, 3
    /*double det = 0;
    double det1 = 1;
    double det2 = 1;
    int nlines = M.size();

    //cout << endl << "nlines = " << nlines << endl;
    for(int i = 0; i < nlines; i++)
    {
        int j = i;
        int k = i;

        vector<Vec>::iterator vecit;
        for(vecit = M.begin(); vecit != M.end(); vecit++)
        {
            det1 = det1 * (*vecit)[j];
            det2 = det2 * (*vecit)[k];
            j++;
            k--;

            if(j > (nlines - 1))
            {
                j = 0;
            }

            if(k < 0)
            {
                k = nlines - 1;
            }
        }

        det = det + det1 - det2;
        det1 = 1;
        det2 = 1;
        //cout << endl << det << endl;
    }

    if(nlines == 1)
    {
        det = (*M.begin())[0];
    }

    if(nlines == 2)
    {
        vector<Vec>::iterator vecit1, vecit2;
        vecit1 = M.begin();
        vecit2 = vecit1++;
        det = ((*vecit1)[0])*((*vecit2)[1]) - ((*vecit1)[1])*((*vecit2)[0]);
    }*/

    //Eliminação de Gauss
    FCmatrixFull M0(M);
    int l;
    int init;
    int col = 0;
    int GLOBAL = 0;
    int teste = 0;
    int count = 0;
    double CO;

    while(1)
    {
        teste = 0;

        //Ordenação das linhas
        
        init = GLOBAL;
        l = GLOBAL;
        for(int k = GLOBAL; k < M0.GetFN(); k++)
        {
            for(int i = init; i < M0.GetFM(); i++)
            {
                if(M0[i][k] != 0)
                {
                    //troca as linhas da matriz
                    M0.swapRows(l, i);

                    if(l != i)
                    {
                        count++;
                    }

                    l++;
                }
            }

            init += l+1;
            l = init - 1;

            if(init >= M0.GetFM())
            {
                break;
            }
        }

        //cout << endl << "Ordenação das Linhas!";
        //M0.Print();
        //Fim da Ordenação das linhas

        //Eliminação de Gauss

        init = GLOBAL;
        for(int k = col; k < M0.GetFN(); k++)
        {
            if( M0[GLOBAL][col] != 0 )
            {
                for(int i = init + 1; i < M0.GetFM(); i++)
                {
                    if(M0[i][k] != 0)
                    {
                        CO =  M0[i][col] / M0[GLOBAL][col];
                        //subtração de linhas da matriz A
                        M0[i] = M0.GetRow(i) - M0.GetRow(init) * CO;
                    }
                }
                GLOBAL ++;
                teste = 1;
                break;    
            }
            col ++;
        }

        //cout << endl << "Eliminação de Gauss";
        //M0.Print();
        //Fim da Eliminação de Gauss

        if(init >= M0.GetFM()-1 || col >= M0.GetFN()-1)
        {
            teste = 0;
        }

        if(teste == 0)
        {
            break;
        }
    }    

    //cálculo do determinante
    double det = 1;
    for(int i = 0; i < (*this).GetFM(); i++)
    {
        det *= M0[i][i];
    }
    
    if(count % 2 == 1)
    {
        det = -det;
    }

    BENCHMARK_END(0);
    ECHOF;
    
    return det;
}

Vec& FCmatrixFull::operator[] (int n) {

    return M[rowindices[n]];
}

FCmatrixFull FCmatrixFull::operator=(const FCmatrix& Mx) {

    ECHOS;
    if(&Mx != this)
    {
        M = Mx.GetM();

        int fm = Mx.GetFM();
        int fn = Mx.GetFN();

        rowindices = new int [fm];
        colindices = new int [fn];

        for(int i = 0; i < fm; i++)
        {
            rowindices[i] = i;
        }

        for(int i = 0; i < fn; i++)
        {
                colindices[i] = i;
        }            
    }
    ECHOF;
    return (*this);
}

FCmatrixFull FCmatrixFull::operator+(const FCmatrix& Mx) {
    ECHOS;
    int SIZE = (*M.begin()).size();
    double x;
    double a[SIZE];
    Vec V0;
    vector<Vec> M0;

    if((*this).GetFM() != Mx.GetFM() || (*this).GetFN() != Mx.GetFN())
    {
        cout << "As matrizes não têm dimensões apropriadas para a soma!" << endl;
        FCmatrixFull NEW(M0);
        return NEW;
    }

    for(int i = 0; i < M.size(); i++)
    {
        for(int j = 0; j < SIZE; j++)
        {
            x = (M[i])[j];
            x += Mx.GetCol(i)[j];
            a[j] = x;
        }
        V0.SetEntries(SIZE, a);
        M0.push_back(V0);
    }

    FCmatrixFull NEW(M0);

    ECHOF;
    return NEW;
}

FCmatrixFull FCmatrixFull::operator-(const FCmatrix& Mx) {
    ECHOS;
    int SIZE = (*M.begin()).size();
    double x;
    double a[SIZE];
    Vec V0;
    vector<Vec> M0;

    if((*this).GetFM() != Mx.GetFM() || (*this).GetFN() != Mx.GetFN())
    {
        cout << "As matrizes não têm dimensões apropriadas para a subtração!" << endl;
        FCmatrixFull NEW(M0);
        return NEW;
    }

    for(int i = 0; i < M.size(); i++)
    {
        for(int j = 0; j < SIZE; j++)
        {
            x = (M[i])[j];
            x -= Mx.GetCol(i)[j];
            a[j] = x;
        }
        V0.SetEntries(SIZE, a);
        M0.push_back(V0);
    }

    FCmatrixFull NEW(M0);

    ECHOF;
    return NEW;
}

FCmatrixFull FCmatrixFull::operator*(const FCmatrix& Mx) {


    ECHOS;
    vector<Vec> M0;
    int fm1 = (*this).GetFM();
    int fn1 = (*this).GetFN();
    int fm2 = Mx.GetFM();
    int fn2 = Mx.GetFN();
    double mul;
    double ent[fn2];
    Vec V0;
    
    //cout << "fm1 = " << fm1 << " fn1 = " << fn1 << " fm2 = " << fm2 << " fn2 = " << fn2;

    if(fn1 != fm2)
    {
        cout << "As matrizes não têm dimensões que permitam a multiplicação!" << endl;
	#ifdef DEBUG
	printf("M1 dim %dx%d, M2 dim %dx%d\n", fm1, fn1, fm2, fn2);
	printf("\n[FCmatrixFull::operator* FCmatrix] END\n\n");
	#endif
        FCmatrixFull NEW(M0);
        return NEW;
    }
    else
    {
        for(int i = 0; i < fm1; i++)
        {
            for(int j = 0; j < fn2; j++)
            {
                mul = (*this).GetRow(i).dot(Mx.GetCol(j));
                ent[j] = mul;
            }
            V0.SetEntries(fn2, ent);
            M0.push_back(V0);
        }
    }

    //cout << endl << "Hello" << endl;
    FCmatrixFull NEW(M0);

    ECHOF;
    return NEW;
}

FCmatrixFull FCmatrixFull::operator*(const Vec& V) {

    ECHOS;
    
    double *mul;
    int fm = M.size();
    int fn = (*M.begin()).size();
    vector<Vec> M0;
    Vec V0[fm];

    if (V.size() != fn)	{
        FCmatrixFull NEW(M0);
        
        cout << "A matriz e o vetor não têm dimensões que permitam a multiplicação!" << endl;
	
	ECHOF;
        return NEW;
    }

    for(int i = 0; i < fm; i++)
    {
        *mul = 0;

        for(int j = 0; j < fn; j++)
        {
            *mul += (M[i])[j] * V[j];
        }

        V0[i].SetEntries(1, mul);   
        M0.push_back(V0[i]);
    }

    FCmatrixFull NEW(M0);

    ECHOF;
    return NEW;
}

FCmatrixFull FCmatrixFull::operator*(double lambda) {

    #ifdef DEBUG
    printf("[FCmatrixFull::operator* double] START\n");
    #endif
    int fm = M.size();
    int fn = (*this).GetFN();
    vector<Vec> M0;
    double mul[fn];
    Vec V0;

    for(int i = 0; i < fm; i++)
    {
        for(int j = 0; j < fn; j++)
        {
            mul[j] = (M[i])[j]*lambda;
        }

        V0.SetEntries(fn, mul);   
        M0.push_back(V0);
    }

    FCmatrixFull NEW(M0);

    #ifdef DEBUG
    printf("[FCmatrixFull::operator* double] END\n");
    #endif
    
    return NEW;
}


Vec FCmatrixFull::dot(const Vec& v)
{
    int N = M.size();
    
    Vec temp(N);

    for(int i = 0; i<N; i++){
	temp[i] = GetRow(i).dot(v);
    }

    return temp;
}

string FCmatrixFull::PrintToString(string format) const {
    #ifdef DEBUG
    printf("\n[FCmatrixFull::PrintToString] START\n");
    #endif
    string temp;

    for(int i = 0; i<GetFM(); i++){
	temp += GetRow(i).PrintToString(format) + "\n";
    }

    #ifdef DEBUG
    printf("[FCmatrixFull::PrintToString] END\n");
    #endif
    return temp;
}

void FCmatrixFull::Print(string format) const {
    cout << '\n' << PrintToString(format);
}

void FCmatrixFull::PrintIndices() {

    int fm = (*this).GetFM();
    int fn = (*this).GetFN();
    
    cout << endl;
    for(int i = 0; i < fm; i++)
    {
        cout << (*this).rowindices[i] << "  ";
    }
    cout << endl;
    for(int i = 0; i < fn; i++)
    {
        cout << (*this).colindices[i] << "  ";
    }
    cout << endl;
}

int FCmatrixFull::Getrowindices(int i){
    
    return rowindices[i];
}

int* FCmatrixFull::Getcolindices(){
    
    return colindices;
}

void FCmatrixFull::swapRows(int i,int j) {

    if (i != j)
    {
        int l = rowindices[i];
        rowindices[i] = rowindices[j];
        rowindices[j] = l;
        //cout << endl << "swap rows " << i << " e " << j << endl;
    }
}

void FCmatrixFull::swapCols(int i,int j) {

    if (i != j)
    {
        int l = colindices[i];
        colindices[i] = colindices[j];
        colindices[j] = l;
        //cout << endl << "swap cols " << colindices[i] << " e " << rowindices[j] << endl;
    }
}

FCmatrixFull FCmatrixFull::Transpose() {

    int fn = M[0].size();
    vector<Vec> M0;
      
    for(int i = 0; i < fn; i++)
    {
        M0.push_back((*this).GetCol(i));
    }

    int * ri = rowindices;
    int * ci = colindices;

    FCmatrixFull NEW(M0);

    NEW.rowindices = ci;
    NEW.colindices = ri;

    return NEW;
}

FCmatrixFull FCmatrixFull::EliminateRow(int i) {

    vector<Vec> M0;
    int * rowindices0 = new int [(*this).GetFM() - 1];

    int k = 0;
    for(int j = 0; j < (*this).GetFM(); j++)
    {
        if(j != i)
        {
            M0.push_back((*this).GetRow(j));
            rowindices0[k] = k;
            k++;
        }
    }

    delete[] rowindices;
    rowindices = rowindices0;
    M = M0;
    return *this;
}

FCmatrixFull FCmatrixFull::EliminateCol(int i) {

    vector<Vec> M0;
    Vec V0;
    double *a = new double [(*this).GetFM()];

    int * colindices0 = new int [(*this).GetFN() - 1];

    int k = 0;
    for(int l = 0; l < (*this).GetFM(); l++)
    {
        for(int j = 0; j < (*this).GetFN(); j++)
        {
            if(j != i)
            {
                a[k] = M[rowindices[l]][k];
                colindices0[k] = k;
                k++;
            }
        }
        V0.SetEntries((*this).GetFN()-1, a);
        M0.push_back(V0);
        k = 0;
    }
    
    delete[] colindices;
    delete[] a;
    colindices = colindices0;
    M = M0;
    return *this;
}

FCmatrixFull::~FCmatrixFull() {
    #ifdef DEBUG
    printf("\n[FCmatrixFull::~FCmatrixFull] START\n\n");
    #endif
    
    delete[] rowindices;
    delete[] colindices;

    #ifdef DEBUG
    printf("[FCmatrixFull::~FCmatrixFull] END\n\n");
    #endif
    //cout << endl << "Destructor de FCmatrixFull" << endl;
}
