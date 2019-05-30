#include "FCmatrix.h"
#include <vector>
#include "Vec.h"
using namespace std;

/*FCmatrix::FCmatrix() {

    vector<Vec> M0;

    classname = "";
    M = M0;
}*/

FCmatrix::FCmatrix(double** fM, int fm, int fn) {

    /*cout << "fm = " << fm << endl;
    cout << "fn = " << fn << endl;*/

    for(int i = 0; i < fm; i++)
    {
        Vec V(fn, fM[i]);
        M.push_back(V);
        
        //Print para confirmar
        /*for(int j = 0; j < fn; j++)
        {
            cout << V[j] << "  ";
        }*/
    }
    

    //cout << endl << "Constructor type 1 called!" << endl;
}

FCmatrix::FCmatrix(double* fM, int fm, int fn) {

    Vec V(fm, fM);
    for(int i = 0; i < fm; i++)
    {
        M.push_back(V);
    }
    
    //cout << endl << "Constructor type 2 called!" << endl;
}

FCmatrix::FCmatrix(const vector<Vec> V) {

    for(int i = 0; i < V.size(); i++)
    {
        M.push_back(V[i]);
    }
    
    //cout << endl << "Copy Constructor!" << endl;
}

int FCmatrix::GetFM() const {

    return M.size();
}

int FCmatrix::GetFN() const{

    return M[0].size();
}

string FCmatrix::Getclassname() const {

    return classname;
}

Vec& FCmatrix::operator[] (int n) {

    return M[n];
} 

Vec FCmatrix::GetRow(int i) const {

    return M[i];
}

Vec FCmatrix::GetCol(int i) const {

    Vec V;
    double a[M.size()];
    vector<Vec>::iterator vecit;

    for(int j = 0; j < M.size(); j++)
    {
        a[j] = (M[j])[i];
    }

    V.SetEntries(M.size(), a);

    return V;
}

vector<Vec> FCmatrix::GetM() const {
    return M;
}

vector<Vec> FCmatrix::GetRawM() const{
    return M;
}

double FCmatrix::Determinant() {

    int det = 0;
    int det1 = 1;
    int det2 = 1;
    int nlines = M.size();

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
    }
    
    return det;
}

void FCmatrix::Print(string format) const {
    
    cout << '\n' << PrintToString(format);
}    
    
double FCmatrix::GetRowMax(int i=0) {

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

double FCmatrix::GetColMax(int j=0) {

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


FCmatrix::~FCmatrix() {
}
