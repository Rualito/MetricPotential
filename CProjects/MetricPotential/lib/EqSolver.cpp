#include "FCmatrix.h"
#include "FCmatrixFull.h"
#include "FCmatrixBanded.h"
#include "Vec.h"
#include "EqSolver.h"
#include <math.h>

// soma int
#define SOMA(x) ((x)*((x)+1)/2)

//#define DEBUG
//#define COUNT
#include "Debug.hpp"

//#define BENCHMARK
#include "BenchMark.hpp"

EqSolver::EqSolver() {

    M = nullptr;
}

EqSolver::EqSolver(const FCmatrix& Mx, const Vec& b0) {
    ECHOS;
    
    if(Mx.Getclassname() == "FCmatrixFull")
    {
	M = new FCmatrixFull (Mx.GetM());
    }
    else if (Mx.Getclassname() == "FCmatrixBanded")
    {
	M = new FCmatrixBanded(Mx);
    }

    b = b0;
    ECHOF;
}

void EqSolver::SetConstants(const Vec& b0) {

    b = b0;
}

void EqSolver::SetMatrix(const FCmatrix& Mx) {

    *M = Mx;
}

void EqSolver::GaussElimination(FCmatrixFull &M0, Vec& b0) {
    ECHOS;
    //cout << endl << M0.GetFM() << endl << M0.GetFN() << endl;

    int l;
    int init;
    int col = 0;
    int GLOBAL = 0;
    int teste = 0;
    double CO;
    double c;
    //M0.Print();

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
                    //troca os elementos do Vec
                    c = b0[l];
                    b0.Seti(l, b0[i]);
                    b0.Seti(i, c);

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
        //b0.Print();
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
                        CO = M0[i][col] / M0[GLOBAL][col];
                        
                        //subtração de linhas da matriz A
                        M0[i] = M0.GetRow(i) - M0.GetRow(init) * CO;
                        //subtração de elementos do Vec
                        b0.Seti(i, b0[i] - b0[init] * CO );
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
        //b0.Print();
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
    //cout << endl << "acabou" << endl;
    ECHOF;
}

Vec EqSolver::GaussEliminationSolver() {
    ECHOS;
    if(M->GetFN() != b.size())
    {
        cout << endl << "A matriz A e o vetor b não têm dimensões compatíveis!" << endl;

        throw runtime_error("Wrong matrix dimensions");
	    return Vec();
    }

    FCmatrixFull M1(*M);
    Vec b0(b);
    
    //Eliminação de Gauss
    GaussElimination(M1, b0);

    //(M1).Print();

    //Eliminação de linhas nulas
    int count = 0;
    int teste0 = 0;
    for(int i = M1.GetFM() - 1; i >= 0 ; i--)
    {
        for(int j = 0; j < M1.GetFN(); j++)
        {
            if(M1[i][j] == 0)
            {
                teste0 = 1;
            }
            else
            {
                teste0 = 0;
                break;
            }
        }

        if(teste0 != 0)
        {
            M1.EliminateRow(i);
            teste0 = 0;

            if(b0[i] != 0 && (M1.GetFM() + 1) == b0.size())
            {
                cout << endl << "O sistema é impossível  (caso 0 = constante)!" << endl;
                Vec SOLUCAO;
                return SOLUCAO;
            }
        }
        else
        {
            break;
        }
    }
    //M1.Print();
    int FM = M1.GetFM();

    //Eliminação de colunas nulas
    teste0 = 0;
    for(int i = M1.GetFN() - 1; i >= 0 ; i--)
    {
        for(int j = 0; j < M1.GetFM(); j++)
        {
           if(M1[j][i] == 0)
            {
                teste0 = 1;
            }
            else
            {
                teste0 = 0;
                break;
            }
        }

        if(teste0 == 1)
        {
            M1.EliminateCol(i);
            count ++;
        }
        else
        {
            break;  
        }
    }
    //M1.Print();
    int FN = M1.GetFN();

    //Eliminação de Gauss-Jordan
    //Determinação da solução de Ax = b (Ux = c)
    double *a = new double [M1.GetFM()]();

    teste0 = 0;
    if(FM == FN)
    {
        for(int i = M1.GetFM() - 1; i >= 0 ; i--)
        {
            for(int j = i + 1; j < M1.GetFN(); j++)
            {
                a[i] += - a[j] * M1[i][j];
            }
            a[i] += b0[i];
            a[i] = a[i] / (M1[i][i]);

            //cout << "x" << i << " = " << a[i] << endl;
        }
        
        if(count > 0)
        {
            cout << endl << "A solução deveria ter " << (*M).GetFM() << " elementos!" << endl;
            cout << "Faltam elementos da solução porque o sistema é indeterminado!" << endl;
            cout << "Os elementos em falta são a variáveis independentes!" << endl;
        }

        Vec SOLUCAO(M1.GetFM(), a);
        delete[] a;

        //imprime a solução
        /*cout << endl << "Solução do Sistema:" << endl;
        for(int i = 0; i < M1.GetFM(); i++)
        {
            cout << SOLUCAO[i] << "  ";
        }
        cout << endl;*/
	ECHOF;
        return SOLUCAO;
    }
    else
    {
        cout << endl << "O sistema não é possível e determinado!" << endl;
        Vec SOLUCAO;
        delete[] a;   
        return SOLUCAO;
    }
    
}

void EqSolver::LUdecomposition(FCmatrixFull& M0, vector<int> index) {
    ECHOS;
    //Decompõe a matriz A numa matriz L (lower) e U (upper)

    double c;

    //constrói matriz justaposta [L \ U]
    for(int i = 0; i < M0.GetFM(); i++)
    {
        for(int j = 0; j < M0.GetFN(); j++)
        {
            if(j >= i)
            {                 
                //preenche entradas da matriz U (j >= i)
                c = 0;
                for(int k = 0; k < i; k++)
                {
                    c += - M0[index[k]][j] * M0[index[i]][k];
                }
                c += M0[index[i]][j];                   
                M0[index[i]].Seti(j, c);
            }
            else
            {
                //preenche entradas não diagonais da matriz L (i < j)
                c = 0;
                for(int k = 0; k < j; k++)
                {
                    c += - M0[index[k]][j] * M0[index[i]][k];
                }
                c += M0[index[i]][j];   
                c = c/M0[index[j]][j];
                M0[index[i]].Seti(j, c);       
            }
        }
    }

    /*cout << endl << "Matriz L e U" << endl;
    M0.Print();*/
    ECHOF;
}

void EqSolver::LUdecomposition(FCmatrixBanded& M0) {
    ECHOS;
    //Decompõe a matriz A numa matriz L (lower) e U (upper)

    double c;

    //constrói matriz justaposta [L \ U]
    for(int i = 1; i < M0.GetFM(); i++)
    {
        if(M0.GetVec(1)[i-1] == 0)
        {
            cout << endl << "A decomposição LU da matriz é impossível sem trocas de linhas!" << endl;
            exit(1);
        }
        M0.Seti(1, i, M0.GetVec(1)[i] - M0.GetVec(0)[i-1]*M0.GetVec(2)[i-1]/M0.GetVec(1)[i-1]);
        M0.Seti(2, i-1, M0.GetVec(2)[i-1]/M0.GetVec(1)[i-1]);
    }
    //M0.PrintFull();
    ECHOF;
}

Vec EqSolver::LUdecompositionSolver() {
    ECHOS;
    //retorna a solução pelo método LU
    
    if((*M).GetFN() != b.size())
    {
        cout << endl << "A matriz A e o vetor b não têm dimensões compatíveis!" << endl;
	    throw runtime_error("Wrong matrix dimensions");
	
        Vec SOLUCAO;
        return SOLUCAO;
    }

    //Banded
    if((*M).Getclassname() == "FCmatrixBanded")
    {
        FCmatrixBanded M1(*M);

        //A matriz *M é transformada em [L \ U]
        LUdecomposition(M1);
        //M1.PrintFull();

        //Obtenção da solução de Ly = b
        int size1 = M1.GetFM();
        double *a = new double [size1]();

        a[0] = b[0];
        for(int i = 1; i < size1; i++)
        {
            a[i] += b[i] - a[i-1]*M1[i][i-1];
        }
        //(*M).Print();

        Vec y (size1, a);
        //y.Print();

        //Obtenção da solução de Ux = y
        double *A = new double [size1]();

        A[size1-1] = y[size1-1]/M1[size1-1][size1-1];
        for(int i = size1 - 2; i >= 0 ; i--)
        {
            A[i] = (y[i]-M1[i][i+1]*A[i+1])/M1[i][i];
        }

        Vec SOLUCAO (size1, A);

        //SOLUCAO.Print();
        delete[] a;
        delete[] A;
	ECHOF;
	return SOLUCAO;        
    }

    //Full
    FCmatrixFull M1(*M);
    //(*M).Print();
    
    //Obtenção da matriz P (permutação) para a decomposição PA = PLU
    Vec b0(b);
    FCmatrixFull M0(*M);
    GaussElimination(M0, b0);
    vector<int> vri;
    for(int i = 0; i < M0.GetFM(); i++)
    {
        vri.push_back(M0.Getrowindices(i));
    }

    //Decomposição LU (apenas matrizes quadradas)
    if(M->GetFM() == M->GetFN())
    {
        int size1 = (*M).GetFM();
        int count = 0;
        int teste0 = 0;

        //Eliminação de linhas nulas
        for(int i = M1.GetFM() - 1; i >= 0 ; i--)
        {
            for(int j = i; j < M1.GetFN(); j++)
            {
                if(M1[i][j] == 0)
                {
                    teste0 = 1;
                }
                else
                {
                    teste0 = 0;
                    break;
                }
            }

            if(teste0 != 0)
            {
                M1.EliminateRow(i);
                count++;
                teste0 = 0;

                if(b[i] != 0)
                {
                    cout << endl << "O sistema é impossível  (caso 0 = constante)!" << endl;
		    ECHOF;
                    Vec SOLUCAO;
                    return SOLUCAO;
                }
            }
            else
            {
                break;
            }
        }
        //(*M).Print();

        //Eliminação de colunas nulas
        teste0 = 0;
        for(int i = M1.GetFN() - 1; i >= 0 ; i--)
        {
            for(int j = 0; j < M1.GetFM(); j++)
            {
                if((*M)[j][i] == 0)
                {
                    teste0 = 1;
                }
                else
                {
                    teste0 = 0;
                    break;
                }
            }

            if(teste0 == 1)
            {
                M1.EliminateCol(i);
                count ++;
            }
            else
            {
                break;  
            }
        }
        //M1.Print();

        int size2 = M1.GetFM();
        int size3 = M1.GetFN();

        if(size1 != size2 || size1 != size3)
        {
            cout << endl << "O sistema não é possível e determinado!" << endl;
            Vec SOLUCAO;
            return SOLUCAO;
        }
  
        //A matriz *M é transformada em [L \ U]
        LUdecomposition(M1, vri);
        //M1.Print();

        //Obtenção da solução de Ly = b
        double *a = new double [size1]();

        for(int i = 0; i < size1; i++)
        {
            for(int j = 0; j < i; j++)
            {
                a[vri[i]] -= a[vri[j]]*M1[vri[i]][j];
            }

            a[vri[i]] += b[vri[i]];
        }
        //(*M).Print();

        Vec y (size1, a);
        //y.Print();

        //Obtenção da solução de Ux = y
        double *A = new double [size1]();

        for(int i = size1 - 1; i >= 0 ; i--)
        {
            for(int j = i + 1; j < size1; j++)
            {
                A[vri[i]] += - A[vri[j]] * M1[vri[i]][j];
            }
            A[vri[i]] += y[vri[i]];
            A[vri[i]] = A[vri[i]] / M1[vri[i]][i];
        }

        Vec SOLUCAO0(size1, A);
        Vec SOLUCAO(size1, A);

        for(int i = 0; i < size1; i++)
        {
            SOLUCAO.Seti(i, SOLUCAO0[vri[i]]);
        }

        //SOLUCAO.Print();
        delete[] a;
        delete[] A;
        vri.clear();
	ECHOF;
	return SOLUCAO;
    }
    else
    {
        cout << endl << "A matriz não é quadrada, logo a decomposição LU não se aplica!" << endl;
        Vec SOLUCAO;
        return SOLUCAO;
    }
}

Vec EqSolver::JacobiIterator(double tol) {

    if(M->GetFN() != b.size())
    {
        cout << endl << "A matriz A e o vetor b não têm dimensões compatíveis!" << endl;
	throw runtime_error("Wrong matrix dimensions");
	
        Vec SOLUCAO;
        return SOLUCAO;
    }
    
    // linear system of m unknowns
    int m = b.size();
    Vec x(m, 0.); 
    Vec x_aux(m, 0.); //zero’s
    bool btol = false;
    int it = 0.; 
    double eps = tol;
    
    while (!btol && (it++ < 1000)) 
    {
        x_aux = x;
        
        for (int i=0; i<m; i++)
        {
            x.Seti(i, 0.);
            
            for (int j=0; j<m; j++)
            {
                if (i != j) 
                {
                    x.Seti(i, x[i] -(*M)[i][j]*x_aux[j]);
                }
            }

            x.Seti(i, x[i] + b[i]);
            x.Seti(i, x[i] / (*M)[i][i]);
            //guarantee that all vector entries are converging equally
            if (fabs(x[i]-x_aux[i]) < eps) btol = true;
            else btol = false;

            //só fica true se para todas as componentes de x o erro for < eps
        }
    }

    return x;
}

Vec EqSolver::GaussSeidelIterator(double tol) {

    if(M->GetFN() != b.size())
    {
        cout << endl << "A matriz A e o vetor b não têm dimensões compatíveis!" << endl;
	throw runtime_error("Wrong matrix dimensions");
	
        Vec SOLUCAO;
        return SOLUCAO;
    }
    
    // linear system of m unknowns
    int m = b.size();
    Vec x(m, 0.); 
    Vec x_aux(m, 0.); //zero’s
    bool btol = false;
    int it = 0.; 
    double eps = tol;
    
    while (!btol && (it++ < 1000)) 
    {
        x_aux = x;
        
        for (int i=0; i<m; i++)
        {
            x.Seti(i, 0.);
            
            for (int j=0; j<m; j++)
            {
                if (i != j) 
                {
                    x.Seti(i, x[i] -(*M)[i][j]*x[j]);
                }
            }
            
            x.Seti(i, x[i] + b[i]);
            x.Seti(i, x[i] / (*M)[i][i]);
            //guarantee that all vector entries are converging equally
            if (fabs(x[i]-x_aux[i]) < eps) btol = true;
            else btol = false;

            //só fica true se para todas as componentes de x o erro for < eps
        }
    }

    return x;
}

Vec EqSolver::TridiagonalSolver(){
    ECHOS;
    BENCHMARK_START(1);
    FCmatrixBanded MB(*M);

    Vec gamma(MB.GetVec(0)); // starts with i=0, ends at N-2
    Vec beta(MB.GetVec(1)); // starts with i=0, ends at N-1
    Vec alpha(MB.GetVec(2)); // starts with i=1, ends at N-1

    int N = beta.size();
    
    Vec g(N-1);
    Vec h(N-1);
    
    g[N-2] = -alpha[N-2]/beta[N-1];
    h[N-2] = b[N-1]/beta[N-1];
    
    for(int i = N-2; i>0; i--){
	g[i-1] = -alpha[i-1]/(beta[i]+gamma[i]*g[i]);
	h[i-1] = (b[i] - gamma[i]*h[i])/(beta[i] + gamma[i]*g[i]);
    }

    Vec x(N);

    x[0] = (b[0] - gamma[0]*h[0])/(beta[0] + gamma[0]*g[0]);
    
    for(int i = 1; i<N; i++){
	x[i] = g[i-1]*x[i-1] + h[i-1];
    }
    BENCHMARK_END(1);
    ECHOF;
    return x;
}
     
void EqSolver::Print() {

    cout << endl << "Matriz A:" << endl;

    if((*M).Getclassname() == "FCmatrixFull")
    {
        M->Print();
    }
    else
    {
        FCmatrixBanded M0(*M);
        M0.PrintFull();
    }
    

    /*
    for(int i = 0; i < M->GetFM(); i++)
    {
        for(int j = 0; j < M->GetFN(); j++)
        {
            cout << (*M)[i][j] << "  ";
        }

        cout << endl;
	}*/

    cout << endl << "Vector b:" << endl;

    b.Print();
    /*
    for(int i = 0; i < b.size(); i++)
    {
        cout << b[i] << "  ";
	}*/

    cout << endl;

}

EqSolver::~EqSolver() {

    delete M;
}
