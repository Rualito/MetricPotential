#ifndef HEQSOLVERh
#define HEQSOLVERh
#include "Vec.h"
#include "FCmatrix.h"
#include "FCmatrixFull.h"
#include "FCmatrixBanded.h"
#include <vector>
using namespace std;

class EqSolver {
    
    public:
        EqSolver();
        EqSolver(const FCmatrix&, const Vec&);
        // set
        void SetConstants(const Vec&);
        void SetMatrix(const FCmatrix&);
        //eliminação de Gauss:
        //resolução do sistema pelo método de eliminação de Gauss
        Vec GaussEliminationSolver();
        Vec LUdecompositionSolver();
        Vec JacobiIterator(double tol = 1.e-6);
        Vec GaussSeidelIterator(double tol = 1.e-6);
        void Print();
        ~EqSolver();

	Vec TridiagonalSolver();
	
    private:
        //decomposição LU com |L|=1
        void LUdecomposition(FCmatrixFull&, vector<int>);
        void LUdecomposition(FCmatrixBanded&);
        /* return triangular matrix and changed vector of constants */
        void GaussElimination(FCmatrixFull&, Vec&);
        FCmatrix *M; //matriz de coeffs
        Vec b; //vector de constantes
};



#endif
