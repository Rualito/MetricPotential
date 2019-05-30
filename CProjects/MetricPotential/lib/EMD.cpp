
#include "EMD.hpp"
#include "IMF.hpp"

//#define DEBUG
//#define COUNT
#include "Debug.hpp"

//#define BENCHMARK
#include "BenchMark.hpp"
#define ABS(x) (x>0?(x) : -(x) )

EMD::EMD(){;}

EMD::EMD(int n, double *fx, double *fy) : DataPoints(n,fx,fy) {
    ECHOS;
    
    maxMedVar = 0.05;
    maxVar = 0.5;
    medVarRate = 0.05;

    (*this)+=0;
    ECHOF;
}

EMD::~EMD() {
    IMFs.clear();
    ECHOF;
}

int EMD::size() const{
    return IMFs.size();
}


IMF& EMD::operator[](int p){
    return IMFs[p];
}

void EMD::SetVarianceParameters(double maxMed, double maxVAR, double medVRate) 
{
    maxMedVar = maxMed;
    maxVar = maxVAR;
    medVarRate = medVRate;
}

void EMD::UpdateAll(){
    for(int i = 0; i<IMFs.size(); i++){
	IMFs[i].CalculateExtras();
    }
}

bool EMD::operator+=(int steps)
{
    // it returns "did I find the last IMF?"
    // finds "steps" IMFs
    // if its 0, then it generates the first IMF, without sifting
    ECHOS;
        
    int p0 = IMFs.size();
    bool test;
    
    if(p0 == 0)
	IMFs.emplace_back(N, x, y);

    if(!steps)
	return false;
    
    test = SiftingProcess();
    
    if(!test)
	return false;
    
    if(steps>-1){
	for(int i = 0; i<steps && test; i++){
	    IMFs.emplace_back(N, x, IMFs.back().GetMean());
	    test = SiftingProcess();
	}
    } else {
	while(test){
	    IMFs.emplace_back(N, x, IMFs.back().GetMean());
	    test = SiftingProcess();
	}
    }
	    
    ECHOF;
    return test;
}

bool EMD::SiftingProcess()
{
    // returns "should the next IMF be calculated? (does the next IMF exist)"
    ECHOS;
    bool flagContinue; // corresponds to not being the last IMF
    bool flagVar;

    IMF& imf = IMFs.back();
    
    // while its not the last IMF and while the variance isnt checked
    do {
        flagContinue = ++imf;
	flagVar = CheckVariance(imf);
    } while(flagContinue && !flagVar);

    if(!flagContinue){
	IMFs.pop_back();
    }
    
    ECHOF;
    return flagContinue;
}

bool EMD::CheckVariance(const IMF& imf)
{
    ECHOS;
    int count = 0;
    
    double *var = imf.GetVariance();
    for(int i = 0; i<N; i++){
	if(var[i] > maxMedVar){
	    ++count;
	} else if(var[i] > maxVar){
	    ECHOF;
	    return false;
	}
    }
    
    // the ratio of variances that are bigger than the recomended maximum
    // is less than the given ratio
    ECHOF;
    return ( ( ((double)count)/N ) < medVarRate ); 
}
/*

  void EMD::DrawSplines(int index){
    IMFs[i]
*/
void EMD::DrawIMFs(int perPad)
{
    ECHOS;
    int IMFcount = IMFs.size();

    printf(">> imf %d\n", IMFcount );
    
    if(perPad < 1)
	perPad = -1;
    
    int padn = IMFcount/perPad;
    if(padn <= 0)
	padn = 1;
    
    cFCgraphics g[padn];
    int j;
    TPad *pad[padn];
    string padName;
    for(int i = 0; i<IMFcount; i++){
	j = i/perPad;
	
	if(i%perPad == 0){
	    if(i>0){
		g[j].AddObject(pad[j-1]);
	    }
	    
	    g[j].Start();
	    padName = "EMD pad "+j;
	    pad[j] = g[j].CreatePad(padName);
	}
	
	TGraph *tg = new TGraph(N, x, IMFs[i].GetIntrinsic());

	tg->SetMarkerStyle(20);
	tg->SetLineColor(j%10);
	tg->SetLineWidth(1);
	
	if(i == 0){
	    g[j].AddObject(tg, padName, "AL");
	    
	    continue;
	}
	g[j].AddObject(tg, padName, "L");
    }

    for(int i=0; i<padn; i++){
	g[i].Draw();
    }

    // throw runtime_error("Draw");
    ECHOF;
}
