#include "DataPoints.hpp"

//#define COUNT
//#define DEBUG
#include "Debug.hpp"

DataPoints::DataPoints()
{
    ECHOS;
    graphics = false;

    N = 0;
    x = y = NULL;
    
    ECHOF;
}


DataPoints::DataPoints(int n, double *xp, double *yp) : N(n)
{
    ECHOS;
    graphics = false;
    
    if(n && xp && yp){
	x = new double[n];
	y = new double[n];

	for(int i = 0; i<n; i++){
	    x[i] = xp[i];
	    y[i] = yp[i];
	}
    }
    ECHOF;
}

void DataPoints::StartGraphics(){
    if(!graphics){
	graphics = true;
	G.Start();
	
    }
}

DataPoints::~DataPoints()
{
    ECHOS;
    if(x)
	delete[] x;
    if(y)
	delete[] y;

    ECHOF;
}

double* DataPoints::GetArrayX(){
    return x;
}

double* DataPoints::GetArrayY(){
    return y;
}

int DataPoints::size() const{
    return N;
}

void DataPoints::Draw(double step)
{
    ECHOS;

    StartGraphics();
    static int calln = 0;

    COUNT_P(2);
    
    calln++;
    char padName[24];
    
    sprintf(padName, "pad%d", calln);
    
    TPad *pad1 = G.CreatePad(padName);

    TGraph *g1 = new TGraph(N,x,y);

    COUNT_P(2);
    
    g1->SetMarkerStyle(33);
    g1->SetMarkerColor(kRed);
    g1->SetMarkerSize(2);
    
    if(step != 0){
	//note: for this to work the points have to be ordered
	int stepn = (int)((x[N-1]-x[0])/step);
	printf("stepn %d\n",stepn);
	
	double xi[stepn];
	double yi[stepn];

	COUNT_P(2);
    
	for(int i = 0; i<stepn; i++){
	    xi[i] = x[0] + i*step;
	    yi[i] = Interpolate(xi[i]);
	    
	}
	COUNT_P(2);
    
	TGraph *g2 = new TGraph(stepn, xi, yi);
	g2->SetMarkerStyle(20);
	g2->SetMarkerColor(kBlack);
	g2->SetLineWidth(2);
	g2->SetMarkerSize(0.1);
	
	G.AddObject(g2, padName, "APL");
	G.AddObject(g1,padName,"P");
    } else {
	G.AddObject(g1,padName,"APL");
    }
    
    G.AddObject(pad1);
    G.DrawPad(padName);
    
    ECHOF;
}

double DataPoints::Interpolate(double x) const
{
    return 0;
}

void DataPoints::Print(string filename, string format)
{
    ECHOS;
    
    TGraph *g = new TGraph(N,x,y);
    
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kRed);
    g->SetMarkerSize(1.5);
    
    TPad *pad1 = G.CreatePad("pad1");
    
    G.AddObject(g,"pad1","AP");
    G.AddObject(pad1);

    G.Print(filename);
    
    double vec[2*N];
    for(int i = 0; i<N; i++){
	vec[i] = x[i];
	vec[i+N] = y[i];
    }

    FCmatrixFull mat(vec, N, 2);
    
    if(filename == ""){
	mat.Print(format);
    } else {
	FCtools::PrintMatrixToFileT(mat, filename, format);
    }    

    ECHOF;
}
/*
TGraph& DataPoints::GetPointGraph(){

    TGraph *tg= new TGraph(N, x, y);
    
    return *tg;
}

TGraph& DataPoints::GetInterpolatedGraph(double step)
{
    int stepN = (int)((x[N-1]-x[0])/step);
    double xt[stepN];
    double yt[stepN];

    for(int i = 0; i<stepN; i++){
	xt[i] = x[0] + i*step;
	yt[i] = Interpolate(xt[i]);
    }

    TGraph *tg = new TGraph(stepN, xt, yt);

    return *tg;    
    }*/

	
    

