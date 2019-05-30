
//#define DEBUG
//#define COUNT
#include "Debug.hpp"

//#define BENCHMARK
#include "BenchMark.hpp"

#include "IMF.hpp"

#define ABS(x) ((x)>0?(x):-(x))

#define CODE_MAX 4
#define CODE_MIN 2
#define CODE_ZERO 1

int IMF::indexCount = 0;

IMF::IMF() : DataPoints()
{
    ECHOS;
    index = indexCount++;

    queueCalculateExtras = true;
    
    correlation = 0;
    intrinsic = old_intrinsic = variance = mean = frequency = NULL;
    SplineMax = SplineMin = NULL;
	
    ECHOF;
}


IMF::IMF(int n, double* xp, double* yp) : DataPoints(n,xp,yp)
{
    ECHOS;
    
    index = indexCount++;
    frequency = 0;

    if(n && xp && yp){
        intrinsic = new double[n];
	old_intrinsic = new double[n];
	variance = new double[n];
	mean = new double[n];
	frequency = new double[n];
	
	for(int i = 0; i<n; i++){
	    old_intrinsic[i] = intrinsic[i] = yp[i];
	    frequency[i] = 0.;
	}
	
	SplineMax = NULL;
	SplineMin = NULL;
        
	++(*this);
    }
    ECHOF;    
}

IMF::IMF(const IMF& imf) : DataPoints(imf.N, imf.x, imf.y){
    ECHOS;
    
    index = indexCount++;

    intrinsic = new double[N];
    old_intrinsic = new double[N];
    variance = new double[N];
    mean = new double[N];
    frequency = new double[N];
    
    CopyParameters(imf);

    ++(*this);
    ECHOF;
}

const IMF& IMF::operator=(const IMF& imf)
{
    ECHOS;

    if(&imf != this){
	if(imf.N>0){
	    
	    if(N!=imf.N){
		N = imf.N;
		
		if(intrinsic)
		    delete[] intrinsic;
		if(old_intrinsic)
		    delete[] old_intrinsic;
		if(variance)
		    delete[] variance;
		if(frequency)
		    delete[] frequency;
		
		if(SplineMax)
		    delete SplineMax;
		if(SplineMin)
		    delete SplineMin;
		
		intrinsic = new double[N];
		old_intrinsic = new double[N];
		variance = new double[N];
		mean = new double[N];
		frequency = new double[N];
		
		CopyParameters(imf);
	    }
	} else {
	    CopyParameters(IMF());
	}
    }
    
    ECHOF;
    return (*this);
}

IMF::~IMF()
{
    ECHOS;
    
    if(intrinsic)
	delete[] intrinsic;
    
    if(variance)
	delete[] variance;
    
    if(mean)
	delete[] mean;
    
    if(frequency)
	delete[] frequency;
    
    if(SplineMax)
	delete SplineMax;
    
    if(SplineMin)
	delete SplineMin;
    
    ECHOF;
}

void IMF::CopyParameters(const IMF& imf)
{
    ECHOS;

    int i;
    if(imf.intrinsic)
	for(i = 0; i<N; i++)
	    intrinsic[i] = imf.intrinsic[i];
    
    if(imf.old_intrinsic)
	for(i = 0; i<N; i++)
	    old_intrinsic[i] = imf.old_intrinsic[i];

    if(imf.variance)
	for(i = 0; i<N; i++)
	    variance[i] = imf.variance[i];

    if(imf.mean)
	for(i = 0; i<N; i++)
	    mean[i] = imf.mean[i];

    if(imf.frequency)
	for(int i = 0; i<N; i++)
	    frequency[i] = imf.frequency[i];
    
    //printf(">> %d %d\n", (SplineMax?1:0), (SplineMin?1:0));
    SplineMax = new Spline3Interpolator(*imf.SplineMax);
    SplineMin = new Spline3Interpolator(*imf.SplineMin);

    queueCalculateExtras = true;
    
    ECHOF;
}

// ================  Getters Sifting ======================== 
double* IMF::GetIntrinsic() const
{
    return intrinsic;
}

double* IMF::GetVariance() const
{
    return variance;
}

double IMF::operator[](int px) const
{
    return intrinsic[px];
}

// ================  Getters Extra ======================== 

double* IMF::GetFrequency()
{
    CalculateExtras();
    return frequency;
}

double IMF::GetCorrelation(double *ys)  
{
    CalculateCorrelation(ys);
    return correlation;
}

double* IMF::GetMean()
{
    CalculateExtras();
    return mean;
}

int IMF::GetMaxCount()
{
    CalculateExtras();
    return maxCount;
}

int IMF::GetMinCount()
{
    CalculateExtras();
    return minCount;
}

int IMF::GetZeroCount()
{
    CalculateExtras();
    return zeroCount;
}

// ============= Calculators Sifting ============

void IMF::CalculateIntrinsic()
{
    for(int i=0; i<N; ++i){
	old_intrinsic[i] = intrinsic[i];
        intrinsic[i] = intrinsic[i] - SplineMed(x[i]);
    }
}

void IMF::CalculateVariance()
{
    for(int i = 0; i<N; ++i){
	variance[i] = ABS(SplineMed(x[i])/
			  SplineAmp(x[i]));
    }
}

// ================= Calculators Extra ================

void IMF::CalculateMean()
{
    for(int i = 0; i<N; i++){
        mean[i] = y[i] - intrinsic[i];
    }
}

void IMF::CalculateCorrelation(double* ys)
{
    double sum1 = 0;
    double sum2 = 0;
    double sum3 = 0;

    for(int i = 0; i<N; i++){
	sum1 += ys[i] * intrinsic[i];
	sum2 += ys[i] * ys[i];
	sum3 += intrinsic[i] * intrinsic[i];
    }
    
    correlation = sum1/sqrt(sum2*sum3);    
}

void IMF::CountExtrema() 
{
    if(intrinsic){
	maxCount = GetCount(intrinsic, CODE_MAX);
	minCount = GetCount(intrinsic, CODE_MIN);
	zeroCount = GetCount(intrinsic, CODE_ZERO);
    } else {
	maxCount = minCount = zeroCount = -1;
    }
}

void IMF::CalculateFrequency()
{
    ECHOS;
    
    vector<int> zeroesIndx = GetIndexes(intrinsic, CODE_ZERO);
    
    int size = zeroesIndx.size();

    //average point between the zeroes
    double xf[size+2];

    //corresponding frequency
    double pf[size+2];

    // each zero 
    double avgx1;
    double avgx2;
    
    for(int i = 0; i<size-1; i++){
	avgx1 = 0.5*(x[zeroesIndx[i]] + x[zeroesIndx[i]+1]);
	avgx2 = 0.5*(x[zeroesIndx[i+1]] + x[zeroesIndx[i+1]+1]);
	
	xf[i+1] = 0.5*(avgx1 + avgx2);
	pf[i+1] = 0.5/(avgx2-avgx1);
    }
    // boundary conditions
    xf[0] = x[0];
    pf[0] = pf[1];

    xf[size] = x[N-1];
    pf[size] = pf[size-1];
    
    // spline that interpolates those points
    // Spline3Interpolator freqSpline(size-1, xf, pf);

    double a, b;
    
    for(int j=0; j<size+1; ++j){
	//a = (pf[j+1] - pf[j])/(xf[j+1] - xf[j]);
	//b = pf[j]-a*xf[j];

	int i;
	for(i=0; i<N; ++i){
	    if(x[i]>xf[j] && x[i]<xf[j+1])
		frequency[i] = pf[j];
	}
    }
    ECHOF;
}

// ================== Batch Calculators ==================
void IMF::Update() 
{
    CalculateMean();
    CalculateIntrinsic(); // optimize here
    CalculateVariance(); // optimize here
    
    queueCalculateExtras = true;
}

void IMF::CalculateExtras()
{
    if(queueCalculateExtras){

	CountExtrema();
	CalculateFrequency();
	CalculateMean();

	queueCalculateExtras = false;
    }
}

// ================= Spline Interpolators =====================
double IMF::SplineMed(double xp) const{
    if(SplineMax && SplineMin){
	return ((*SplineMax)(xp) + (*SplineMin)(xp))/2.;
    }
    return 0;
}

double IMF::SplineAmp(double xp) const{
    if(SplineMax && SplineMin){
	return ((*SplineMax)(xp) - (*SplineMin)(xp))/2.;
    }
    return 0;
}

bool IMF::CalculateSpline(double *p, uint8_t spl)
{
    // returns "estrema criteria verified?"
    ECHOS;
    
    vector<int> indexList = GetIndexes(p, spl);
    
    int size = indexList.size();

    if(size < 2){ // the number of extrema is 3 or less
	// the number of minima + maxima <= 3
	// -> means that the number of maxima or minima is 1 or less 
	ECHOF;
	return false;
    }
    
    double xt[size+2];
    double yt[size+2];

    for(int i = 1; i<size+1; i++){
	xt[i] = x[indexList[i-1]];
	yt[i] = p[indexList[i-1]];
    }
    // periodicity conditions
    // d = xt[1]-x[0]
    // x[0] - d
    xt[0] = 2*x[0] - xt[1];
    xt[size+1] = 2*x[N-1] - xt[size];
    
    yt[0] = yt[1];
    yt[size+1] = yt[size];

    if(spl & CODE_MAX){
        SplineMax = new Spline3Interpolator(size+2, xt, yt);
    } else {
	SplineMin = new Spline3Interpolator(size+2, xt, yt);
    }
    ECHOF;
    return true;
}

// ================== Main Iterator =====================
bool IMF::operator++()
{   // do one IMF iteration
    // returns "not the last IMF?"
    ECHOS;
    BENCHMARK_START(0);
    // this one is in case that you get to the last IMF, where there arent many extrema 
    if(!CalculateSpline(intrinsic, CODE_MAX) ||
       !CalculateSpline(intrinsic, CODE_MIN) ) {
	ECHOF;
	return false;
    }
    BENCHMARK_END(0);

    BENCHMARK_START(1);
    Update();
    BENCHMARK_END(1);
    
    ECHOF;
    return true;
}


// =============== Draw Functions =====================

void IMF::DrawIntrinsic() {
    ECHOS;
    
    StartGraphics();
    
    string padName = "Intrinsic Pad ";
    padName += index;
    
    TPad *pad = G.CreatePad(padName);

    TGraph *gr = new TGraph(N, x, intrinsic);
    
    gr->SetMarkerStyle(33);
    gr->SetMarkerColor(kBlack);
    gr->SetMarkerSize(0.4);
    gr->SetLineColor(kBlack);
    gr->SetLineWidth(1);
    
    G.AddObject(gr, padName, "APL");
    G.AddObject(pad);

    G.Draw();
    ECHOF;
}

void IMF::DrawSplines(double step) {
    ECHOS;
    StartGraphics();
    
    // to give an unique name
    char padName[32];
    sprintf(padName, "SplinePad %d", index);
    
    TPad *pad = G.CreatePad(padName);
    
    int stepN = (int)((x[N-1]- x[0])/step);

    // these are the points that are going to be drawn
    double xArr[stepN];
    double yArrMin[stepN];
    double yArrMax[stepN];
    double yArrMed[stepN];
    for(int i = 0; i<stepN; i++){
	// the x coordinate is the same for all graphs
	xArr[i] = x[0] + i*step;
	yArrMed[i] = SplineMed(xArr[i]);
	yArrMin[i] = (*SplineMin)(xArr[i]);
	yArrMax[i] = (*SplineMax)(xArr[i]);
    }
    
    TGraph *gSignal = new TGraph(N, x, old_intrinsic);
    gSignal->SetMarkerStyle(33);
    gSignal->SetMarkerColor(kBlack);
    gSignal->SetMarkerSize(0.6);
    gSignal->SetLineColor(kBlack);
    gSignal->SetLineWidth(1);
    gSignal->SetLineStyle(2);
    
    TGraph *gMax = new TGraph(stepN, xArr, yArrMax); 
    gMax->SetMarkerStyle(20);
    gMax->SetMarkerColor(kRed);
    gMax->SetMarkerSize(0.5);
    gMax->SetLineColor(kRed);
    gMax->SetLineWidth(2);
    
    TGraph *gMed = new TGraph(stepN, xArr, yArrMed);
    gMed->SetMarkerStyle(20);
    gMed->SetMarkerColor(kGreen);
    gMed->SetMarkerSize(0.5);
    gMed->SetLineColor(kGreen);
    gMed->SetLineWidth(2);

    TGraph *gMin = new TGraph(stepN, xArr, yArrMin);
    gMin->SetMarkerStyle(20);
    gMin->SetMarkerColor(kBlue);
    gMin->SetMarkerSize(0.5);
    gMin->SetLineColor(kBlue);
    gMin->SetLineWidth(2);

    TGraph *gVar = new TGraph(N, x, variance);
    gVar->SetMarkerStyle(20);
    gVar->SetMarkerColor(6);
    gVar->SetMarkerSize(0.5);
    gVar->SetLineColor(6);
    gVar->SetLineWidth(2);
    
    G.AddObject(gSignal, padName, "AL");
    G.AddObject(gMax, padName, "L");
    G.AddObject(gMed, padName, "L");
    G.AddObject(gMin, padName, "L");
    G.AddObject(gVar, padName, "L");
    
    G.AddObject(pad);
    G.DrawPad(padName);
    
    ECHOF;
}

TGraph* IMF::GetDrawable(int what, int color, double step) const
{
    ECHOS;
    int stepN = N;
    if(step > 0){
	stepN = N/step;
    }
 
    double xp[stepN];
    double yp[stepN];

    if(step>0){
	for(int i = 0; i<stepN; i++){
	    xp[i] = x[0] + step*i;
	}
    }else{
	for(int i = 0; i<N; i++){
	    xp[i] = x[i];
	}
    }
    
    int lineColor = kBlack;
    double lineWidth = 1;
    int lineStyle = 1;
    int markerColor = kBlack;
    double markerSize = 0.5;
    int markerStyle = 33;
     
    switch(what){
    case DRAW_SPLINEMAX:

	for(int i = 0; i<stepN; i++){
	    yp[i] = (*SplineMax)(xp[i]);
	}
	
	lineColor = kRed;
	markerColor = kRed;
	
	break;
    case DRAW_SPLINEMIN:
	for(int i = 0; i<stepN; i++){
	    yp[i] = (*SplineMin)(xp[i]);
	}
	
	lineColor = kBlue;
	markerColor = kBlue;
	
	break;
    case DRAW_SPLINEMED:
	for(int i = 0; i<stepN; i++){
	    yp[i] = SplineMed(xp[i]);
	}

	lineColor = kGreen;
	markerColor = kGreen;
	
	break;
    case DRAW_SPLINEAMP:
	for(int i = 0; i<stepN; i++){
	    yp[i] = SplineAmp(xp[i]);
	}
	
	lineColor = kMagenta;
	markerColor = kMagenta;
	
	break;
    case DRAW_FREQUENCY:
	for(int i = 0; i<N; i++){
	    yp[i] = frequency[i];
	}
	
	lineColor = kYellow-3;
	markerColor = kYellow-3;
	
	break;
    case DRAW_VARIANCE:
	for(int i = 0; i<N; i++){
	    yp[i] = variance[i];
	}

	lineColor = kOrange;
	markerSize = 0.3;
	markerColor = kOrange;
	
	break;
    case DRAW_INTRINSIC:
	for(int i = 0; i<N; i++){
	    yp[i] = intrinsic[i];
	}
	markerSize = 0.3;
	
	break;
    case DRAW_OLDINTRINSIC:
	for(int i = 0; i<N; i++){
	    yp[i] = old_intrinsic[i];
	}
	markerSize = 0.3;
	break;
    case DRAW_MEAN:
	for(int i = 0; i<N; i++){
	    yp[i] =  mean[i];
	    //printf("%d %f %f %f \n", i, y[i], intrinsic[i], mean[i]);
	}

        markerColor = kGreen;
	lineColor = kGreen;
	break;
    default: // case DRAW_DATA
	for(int i = 0; i<N; i++){
	    yp[i] = y[i];
	}

	markerSize = 0.4;
    }
    
    if(color>-1){
	lineColor = color;
	markerColor = color;
    }
    
    TGraph *temp =  new TGraph(stepN, xp, yp);

    temp->SetLineColor(lineColor);
    temp->SetLineWidth(lineWidth);
    temp->SetLineStyle(lineStyle);
    temp->SetMarkerColor(markerColor);
    temp->SetMarkerSize(markerSize);
    temp->SetMarkerStyle(markerStyle);

    ECHOF;
    return temp;
}

// ======================= Point Counters =================
int IMF::GetCount(double* yp, uint8_t t) const
{
    
    int count = 0;

    // if its requested to count the local maxima
    if(t & CODE_MAX){	
	for(int i = 1; i<N-1; i++){
	    // it adds 1 if its true, 0 if false
	    if (yp[i] > yp[i+1] &&
	        yp[i] > yp[i-1]) {
		++count;
	    } else if( yp[i] == yp[i+1] && yp[i] > yp[i-1]){
		int j;
		for(j = 1; yp[i] == yp[i+j+1] && i+j<N-1; j++);
		
		if(i+j>N) break;
		
		if(yp[i] > yp[i+j+1]){
		    ++count;
		    i+=j;
		}
	    }
	}
	// if its requested to count the local minima
    } else if(t & CODE_MIN) {	
	for(int i = 1; i<N-1; i++){
	    // it adds 1 if its true, 0 if false
	    if (yp[i] < yp[i+1] &&
		yp[i] < yp[i-1]) {
		++count;
	    } else if( yp[i] == yp[i+1] && yp[i] < yp[i-1]){
		int j;
		for(j = 1; yp[i] == yp[i+j+1] && i+j<N-1; j++);
		
		if(i+j>N) break;
		
		if(yp[i] < yp[i+j+1]){
		    ++count;
		    i+=j;
		}
	    }
	}
	// if its requested to count the zeroes
    } else if(t & CODE_ZERO){
	for(int i = 0; i<N-1; i++){
	    // it adds 1 if its true, 0 if false
	    if(yp[i]*yp[i+1]<=0){
		if(yp[i+1]==0)
		    ++count;
	    }
	}
    }
    
    return count;
}

vector<int> IMF::GetIndexes(double*yp, uint8_t t) const
{
    vector<int> vi;
    
    // if its requested to count the local maxima
    if(t & CODE_MAX){	
	for(int i = 1; i<N-1; i++){
	    // it adds 1 if its true, 0 if false
	    if (yp[i] > yp[i+1] &&
	        yp[i] > yp[i-1]) {
		vi.push_back(i);
	    } else if( yp[i] == yp[i+1] && yp[i] > yp[i-1]){
		int j;
		for(j = 1; yp[i] == yp[i+j+1] && i+j<N-1; j++);
		
		if(i+j>N) break;
		
		if(yp[i] > yp[i+j+1]){
		    vi.push_back((2*i+j)/2);
		    i+=j;
		}
	    }
	}
	// if its requested to count the local minima
    } else if(t & CODE_MIN) {	
	for(int i = 1; i<N-1; i++){
	    // it adds 1 if its true, 0 if false
	    if (yp[i] < yp[i+1] &&
		yp[i] < yp[i-1]) {
		vi.push_back(i);
	    } else if( yp[i] == yp[i+1] && yp[i] < yp[i-1]){
		int j;
		for(j = 1; yp[i] == yp[i+j+1] && i+j<N-1; j++);
		
		if(i+j>N) break;
		
		if(yp[i] < yp[i+j+1]){
		    vi.push_back((2*i+j)/2);
		    i+=j;
		}
	    }
	}
	// if its requested to count the zeroes
    } else if(t & CODE_ZERO){
	for(int i = 0; i<N-1; i++){
	    // it adds 1 if its true, 0 if false
	    if(yp[i]*yp[i+1]<=0){
		if(yp[i+1] == 0){
		    vi.push_back(i+1);
		    continue;
		}
		vi.push_back(i);
	    } 
	}
    }
    
    return vi;
}
