#include "FCtools.hpp"

//#define DEBUG
#include "Debug.hpp"
vector<string> FCtools::ReadFileToStrings(string filename)
{
    ECHOS;
    
    fstream file;

    file.open(filename, ios::in);
    
    if(!file.is_open()){
        #ifdef DEBUG
        printf("[File read Failed] %s\n", filename.c_str());
	printf("[FCtools::ReadFileToStrings] END\n\n");
        #endif
	
	return vector<string>();
    }

    string line;

    vector<string> svec;
    
    while(getline(file, line)){
        svec.push_back(line);
    }
    
    file.close();
    

    ECHOF;
    
    return svec;
}

vector<string> FCtools::ParseStringToVector(string str)
{
    // parses a string into a vector of strings, ignoring comments
    ECHOS;
    
    int delim = 0;
    bool divider = false;
    vector<string> sreturn;
    //printf(">>%s\n", str.c_str());
    for(int i = 0; i<str.size(); i++){
	// check for spaces
	if((str[i] == ' ') /*&& !divider*/){
	    string temp = str.substr(delim, i-delim);
	    sreturn.push_back(temp);
	    divider = true;
	    delim = i+1;
	}// check for commemnts with '#'
	else if (str[i] == '#' && !divider){
	    string temp = str.substr(delim, i-delim);
	    sreturn.push_back(temp);
	    
	    // delim = i+1;
	    break;
	}// check for comments with //
	else if(i<str.size()-1){
	    if(str[i] == '/' && str[i+1] == '/'){
		string temp = str.substr(delim, i-delim);
		sreturn.push_back(temp);

          	break;
                // delim = i+1;
	    }
	} else if(i == str.size()-1){
	    string temp = str.substr(delim);
	    sreturn.push_back(temp);
	    divider = false;
	    delim = i+1;
	    break;
	}	
    }

    ECHOF;
    return sreturn;
}

vector<Vec> FCtools::ReadFileTToVecs(string filename)
{
    ECHOS;
    vector<string> lines = ReadFileToStrings(filename);
    
    if(lines.empty()){
	#ifdef DEBUG
	printf("[File Empty] %s\n", filename.c_str());
	printf("[FCtools::ReadFileTToVecs] END\n\n");
        #endif
    
	return vector<Vec>();
    }

    vector<Vec> temp;
    
    for(int i = 0; i<lines.size(); i++){
	vector<string> line = ParseStringToVector(lines[i]);

	double *dArray = new double[line.size()];
	
	for(int j = 0; j<line.size(); j++){
	    dArray[j] = atof(line[j].c_str());
	    //printf("<> %s ", line[j].c_str());
	}
	if(line.size()>0)
	    temp.push_back(Vec(line.size(), dArray));

	delete[] dArray;
    }

    ECHOF;
    
    return temp;
}

void FCtools::PrintMatrixToFileT(const FCmatrix& matrix, string filename, string format)
{
    #ifdef DEBUG
    printf("\n[FCtools::PrintMatrixToFile] START\n");
    #endif
    
    FILE *file = fopen(filename.c_str(), "wt");

    if(file == nullptr){
	#ifdef DEBUG
	printf("[File write Failed] %s\n", filename.c_str());
	#endif
    }

    string matstr = matrix.PrintToString(format);
    
    #ifdef DEBUG
    printf("[Matrix to print]\n%s\n", matstr.c_str());
    #endif

    fprintf(file, "%s", matstr.c_str());

    fclose(file);

    #ifdef DEBUG
    printf("\n[FCtools::PrintMatrixToFile] END\n");
    #endif
}

    
	
	    
    


   
