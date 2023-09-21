#include <list> 
#include <sstream>
#include "peaks.h"
#include "stat.h"
#include <cstring>
#include <algorithm>

struct gene
{
    int dim;
    double height;
    double alpha;
    std::vector <double> center;
    std::vector <double> width;
	double fitness;
};

vector <double> findApproxLinewidth(vector <std::vector<double> > & listOfPoints, int numOfPoints, int dim, vector <double> currentMax);
vector <double> findMax(vector < vector <double> >  & listOfPoints, int numOfPoints, int dim);
void initializePopulation(vector <gene> & pop, double maxAlpha, double minAlpha, int dim, int popsize, vector <vector <double > > & field, int numofPoints, double mutsize, std::vector <double> max, std::vector <double> width, vector <double> edges);
void Repair(gene& g);
void findStats4D(double * MaxH, double * MinH, double * SDx, double * SDy, double * SDz, double * SDw, vector <std::vector<double> > & listOfPoints, int numOfPoints, double * minC1, double * maxC1, double * minC2, double * maxC2, double * minC3, double * maxC3, double * minC4, double * maxC4, int dim);
void tokenize(std::string const str, const char delim, std::vector<double> * out);
void ReadData(string input1, vector <std::vector<double> > * listOfPoints);
double findFitness(double scale, int dim, gene a, int numofPoints, vector < vector <double> > & field, bool min, vector <double> edges);
void gene_inital(gene& g, int dim);
void printGenestoFile(vector <gene> &a, string filename, double cutoff);
double findSinglePointInModel(int dim, gene a, int numofPoints, vector < vector <double> > & field, int i);
double peakVolume(gene a);
double vL(gene a);
double vG(gene a);
