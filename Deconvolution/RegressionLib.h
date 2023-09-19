//#include <cmath>
#include <list> 
#include <sstream>
//#include <vector>
//#include <iostream>
//#include <fstream>
#include "peaks.h"
#include "stat.h"
//#include <time.h>
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

struct gene2hills
{
	double fitness;
	struct gene hill1;
	struct gene hill2;
	int dim;
};

vector <double> findApproxLinewidth(vector <std::vector<double> > & listOfPoints, int numOfPoints, int dim, vector <double> currentMax);
//vector <double> findMax(list <std::vector<double> >  listOfPoints, int numOfPoints, int dim);
vector <double> findMax(vector < vector <double> >  & listOfPoints, int numOfPoints, int dim);
void printtoFile(vector <double> a, string filename);
void printtoTwoFile(vector <double> a, vector <double> b, string filename);
double distanceFromTrue(gene a, gene b);
void printGene(gene &a);
//void initializePopulation(vector <gene> & pop, double maxAlpha, double minAlpha, int dim, int popsize, double ** field, int numofPoints, vector <std::vector<double> >  listOfPoints, double mutsize);
void initializePopulation(vector <gene> & pop, double maxAlpha, double minAlpha, int dim, int popsize, vector <vector <double > > & field, int numofPoints, double mutsize, std::vector <double> max, std::vector <double> width, vector <double> edges);
void Repair(gene& g);
void findStats4D(double * MaxH, double * MinH, double * SDx, double * SDy, double * SDz, double * SDw, vector <std::vector<double> > & listOfPoints, int numOfPoints, double * minC1, double * maxC1, double * minC2, double * maxC2, double * minC3, double * maxC3, double * minC4, double * maxC4, int dim);
void findStats4D(double * MaxH, double * MinH, double * SDx, double * SDy, double * SDz, double * SDw, vector <std::vector<double> > & listOfPoints, int numOfPoints, double * minC1, double * maxC1, double * minC2, double * maxC2, double * minC3, double * maxC3, double * minC4, double * maxC4);
void tokenize(std::string const str, const char delim, std::vector<double> * out);
void ReadData(string input1, vector <std::vector<double> > * listOfPoints);
double findFitness(double scale, int dim, gene a, int numofPoints, vector < vector <double> > & field, bool min, vector <double> edges);
void gene_inital(gene& g, int dim);
void printGenetoFile(gene &a, string filename);
void printGenestoFile(vector <gene> &a, string filename, double cutoff);
double findSinglePointInModel(int dim, gene a, int numofPoints, vector < vector <double> > & field, int i);
double peakVolume(gene a);
double vL(gene a);
double vG(gene a);
//2 hill functions
void gene2hills_inital(gene2hills& g, int dim);
double findFitness2hills(double scale, int dim, gene2hills a, int numofPoints, double ** field, bool min);
void initializePopulation2Hills(vector <gene2hills> & pop, double boundHTop, double boundHBot, double maxAlpha, double minAlpha, double maxC1, double minC1, double maxC2, double minC2, double maxC3, double minC3, double maxC4, double minC4,  double maxW1, double minW1, double maxW2, double minW2, double maxW3, double minW3, double maxW4, double minW4, int dim, int popsize, double ** field, int numofPoints, vector <std::vector<double> > & listOfPoints, double mutsize, vector <double> edges);
void print2HillGenestoFile(vector <gene2hills> &a, string filename);
