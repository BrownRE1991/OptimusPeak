#include "RegressionLib.h"

gene DE(std::vector <std::vector <double> > cluster, int numOfPoints, int dim, std::vector <double> max, std::vector <double> width);
void nextGen(gene & a, int i, int r1, int r2, vector <gene> & pop, double scale, int dim, int numofPoints, vector <vector <double> > & field, double mutsize, vector <double> edges);
void DEregenerate(vector <gene> & pop, double scale, int dim, int numofPoints, vector <vector <double> > & field, double mutsize, vector <double> edges);
vector <double> meanFit(vector <gene> population);
vector <double> meanFit(vector <double> population);
