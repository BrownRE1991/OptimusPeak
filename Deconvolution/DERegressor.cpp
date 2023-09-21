#include "RegressionLib.h"

using namespace std;

void nextGen(gene & a, int i, int r1, int r2, vector <gene> & pop, double scale, int dim, int numofPoints, double ** field, double mutsize)
{
    gene_inital(a, dim);
	double mod = Gauss(mutsize);
	a.height = pop[i].height + mod*(pop[r1].height - pop[r2].height);
	
	mod = Gauss(mutsize);
	a.alpha = pop[i].alpha + mod*(pop[r1].alpha - pop[r2].alpha);
	
	for(int x = 0; x < a.dim; x++)
	{
		mod = Gauss(mutsize);
		a.center[x] = pop[i].center[x] + mod*(pop[r1].center[x] - pop[r2].center[x]);
		mod = Gauss(mutsize);
		a.width[x] = pop[i].width[x] + mod*(pop[r1].width[x] - pop[r2].width[x]);
	}
	
	Repair(a);
	a.fitness = findFitness(scale, dim, a, numofPoints, field, true);
}

void DEregenerate(vector <gene> & pop, double scale, int dim, int numofPoints, double ** field, double mutsize)
{
	int r1 = 0;
	int r2 = 0;
	gene v;
	for(int i = 0; i < pop.size(); i++)
	{
		r1 = rand()%pop.size();
		while(r1 == i)
		{
			r1 = rand()%pop.size();
		}
		r2 = rand()%pop.size();
		while(r2 == i || r2 == r1)
		{
			r2 = rand()%pop.size();
		}
		nextGen(v, i, r1, r2, pop, scale, dim, numofPoints, field, mutsize);
		if(v.fitness <= pop[i].fitness)
		{
            pop[i].fitness = v.fitness;
            pop[i].height = v.height;
            pop[i].alpha = v.alpha;
            for(int y = 0; y < dim; y++)
            {
                pop[i].center[y] = v.center[y];
                pop[i].width[y] = v.width[y];
            }
		}
	}
}

vector <double> meanFit(vector <gene> population)
{
	vector <double> stats;
	vector <double> fitnesses;
	double sum1 = 0;
	double sum2 = 0;
	for(int y = 0; y < population.size(); y++)
	{
		fitnesses.push_back(population[y].fitness);
		sum1 = sum1 + population[y].fitness;
	}
	double mean = sum1/fitnesses.size();
	stats.push_back(mean);
	for(int y = 0; y < fitnesses.size(); y++)
	{
		sum2 = sum2 + (fitnesses[y] - mean)*(fitnesses[y] - mean);
	}
	double var = sum2/fitnesses.size();
	stats.push_back(var);
	return(stats);
}

vector <double> meanFit(vector <double> population)
{
	vector <double> stats;
	vector <double> fitnesses;
	double sum1 = 0;
	double sum2 = 0;
	for(int y = 0; y < population.size(); y++)
	{
		fitnesses.push_back(population[y]);
		sum1 = sum1 + population[y];
	}
	double mean = sum1/fitnesses.size();
	stats.push_back(mean);
	for(int y = 0; y < fitnesses.size(); y++)
	{
		sum2 = sum2 + (fitnesses[y] - mean)*(fitnesses[y] - mean);
	}
	double var = sum2/fitnesses.size();
	stats.push_back(var);
	return(stats);
}

int main(int argc, char ** argv)
{
	srand((int)time(NULL));
	int sets = 3;
	int popsize = 10;
	int generations = 300000;
	double mutsize = 0.5;
	string filename;
    if(argc > 1)
    {
        filename = argv[1];
    }
	list <std::vector<double> >  listOfPoints;
    int dim = 0;
    int numOfPoints = 0;
    ReadData(filename, &listOfPoints);
    dim = (int)(listOfPoints.begin()->size() - 1);
	if(dim == -1)
	{
		cout <<"ERROR! FILE EMPTY!!!\n";
		return(1);
	}
    numOfPoints = (int)(listOfPoints.size());
	int fileNum = 0;
	gene trueGene;
	gene_inital(trueGene, dim);
	int numofModels = 1;
	if(argc > 2)
	{
	    {
			if(stof(argv[2]) < 5)
			{
				cout << "Population Size not Large Enough. Try > 5.\n";
				return(1);
			}
	        popsize = stof(argv[2]);
	    }
	    if(argc > 3)
	    {
	        mutsize = stof(argv[3]);
		}
		if(argc > 4)
		{
			numofModels = stof(argv[4]);
			if(numofModels > 2){
				cout << "Error. Can't do more than 2 hills yet.\n";
			}
		}
		if(argc > 5)
		{
			fileNum = stoi(argv[5]);
		}
    }
	if(numOfPoints < 20*numofModels)
	{
		return(0);
	}
	
    double ** data = nullptr;
    data = (double**)malloc(numOfPoints*sizeof(double*));
    if(data != nullptr)
    {
        for(int x = 0; x < numOfPoints; x++)
        {
            data[x] = nullptr;
            data[x] = (double*)malloc((dim+1)*sizeof(double));
            if(data[x] == nullptr)
            {
                cout << "Memory Error\n";
                return(1);
            }
        }
    }
    std::list <std::vector<double> > :: iterator it;
    std::vector<double>:: iterator it2;
    int count1 = 0;
    int count2 = 0;
    for(it = listOfPoints.begin(); it != listOfPoints.end(); it++)
    {
        if(count1 < numOfPoints)
        {
            count2 = 0;
            for (it2 = it->begin(); it2 != it->end(); it2++)
            {
                if(count2 <= dim)
                {
                    if(count1 < numOfPoints)
                    {
                        data[count1][count2] = *it2;
                        count2++;
                    }
                }
            }
        }
        count1++;
     }
	
	double MaxH = 0;
	double MinH = 1000000;
	double SDx = 0;
	double SDy = 0;
	double SDz = 0;
	double SDw = 0;
	double minC1 = 1000000;
	double maxC1 = 0;
	double minC2 = 1000000;
	double maxC2 = 0;
	double minC3 = 1000000;
	double maxC3 = 0;
	double minC4 = 1000000;
	double maxC4 = 0;
	double maxAlpha = 2.0;
	double minAlpha = -1.0;
	findStats4D(&MaxH, &MinH, &SDx, &SDy, &SDz, &SDw, listOfPoints, numOfPoints, &minC1, &maxC1, &minC2, &maxC2, &minC3, &maxC3, &minC4, &maxC4);
	double scale = MaxH - MinH;
	
	double boundHTop = 2*MaxH;
	double boundHBot = 0;
	
	if(MinH >= 0)
	{
		boundHBot = MinH/2;
	}
	if(MinH < 0)
	{
		boundHBot = MinH*2;
	}
	
	double maxW1 = 2*SDx;
	double minW1 = 0.0;
	double maxW2 = 0;
	double minW2 = 0;
	double maxW3 = 0;
	double minW3 = 0;
	double maxW4 = 0;
	double minW4 = 0;
	if(dim > 1)
	{
		maxW2 = 2*SDy;
		minW2 = 0.0;
	}
	if(dim > 2)
	{
		maxW3 = 2*SDz;
		minW3 = 0.0;
	}
	if(dim > 3)
	{
		maxW4 = 2*SDw;
		minW4 = 0.0;
	}
	
	if(numofModels == 1)
	{
		vector <gene> population;
		gene bestGene;
		vector <gene> bestGenes;
		gene_inital(bestGene, dim);
		double dist = 0;
	    //cout << "here\n";
		if(dim == 3)
		{
			generations = generations*2;
			
		}
		vector <double> variances;
		vector <double> variances2;
		vector <double> best;
		double sum1 = 0;
		double sum2 = 0;
		double mean = 0;
		for(int z = 0; z < sets; z++)
		{
			dist = 0;
	        initializePopulation(population, maxAlpha, minAlpha, dim, popsize, data, numOfPoints, listOfPoints, mutsize);
		    for(int i = 0; i < generations; i++)
		    {
				DEregenerate(population, scale, dim, numOfPoints, data, mutsize);
				if(i%1000 == 0)
				{
					for(int x = 0; x < population.size(); x++)
					{
						if(population[x].fitness < bestGene.fitness )
						{
							bestGene.fitness = population[x].fitness;
							bestGene.height = population[x].height;
							bestGene.alpha = population[x].alpha;
							for(int y = 0; y < dim; y++)
							{
								bestGene.center[y] = population[x].center[y];
								bestGene.width[y] = population[x].width[y];
							}
						}
					}
					vector <double> stat = meanFit(population);
					variances2.push_back(stat[0]);
					variances.push_back(stat[1]);
					best.push_back(bestGene.fitness);
				}
		    }
			for(int x = 0; x < population.size(); x++)
			{
				if(population[x].fitness < bestGene.fitness )
				{
					bestGene.fitness = population[x].fitness;
					bestGene.height = population[x].height;
					bestGene.alpha = population[x].alpha;
					for(int y = 0; y < dim; y++)
					{
						bestGene.center[y] = population[x].center[y];
						bestGene.width[y] = population[x].width[y];
					}
				}
			}
			bestGenes.push_back(bestGene);
			gene_inital(bestGene, dim);
		}
		
		vector <double> stat = meanFit(variances);
		sum1 = 0;
		sum2 = 0;
		for(int y = 0; y < variances.size(); y++)
		{
			sum1 = sum1 + variances[y];
		}
		mean = sum1/variances.size();
		
	    fstream aus;  //output file
		filename = "MeanOverTime.txt";
	    aus.open(filename,ios::out);
	    aus.precision(10);
		for(int y = 0; y < variances2.size(); y++)
		{
			aus << variances2[y] << "\n";
		}
		aus.close();
		
		filename = "VarOverTime.txt";
	    aus.open(filename,ios::out);
	    aus.precision(10);
		for(int y = 0; y < variances.size(); y++)
		{
			aus << variances[y] << "\n";
		}
		aus.close();
		filename = "BestOverTime/BestOverTime_1model_" + to_string(fileNum) + ".txt";
	    aus.open(filename,ios::out);
	    aus.precision(10);
		for(int y = 0; y < best.size(); y++)
		{
			aus << best[y] << "\n";
		}
		aus.close();
		
		
		string outfileName2 = "Models/Model_ASRG" + to_string(fileNum) + ".txt";
		printGenestoFile(bestGenes, outfileName2);
		outfileName2 = "../GoodnessOfFitTests/Models/Model_ASRG" + to_string(fileNum) + ".txt";
		printGenestoFile(bestGenes, outfileName2);
	}
	if(numofModels == 2)
	{
		generations = 2*generations;
		vector <double> variances;
		vector <double> variances2;
		vector <double> best;
		double sum1 = 0;
		double sum2 = 0;
		double mean = 0;
		popsize = popsize*2;
		vector <gene2hills> population;
		gene2hills bestGene;
		vector <gene2hills> bestGenes;
		gene2hills_inital(bestGene, dim);
		double dist = 0;
		if(dim == 3)
		{
			generations = generations*2;
		}

		for(int z = 0; z < sets; z++)
 		{
 	        initializePopulation2Hills(population, boundHTop, boundHBot, maxAlpha, minAlpha, maxC1, minC1, maxC2, minC2, maxC3, minC3, maxC4, minC4,  maxW1, minW1, maxW2, minW2, maxW3, minW3, maxW4, minW4, dim, popsize, data, numOfPoints, listOfPoints, mutsize);
 		    for(int i = 0; i < generations; i++)
 		    {
 				DEregenerate2Hills(population, scale, dim, numOfPoints, data, mutsize);
				if(i%1000 == 0)
				{
					for(int x = 0; x < population.size(); x++)
					{
						if(population[x].fitness < bestGene.fitness )
						{
							bestGene.fitness = population[x].fitness;
							bestGene.hill1.height = population[x].hill1.height;
							bestGene.hill2.height = population[x].hill2.height;
							bestGene.hill1.alpha = population[x].hill1.alpha;
							bestGene.hill2.alpha = population[x].hill2.alpha;
							for(int y = 0; y < dim; y++)
							{
								bestGene.hill1.center[y] = population[x].hill1.center[y];
								bestGene.hill1.width[y] = population[x].hill1.width[y];
								bestGene.hill2.center[y] = population[x].hill2.center[y];
								bestGene.hill2.width[y] = population[x].hill2.width[y];
							}
						}
					}
					vector <double> stat = meanFit(population);
					variances2.push_back(stat[0]);
					variances.push_back(stat[1]);
					best.push_back(bestGene.fitness);
				}
 		    }
			bestGenes.push_back(bestGene);
			gene2hills_inital(bestGene, dim);
 		}	
		
		vector <double> stat = meanFit(variances);
		sum1 = 0;
		sum2 = 0;
		for(int y = 0; y < variances.size(); y++)
		{
			sum1 = sum1 + variances[y];
		}
		mean = sum1/variances.size();
		
	    fstream aus;  //output file
		filename = "MeanOverTime2.txt";
	    aus.open(filename,ios::out);
	    aus.precision(10);
		for(int y = 0; y < variances2.size(); y++)
		{
			aus << variances2[y] << "\n";
		}
		aus.close();
		
		filename = "VarOverTime2.txt";
	    aus.open(filename,ios::out);
	    aus.precision(10);
		for(int y = 0; y < variances.size(); y++)
		{
			aus << variances[y] << "\n";
		}
		aus.close();
		filename = "BestOverTime/BestOverTime_2models_" + to_string(fileNum) + ".txt";
	    aus.open(filename,ios::out);
	    aus.precision(10);
		for(int y = 0; y < best.size(); y++)
		{
			aus << best[y] << "\n";
		}
		aus.close();
		
		string outfileName2 = "Models_2/Models_2ASRG" + to_string(fileNum) + ".txt";
		print2HillGenestoFile(bestGenes, outfileName2);
		outfileName2 = "../GoodnessOfFitTests/Models_2/Models_2ASRG" + to_string(fileNum) + ".txt";
		print2HillGenestoFile(bestGenes, outfileName2);
	}
	
    for(int x = 0; x < numOfPoints; x++)
    {
        free(data[x]);
    }
    free(data);
    return(0);
}
