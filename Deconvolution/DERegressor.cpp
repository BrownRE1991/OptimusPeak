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
	//return(a);
}

void nextGen2Hills(gene2hills & a, int i, int r1, int r2, vector <gene2hills> & pop, double scale, int dim, int numofPoints, double ** field, double mutsize)
{
    gene2hills_inital(a, dim);
	double mod = Gauss(mutsize);
	a.hill1.height = pop[i].hill1.height + mod*(pop[r1].hill1.height - pop[r2].hill1.height);
	mod = Gauss(mutsize);
	a.hill2.height = pop[i].hill2.height + mod*(pop[r1].hill2.height - pop[r2].hill2.height);
	
	mod = Gauss(mutsize);
	a.hill1.alpha = pop[i].hill1.alpha + mod*(pop[r1].hill1.alpha - pop[r2].hill1.alpha);
	mod = Gauss(mutsize);
	a.hill2.alpha = pop[i].hill2.alpha + mod*(pop[r1].hill2.alpha - pop[r2].hill2.alpha);
	
	for(int x = 0; x < a.dim; x++)
	{
		mod = Gauss(mutsize);
		a.hill1.center[x] = pop[i].hill1.center[x] + mod*(pop[r1].hill1.center[x] - pop[r2].hill1.center[x]);
		mod = Gauss(mutsize);
		a.hill2.center[x] = pop[i].hill2.center[x] + mod*(pop[r1].hill2.center[x] - pop[r2].hill2.center[x]);
		mod = Gauss(mutsize);
		a.hill1.width[x] = pop[i].hill1.width[x] + mod*(pop[r1].hill1.width[x] - pop[r2].hill1.width[x]);
		mod = Gauss(mutsize);
		a.hill2.width[x] = pop[i].hill2.width[x] + mod*(pop[r1].hill2.width[x] - pop[r2].hill2.width[x]);
	}
	
	Repair(a.hill1);
	Repair(a.hill2);
	//cout << "\n\nNext Gene:\n";
	a.fitness = findFitness2hills(scale, dim, a, numofPoints, field, true);
	a.hill1.fitness = a.fitness;
	a.hill2.fitness = a.fitness;
	//printGene(a.hill1);
	//printGene(a.hill2);
	//return(a);
}

void DEregenerate(vector <gene> & pop, double scale, int dim, int numofPoints, double ** field, double mutsize)
{
	int r1 = 0;
	int r2 = 0;
	gene v;
    //cout << "in here! popsize: " << pop.size() << "\n";
	for(int i = 0; i < pop.size(); i++)
	{
        //cout << "more in here!\n";
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
        //cout << "more in here1!\n";
		nextGen(v, i, r1, r2, pop, scale, dim, numofPoints, field, mutsize);
        //cout << "more in here2!\n";
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
        //cout << "more in here3!\n";
	}
//    for(int i = 0; i < dim; i++)
//    {
//        v.center.pop_back();
//        v.width.pop_back();
//    }
    //cout << "more in here4!\n";
}

void DEregenerate2Hills(vector <gene2hills> & pop, double scale, int dim, int numofPoints, double ** field, double mutsize)
{
	int r1 = 0;
	int r2 = 0;
	gene2hills v;
	gene2hills_inital(v, pop[0].dim);
    //cout << "in here! popsize: " << pop.size() << "\n";
	for(int i = 0; i < pop.size(); i++)
	{
        //cout << "more in here!\n";
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
        //cout << "more in here1!\n";
		nextGen2Hills(v, i, r1, r2, pop, scale, dim, numofPoints, field, mutsize);
        //cout << "more in here2!\n";
		//cout << "Is fitness better?\n";
		//cout << v.fitness << " " << pop[i].fitness << " " << (v.fitness <= pop[i].fitness) << "\n";
		if(v.fitness <= pop[i].fitness)
		{
			pop[i].fitness = v.fitness;
            pop[i].hill1.fitness = v.hill1.fitness;
			pop[i].hill2.fitness = v.hill2.fitness;
            pop[i].hill1.height = v.hill1.height;
			pop[i].hill2.height = v.hill2.height;
            pop[i].hill1.alpha = v.hill1.alpha;
			pop[i].hill2.alpha = v.hill2.alpha;
            for(int y = 0; y < dim; y++)
            {
                pop[i].hill1.center[y] = v.hill1.center[y];
                pop[i].hill1.width[y] = v.hill1.width[y];
                pop[i].hill2.center[y] = v.hill2.center[y];
                pop[i].hill2.width[y] = v.hill2.width[y];
            }
		}
		//cout << "Next Gen Fitness: " << pop[i].fitness << "\n";
        //cout << "more in here3!\n";
	}
   // for(int i = 0; i < dim; i++)
//    {
//        v.center.pop_back();
//        v.width.pop_back();
//    }
    //cout << "more in here4!\n";
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
	
vector <double> meanFit(vector <gene2hills> population)
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
	//int generations = 20000;
	int generations = 300000;
	//int generations = 100;
	double mutsize = 0.5;
	string filename;
    if(argc > 1)
    {
        filename = argv[1];
    }
    //cout << generations << "\n";
	list <std::vector<double> >  listOfPoints;
    int dim = 0;
    int numOfPoints = 0;
    //vector <double> a1Distances(sets);
    ReadData(filename, &listOfPoints);
    dim = (int)(listOfPoints.begin()->size() - 1);
	//cout << dim << "\n";
	if(dim == -1)
	{
		cout <<"ERROR! FILE EMPTY!!!\n";
		return(1);
	}
    numOfPoints = (int)(listOfPoints.size());
	int fileNum = 0;
	// if(numOfPoints < 20)
// 	{
// 		cout << "Not enough points.\n";
// 		return(0);
// 	}
	//cout << numOfPoints << "\n";
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
			//cout << "popsize:"  << popsize << "\n";
	    }
	    if(argc > 3)
	    {
	        mutsize = stof(argv[3]);
			//cout << "mutsize:"  << mutsize << "\n";
		}
		if(argc > 4)
		{
			numofModels = stof(argv[4]);
			if(numofModels > 2){
				cout << "Error. Can't do more than 2 hills yet.\n";
			}
			//cout << "numofModels:"  << numofModels << "\n";
		}
		if(argc > 5)
		{
			fileNum = stoi(argv[5]);
			//cout << fileNum << "\n";
		}
    }
	if(numOfPoints < 20*numofModels)
	{
		//this print statement will be commented out by default. It is good for debugging, not so good for concurrent running.
		//cout << "Not enough points.\n";
		return(0);
	}
    // if(argc > 2)
//     {
//         trueGene.alpha = stof(argv[2]);
//     }
//     if(argc > 3)
//     {
//         trueGene.height = stof(argv[3]);
//     }
// 	if(dim == 2)
// 	{
// 	    if(argc > 4)
// 	    {
// 	        trueGene.center[0] = stof(argv[4]);
// 	    }
// 	    if(argc > 5)
// 	    {
// 	        trueGene.center[1] = stof(argv[5]);
// 	    }
// 	    if(argc > 6)
// 	    {
// 	        trueGene.width[0] = stof(argv[6]);
// 	    }
// 	    if(argc > 7)
// 	    {
// 	        trueGene.width[1] = stof(argv[7]);
// 	    }
// 		if(argc > 9)
// 		{
// 			popsize = stoi(argv[9]);
// 			cout << "popsize: " << popsize << "\n";
// 		}
// 		if(argc > 10)
// 		{
// 			mutsize = stof(argv[10]);
// 			cout << "mutsize: " << mutsize << "\n";
// 		}
// 	}
	//cout << "HERE\n";
	
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
	//cout << "HERE2\n";
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
	 //cout << "HERE3\n";
	
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
	//cout << scale << "\n";
	//cout << "HERE4\n";
	
 	//trueGene.fitness = findFitness(scale, dim, trueGene, numOfPoints, data, true);
 	//cout << "True Solution:\n";
 	//printGene(trueGene);
	
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
	
	//cout << "Boundries: \n";
	//cout << boundHTop << " " << boundHBot << "\n";
	//cout << maxAlpha << " " << minAlpha << "\n";
	//cout << maxC1 << " " << minC1 << "\n";
	//cout << maxC2 << " " << minC2 << "\n";
	
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
	//cout << maxW1 << " " << minW1 << "\n";
	//cout << maxW2 << " " << minW2 << "\n";
    
    //sum1 = 0;
    //sum2 = 0;
	
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
	        //cout << "here0\n";
			dist = 0;
	        initializePopulation(population, maxAlpha, minAlpha, dim, popsize, data, numOfPoints, listOfPoints, mutsize);
	        //cout << "here1\n";
		    for(int i = 0; i < generations; i++)
		    {
	            //cout << "here1 gen: " << i << "\n";
				DEregenerate(population, scale, dim, numOfPoints, data, mutsize);
	            //cout << "here1.5\n";
				// for(int x = 0; x < population.size(); x++)
// 				{
// 					if(population[x].fitness < bestGene.fitness )
// 					{
// 						bestGene.fitness = population[x].fitness;
// 						bestGene.height = population[x].height;
// 						bestGene.alpha = population[x].alpha;
// 						for(int y = 0; y < dim; y++)
// 						{
// 							bestGene.center[y] = population[x].center[y];
// 							bestGene.width[y] = population[x].width[y];
// 						}
// 					}
// 					//cout << "Best Fitness: " << bestGene.fitness << "\n";
// 				}
				//cout << "Best Fitness: " << bestGene.fitness << "\n";
				//printGene(bestGene);
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
						//cout << "Best Fitness: " << bestGene.fitness << "\n";
					}
					vector <double> stat = meanFit(population);
					//cout << stat[0] << "\n";
					variances2.push_back(stat[0]);
					variances.push_back(stat[1]);
					// if((best.size() > 20) & (i%20000 == 0))
// 					{
// 						if(abs(bestGene.fitness - best[best.size()-21]) > 0.00001)
// 						{
// 							cout << i << "\n";
// 							generations = generations + 30000;
// 						}
// 					}
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
				//cout << "Best Fitness: " << bestGene.fitness << "\n";
			}
			//cout << "Done Run\n";
			// vector <double> stat = meanFit(population);
// 			//cout << "var: " << var << "\n";
// 			variances.push_back(stat[1]);
			
			// vector <double> fitnesses2;
// 			sum1 = 0;
// 			sum2 = 0;
// 			for(int y = 0; y < population.size(); y++)
// 			{
// 				fitnesses2.push_back(population[y].center[1]);
// 				sum1 = sum1 + population[y].center[1];
// 			}
// 			mean = sum1/fitnesses2.size();
// 			//cout << "Mean: " << mean << "\n";
// 			for(int y = 0; y < fitnesses2.size(); y++)
// 			{
// 				sum2 = sum2 + (fitnesses2[y] - mean)*(fitnesses2[y] - mean);
// 			}
// 			var = sum2/fitnesses2.size();
// 			//cout << "var: " << var << "\n";
// 			variances2.push_back(var);
			//cout << "Best Fitness: " << bestGene.fitness << "\n";
		    //dist = distanceFromTrue(trueGene, bestGene);
			//a1Distances[z] = dist;
			bestGenes.push_back(bestGene);
			gene_inital(bestGene, dim);
		}
		
		//cout << "got to here\n";
		vector <double> stat = meanFit(variances);
		sum1 = 0;
		sum2 = 0;
		for(int y = 0; y < variances.size(); y++)
		{
			sum1 = sum1 + variances[y];
		}
		mean = sum1/variances.size();
		// sum1 = 0;
// 		sum2 = 0;
// 		for(int y = 0; y < variances2.size(); y++)
// 		{
// 			sum1 = sum1 + variances2[y];
// 		}
		//double mean2 = sum1/variances2.size();
		//cout << "var mean: " << stat[0] << "\n";
		//cout << "var mean2: " << mean2 << "\n";
		
	    fstream aus;  //output file
		filename = "MeanOverTime.txt";
	    aus.open(filename,ios::out);
	    aus.precision(10);
		for(int y = 0; y < variances2.size(); y++)
		{
			aus << variances2[y] << "\n";
		}
		aus.close();
		
	    //fstream aus;  //output file
		filename = "VarOverTime.txt";
	    aus.open(filename,ios::out);
	    aus.precision(10);
		for(int y = 0; y < variances.size(); y++)
		{
			aus << variances[y] << "\n";
		}
		aus.close();
		filename = "BestOverTime/BestOverTime_1model_" + to_string(fileNum) + ".txt";
		// if(numofModels == 1)
// 		{
// 			filename = "BestOverTime/BestOverTime_1model_" + str(fileNum) + ".txt";
// 		}
		// elseif(numofModels == 2)
// 		{
// 			filename = "BestOverTime/BestOverTime_2models_" + str(fileNum) + ".txt";
// 		}
		//cout << filename << "\n";
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
		//cout << "num of models = 2\n";
		//generations = generations*2;
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
	    //cout << "here\n";
		for(int z = 0; z < sets; z++)
 		{
 	        //cout << "here0\n";
			//dist = 0;
 	        initializePopulation2Hills(population, boundHTop, boundHBot, maxAlpha, minAlpha, maxC1, minC1, maxC2, minC2, maxC3, minC3, maxC4, minC4,  maxW1, minW1, maxW2, minW2, maxW3, minW3, maxW4, minW4, dim, popsize, data, numOfPoints, listOfPoints, mutsize);
	        //cout << "here1\n";
			// for(int x = 0; x < population.size(); x++)
// 			{
// 				cout << "\nHill1:\n";
// 				printGene(population[x].hill1);
// 				cout << "\nHill2:\n";
// 				printGene(population[x].hill2);
// 			}
 		    for(int i = 0; i < generations; i++)
 		    {
 	            //cout << "here1 gen: " << i << "\n";
 				DEregenerate2Hills(population, scale, dim, numOfPoints, data, mutsize);
	            //cout << "here1.5\n";
				// for(int x = 0; x < population.size(); x++)
// 				{
// 					if(population[x].fitness < bestGene.fitness )
// 					{
// 						bestGene.fitness = population[x].fitness;
// 						bestGene.hill1.height = population[x].hill1.height;
// 						bestGene.hill2.height = population[x].hill2.height;
// 						bestGene.hill1.alpha = population[x].hill1.alpha;
// 						bestGene.hill2.alpha = population[x].hill2.alpha;
// 						for(int y = 0; y < dim; y++)
// 						{
// 							bestGene.hill1.center[y] = population[x].hill1.center[y];
// 							bestGene.hill1.width[y] = population[x].hill1.width[y];
// 							bestGene.hill2.center[y] = population[x].hill2.center[y];
// 							bestGene.hill2.width[y] = population[x].hill2.width[y];
// 						}
// 					}
// 					// cout << "Best Genes:\n";
// // 					printGene(bestGene.hill1);
// // 					printGene(bestGene.hill2);
// 				}
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
						// cout << "Best Genes:\n";
	// 					printGene(bestGene.hill1);
	// 					printGene(bestGene.hill2);
					}
					vector <double> stat = meanFit(population);
					//cout << stat[0] << "\n";
					variances2.push_back(stat[0]);
					variances.push_back(stat[1]);
					// if((best.size() > 20) & (i%20000 == 0))
// 					{
// 						if(abs(bestGene.fitness - best[best.size()-21]) > 0.00001)
// 						{
// 							cout << i << "\n";
// 							generations = generations + 30000;
// 						}
// 					}
					best.push_back(bestGene.fitness);
				}
				//cout << "Best Fitness: " << bestGene.fitness << "\n";
				// cout << "Best Genes:\n";
				// printGene(bestGene.hill1);
				// printGene(bestGene.hill2);
 		    }
			//cout << "Best Fitness: " << bestGene.fitness << "\n";
// 		    //dist = distanceFromTrue(trueGene, bestGene);
// 			//a1Distances[z] = dist;
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
		// sum1 = 0;
// 		sum2 = 0;
// 		for(int y = 0; y < variances2.size(); y++)
// 		{
// 			sum1 = sum1 + variances2[y];
// 		}
		//double mean2 = sum1/variances2.size();
		//cout << "var mean: " << stat[0] << "\n";
		//cout << "var mean2: " << mean2 << "\n";
		
	    fstream aus;  //output file
		filename = "MeanOverTime2.txt";
	    aus.open(filename,ios::out);
	    aus.precision(10);
		for(int y = 0; y < variances2.size(); y++)
		{
			aus << variances2[y] << "\n";
		}
		aus.close();
		
	    //fstream aus;  //output file
		filename = "VarOverTime2.txt";
	    aus.open(filename,ios::out);
	    aus.precision(10);
		for(int y = 0; y < variances.size(); y++)
		{
			aus << variances[y] << "\n";
		}
		aus.close();
		filename = "BestOverTime/BestOverTime_2models_" + to_string(fileNum) + ".txt";
		//cout << filename << "\n";
		// filename = "BestOverTime2.txt";
	    aus.open(filename,ios::out);
	    aus.precision(10);
		for(int y = 0; y < best.size(); y++)
		{
			aus << best[y] << "\n";
		}
		aus.close();
		
		// string outfileName1 = "";
	// 	if(filename[0] != '/')
	// 	{
	// 		outfileName1 = "Output/DEDistances" + filename;
	// 	}
	// 	else
	// 	{
	// 		outfileName1 += "/Users/Pixie/Desktop/OptmizationExperiments/TestingonHills/DE/Output/DEDistances";
	// 		outfileName1 += argv[8];
	// 	}
	// 	printtoFile(a1Distances, outfileName1);
	
		//string outfileName2 = "Model" + string(argv[5]) + ".txt";
		
		string outfileName2 = "Models_2/Models_2ASRG" + to_string(fileNum) + ".txt";
		print2HillGenestoFile(bestGenes, outfileName2);
		outfileName2 = "../GoodnessOfFitTests/Models_2/Models_2ASRG" + to_string(fileNum) + ".txt";
		print2HillGenestoFile(bestGenes, outfileName2);
		// for(int x = 0; x < bestGenes.size(); x++)
	// 	{
			//printGene(bestGene);
			//printGene(bestGenes[x]);
			// if(bestGenes[x].fitness < bestGene.fitness )
	// 		{
	// 			bestGene.fitness = bestGenes[x].fitness;
	// 			bestGene.height = bestGenes[x].height;
	// 			bestGene.alpha = bestGenes[x].alpha;
	// 			for(int y = 0; y < dim; y++)
	// 			{
	// 				bestGene.center[y] = bestGenes[x].center[y];
	// 				bestGene.width[y] = bestGenes[x].width[y];
	// 			}
	// 		}
	//	}
		//print2HillGenestoFile(bestGenes, outfileName2);
	}
	
    for(int x = 0; x < numOfPoints; x++)
    {
        free(data[x]);
    }
    free(data);
    return(0);
}
