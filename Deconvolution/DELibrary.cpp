#include "DELibrary.h"

gene DE(std::vector <std::vector <double> > cluster, int numOfPoints, int dim, std::vector <double> max, std::vector <double> width)
{
	//cout << "in DE\n";
	/*for(int x = 0; x < max.size(); x++)
	{
		cout << max[x] << " ";
	}
	cout << "\n";
	for(int x = 0; x < dim; x++)
	{
		cout << width[x] << " ";
	}
	cout << "\n";*/
	
	/*if(max[0] > 120)
	{
		cout << max[0] << "\n";
		
	}*/
	//cout << "\n";
	gene a;
	gene_inital(a, dim);
	srand((int)time(NULL));
	srand48((int)time(NULL));
	int popsize = 15;
	int generations = 700000;
	//double mutsize = 0.25;
	double mutsize1 = 0.01;
	double mutsize2 = 2;
	if(numOfPoints < 18)
	{
		return(a);
	}
	
	//vector < vector <double> > data;
	//vector <double> temp;
    // double ** data = nullptr;
//     data = (double**)malloc(numOfPoints*sizeof(double*));
//     if(data != nullptr)
//     {
//         for(int x = 0; x < numOfPoints; x++)
//         {
//             data[x] = nullptr;
//             data[x] = (double*)malloc((dim+1)*sizeof(double));
//             if(data[x] == nullptr)
//             {
//                 cout << "Memory Error\n";
//                 return(a);
//             }
//         }
//     }
	// for(int x = 0; x <= dim; x++)
// 	{
// 		temp.push_back(0.0)
// 	}
	// for(int x = 0; x < cluster.size(); x++)
// 	{
// 		vector <double> temp;
// 		for(int x = 0; x <= dim; x++)
// 		{
// 			temp.push_back(0.0);
// 		}
// 		for(int y = 0; y < cluster[x].size(); y++)
// 		{
// 			temp[y] = cluster[x][y];
// 			data[x][y] = cluster[x][y];
// 		}
// 		data.push_back(temp);
// 	}
	//cout << "HERE\n";
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
	findStats4D(&MaxH, &MinH, &SDx, &SDy, &SDz, &SDw, cluster, numOfPoints, &minC1, &maxC1, &minC2, &maxC2, &minC3, &maxC3, &minC4, &maxC4, dim);
	double scale = MaxH - MinH;
	//cout << scale << "\n";
	vector <vector < gene > > states;
	
	//cout << minC1 << " " << maxC1 << "\n";
	//cout << minC2 << " " << maxC2 << "\n";
	//cout << minC3 << " " << maxC3 << "\n";
	//cout << minC4 << " " << maxC4 << "\n";
	
	vector <double> edges;
	edges.push_back(minC1);
	edges.push_back(maxC1);
	edges.push_back(minC2);
	edges.push_back(maxC2);
	if(dim > 2)
	{
		edges.push_back(minC3);
		edges.push_back(maxC3);
	}
	if(dim > 3)
	{
		edges.push_back(minC4);
		edges.push_back(maxC4);
	}
	
	//cout << edges.size() << "\n";
	int sets = 1;
		vector <gene> population;
		gene bestGene;
		vector <gene> bestGenes;
		gene_inital(bestGene, dim);
		//double dist = 0;
		if(dim == 3)
		{
			generations = generations*2;
			
		}
		//vector <double> variances;
		//vector <double> variances2;
		//vector <double> best;
		//double sum1 = 0;
		//double sum2 = 0;
		//double mean = 0;
		//double best_Linewidth = 0;
		for(int z = 0; z < sets; z++)
		{
			//dist = 0;
	        initializePopulation(population, maxAlpha, minAlpha, dim, popsize, cluster, numOfPoints, mutsize1, max, width, edges);
		    for(int i = 0; i < generations; i++)
		    {
				//cout << i << "\n";
				if(bestGene.width[0] < 100)
				{
					DEregenerate(population, scale, dim, numOfPoints, cluster, mutsize2, edges);
					if(i%10000 == 0)
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
					}
					//vector <double> stat = meanFit(population);
					//cout << stat[0] << "\n";
					//variances2.push_back(stat[0]);
					//variances.push_back(stat[1]);
					// if((best.size() > 20) & (i%20000 == 0))
// 					{
// 						if(abs(bestGene.fitness - best[best.size()-21]) > 0.00001)
// 						{
// 							cout << i << "\n";
// 							generations = generations + 30000;
// 						}
// 					}
					//best.push_back(bestGene.fitness);
				}
				//cout << i << " ";
				//cout << generations << " ";
				//cout << generations%100000 << "\n";
				//if(i%100000 == 0)
				//{
					//cout << i << "\n";
					//if(max[0] > 120)
					//{
						/*for(int x = 0; x < population.size(); x++)
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
						}*/
						//cout << max[0] << "\n";
						//vector <gene> holder;
						//gene holder2;
						//gene_inital(holder2, dim);
						//for(int x = 0; x < population.size(); x++)
						//{
							//gene_inital(holder2, dim);
							//holder2.fitness = population[x].fitness;
							//holder2.height = population[x].height;
							//holder2.alpha = population[x].alpha;
							//for(int y = 0; y < dim; y++)
							//{
								//holder2.center[y] = population[x].center[y];
								//holder2.width[y] = population[x].width[y];
							//}
							//holder.push_back(holder2);
						//}
						//states.push_back(holder);
					//}
				//}
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
		
		//vector <double> stat = meanFit(variances);
//		sum1 = 0;
//		sum2 = 0;
//		for(int y = 0; y < variances.size(); y++)
//		{
//			sum1 = sum1 + variances[y];
//		}
//		mean = sum1/variances.size();
		
		// string filename;
// 	    fstream aus;  //output file
// 		filename = "MeanOverTime.txt";
// 	    aus.open(filename,ios::out);
// 	    aus.precision(10);
// 		for(int y = 0; y < variances2.size(); y++)
// 		{
// 			aus << variances2[y] << "\n";
// 		}
// 		aus.close();
//
// 		filename = "VarOverTime.txt";
// 	    aus.open(filename,ios::out);
// 	    aus.precision(10);
// 		for(int y = 0; y < variances.size(); y++)
// 		{
// 			aus << variances[y] << "\n";
// 		}
// 		aus.close();
// 		filename = "BestOverTime_1model_" + to_string(1) + ".txt";
// 		cout << filename << "\n";
// 	    aus.open(filename,ios::out);
// 	    aus.precision(10);
// 		for(int y = 0; y < best.size(); y++)
// 		{
// 			aus << best[y] << "\n";
// 		}
// 		aus.close();
		
		//string outfileName2 = "Decon_Model_ASRG" + to_string(1) + ".txt";
		//if(max[0] > 120)
		//{
			//string outfileName2 = "GenesOverTime";
			//for(int x = 0; x < states.size(); x++)
			//{
				//outfileName2 = "GenesOverTime" + to_string(x) + ".txt";
				//printGenestoFile(states[x], outfileName2, 1);
			//}
		//}
		//outfileName2 = "../GoodnessOfFitTests/Models/Model_ASRG" + to_string(1) + ".txt";
		//printGenestoFile(bestGenes, outfileName2);
	
    // for(int x = 0; x < numOfPoints; x++)
//     {
//         free(data[x]);
//     }
//     free(data);
	//cout << "Got to here\n";
	return(bestGenes[0]);
}

//might be worth putting only one scalar instead of one for each member of the vector
void nextGen(gene & a, int i, int r1, int r2, vector <gene> & pop, double scale, int dim, int numofPoints, vector <vector <double> > & field, double mutsize, vector <double> edges)
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
	a.fitness = findFitness(scale, dim, a, numofPoints, field, true, edges);
	//cout << a.fitness << "\n";
	//cout << "NextGen complete\n";
}

void DEregenerate(vector <gene> & pop, double scale, int dim, int numofPoints, vector <vector <double> > & field, double mutsize, vector <double> edges)
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
		nextGen(v, i, r1, r2, pop, scale, dim, numofPoints, field, mutsize, edges);
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

void outputVector(vector < vector <double> > data, string filename, int size)
{
    fstream aus;  //output file
    aus.open(filename,ios::out);
	aus.precision(7);
	string outstring;
	for(int x = 0; x < size; x++)
	{
		outstring = "";
		for(int y = 0; y < data[x].size(); y++)
		{
			outstring = outstring + to_string(data[x][y]) + " ";
		}
		outstring = outstring + "\n";
		aus << outstring;
	}
	aus.close();
}
