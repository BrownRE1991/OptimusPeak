//#include <cmath>
//#include <list> 
//#include <sstream>
//#include <vector>
//#include <iostream>
//#include <fstream>
//#include "peaks.h"
//#include "stat.h"
//#include <time.h>
//#include <cstring>
#include <thread>
#include "DELibrary.h"
//#include "RegressionLib.h"

void tokenize(std::string const str, const char delim, std::vector<double> * out);
void ReadData(string input1, list <std::vector<double> > * listOfPoints);
//void outputDataset(double ** data, int numOfPoints, int dim, string filename);
vector <double> findMaxSecond(vector <std::vector<double> > & listOfPoints, int numOfPoints, int dim, vector <double> currentMax, vector <double> widths);
double dist(vector <double> vec1, vector <double> vec2);
vector <double> findApproxLinewidth(list <std::vector<double> >  listOfPoints, int numOfPoints, int dim, vector <double> currentMax);
double findMin(vector <std::vector<double> > & listOfPoints, int numOfPoints, int dim);
vector <double> findMax(vector <std::vector<double> > & listOfPoints, int numOfPoints, int dim);
vector <double> findMaxOther(vector <std::vector<double> > & listOfPoints, int numOfPoints, int dim, vector < vector <double> > currentMax, vector < vector <double> > widths);
std::vector <gene> kmodels(vector <std::vector<double> > & data, int numOfPoints, int dim, int k, vector < vector <double> > maxs, vector < vector <double> > widths, double cutoff, double Thresh);
double dist(vector <double> vec1, double * vec2);
void outputVector(vector < vector <double> > data, string filename, int size);
void outputVector2(vector < vector <double> > & data, string filename, int size, int dim);
void outputDatasetClustersV1(vector < vector <double> > & data, int numOfPoints, int dim, std::vector <std::vector <std::vector <double> > > clusters);
void outputDatasetClustersV2(vector < vector <double> > & data, int numOfPoints, int dim, std::vector <std::vector <std::vector <double> > > clusters);
void outputVector3(vector < vector <double> > & data, string filename, int size, int dim);
void removeElement(std::vector<gene>& vec, int index);
void removeElement(vector < vector <double> >& vec, int index);
void outputDataset(vector <vector <double> > & data, int numOfPoints, int dim, string filename);

int main(int argc, const char * argv[])
{
	double cutoff = 0.2;
	srand((int)time(NULL));
	vector <std::vector<double> >  listOfPoints;
	int numOfPoints = 0;
	int numOfPoints2 = 0;
	int dim = 0;
	int file = 0;
	
	string filename;
    if(argc > 1)
    {
        filename = argv[1];
    }
    if(argc > 2)
    {
        file = stoi(argv[2]);
    }
    //cout << filename << endl;
    ReadData(filename, &listOfPoints);
    dim = (int)(listOfPoints.begin()->size() - 1);
	if(dim == -1)
	{
		cout <<"ERROR! FILE EMPTY!!!\n";
		return(1);
	}
    numOfPoints = (int)(listOfPoints.size());
	if(numOfPoints < 20)
	{
		cout << "Not enough points.\n";
		return(0);
	}
	cout << "numOfPoints: " << numOfPoints << "\n";
	cout << "dim: " << dim << "\n";
	
	vector <std::vector<double> >  listOfPoints2;
	
	std::vector <std::vector<double> > :: iterator it;
	double  max = 0;
    for(it = listOfPoints.begin(); it != listOfPoints.end(); it++)
    {
        if((*it)[dim] > max)
		{
			max = (*it)[dim];
		}
    }
	
    std::vector<double>:: iterator it2;
    for(it = listOfPoints.begin(); it != listOfPoints.end(); it++)
    {
        if(abs((*it)[dim]) > abs(max*(1.0/30.0)))
        {
            listOfPoints2.push_back((*it));
        }
     }
	 numOfPoints2 = int(listOfPoints2.size());
	 //cout << listOfPoints2.size() << "\n";

    // double ** data = nullptr;
//     data = (double**)malloc((listOfPoints.size())*sizeof(double*));
//     if(data != nullptr)
//     {
//         for(int x = 0; x < listOfPoints.size(); x++)
//         {
//             data[x] = nullptr;
//             data[x] = (double*)malloc((dim+1)*sizeof(double));
//             if(data[x] == nullptr)
//             {
//                 cout << "Memory Error\n";
//                 return(1);
//             }
//         }
//     }
	
    int count1 = 0;
    int count2 = 0;
	double minIntensity = 100000000;
    for(it = listOfPoints.begin(); it != listOfPoints.end(); it++)
    {
		if(abs(it->at(dim)) < abs(minIntensity))
		{
		 	minIntensity = abs(it->at(dim));
		}
       // if((count1 < numOfPoints))
//        {
//            count2 = 0;
//            for (it2 = it->begin(); it2 != it->end(); it2++)
//            {
//                if(count2 <= dim)
//                {
//                    if(count1 < numOfPoints)
//                    {
//                        data[count1][count2] = *it2;
//                        count2++;
//                    }
//                }
//            }
//        }
//        count1++;
    }
	//cout << "Here " << minIntensity << "\n";
	 
	outputDataset(listOfPoints, numOfPoints2, dim, "CurrentData.txt");
	 
	vector < vector <double> > maxs;
	vector < vector <double> > widths;
	maxs.push_back(findMax(listOfPoints2, numOfPoints2, dim));
	//cout << "HERE4\n";	
	while((maxs[maxs.size()-1][0] > -10000000.0) && (maxs[maxs.size()-1][1] > -10000000.0))
	{
		if(dim == 3)
		{
			if(maxs[maxs.size()-1][2] > -10000000.0)
			{
				widths.push_back(findApproxLinewidth(listOfPoints, numOfPoints, dim, maxs[maxs.size()-1]));
				maxs.push_back(findMaxOther(listOfPoints2, numOfPoints2, dim, maxs, widths));
			}
		}
		else
		{
			widths.push_back(findApproxLinewidth(listOfPoints, numOfPoints, dim, maxs[maxs.size()-1]));
			maxs.push_back(findMaxOther(listOfPoints2, numOfPoints2, dim, maxs, widths));
		}
	}
	//cout << "HERE6\n";
	
	vector < vector <double> > widthsOutput;
	vector <double> temp;
	for(int y = 0; y <  dim*4 + 2; y++)
	{
		temp.push_back(0.0);
	}
	for(int y = 0; y <  widths.size(); y++)
	{
		temp[0] = maxs[y][0] - (widths[y][0]/2);
		temp[1] = maxs[y][1];
		temp[2] = maxs[y][0] + (widths[y][0]/2);
		temp[3] = maxs[y][1];
		temp[4] = widths[y][2];
		temp[5] = maxs[y][0];
		temp[6] = maxs[y][1] - (widths[y][1]/2);
		temp[7] = maxs[y][0];
		temp[8] = maxs[y][1] + (widths[y][1]/2);
		temp[9] = widths[y][2];
		widthsOutput.push_back(temp);
	}
	
	//outputVector(maxs, "maxs.txt", int(widths.size()));
	//outputVector2(widthsOutput, "widths.txt", int(widths.size()), dim);
	//outputVector3(widthsOutput, "widths2.txt", int(widths.size()), dim);
	
	int startingK =  int(widths.size());
	//cout << "startingK: " << startingK << "\n";
	
	//Put k-models here using the starting K, and then remove those models from the dataset and repeat.
	std::vector <gene> results;
	
	//Maybe do something else: Partition, fit give only a fraction of each point to each model, refit, repeat until convergence
	results = kmodels(listOfPoints, numOfPoints, dim, startingK, maxs, widths, cutoff, minIntensity);
	string outfile = "genes" + to_string(file) + ".txt";
	//cout << outfile << endl;
	printGenestoFile(results, outfile, 0.2);
	
    // for(int x = 0; x < numOfPoints2; x++)
//     {
//         free(data[x]);
//     }
//     free(data);
	return(0);
}

//set up starting partition based on distance to peak center
std::vector <std::vector <std::vector <double> > > startingPartition(vector < vector <double> > & data, int numOfPoints, int dim, int k, vector < vector <double> > maxs, vector < vector <double> > widths)
{
	std::vector <std::vector <std::vector <double> > > results;
	std::vector <std::vector <double> > temp;
	std::vector <double> temp1;
	//cout << "In startingPartition\n";
	
	for(int x = 0; x < dim+1; x++)
	{
		temp1.push_back(0.0);
	}
	double sizeIncrease = 1.0;
	for(int x = 0; x < k; x++)
	{
		results.push_back(temp);
	}
	vector <double> distances;
	for(int x = 0; x < k; x++)
	{
		distances.push_back(0.0);
	}
	
	//cout << "In startingPartition Initialization complete\n";
	for(int currentPoint = 0; currentPoint < numOfPoints-1; currentPoint++)
	{
		//cout << currentPoint << "\n";
		for(int x = 0; x < k; x++)
		{
			distances[x] = 0.0;
		}
		//cout << "Got to here0\n";
		for(int x = 0; x < k; x++)
		{
			//cout << x << "\n";
			if(dim == 2)
			{
				distances[x] = ((data[currentPoint][0] - maxs[x][0])/(sizeIncrease*widths[x][0]))*((data[currentPoint][0] - maxs[x][0])/(sizeIncrease*widths[x][0])) + ((data[currentPoint][1] - maxs[x][1])/(sizeIncrease*widths[x][1]))*((data[currentPoint][1] - maxs[x][1])/(sizeIncrease*widths[x][1]));
			}
			if(dim == 3)
			{
				//cout << "havent done that yet\n";
				distances[x] = ((data[currentPoint][0] - maxs[x][0])/(sizeIncrease*widths[x][0]))*((data[currentPoint][0] - maxs[x][0])/(sizeIncrease*widths[x][0])) + ((data[currentPoint][1] - maxs[x][1])/(sizeIncrease*widths[x][1]))*((data[currentPoint][1] - maxs[x][1])/(sizeIncrease*widths[x][1])) + ((data[currentPoint][2] - maxs[x][2])/(sizeIncrease*widths[x][2]))*((data[currentPoint][2] - maxs[x][2])/(sizeIncrease*widths[x][2]));
			}
		}
		//cout << "Got to here1\n";
		std::vector<double>::iterator minDist = min_element(distances.begin(), distances.end());
		int ind =  (int)std::distance(distances.begin(), minDist);
		for(int x = 0; x <dim+1;x++)
		{
			temp1[x] = data[currentPoint][x];
		}
		results[ind].push_back(temp1);
	}
	return(results);
}

std::vector <std::vector <std::vector <double> > > repartition(vector < vector <double> > & data, std::vector <gene> genes, int dim, int numOfPoints, int k, double Thresh)
{
	std::vector <std::vector <std::vector <double> > > results;
	std::vector <std::vector <double> > temp;
	std::vector <double> sections;
	std::vector <double> temp1;
	double height = 0;
	
	//setup data structure peices 
	for(int x = 0; x < dim+1; x++)
	{
		temp1.push_back(0.0);
	}
	for(int x = 0; x < k; x++)
	{
		results.push_back(temp);
	}
	
	//For each point
	//int x = 0;
	for(int x = 0; x < numOfPoints; x++)
	{
		while(sections.size()>0)
		{
			sections.pop_back();
		}
		height = data[x][dim];
		//give only part of the point to each cluster, depending on proportion of height
		for(int w = 0; w <dim; w++)
		{
			temp1[w] = data[x][w];
		}
		double val1 = 0;
		//find hieght of each model at the point
		for(int y = 0; y < k; y++)
		{
			val1 = findSinglePointInModel(dim, genes[y], numOfPoints, data, x);
			sections.push_back(val1);
		}
	
		//find the sum of all the positive and negative models
		int N1 = 0;
		int N2 = 0;
		int sum = 0;
		int sum_neg = 0;
		for(int w = 0; w <sections.size();w++)
		{
			if(sections[w] >= 0)
			{
				sum = sum + sections[w];
				N1 = N1 + 1;
			}
			else
			{
				sum_neg = sum_neg + sections[w];
				N2 = N2 + 1;
			}
		}
		//Determine the proportion of the model that each point will be divided up by. Then use that to set up the different datasets.
		//only add if above threshold
		//yi = hieght of model at position i
		//z = height of data
		//n = number of models in that partition
		double p_neg = 0;
		double p_pos = 0;
		for(int w = 0; w < k; w++)
		{
			if(height < sum_neg)
			{
				p_neg = height;
				p_pos = 0.0;
			}
			if((height >= sum_neg - sum) && (height <= sum - sum_neg))
			{
				//p_neg = sum_neg + ((height - (sum_neg + sum))/N2);
				//p_pos = sum + ((height - (sum_neg + sum))/N1);
				p_neg = sum_neg + ((height - (sum_neg + sum))/2);
				p_pos = sum + ((height - (sum_neg + sum))/2);
			}
			if(height > sum - sum_neg)
			{
				p_pos = height;
				p_neg = 0.0;
			}
			if(sections[w] <= 0.0)
			{
				temp1[dim] = p_neg*(sections[w]/sum_neg);
			}
			if(sections[w] > 0.0)
			{
				temp1[dim] = p_pos*(sections[w]/sum);
			}
			if(abs(temp1[dim]) > abs(Thresh))
			{
				results[w].push_back(temp1);
			}
		}
	}
	//return new clusters
	return(results);
}

//If the genes are right on top of one another, add the heights, keeps A's linewidth, center location and alpha, and return the new gene. Meh just keep a's height.
gene combineGenes(gene a, gene b, std::vector <std::vector <std::vector <double> > > results, int aNum, int bNum, double cut)
{
	//cout << "in combine\n";
	gene c;
	gene_inital(c, a.dim);
	int flag1 = 0;
	flag1 = 0;
	for(int q = 0; q < a.dim; q++)
	{
		if(abs(b.center[q] - a.center[q]) < cut*a.width[q])
		{
			flag1++;
		}
	}
	if(flag1 == a.dim)
	{
		c.height = a.height;
		c.alpha = a.alpha;
		for(int x = 0; x < a.dim; x++)
		{
			c.center[x] = a.center[x];
			c.width[x] = a.width[x];
		}
	}
	//std::vector <std::vector <double> > results2;
	//std::vector <double> temp;
	//for(int x = 0; x < c.dim+1; x++)
	//{
		//temp.push_back(0.0);
	//}
	//cout << "got to here\n";
	//if(results.size() > 0)
	//{
		//for(int x = 0; x < results[aNum].size(); x++)
		//{
			//results2.push_back(temp);
		//}
		//for(int x = 0; x < results[bNum].size(); x++)
		//{
			//results2.push_back(temp);
		//}
	//}
	c.fitness = a.fitness;
	//cout << "finished combine\n";
	return(c);
}


int removeRunaways(std::vector <gene> & genes, vector < vector <double> > & max, vector < vector <double> > & width, vector <double> & mins, vector <double> & Uplim, int & k, int dim)
{
	//int widthFlag = 0;
	int goAgainFlag = 1;
	int flag3 = 0;
	int worstInd = -1;
	//cout << "Removing runaways\n";
	while(goAgainFlag == 1)
	{
		goAgainFlag = 0;
		flag3 = 0;
		for(int x = (k-1); x >= 0; x--)
		{
			for(int y = 0; y < dim; y++)
			{
				if(genes[x].center[y] < mins[y])
				{
					flag3 = 1;
					worstInd = x;
				}
				if(genes[x].center[y] > Uplim[y])
				{
					flag3 = 1;
					worstInd = x;
				}
			}
			if((flag3 == 1) & (k > 0) & (worstInd > -1))
			{
				//cout << mins[0] << " " << mins[1] << " " << Uplim[0] << " " << Uplim[1] << "\n";
				//printGene(genes[worstInd]);
				removeElement(genes, worstInd);
				removeElement(max, worstInd);
				removeElement(width, worstInd);
				k = k - 1;
				//cout << "Redoing with k-1 = " << k << "\n";
				goAgainFlag = 1;
				flag3 = 0;
			}
		}
		// if((flag3 == 1) & (k > 0))
		// {
		// 	cout << mins[0] << " " << mins[1] << " " << Uplim[0] << " " << Uplim[1] << "\n";
		// 	printGene(genes[worstInd]);
		// 	removeElement(genes, worstInd);
		// 	removeElement(max, worstInd);
		// 	removeElement(width, worstInd);
		// 	k = k - 1;
		// 	cout << "Redoing with k-1 = " << k << "\n";
		// 	goAgainFlag = 1;
		// }
	}
	//cout << "Completed removeRunaways\n";
	return(flag3);
}

int checkLineWidths(std::vector <gene> & genes, vector < vector <double> > & max, vector < vector <double> > & width, int & k, int dim)
{
	int flag2 = 0;	
	int widthFlag = 0;
	//int goAgainFlag = 1;
	int worstInd = -1;
	//cout << "Removing too thin and too thick\n";
	// while(goAgainFlag == 1)
	// {
		//goAgainFlag = 0;
		for(int x = (k-1); x >= 0; x--)
		{
			widthFlag = 0;
			for(int y = 0; y < dim; y++)
			{
				//cout << x << "\n";
				if((genes[x].width[y] < 0.0001) | (genes[x].width[y] > 10))
				{
					flag2 = 1;
					widthFlag = 1;
					worstInd = x;
					//printGene(genes[x]);
				}
			}
			if(widthFlag == 1)
			{
				if(k > 0)
				{
					removeElement(genes, worstInd);
					removeElement(max, worstInd);
					removeElement(width, worstInd);
					k = k - 1;
					//goAgainFlag = 1;
				}
			}
		}
		// if(widthFlag == 1)
		// {
		// 	if(k > 0)
		// 	{
		// 		removeElement(genes, worstInd);
		// 		removeElement(max, worstInd);
		// 		removeElement(width, worstInd);
		// 		k = k - 1;
		// 		goAgainFlag = 1;
		// 	}
		// }
		//}
	//cout << "Completed checkLineWidths\n";
	return(flag2);
}

int runCombine(int & k, std::vector <gene> & genes, int dim, vector < vector <double> > & max, vector < vector <double> > & width, std::vector <std::vector <std::vector <double> > > & results, double cut)
{
	gene temp;
	int flag1 = 0;
	//double cut = 0.5;
	//double cut = 1.0;
	//maybe try 0.75 next? or 0.8?
	if(k > 1)
	{
		//cout << "Here6\n";
		for(int x = 0; x < k; x++)
		{
			//cout << x << " Here7\n";
			for(int y = (k-1); y >= 0; y--)
			{
				//cout << y << " Here8\n";
				if(y != x)
				{
					//cout << y << " Here9\n";
					flag1 = 0;
					if(genes.size() > 0)
					{
						//cout << y << " Here10\n";
						for(int q = 0; q < dim; q++)
						{
							//cout << y << " Here10.1\n";
							//cout << x << "\n";
							if(y < k && x < k)
							{
								//cout << "good\n";
								if(abs(genes[y].center[q] - genes[x].center[q]) < cut*(genes[x].width[q]))
								{
									//cout << y << " Here10.2\n";
									//printGene(genes[y]);
									//printGene(genes[x]);
									//flag1++;
									if(genes[x].height > genes[y].height)
									{
										flag1++;
									}
									// if(genes[x].width[q] < genes[y].width[q]*2)
		// 							{
		// 								flag1++;
		// 							}
									//cout << y << " Here10.3\n";
								}
							}
							/*if(abs(genes[y].center[q] - genes[x].center[q]) < 0.5*(genes[x].width[q]))
							{
								cout << y << " Here10.2\n";
								//printGene(genes[y]);
								//printGene(genes[x]);
								//flag1++;
								if(genes[x].height > genes[y].height)
								{
									flag1++;
								}
								// if(genes[x].width[q] < genes[y].width[q]*2)
	// 							{
	// 								flag1++;
	// 							}
								cout << y << " Here10.3\n";
							}*/
							//cout << y << " Here10.4\n";
						}
						//cout << y << " Here11\n";
						if(flag1 == dim)
						{
							temp = combineGenes(genes[x], genes[y], results, x, y, cut);
							//cout << "Got to here1\n";
							//if the returned gene has a fitness below 10^18, then the combine was a success.
							if(temp.fitness < 1e+18)
							{
								genes[x] = temp;
								//cout << "Got to here2.1\n";
								removeElement(genes, y);
								//cout << "Got to here2.2\n";
								removeElement(max, y);
								//cout << "Got to here2.3\n";
								removeElement(width, y);
								k = k-1;
								//flag1 = 1;
							}
							//cout << "Got to here2\n";
						
						}
						//cout << y << " Here12\n";
					}
					if(genes.size() == 0)
					{
						for(int q = 0; q < dim; q++)
						{
							if(abs(max[y][q] - max[x][q]) < cut*(width[x][q]))
							{
								if(max[x][dim] > max[y][dim])
								{
									flag1++;
								}
							}
						}
						if(flag1 == dim)
						{
							removeElement(max, y);
							removeElement(width, y);
							k = k-1;
						}
					}
					//cout << "Got to here3.1\n";
				}
				//cout << "Got to here3.2\n";
			}
			//cout << "Got to here3.3\n";
		}
		//cout << "Got to here3.4\n";
	}
	//cout << "Got to here3\n";
	return(k);
}


void Multithread_DE(std::vector <std::vector <double> > cluster, int numOfPoints, int dim, std::vector <double> & max, std::vector <double> & width, gene & best)
{
	best = DE(cluster, numOfPoints, dim, max, width);
	for(int y = 0; y < dim; y++)
	{
		max[y] = best.center[y];
		width[y] = best.width[y];
	}
	max[dim] = best.height;
}


/**
What you need:
set of S points in Rn
a value for k as the number of clusters
choice of statistical model
Output: category function
Steps:
1) partition data randomly into k clusters 
repeat:
	2) fit a model to each cluster
	3) assign points to clusters based on RMSE to closest model
until: no points change clusters
Report final cluster assignment and statistical model

for us:
cut high (1/5 of max maybe), and initialize to the number of likely cluster model centers.
Once these models have been found subtract the model from the data, and restart.
*/
//Maybe: Partition, fit give only a fraction of each point to each model, refit, repeat until convergence
std::vector <gene> kmodels(vector <std::vector<double> > & data, int numOfPoints, int dim, int k, vector < vector <double> > maxs, vector < vector <double> > widths, double cutoff, double Thresh)
{
	cout << "Inside K-models\n";
	vector <double> mins;
	vector <double> Uplim;
	int max_Proc = std::thread::hardware_concurrency();
	//int max_Proc = 1;
	cout << max_Proc << "\n";
	
	//double temp1 = 0;
	//double temp2 = 1000000000;
	for(int x = 0; x < dim; x++)
	{
		mins.push_back(10000000.0);
		Uplim.push_back(-10000000.0);
	}
	for(int x = 0; x < numOfPoints; x++)
	{
		for(int y = 0; y < dim; y++)
		{
			if(data[x][y] < mins[y])
			{
				mins[y] = data[x][y];
			}
			if(data[x][y] > Uplim[y])
			{
				Uplim[y] = data[x][y];
			}
		}
	}
	for(int x = 0; x < dim; x++)
	{
		cout << mins[x] << " "; 
	}
	cout << "\n";
	for(int x = 0; x < dim; x++)
	{
		cout << Uplim[x] << " "; 
	}
	cout << "\n";
	
	std::vector <std::vector <std::vector <double> > > results;
	std::vector <gene> genes;
	gene a;
	
	int flag = 1;
	double worst = 0.0;
	int worstInd = -1;
	vector < vector <double> > max;
	vector < vector <double> > width;
	vector <double> max2;
	vector <double> width2;
	for(int y = 0; y < dim+1; y++)
	{
		max2.push_back(0.0);
	}
	for(int y = 0; y < dim; y++)
	{
		width2.push_back(0.0);
	}
	for(int x = 0; x < k; x++)
	{
		for(int y = 0; y < dim+1; y++)
		{
			max2[y] = maxs[x][y];
		}
		max.push_back(max2);
		for(int y = 0; y < dim; y++)
		{
			width2[y] = widths[x][y];
		}
		width.push_back(width2);
	}
	//cout << "Here\n";
	
	int count = 0;
	
	//partition
	results = startingPartition(data, numOfPoints, dim, k, max, width);
	//outputDatasetClustersV1(data, numOfPoints, dim, results);
	//cout << "Here2\n";
	//outputVector(max, "maxs1.txt", int(widths.size()));
	
	//if any partitions have too few points, remove the shortest model, and do starting partition again
	
	//int flag_tooSmall = 0;
	int flag_goAgain = 1;
	int toRemove = -1;
	int smallest = 18;
	int nextSmallest = 18;
	double smallest_H = 1000000000;
	while(flag_goAgain == 1)
	{
		smallest = 18;
		nextSmallest = 18;
		flag_goAgain = 0;
		toRemove = -1;
		for(int x = 0; x < k; x++)
		{
			if((results[x].size() < smallest) & (k > 1))
			{
				smallest = (int)results[x].size();
				toRemove = x;
				flag_goAgain = 1;
			}
		}
		if(flag_goAgain == 1)
		{
			//double smallest_H = 1000000000;
			smallest_H = 1000000000;
			flag_goAgain = 0;
			for(int x = 0; x < k; x++)
			{
				if((abs(max[x][dim]) < abs(smallest_H)) & (k > 1))
				{
					if(x != 0)
					{
						smallest_H = max[x][dim];
						toRemove = x;
						flag_goAgain = 1;
					}
				}
			}
			//assuming you find somethiing to remove that is not the very first model
			if(toRemove > 0)
			{
				removeElement(max, toRemove);
				removeElement(width, toRemove);
				k = k-1;
				cout << "k-1 = " << k << "\n";
				results = startingPartition(data, numOfPoints, dim, k, max, width);
			}
		}
	}
	
	//outputDatasetClustersV1(data, numOfPoints, dim, results);
	//outputVector(max, "maxs.txt", int(max.size()));
	
	for(int x = 0; x < k; x++)
	{
		for(int y = 0; y < dim+1; y++)
		{
			cout << max[x][y] << " ";
		}
		cout << "\n";
	}
	for(int x = 0; x < k; x++)
	{
		for(int y = 0; y < dim; y++)
		{
			cout << width[x][y] << " ";
		}
		cout << "\n";
	}
	count++;

	for(int x = 0; x < k; x++)
	{
		genes.push_back(a);
	}
	
	cout << "Here4\n";
	//fit
	cout << k << "\n";
	if(k == 1)
	{
		genes[0] = DE(results[0], int(results[0].size()), dim, max[0], width[0]);
	}
	cout << k << "\n";
	if(k > 1 & k <= max_Proc)
	{
		cout << "In here!\n";
		vector<std::thread> threads(k);
		for(int x = 0; x < k; x++)
		{
			//threads.push_back(std::thread(&Multithread_DE, results[x], int(results[x].size()), dim, std::ref(max[x]), std::ref(width[x]), std::ref(genes[x])));
			threads[x] = std::thread(&Multithread_DE, results[x], int(results[x].size()), dim, std::ref(max[x]), std::ref(width[x]), std::ref(genes[x]));
		}
		for(int x = 0; x < k; x++)
		{
			threads[x].join();
		}
	}
	if(k > 1 & k > max_Proc)
	{
		int count = 0;
		cout << "In here2!\n";
		vector<std::thread> threads(k);
		while(count < k)
		{	
			cout << count << "\n";
			for(int x = count; x < max_Proc + count; x++)
			{
				//threads.push_back(std::thread(&Multithread_DE, results[x], int(results[x].size()), dim, std::ref(max[x]), std::ref(width[x]), std::ref(genes[x])));
				if(x < k)
				{
					threads[x] = std::thread(&Multithread_DE, results[x], int(results[x].size()), dim, std::ref(max[x]), std::ref(width[x]), std::ref(genes[x]));
				}
			}
			for(int x = count; x < max_Proc+ count; x++)
			{
				if(x<k)
				{
					threads[x].join();
				}
			}
			count = count + max_Proc;
		}
	}
	cout << k << "\n";
	for(int x = 0; x < k; x++)
	{
		printGene(genes[x]);
	}
	/*for(int x = 0; x < k; x++)
	{
		genes[x] = DE(results[x], int(results[x].size()), dim, max[x], width[x]);
		printGene(genes[x]);
		for(int y = 0; y < dim; y++)
		{
			max[x][y] = genes[x].center[y];
			width[x][y] = genes[x].width[y];
		}
		max[x][dim] = genes[x].height;
	}*/
	
	//cout << "Here5\n";
	//outputVector(max, "maxs1.txt", int(max.size()));
	

	int flag3 = 0;
	// removeRunaways(genes, max, width, mins, Uplim, k, dim);
// 	cout << k << "\n";
	
	//flag = 0;s
	//cout << "                   HERE21\n";
	int flag2 = 0;	
	flag2 = checkLineWidths(genes, max, width, k, dim);
	flag = 1;
	//outputVector(max, "maxs3.txt", int(max.size()));
	//cout << "                   HERE22\n";
	
	//combine
	k = runCombine(k, genes, dim, max, width, results, 1.0);
	//outputVector(max, "maxs2.txt", int(max.size()));
	//cout << "Here5.5\n";
	gene temp;
	int flag1 = 0;
	//cout << "Here6\n";
	
	// for(int x = 0; x < k; x++)
// 	{
// 		for(int y = 0; y < dim; y++)
// 		{
// 			max[x][y] = genes[x].center[y];
// 			width[x][y] = genes[x].width[y];
// 		}
// 		max[x][dim] = genes[x].height;
// 	}

	//formerlly, the initial parition is done for every k decrease. Might want to reduce that. Remove the worst model, or the model with erroneous parameters, and then keep the fittings from before rather than redo them.
	//what if we removed all erroneous models at once, rather than one at a time? With the double repartition, it takes too long for cluster 16 (k=12 to start, ends with k=11?)
	flag2 = 0;	
	while((flag == 1) & (k > 0))
	{
		cout << "Repeat " << count << "\n";
		cout << "k: " << k << "\n";
		count++;
		
		//repartition
		results = repartition(data, genes, dim, numOfPoints, k, Thresh);
		if(k == 1)
		{
			genes[0] = DE(results[0], int(results[0].size()), dim, max[0], width[0]);
		}
		cout << k << "\n";
		if(k > 1 & k <= max_Proc)
		{
			cout << "In here!\n";
			vector<std::thread> threads(k);
			for(int x = 0; x < k; x++)
			{
				//threads.push_back(std::thread(&Multithread_DE, results[x], int(results[x].size()), dim, std::ref(max[x]), std::ref(width[x]), std::ref(genes[x])));
				threads[x] = std::thread(&Multithread_DE, results[x], int(results[x].size()), dim, std::ref(max[x]), std::ref(width[x]), std::ref(genes[x]));
			}
			for(int x = 0; x < k; x++)
			{
				threads[x].join();
			}
		}
		if(k > 1 & k > max_Proc)
		{
			int count = 0;
			cout << "In here2!\n";
			vector<std::thread> threads(k);
			while(count < k)
			{	
				cout << count << "\n";
				for(int x = count; x < max_Proc + count; x++)
				{
					//threads.push_back(std::thread(&Multithread_DE, results[x], int(results[x].size()), dim, std::ref(max[x]), std::ref(width[x]), std::ref(genes[x])));
					if(x < k)
					{
						threads[x] = std::thread(&Multithread_DE, results[x], int(results[x].size()), dim, std::ref(max[x]), std::ref(width[x]), std::ref(genes[x]));
					}
				}
				for(int x = count; x < max_Proc+ count; x++)
				{
					if(x<k)
					{
						threads[x].join();
					}
				}
				count = count + max_Proc;
			}
		}
		/*for(int x = 0; x < k; x++)
		{
			//genes[x] = DE(results[x], int(results[x].size()), dim, max[x], width[x]);
			cout << "Working on it\n";
			std::thread t1(&Multithread_DE, results[x], int(results[x].size()), dim, std::ref(max[x]), std::ref(width[x]), std::ref(genes[x]));
			cout << "got through it!\n";
			t1.join();
			cout << "done!\n";
			printGene(genes[x]);
			/*for(int y = 0; y < dim; y++)
			{
				max[x][y] = genes[x].center[y];
				width[x][y] = genes[x].width[y];
			}
			max[x][dim] = genes[x].height;
		}*/
		
		if(k > 1)
		{
			//repartition
			results = repartition(data, genes, dim, numOfPoints, k, Thresh);
			//outputDatasetClustersV2(data, numOfPoints, dim, results);
			if(k == 1)
			{
				genes[0] = DE(results[0], int(results[0].size()), dim, max[0], width[0]);
			}
			cout << k << "\n";
			if(k > 1 & k <= max_Proc)
			{
				cout << "In here!\n";
				vector<std::thread> threads(k);
				for(int x = 0; x < k; x++)
				{
					//threads.push_back(std::thread(&Multithread_DE, results[x], int(results[x].size()), dim, std::ref(max[x]), std::ref(width[x]), std::ref(genes[x])));
					threads[x] = std::thread(&Multithread_DE, results[x], int(results[x].size()), dim, std::ref(max[x]), std::ref(width[x]), std::ref(genes[x]));
				}
				for(int x = 0; x < k; x++)
				{
					threads[x].join();
				}
			}
			if(k > 1 & k > max_Proc)
			{
				int count = 0;
				cout << "In here2!\n";
				vector<std::thread> threads(k);
				while(count < k)
				{	
					cout << count << "\n";
					for(int x = count; x < max_Proc + count; x++)
					{
						//threads.push_back(std::thread(&Multithread_DE, results[x], int(results[x].size()), dim, std::ref(max[x]), std::ref(width[x]), std::ref(genes[x])));
						if(x < k)
						{
							threads[x] = std::thread(&Multithread_DE, results[x], int(results[x].size()), dim, std::ref(max[x]), std::ref(width[x]), std::ref(genes[x]));
						}
					}
					for(int x = count; x < max_Proc+ count; x++)
					{
						if(x<k)
						{
							threads[x].join();
						}
					}
					count = count + max_Proc;
				}
			}
			//Refit
			/*for(int x = 0; x < k; x++)
			{
				genes[x] = DE(results[x], int(results[x].size()), dim, max[x], width[x]);
				printGene(genes[x]);
				for(int y = 0; y < dim; y++)
				{
					max[x][y] = genes[x].center[y];
					width[x][y] = genes[x].width[y];
				}
				max[x][dim] = genes[x].height;
			}*/
		}

		
		cout << "Checking linewidths\n";
		//if any have weird linewidths, remove them
		flag2 = checkLineWidths(genes, max, width, k, dim);
		flag = 0;
		worst = 0.0;
		worstInd = -1;
		//flag2 = 0;
		//widthFlag = 0;
		//goAgainFlag = 1;
		//if a gene is not within 4 linewidths of another hill, remove it. this may not work
		//maybe within range of data? 
		cout << "Checking Range\n";
		flag3 = removeRunaways(genes, max, width, mins, Uplim, k, dim);
		
		cout << "Checking fitness\n";
		//if any fit worse than 2.5% (cutoff) normalized error, remove the worst fitting model, and try again with one less K.
		for(int x = 0; x < genes.size(); x++)
		{
			if(genes[x].fitness > cutoff)
			{
				flag = 1;
				if(genes[x].fitness > worst)
				{
					worst = genes[x].fitness;
					worstInd = x;
				}
			}
		}
		//removes worst hill
		if(flag == 1 & k > 0)
		{
			printGene(genes[worstInd]);
			removeElement(genes, worstInd);
			removeElement(max, worstInd);
			removeElement(width, worstInd);
			k = k - 1;
			cout << "Redoing with k-1 = " << k << "\n";
		}
		
		flag1 = runCombine(k, genes, dim, max, width, results, 1.0);

		if(flag1 < k)
		{
			flag = 1;
		         k = flag1;
		}
		if(flag2 == 1)
		{
			flag = 1;
		}
		if(flag3 == 1)
		{
			flag = 1;
		}
	}
	if(k == 0)
	{
		cout << "Finishing with 1 hill\n";
		k = 1;
		for(int x = 0; x < k; x++)
		{
			genes.push_back(a);
		}
		for(int y = 0; y < dim+1; y++)
		{
			max2.push_back(0.0);
		}
		for(int y = 0; y < dim; y++)
		{
			width2.push_back(0.0);
		}
		for(int x = 0; x < k; x++)
		{
			for(int y = 0; y < dim+1; y++)
			{
				max2[y] = maxs[x][y];
			}
			max.push_back(max2);
			for(int y = 0; y < dim; y++)
			{
				width2[y] = widths[x][y];
			}
			width.push_back(width2);
		}
		//results = repartition(data, genes, dim, numOfPoints, k, Thresh);
		results = startingPartition(data, numOfPoints, dim, k, maxs, widths);
		//outputDatasetClustersV1(data, numOfPoints, dim, results);
		//fit
		for(int x = 0; x < k; x++)
		{
			genes[x] = DE(results[x], int(results[x].size()), dim, maxs[x], widths[x]);
			printGene(genes[x]);
			for(int y = 0; y < dim; y++)
			{
				max[x][y] = genes[x].center[y];
				width[x][y] = genes[x].width[y];
			}
			max[x][dim] = genes[x].height;
		}
		
		//add a bit here that removes the last peak if it is bad.
		//flag2 = 0;
		flag2 = checkLineWidths(genes, max, width, k, dim);

		cout << "Checking Range\n";
		//flag3 = 0;
		flag3 = removeRunaways(genes, max, width, mins, Uplim, k, dim);
		cout << k << "\n";
		cout << genes.size() << "\n";
	}
	//runCombine(k, genes, dim, max, width, results, 1.0);
	//outputVector(max, "maxs_final.txt", int(max.size()));
	cout << "All done\n";
	//report
	cout << k << "\n";
	cout << genes.size() << "\n";
	//outputDatasetClustersV2(data, numOfPoints, dim, results);
	return(genes);
}

vector <double> findMaxOther(vector <std::vector<double> > & listOfPoints, int numOfPoints, int dim, vector < vector <double> > currentMax, vector < vector <double> > widths)
{
	double MaxH = 0;
	vector <double> max2;
	for(int x = 0; x < dim+1; x++)
	{
		max2.push_back(-10000000.0);
	}
	
	//double w = 0;
	
	//w = (*min_element(widths.begin(), widths.end()))/2;
	//cout << w << "\n";
	double elipse = 0;
	int counter = 0;
	double sizeIncrease = 1.0;
	
	vector <std::vector<double> > :: iterator it;
	for(it = listOfPoints.begin(); it != listOfPoints.end(); it++)
	{
		//cout << currentMax[0] << " " << currentMax[1] << " " << currentMax[2] << "\n";
		//cout << dist((*it), currentMax) << "\n";
		elipse = 0;
		counter = 0;
		for(int x = 0; x < currentMax.size(); x++)
		{
			if(dim == 2)
			{
				elipse = (((*it)[0] - currentMax[x][0])/(sizeIncrease*widths[x][0]))*(((*it)[0] - currentMax[x][0])/(sizeIncrease*widths[x][0])) + (((*it)[1] - currentMax[x][1])/(sizeIncrease*widths[x][1]))*(((*it)[1] - currentMax[x][1])/(sizeIncrease*widths[x][1]));
			}
			if(dim == 3)
			{
				elipse = (((*it)[0] - currentMax[x][0])/(sizeIncrease*widths[x][0]))*(((*it)[0] - currentMax[x][0])/(sizeIncrease*widths[x][0])) + (((*it)[1] - currentMax[x][1])/(sizeIncrease*widths[x][1]))*(((*it)[1] - currentMax[x][1])/(sizeIncrease*widths[x][1])) + (((*it)[2] - currentMax[x][2])/(sizeIncrease*widths[x][2]))*(((*it)[2] - currentMax[x][2])/(sizeIncrease*widths[x][2]));
			}
			if(dim == 4)
			{
				elipse = (((*it)[0] - currentMax[x][0])/(sizeIncrease*widths[x][0]))*(((*it)[0] - currentMax[x][0])/(sizeIncrease*widths[x][0])) + (((*it)[1] - currentMax[x][1])/(sizeIncrease*widths[x][1]))*(((*it)[1] - currentMax[x][1])/(sizeIncrease*widths[x][1])) + (((*it)[2] - currentMax[x][2])/(sizeIncrease*widths[x][2]))*(((*it)[2] - currentMax[x][2])/(sizeIncrease*widths[x][2])) + (((*it)[3] - currentMax[x][3])/(sizeIncrease*widths[x][3]))*(((*it)[3] - currentMax[x][3])/(sizeIncrease*widths[x][3]));
			}
			if(elipse > 1)
			{
				counter++;
			}
		}
		//cout << counter << "\n";
		if(counter == currentMax.size())
		{
			if(dim == 2)
			{
				if(abs((*it).at(2)) > abs(MaxH))
				{
					(MaxH) = (*it).at(2);
					max2[0] = (*it).at(0);
					max2[1] = (*it).at(1);
					max2[2] = (*it).at(2);
				}
			}
			if(dim == 3)
			{
				if(abs((*it).at(3)) > abs(MaxH))
				{
					(MaxH) = (*it).at(3);
					max2[0] = (*it).at(0);
					max2[1] = (*it).at(1);
					max2[2] = (*it).at(2);
					max2[3] = (*it).at(3);
				}
			}
	 		if(dim == 4)
	 		{
	 			if(abs((*it).at(4)) > abs(MaxH))
	 			{
					(MaxH) = (*it).at(4);
					max2[0] = (*it).at(0);
					max2[1] = (*it).at(1);
					max2[2] = (*it).at(2);
					max2[3] = (*it).at(3);
					max2[4] = (*it).at(4);
	 			}
	 		}
		}
	}
	return(max2);
}


// void ReadData(string input1, list <std::vector<double> > * listOfPoints)
// {//read in the data
//
// 	fstream input;
// 	char buf1[256];
// 	string buf2;
// 	//list <string> lines;
//     string line;
//
// 	string inputFile = input1;
//
// 	//cout << inputFile << "\n";
//     std::vector<double> out;
// 	input.open(inputFile,ios::in);
// 	input.getline(buf1, 255); //Get first line
//     line = buf1;
// 	//cout << line;
//     tokenize(line, ' ', &out);
//     std::vector<double> buffer = std::vector<double> (out);
//     out.erase(out.begin(), out.end());
//     listOfPoints->push_back(buffer);
//     //char flag = 0;
//
// 	while((strlen(buf1)>0))//while we have data
// 	{
//         input.getline(buf1, 255);
//         if(input.eof())
//         {
//             input.close();
//             return;
//         }
//         line = buf1;
//         tokenize(line, ' ', &out);
//         std::vector<double> buffer = std::vector<double> (out);
//         out.erase(out.begin(), out.end());
//         listOfPoints->push_back(buffer);
// 		//lines.push_back(buf1);
// 	}
// 	input.close();
// }
//
// void tokenize(std::string const str, const char delim, std::vector<double> * out)
// {
// 	std::stringstream ss(str);
// 	std::string s;
// 	std::vector<string> temp;
// 	while(std::getline(ss, s, delim))
// 	{
// 		temp.push_back(s);
// 	}
//
// 	std::vector<string>:: iterator it;
//     for(it = temp.begin(); it != temp.end(); ++it)
// 	{
// 		out->push_back(stod(*it));
// 	}
// }

void outputVector2(vector < vector <double> > & data, string filename, int size, int dim)
{
	//cout << "In outputVector\n";
    fstream aus;  //output file
    aus.open(filename,ios::out);
	aus.precision(7);
	string outstring;
	for(int x = 0; x < size; x++)
	{
		outstring = "";
		for(int y = 0; y < dim; y++)
		{
			outstring = outstring + to_string(data[x][y]) + " ";
		}
		outstring = outstring + to_string(data[x][data[x].size()-1]);
		outstring = outstring + "\n";
		//outstring = outstring + "\n";
		//cout << outstring;
		for(int y = dim; y < dim*2; y++)
		{
			outstring = outstring + to_string(data[x][y]) + " ";
		}
		outstring = outstring + to_string(data[x][data[x].size()-1]);
		outstring = outstring + "\n";
		outstring = outstring + "\n";
		//outstring = outstring + "\n";
		//cout << outstring;
		for(int y = dim*2+1; y < dim*3+1; y++)
		{
			outstring = outstring + to_string(data[x][y]) + " ";
		}
		outstring = outstring + to_string(data[x][data[x].size()-1]);
		outstring = outstring + "\n";
		//outstring = outstring + "\n";
		//cout << outstring;
		for(int y = dim*3+1; y < dim*4+1; y++)
		{
			outstring = outstring + to_string(data[x][y]) + " ";
		}
		outstring = outstring + to_string(data[x][data[x].size()-1]);
		outstring = outstring + "\n";
		outstring = outstring + "\n";
		//outstring = outstring + "\n";
		//cout << outstring;
		aus << outstring;
	}
	aus.close();
}

void outputVector3(vector < vector <double> > & data, string filename, int size, int dim)
{
	//cout << "In outputVector\n";
    fstream aus;  //output file
    aus.open(filename,ios::out);
	aus.precision(7);
	string outstring;
	for(int x = 0; x < size; x++)
	{
		outstring = "";
		for(int y = 0; y < dim; y++)
		{
			outstring = outstring + to_string(data[x][y]) + " ";
		}
		outstring = outstring + to_string(data[x][data[x].size()-1]);
		outstring = outstring + "\n";
		outstring = outstring + "\n";
		//cout << outstring;
		for(int y = dim; y < dim*2; y++)
		{
			outstring = outstring + to_string(data[x][y]) + " ";
		}
		outstring = outstring + to_string(data[x][data[x].size()-1]);
		outstring = outstring + "\n";
		outstring = outstring + "\n";
		outstring = outstring + "\n";
		//cout << outstring;
		for(int y = dim*2+1; y < dim*3+1; y++)
		{
			outstring = outstring + to_string(data[x][y]) + " ";
		}
		outstring = outstring + to_string(data[x][data[x].size()-1]);
		outstring = outstring + "\n";
		outstring = outstring + "\n";
		//cout << outstring;
		for(int y = dim*3+1; y < dim*4+1; y++)
		{
			outstring = outstring + to_string(data[x][y]) + " ";
		}
		outstring = outstring + to_string(data[x][data[x].size()-1]);
		outstring = outstring + "\n";
		outstring = outstring + "\n";
		outstring = outstring + "\n";
		//cout << outstring;
		aus << outstring;
	}
	aus.close();
}

void outputVector(vector < vector <double> > & data, string filename, int size)
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

void outputDataset(vector <vector <double> > & data, int numOfPoints, int dim, string filename)
{
    fstream aus;  //output file
    aus.open(filename,ios::out);
	aus.precision(7);
	string outstring;
	for(int x = 0; x < numOfPoints; x++)
	{
		outstring = "";
		for(int y = 0; y < dim+1; y++)
		{
			outstring = outstring + to_string(data[x][y]) + " ";
		}
		outstring = outstring + "\n";
		aus << outstring;
	}
	aus.close();
}

void outputDatasetClustersV1(vector < vector <double> > & data, int numOfPoints, int dim, std::vector <std::vector <std::vector <double> > > clusters)
{
    fstream aus;  //output file
	string filename;
	aus.precision(7);
	string outstring;
	for(int x = 0; x < clusters.size(); x++)
	{
		filename = "cluster" + to_string(x) + ".txt";
    	aus.open(filename,ios::out);
		aus.precision(7);
		for(int z = 0; z < clusters[x].size(); z++)
		{
			outstring = "";
			for(int w = 0; w < clusters[x][z].size(); w++)
			{
				outstring = outstring + to_string(clusters[x][z][w]) + " ";
			}
			outstring = outstring + "\n";
			aus << outstring;
		}
		aus.close();
	}
}

void outputDatasetClustersV2(vector < vector <double> > & data, int numOfPoints, int dim, std::vector <std::vector <std::vector <double> > > clusters)
{
    fstream aus;  //output file
	string filename;
	aus.precision(7);
	string outstring;
	for(int x = 0; x < clusters.size(); x++)
	{
		filename = "Round2_cluster" + to_string(x) + ".txt";
    	aus.open(filename,ios::out);
		aus.precision(7);
		for(int z = 0; z < clusters[x].size(); z++)
		{
			outstring = "";
			for(int w = 0; w < clusters[x][z].size(); w++)
			{
				outstring = outstring + to_string(clusters[x][z][w]) + " ";
			}
			outstring = outstring + "\n";
			aus << outstring;
		}
		aus.close();
	}
}

void removeElement(std::vector<gene>& vec, int index)
{
    std::vector<gene>::iterator it = vec.begin();
    std::advance(it, index);
    vec.erase(it);
}

void removeElement(vector < vector <double> >& vec, int index)
{
    vector < vector <double> >::iterator it = vec.begin();
    std::advance(it, index);
    vec.erase(it);
}
