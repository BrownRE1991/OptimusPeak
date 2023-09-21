#include <thread>
#include "DELibrary.h"

void ReadData(string input1, list <std::vector<double> > * listOfPoints);
vector <double> findApproxLinewidth(list <std::vector<double> >  listOfPoints, int numOfPoints, int dim, vector <double> currentMax);
double findMin(vector <std::vector<double> > & listOfPoints, int numOfPoints, int dim);
vector <double> findMaxOther(vector <std::vector<double> > & listOfPoints, int numOfPoints, int dim, vector < vector <double> > currentMax, vector < vector <double> > widths);
std::vector <gene> kmodels(vector <std::vector<double> > & data, int numOfPoints, int dim, int k, vector < vector <double> > maxs, vector < vector <double> > widths, double cutoff, double Thresh);
void removeElement(std::vector<gene>& vec, int index);
void removeElement(vector < vector <double> >& vec, int index);

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
	
    int count1 = 0;
    int count2 = 0;
	double minIntensity = 100000000;
    for(it = listOfPoints.begin(); it != listOfPoints.end(); it++)
    {
		if(abs(it->at(dim)) < abs(minIntensity))
		{
		 	minIntensity = abs(it->at(dim));
		}
    }
	 
	vector < vector <double> > maxs;
	vector < vector <double> > widths;
	maxs.push_back(findMax(listOfPoints2, numOfPoints2, dim));	
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
	
	int startingK =  int(widths.size());
	
	//Put k-models here using the starting K, and then remove those models from the dataset and repeat.
	std::vector <gene> results;
	
	//Maybe do something else: Partition, fit give only a fraction of each point to each model, refit, repeat until convergence
	results = kmodels(listOfPoints, numOfPoints, dim, startingK, maxs, widths, cutoff, minIntensity);
	string outfile = "genes" + to_string(file) + ".txt";
	printGenestoFile(results, outfile, 0.2);
	
	return(0);
}

//set up starting partition based on distance to peak center
std::vector <std::vector <std::vector <double> > > startingPartition(vector < vector <double> > & data, int numOfPoints, int dim, int k, vector < vector <double> > maxs, vector < vector <double> > widths)
{
	std::vector <std::vector <std::vector <double> > > results;
	std::vector <std::vector <double> > temp;
	std::vector <double> temp1;
	
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
	
	for(int currentPoint = 0; currentPoint < numOfPoints-1; currentPoint++)
	{
		for(int x = 0; x < k; x++)
		{
			distances[x] = 0.0;
		}
		for(int x = 0; x < k; x++)
		{
			if(dim == 2)
			{
				distances[x] = ((data[currentPoint][0] - maxs[x][0])/(sizeIncrease*widths[x][0]))*((data[currentPoint][0] - maxs[x][0])/(sizeIncrease*widths[x][0])) + ((data[currentPoint][1] - maxs[x][1])/(sizeIncrease*widths[x][1]))*((data[currentPoint][1] - maxs[x][1])/(sizeIncrease*widths[x][1]));
			}
			if(dim == 3)
			{
				distances[x] = ((data[currentPoint][0] - maxs[x][0])/(sizeIncrease*widths[x][0]))*((data[currentPoint][0] - maxs[x][0])/(sizeIncrease*widths[x][0])) + ((data[currentPoint][1] - maxs[x][1])/(sizeIncrease*widths[x][1]))*((data[currentPoint][1] - maxs[x][1])/(sizeIncrease*widths[x][1])) + ((data[currentPoint][2] - maxs[x][2])/(sizeIncrease*widths[x][2]))*((data[currentPoint][2] - maxs[x][2])/(sizeIncrease*widths[x][2]));
			}
		}
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
	c.fitness = a.fitness;
	return(c);
}


int removeRunaways(std::vector <gene> & genes, vector < vector <double> > & max, vector < vector <double> > & width, vector <double> & mins, vector <double> & Uplim, int & k, int dim)
{
	int goAgainFlag = 1;
	int flag3 = 0;
	int worstInd = -1;
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
				removeElement(genes, worstInd);
				removeElement(max, worstInd);
				removeElement(width, worstInd);
				k = k - 1;
				goAgainFlag = 1;
				flag3 = 0;
			}
		}
	}
	return(flag3);
}

int checkLineWidths(std::vector <gene> & genes, vector < vector <double> > & max, vector < vector <double> > & width, int & k, int dim)
{
	int flag2 = 0;	
	int widthFlag = 0;
	int worstInd = -1;
		for(int x = (k-1); x >= 0; x--)
		{
			widthFlag = 0;
			for(int y = 0; y < dim; y++)
			{
				if((genes[x].width[y] < 0.0001) | (genes[x].width[y] > 10))
				{
					flag2 = 1;
					widthFlag = 1;
					worstInd = x;
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
				}
			}
		}
	return(flag2);
}

int runCombine(int & k, std::vector <gene> & genes, int dim, vector < vector <double> > & max, vector < vector <double> > & width, std::vector <std::vector <std::vector <double> > > & results, double cut)
{
	gene temp;
	int flag1 = 0;
	if(k > 1)
	{
		for(int x = 0; x < k; x++)
		{
			for(int y = (k-1); y >= 0; y--)
			{
				if(y != x)
				{
					flag1 = 0;
					if(genes.size() > 0)
					{
						for(int q = 0; q < dim; q++)
						{
							if(y < k && x < k)
							{
								if(abs(genes[y].center[q] - genes[x].center[q]) < cut*(genes[x].width[q]))
								{
									if(genes[x].height > genes[y].height)
									{
										flag1++;
									}
								}
							}
						}
						if(flag1 == dim)
						{
							temp = combineGenes(genes[x], genes[y], results, x, y, cut);
							//if the returned gene has a fitness below 10^18, then the combine was a success.
							if(temp.fitness < 1e+18)
							{
								genes[x] = temp;
								removeElement(genes, y);
								removeElement(max, y);
								removeElement(width, y);
								k = k-1;
							}
						}
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
				}
			}
		}
	}
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
cut high (1/30 of max maybe), and initialize to the number of likely cluster model centers.
Once these models have been found subtract the model from the data, and restart.
*/
std::vector <gene> kmodels(vector <std::vector<double> > & data, int numOfPoints, int dim, int k, vector < vector <double> > maxs, vector < vector <double> > widths, double cutoff, double Thresh)
{
	cout << "Inside K-models\n";
	vector <double> mins;
	vector <double> Uplim;
	int max_Proc = std::thread::hardware_concurrency();
	cout << max_Proc << "\n";
	
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
	
	int count = 0;
	
	//partition
	results = startingPartition(data, numOfPoints, dim, k, max, width);
	
	//if any partitions have too few points, remove the shortest model, and do starting partition again
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

	int flag3 = 0;
	int flag2 = 0;	
	flag2 = checkLineWidths(genes, max, width, k, dim);
	flag = 1;
	
	//combine
	k = runCombine(k, genes, dim, max, width, results, 1.0);
	gene temp;
	int flag1 = 0;
	

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
		
		if(k > 1)
		{
			//repartition
			results = repartition(data, genes, dim, numOfPoints, k, Thresh);
			if(k == 1)
			{
				genes[0] = DE(results[0], int(results[0].size()), dim, max[0], width[0]);
			}
			cout << k << "\n";
			if(k > 1 & k <= max_Proc)
			{
				vector<std::thread> threads(k);
				for(int x = 0; x < k; x++)
				{
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
				vector<std::thread> threads(k);
				while(count < k)
				{	
					cout << count << "\n";
					for(int x = count; x < max_Proc + count; x++)
					{
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
		}

		
		cout << "Checking linewidths\n";
		//if any have weird linewidths, remove them
		flag2 = checkLineWidths(genes, max, width, k, dim);
		flag = 0;
		worst = 0.0;
		worstInd = -1;
		//cout << "Checking Range\n";
		flag3 = removeRunaways(genes, max, width, mins, Uplim, k, dim);
		
		//cout << "Checking fitness\n";
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
			removeElement(genes, worstInd);
			removeElement(max, worstInd);
			removeElement(width, worstInd);
			k = k - 1;
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
		results = startingPartition(data, numOfPoints, dim, k, maxs, widths);
		//fit
		for(int x = 0; x < k; x++)
		{
			genes[x] = DE(results[x], int(results[x].size()), dim, maxs[x], widths[x]);
			for(int y = 0; y < dim; y++)
			{
				max[x][y] = genes[x].center[y];
				width[x][y] = genes[x].width[y];
			}
			max[x][dim] = genes[x].height;
		}
		flag2 = checkLineWidths(genes, max, width, k, dim);

		flag3 = removeRunaways(genes, max, width, mins, Uplim, k, dim);
	}
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
	
	double elipse = 0;
	int counter = 0;
	double sizeIncrease = 1.0;
	
	vector <std::vector<double> > :: iterator it;
	for(it = listOfPoints.begin(); it != listOfPoints.end(); it++)
	{
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
