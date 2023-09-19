#include <iostream>
#include <fstream>
#include <cmath>
#include <list> 
#include <sstream>
#include <vector>
//#include <ctime>
#include <algorithm>
//#include <pthread.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <cstring>
//#include <iterator>
#include <string.h>


using namespace std;
pthread_mutex_t m;
static int ** listOfRowsGraph;

void ReadData_bin(string filename, std::vector <std::vector<double> > * listOfPoints);

//from K-cluster
void ReadData(int argc, char ** argv, std::vector <std::vector<double> > * listOfPoints, int position);
void tokenize(std::string const str, const char delim, std::vector<double> * out);
std::vector<int> GenZeroVec(int dim);
double ** grabDimensions(int dim, double ** data, int numOfPoints);
int ** makeGraph(int dim);
double ** generate(int clusters, int dim, double ** dimensionVariables);
double distance2D(double pointx1, double pointy1, double pointx2, double pointy2, int max1, int max2, int min1, int min2);
double distance3D(double pointx1, double pointy1, double pointz1, double pointx2, double pointy2, double pointz2);
int * sortAndReCalc(int numOfPoints, int clusters, double ** clustercenters, int dim, double ** data);
void updateGraph(int clusters,double ** clustercenters, int ** listOfRowsGraph, int run, int numOfPoints, double ** data, int dim);
void * wrapper(void * pass1);
void Kcluster(double ** data, int dim, int numOfPoints, double scale1, double scale2, double ** dimensionVariables);
double ** normalizeSet(int numOfPoints, int dim, std::vector <std::vector<double> > * listOfPoints, double *** dimensionVariables);

//from graphcut
int numOfClusters(std::vector <std::vector<double> > * Graph, std::vector <std::vector<double> > * listOfPoints, double cut, int flag, std::vector <int> * finalList, int dim, int totalRuns);
std::vector <int> GraphCut(std::vector <std::vector<double> > Graph, std::vector <std::vector<double> > listOfPoints, int limit, int dim, int totalRuns);
void normalize(std::vector <std::vector <double> > * listOfRowsGraph, int runs);
void * wrapper2(void * pass1);

//from ASRG
void ASRG(std::list <std::vector <double> > Seeds, std::vector <std::vector<double> >  listOfPoints, double ** dimensionVariables, struct imageStats imstat, std::vector <std::vector <std::vector <int> > > dataGraph, int numOfPoints, int dim);
void ASRGsingleCluster2D(int start, std::vector <int> * cluster, int currentCluster, int dim, std::vector <std::vector<double> > * listOfPoints, double ** dimensionVariables, struct imageStats imstat, std::vector <std::vector <std::vector <int> > > * dataGraph, int numOfPoints);
void writeFile2(std::vector <int> finalList, std::vector <std::vector<double> > listOfPoints, int numOfClusters, int numOfPoints, int dim, int flag);
int findListIndex(std::vector <double> input, std::vector <std::vector<double> >  listOfPoints, int numOfPoints, int dim);
std::vector <std::vector <std::vector <int> > > generateImageGraph(double ** dimensionVariables, double ** data, int dim, int numOfPoints, struct imageStats * imstat, vector < vector <int> > * transistionList);
std::vector <int> convertToIndex(vector<double> input, double ** dimensionVariables, struct imageStats imstat, int dim);
std::vector <double> convertToFrequancy(vector<int> input, double ** dimensionVariables, struct imageStats imstat, int dim);
	
int maxC = 100;
int minC = 10;

struct pass
{
    double ** data;
    int dim;
    int numOfPoints;
    double scale1;
    double scale2;
    double ** dimesnionVariables;
};

struct pass2
{
    std::vector <std::vector<double> > * Graph;
    std::vector <std::vector<double> > * listOfPoints;
    double cut;
    int flag;
    std::vector<double> * vectorOfClusters;
    int index;
    int dim;
    int totalRuns;
};

struct imageStats
{
    double smallestDistx;
    double smallestDisty;
    double smallestDistz;
    double smallestDistw;
    double xintervals;
    double yintervals;
    double zintervals;
    double wintervals;
};

void ReadData(int argc, char ** argv, std::vector <std::vector<double> > * listOfPoints, int position)
{//read in the data

    fstream input;
    char buf1[255];
    string buf2;
    list <string> lines;
    std::vector<double> out;
    std::vector<double> out2;
    out.erase(out.begin(), out.end());
    //int totalLength = length*length;
    
    string inputFile = argv[position];
    cout << inputFile << "\n";
    
    if(position == 1)
    {
        input.open(inputFile,ios::in);
    
        input.getline(buf1, 255); //Get first line
    
        while(strlen(buf1)>3)//while we have data
        {
            lines.push_back(buf1);
            input.getline(buf1, 255);
        }
        input.close();
    
        list <string> :: iterator it;
        for(it = lines.begin(); it != lines.end(); ++it)
        {
            buf2 = *it;
            tokenize(buf2, ' ', &out);
            std::vector<double> buffer = std::vector<double> (out);
            out.erase(out.begin(), out.end());
            listOfPoints->push_back(buffer);
        }
    }
}

void tokenize(std::string const str, const char delim, std::vector<double> * out)
{
	std::stringstream ss(str);
	std::string s;
	std::vector<string> temp;
	while(std::getline(ss, s, delim))
	{
		temp.push_back(s);
	}
	
	std::vector<string>:: iterator it;  
    for(it = temp.begin(); it != temp.end(); ++it) 
	{
		out->push_back(stod(*it));
	}
}

std::vector<int> GenZeroVec(int dim)
{
	std::vector<int> out;
	for(int x = 0; x < dim; x++)
	{
		out.push_back(0);
	}
	return(out);
}

//std::vector <std::vector<int> > makeGraph(int dim)
//{
//    std::vector <std::vector<int> > listOfRowsGraph;
//	for(int x = 0; x < dim; x++)
//	{
//		std::vector<int> out;
//		out = GenZeroVec(dim);
//		std::vector<int> buffer = std::vector<int>(out);
//		out.erase(out.begin(), out.end());
//		listOfRowsGraph.push_back(buffer);
//	}
//    return(listOfRowsGraph);
//}

int ** makeGraph(int dim)
{
    int ** listOfRowsGraph;
    listOfRowsGraph = (int**)malloc(dim*sizeof(int*));
    for(int x = 0; x < dim; x++)
    {
        listOfRowsGraph[x] = (int*)malloc(dim*sizeof(int));
        for(int y = 0; y < dim; y++)
        {
            listOfRowsGraph[x][y] = 0;
        }
        //std::vector<int> buffer(dim);
        //listOfRowsGraph.push_back(buffer);
    }
    return(listOfRowsGraph);
}

double ** generate(int clusters, int dim, double ** dimensionVariables)
{
	double ** result = NULL;
	double maxX = dimensionVariables[0][1];
	double minX = dimensionVariables[0][0];
	double maxY = 0;
	double minY = 0;
	double maxZ = 0;
	double minZ = 0;
	double maxW = 0;
	double minW = 0;
	if(dim > 1)
	{
		maxY = dimensionVariables[1][1];
		minY = dimensionVariables[1][0];
	}
	if(dim > 2)
	{
		maxZ = dimensionVariables[2][1];
		minZ = dimensionVariables[2][0];
	}
	if(dim > 3)
	{
		maxW = dimensionVariables[3][1];
		minW = dimensionVariables[3][0];
	}
	result = (double**)(malloc(clusters*sizeof(double*)));
	for(int x = 0; x < clusters; x++)
	{
		result[x] = (double*)(malloc(dim*sizeof(double)));
		for(int y = 0; y < dim; y++)
 		{
			result[x][0] = double(((rand()%100)/100.0)*(maxX - minX) + minX);
			if(dim > 1)
			{
				result[x][1] = double(((rand()%100)/100.0)*(maxY - minY) + minY);
			}
			if(dim > 2)
			{
				result[x][2] = double(((rand()%100)/100.0)*(maxZ - minZ) + minZ);
			}
			if(dim > 3)
			{
				result[x][3] = double(((rand()%100)/100.0)*(maxW - minW) + minW);
			}
		}
	}
	return(result);
}

double distance2D(double pointx1, double pointy1, double pointx2, double pointy2, int max1, int max2, int min1, int min2)
{
	double x = abs(pointx1 - pointx2);
	double y = abs(pointy1 - pointy2);
	x = x*x;
	y = y*y;
	double result = sqrt(x+y);
	return(result);
}

double distance3D(double pointx1, double pointy1, double pointz1, double pointx2, double pointy2, double pointz2)
{
	double x = abs(pointx1 - pointx2);
	double y = abs(pointy1 - pointy2);
	double z = abs(pointz1 - pointz2);
	x = x*x;
	y = y*y;
	z = z*z;
	double result = sqrt(x+y+z);
	return(result);
}

int * sortAndReCalc(int numOfPoints, int clusters, double ** clustercenters, int dim, double ** data)
{
	list <int> SortedCluster[clusters];
	list <list <int> > :: iterator it1;
    double dist = 0;
	for(int x = 0; x < numOfPoints; x++)
	{
		double min = 10000;
		int whichCluster = -1;
        double pointx1;
        double pointy1;
        double pointz1;
        double pointx2;
        double pointy2;
        double pointz2;
		for(int y = 0; y < clusters; y++)
		{
            if(dim == 2)
            {
                pointx1 = data[x][0];
                pointy1 = data[x][1];
                pointx2 = clustercenters[y][0];
                pointy2 = clustercenters[y][1];
                dist = distance2D(pointx1, pointy1, pointx2, pointy2, 1, 1, 0, 0);
            }
            if(dim == 3)
            {
                pointx1 = data[x][0];
                pointy1 = data[x][1];
                pointz1 = data[x][2];
                pointx2 = clustercenters[y][0];
                pointy2 = clustercenters[y][1];
                pointz2 = clustercenters[y][2];
                dist = distance3D(pointx1, pointy1, pointz1, pointx2, pointy2, pointz2);
            }
			if(dist < min)
			{
				min = dist;
				whichCluster = y;
			}
		}
		SortedCluster[whichCluster].push_back(x);
	}
	
	list <int> :: iterator it3;
	
	int * ClusterCounts = (int*)malloc(clusters*sizeof(int));
	//cout << "ClusterCounts\n";
	int sum = 0;
	for(int w = 0; w < clusters; w++)
	{
		ClusterCounts[w] = (int)SortedCluster[w].size();
		//cout << ClusterCounts[w] << "\n";
		sum = sum + ClusterCounts[w];
	}
	if(sum != numOfPoints)
	{
		cout << "ERROR\n";
	}
	
	//calcNewCenters
	
	for(int x = 0; x < clusters; x++)
	{
		double sum1 = 0;
		double sum2 = 0;
		double sum3 = 0;
		double sum4 = 0;
		int length = (int)SortedCluster[x].size();
		if(length > 0)
		{
			for(it3 = SortedCluster[x].begin(); it3 != SortedCluster[x].end(); it3++)
			{
				sum1 = sum1 + data[(*it3)][0];
				if(dim > 1)
				{
					sum2 = sum2 + data[(*it3)][1];
				}
				if(dim > 2)
				{
					sum3 = sum3 + data[(*it3)][2];
				}
				if(dim > 3)
				{
					sum4 = sum4 + data[(*it3)][3];
				}
			}
			clustercenters[x][0] = sum1/length;
			if(dim > 1)
			{
				clustercenters[x][1] = sum2/length;
			}
			if(dim > 2)
			{
				clustercenters[x][2] = sum3/length;
			}
			if(dim > 3)
			{
				clustercenters[x][3] = sum4/length;
			}
		}
	}
	return(ClusterCounts);
}

void updateGraph(int clusters, double ** clustercenters, int ** listOfRowsGraph2, int run, int numOfPoints, double ** data, int dim)
{
	list <int> SortedCluster[clusters];
	list <list <int> > :: iterator it1;
    double dist = 0;
    for(int x = 0; x < numOfPoints; x++)
    {
        double min = 10000;
        int whichCluster = -1;
        double pointx1;
        double pointy1;
        double pointz1;
        double pointx2;
        double pointy2;
        double pointz2;
        for(int y = 0; y < clusters; y++)
        {
            if(dim == 2)
            {
                pointx1 = data[x][0];
                pointy1 = data[x][1];
                pointx2 = clustercenters[y][0];
                pointy2 = clustercenters[y][1];
                dist = distance2D(pointx1, pointy1, pointx2, pointy2, 1, 1, 0, 0);
            }
            if(dim == 3)
            {
                pointx1 = data[x][0];
                pointy1 = data[x][1];
                pointz1 = data[x][2];
                pointx2 = clustercenters[y][0];
                pointy2 = clustercenters[y][1];
                pointz2 = clustercenters[y][2];
                dist = distance3D(pointx1, pointy1, pointz1, pointx2, pointy2, pointz2);
            }
            if(dist < min)
            {
                min = dist;
                whichCluster = y;
            }
        }
        SortedCluster[whichCluster].push_back(x);
    }
	
	list <int> :: iterator it3;
	list <int> :: iterator it4;
	for(int y = 0; y < clusters; y++)
	{
        if(SortedCluster[y].size() > 0)
        {
            for(it3 = SortedCluster[y].begin(); it3 != SortedCluster[y].end(); it3++)
            {
                for(it4 = SortedCluster[y].begin(); it4 != SortedCluster[y].end(); it4++)
                {
                    (listOfRowsGraph2)[(*it3)][(*it4)]++;
                }
            }
        }
	}
}

double ** grabDimensions(int dim, double ** data, int numOfPoints)
{
    double ** dimensionVariables = NULL;
    dimensionVariables = (double**)malloc(dim*sizeof(double*));
    if(dimensionVariables == NULL)
    {
        cout << "ERROR\n";
        //return(1);
    }
    //Dimension variables holds the min and max for each dimension to be used for initialization of cluster centers
    for(int x = 0; x < dim; x++)
    {
        dimensionVariables[x] = (double*)malloc(2*sizeof(double));
        if(dimensionVariables[x] == NULL)
        {
            cout << "ERROR\n";
        }
        dimensionVariables[x][0] = 10000.0;
        dimensionVariables[x][1] = 0.0;
    }
    for(int x = 0; x< numOfPoints; x++)
    {
        if(data[x][0] < dimensionVariables[0][0])
        {
            dimensionVariables[0][0] = data[x][0];
        }
        if(data[x][0] > dimensionVariables[0][1])
        {
            dimensionVariables[0][1] = data[x][0];
        }
        if(dim > 1)
        {
            if(data[x][1] < dimensionVariables[1][0])
            {
                dimensionVariables[1][0] = data[x][1];
               }
                if(data[x][1] > dimensionVariables[1][1])
            {
                dimensionVariables[1][1] = data[x][1];
            }
        }
        if(dim > 2)
        {
            if(data[x][2] < dimensionVariables[2][0])
            {
                dimensionVariables[2][0] = data[x][2];
            }
            if(data[x][2] > dimensionVariables[2][1])
            {
                dimensionVariables[2][1] = data[x][2];
            }
        }
    }
    return(dimensionVariables);
}

double ** normalizeSet(int numOfPoints, int dim, std::vector <std::vector<double> > * listOfPoints, double *** dimensionVariables)
{
    std:vector <std::vector<double> > :: iterator it;
    std::vector<double>:: iterator it2;
    double ** data = (double**)malloc(numOfPoints*sizeof(double*)); //this will hold the normalized dataset
    if(data != NULL)
    {
        for(int x = 0; x < numOfPoints; x++)
        {
            data[x] = (double*)malloc(dim*sizeof(double));
        }
    }
    //double ** dimensionVariables = NULL;
    int count1 = 0;
    int count2 = 0;
    //find range of both dimensions
    for(it = listOfPoints->begin(); it != listOfPoints->end(); it++)
    {
        if(count1 < numOfPoints)
        {
            count2 = 0;
            for (it2 = it->begin(); it2 != it->end(); it2++)
            {
                if(count2 < dim)
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
    (*dimensionVariables) = grabDimensions(dim, data, numOfPoints);
    
    //normalize the dimensions
    for(int x = 0; x < numOfPoints; x++)
    {
        //count2 = 0;
        for(int y = 0; y < dim; y++)
        {
            if(x < numOfPoints)
            {
                if(y == 0)
                {
                    data[x][y] = (data[x][y] - (*dimensionVariables)[0][0])/((*dimensionVariables)[0][1] - (*dimensionVariables)[0][0]);
                }
                if(y == 1)
                {
                    data[x][y] = (data[x][y] - (*dimensionVariables)[1][0])/((*dimensionVariables)[1][1] - (*dimensionVariables)[1][0]);
                }
                if(y == 2)
                {
                    data[x][y] = (data[x][y] - (*dimensionVariables)[2][0])/((*dimensionVariables)[2][1] - (*dimensionVariables)[2][0]);
                }
                if(y == 3)
                {
                    data[x][y] = (data[x][y] - (*dimensionVariables)[3][0])/((*dimensionVariables)[3][1] - (*dimensionVariables)[3][0]);
                }
            }
        }
    }
    for(int x = 0; x < dim; x++)
    {
        free((*dimensionVariables)[x]);
    }
    free((*dimensionVariables));
    (*dimensionVariables) = grabDimensions(dim, data, numOfPoints);
    return(data);
}

void Kcluster(double ** data, int dim, int numOfPoints, double scale1, double scale2, double ** dimensionVariables)
{
    int maxC = numOfPoints/scale1;
    int minC = numOfPoints/scale2;
    
    int run = 0;
    int totalRuns = 25;
    int clusters = 10;
    double ** clustercenters = nullptr;
    int * OldClusterCounts = nullptr;
    int * ClusterCounts = nullptr;
    int differentFlag = 0;
    while(run < totalRuns)
    {
        clusters = int(((rand()%100)/100.0)*(maxC - minC) + minC);
        if(clustercenters != nullptr)
        {
            free(clustercenters);
        }
        clustercenters = generate(clusters, dim, dimensionVariables);
    
        int continueFlag = 0;
        int startFlag = 0;
    
        while(continueFlag < 2)
        {

            //sort, record changes, calcNewCenter
            if(ClusterCounts != nullptr)
            {
                free(ClusterCounts);
                ClusterCounts = nullptr;
            }
            ClusterCounts = sortAndReCalc(numOfPoints, clusters, clustercenters, dim, data);
            if(startFlag == 0)
            {
                OldClusterCounts = (int*)malloc(clusters*(sizeof(int)));
                for(int x = 0; x < clusters; x++)
                {
                    OldClusterCounts[x] = ClusterCounts[x];
                }
                free(ClusterCounts);
                ClusterCounts = nullptr;
                startFlag = 1;
            }
            else
            {
                differentFlag = 0;
                for(int x = 0; x < clusters; x++)
                {
                    if(OldClusterCounts[x] != ClusterCounts[x])
                    {
                        differentFlag = 1;
                    }
                }
                if(differentFlag == 1)
                {
                    for(int x = 0; x < clusters; x++)
                    {
                        OldClusterCounts[x] = ClusterCounts[x];
                    }
                    free(ClusterCounts);
                    ClusterCounts = nullptr;
                }
                if(differentFlag == 0)
                {
                    continueFlag++;
                    if(continueFlag < 10)
                    {
                        free(ClusterCounts);
                        ClusterCounts = nullptr;
                    }
                }
            }
        }//if changes go back to sort
        //if no changes report clusters
        if(OldClusterCounts != nullptr)
        {
            free(OldClusterCounts);
            OldClusterCounts = nullptr;
        }

        //then update graph for graphcut
        pthread_mutex_lock(&m);
        updateGraph(clusters, clustercenters, listOfRowsGraph, run, numOfPoints, data, dim);
        pthread_mutex_unlock(&m);
        run++;
        for(int x = 0; x < clusters; x++)
        {
            free(clustercenters[x]);
        }
        free(clustercenters);
        clustercenters = nullptr;
        free(ClusterCounts);
        ClusterCounts = nullptr;
        if(OldClusterCounts != nullptr)
        {
            free(OldClusterCounts);
            OldClusterCounts = nullptr;
        }
    }
    //do this a bunch of times with random cluster centers
}

void * wrapper(void * pass1)
{
    struct pass * temp = (struct pass *)pass1;
    Kcluster((temp->data), temp->dim, temp->numOfPoints, temp->scale1, temp->scale2, temp->dimesnionVariables);
    return(nullptr);
}

int numOfClusters(std::vector <std::vector<double> > * Graph, std::vector <std::vector<double> > * listOfPoints, double cut, int flag, std::vector <int> * finalList, int dim, int totalRuns)
{
    int length = (int)listOfPoints->size();
    int currentCluster = 0;
    for(int x = 0; x < length; x++)
    {
        finalList->push_back(0);
    }
    std::vector <int> temp;
    int counter = 0;
    int nextStart = 0;
    std::list <std::vector<int> >:: iterator it2;
    int inFlag = 0;
    while(nextStart < length)
    {
        currentCluster++;
        counter = 0;
        for(int x = 0; x < length; x++)
        {
            if((*Graph)[nextStart][x] > cut)
            {
                if((*finalList)[x] < 1)
                {
                    (*finalList)[x] = currentCluster;
                }
            }
        }
        inFlag = 0;
        for(int y = 0; y < counter; y++)
        {
            for(int x = 0; x < length-x; x++)
            {
                if((*Graph)[y][x] > cut)
                 {
                     if((*finalList)[x] < 1)
                     {
                         (*finalList)[x] = currentCluster;
                     }
                }
            }
        }
        for(int y = 0; y < counter; y++)
        {
            for(int x = 0; x < length; x++)
            {
                if((*Graph)[y][x] > cut)
                 {
                     if((*finalList)[x] < 1)
                     {
                         (*finalList)[x] = currentCluster;
                     }
                }
            }
        }
        int w = 0;
        while((*finalList)[w]>0)
        {
            w++;
        }
        
        nextStart = w;
    }
    if(flag == 1)
    {
        //set up output for the clusters
        writeFile2((*finalList), (*listOfPoints), currentCluster, length, dim, 0);
    }
    return(currentCluster);
}

void * wrapper2(void * pass1)
{
    int result = 0;
    struct pass2 * temp = (struct pass2 *)pass1;
    std::vector <int> finalList2;
    result = numOfClusters((temp->Graph),(temp->listOfPoints), temp->cut, temp->flag, &finalList2, temp->dim, temp->totalRuns);
    pthread_mutex_lock(&m);
    //cout << "CUT: " << temp->cut << " Clusters: " << result << "\n";
    (*temp->vectorOfClusters)[temp->index] = result;
    pthread_mutex_unlock(&m);
    return(nullptr);
}

std::vector <int> GraphCut(std::vector <std::vector<double> > Graph, std::vector <std::vector<double> > listOfPoints, int limit, int dim, int totalRuns)
{
    std::vector <int> finalList;
    int length = 0;
    float steps = 0.1;
    length = (int)listOfPoints.size();
    if(length > 100000)
    {
        steps = 0.1;
    }
    std::vector<struct pass2> forMulti;
    std::vector<double> vectorOfClusters;
    struct pass2 temp;
    for(float x = 0; x < 1; x = x + steps)
    {
        vectorOfClusters.push_back(0); //initalize vectorOfClusters
        temp.Graph = &Graph;   //set up variables to pass to the threads
        temp.listOfPoints = &listOfPoints;
        temp.cut = x;
        temp.flag = 0;
        temp.index = x/steps;
        temp.vectorOfClusters = &vectorOfClusters;
        temp.dim = dim;
        temp.totalRuns = totalRuns;
        forMulti.push_back(temp);
    }//this works
    
    //determine a safe number of threads that does not overload the memory
    int Numthreads = 3;
    //3188 very good at 5 threads, more than that just seems to swap too much
    if(length < 5000)
    {
        Numthreads = 3;
    }
    if(length > 10000)
    {
        Numthreads = 3;
    }
    if(length > 100000)
    {
        Numthreads = 2;
    }
    if(length > 1000000)
    {
        Numthreads = 1;
    }
    pthread_t threads[Numthreads];
    int rc;
    int i;
    void * out_void = nullptr;
    
    int counter = 0;
    //run numOfClusters for each cut as a seperate thread
    while(counter < 1/steps)
    {
        if(abs(counter - (1/steps)) < (Numthreads))
        {
            Numthreads = abs(counter - (1/steps));
        }
        for( i = 0; i < Numthreads; i++ )
        {
            rc = pthread_create(&threads[i], NULL, &wrapper2, (void *)&forMulti[counter+i]);
            if (rc) {
                cout << "Error:unable to create thread," << rc << endl;
                exit(-1);
            }
        }
        //cout << "done this bit\n";
        out_void = nullptr;
        for( i = 0; i < Numthreads; i++ )
        {
            pthread_join(threads[i], &out_void);
            counter++;
        }
        cout << "COUNTER: " << counter << "\n";
    }

    //output clusters
    for(float x = 0; x < 1; x = x + steps)
    {
        cout << "CUT: " << x << " Clusters: " << (vectorOfClusters[x/steps]) << "\n";
        //cout << x*steps << " " << (vectorOfClusters[x/steps]) << "\n";
    }//havent tried this yet
    
    
    //this section finds an optimal cut location (still working on it)
    //std::vector<double> :: iterator it2;
    std::vector<double> :: iterator it3;
    int mostCommon = 0;
    int holder;
    int most = 0;
    int count = 0;
    int first1 = 100;
    int last1 = 0;
    int first2 = 100;
    int last2 = 0;
    //try triangle?
    it3 = vectorOfClusters.begin();
    //double firstPoint = (*it3);
    //double len = vectorOfClusters.size();
    //it3 = next(it3,len-1);
    //double lastPoint = (*it3);
    //double slope = (lastPoint-firstPoint)/(0.95);
    //double b = lastPoint - (slope*0.95);
    //cout << firstPoint << " " << lastPoint << " " << slope << " " << b << "\n";
    //maybe use triangle to find good cut rather than most stable (too many possibilities to get most stable)
    //Maybe just pick first large jump? Want to be more conservative.
    
    
    //find the longest stretch that has a difference with the middle value of less than limit
    for(it3 = vectorOfClusters.begin(); it3 != vectorOfClusters.end(); it3++)
    {
        count = 0;
        holder = *it3;
        //cout << "First index: " << (*it3) << "\n";
        for(int x = 0; x < 1/steps; x++)
        {
            //cout << vectorOfClusters[x] << " " << (holder- vectorOfClusters[x]) << "\n";
            if(abs(holder - vectorOfClusters[x]) < limit)
            {
                count++;
                if(x < first1)
                {
                    first1 = x;
                }
                if(x > last1)
                {
                    last1 = x;
                }
            }
        }
        //cout << "Count: " << count << "\n";
        //cout << "Most: " << most << "\n";
        //cout << "MostCommon: " << mostCommon << "\n";
        //cout << "Starts at " << first1 << " Ends at " << last1 << "\n";
        if(count > most)
        {
            mostCommon = holder;
            most = count;
            first2 = first1;
            last2 = last1;
            //first1 = 100;
            //last1 = 0;
        }
        first1 = 100;
        last1 = 0;
        //try this later to restart first and last every time?
    }
    
    //use the stretch found in the preivous section to determine cut point
    //cout << first2 << " " << last2 << "\n";
    float middle = (first2 + last2)/2.0;
    
    //cout <<  middle << " " << middle*steps << "\n";
    
    numOfClusters(&Graph, &listOfPoints, middle*steps, 1, &finalList, dim, totalRuns);
    return(finalList);
}

void normalize(std::vector <std::vector <double> > * listOfRowsGraph, int runs)
{
    std::vector <std::vector<double> >:: iterator it;
    std::vector<double>:: iterator it3;
    for(it = listOfRowsGraph->begin(); it != listOfRowsGraph->end(); it++)
    {
        for(it3 = (*it).begin(); it3 != (*it).end(); it3++)
        {
            (*it3) = (*it3)/runs;
        }
    }
}

std::vector <double> convertToFrequancy(vector<int> input, double ** dimensionVariables, struct imageStats imstat, int dim)
{
    std::vector <double> temp;
    double xindex = input[0]*imstat.smallestDistx + dimensionVariables[0][0];
    temp.push_back(xindex);
    double yindex = input[1]*imstat.smallestDisty + dimensionVariables[1][0];
    temp.push_back(yindex);
    double zindex = 0;
    double windex = 0;
    if(dim > 2)
    {
        zindex = input[2]*imstat.smallestDistz + dimensionVariables[2][0];
        temp.push_back(zindex);
    }
    if(dim > 3)
    {
        windex = input[3]*imstat.smallestDistw + dimensionVariables[3][0];
        temp.push_back(windex);
    }
    return(temp);
}
std::vector <int> convertToIndex(vector<double> input, double ** dimensionVariables, struct imageStats imstat, int dim)
{
    std::vector <int> temp;
    int xindex = int((input[0] - dimensionVariables[0][0])/imstat.smallestDistx);
    temp.push_back(xindex);
    int yindex = int((input[1] - dimensionVariables[1][0])/imstat.smallestDisty);
    temp.push_back(yindex);
    int zindex = 0;
    int windex = 0;
    if(dim > 2)
    {
        zindex = int((input[2] - dimensionVariables[2][0])/imstat.smallestDistz);
        temp.push_back(zindex);
    }
    if(dim > 3)
    {
        windex = int((input[3] - dimensionVariables[3][0])/imstat.smallestDistw);
        temp.push_back(windex);
    }
    return(temp);
}

//*Maybe find a way to save a list of each datapoint and its coordinates in the graph so we don't have to keep looking?*//
std::vector <std::vector <std::vector <int> > > generateImageGraph(double ** dimensionVariables, double ** data, int dim, int numOfPoints, struct imageStats * imstat, vector < vector <int> > * transistionList)
{
    std::vector <std::vector <std::vector <int> > > returnable;
    std::vector <double> xcoordinates;
    std::vector <double> ycoordinates;
    std::vector <double> zcoordinates;
    std::vector <double> wcoordinates;
    
    for(int x = 0; x < numOfPoints; x++)
    {
        int inflag = 0;
        if(xcoordinates.size() == 0)
        {
            xcoordinates.push_back(data[x][0]);
        }
        if(xcoordinates.size() != 0)
        {
            for(int y = 0; y < xcoordinates.size(); y++)
            {
                if(xcoordinates[y] == data[x][0])
                {
                    inflag = 1;
                }
            }
            if(inflag == 0)
            {
                xcoordinates.push_back(data[x][0]);
            }
        }
        if(ycoordinates.size() == 0)
        {
            ycoordinates.push_back(data[x][1]);
        }
        if(ycoordinates.size() != 0)
        {
            for(int y = 0; y < ycoordinates.size(); y++)
            {
                if(ycoordinates[y] == data[x][1])
                {
                    inflag = 1;
                }
            }
            if(inflag == 0)
            {
                ycoordinates.push_back(data[x][1]);
            }
        }
        if(dim > 2)
        {
            if(zcoordinates.size() == 0)
            {
                zcoordinates.push_back(data[x][2]);
            }
            if(zcoordinates.size() != 0)
            {
                for(int y = 0; y < zcoordinates.size(); y++)
                {
                    if(zcoordinates[y] == data[x][2])
                    {
                        inflag = 1;
                    }
                }
                if(inflag == 0)
                {
                    zcoordinates.push_back(data[x][2]);
                }
            }
        }
        if(dim > 3)
        {
            if(wcoordinates.size() == 0)
            {
                wcoordinates.push_back(data[x][3]);
            }
            if(wcoordinates.size() != 0)
            {
                for(int y = 0; y < wcoordinates.size(); y++)
                {
                    if(wcoordinates[y] == data[x][3])
                    {
                        inflag = 1;
                    }
                }
                if(inflag == 0)
                {
                    wcoordinates.push_back(data[x][3]);
                }
            }
        }
    }
    imstat->smallestDistx = 10000;
    for(int y = 0; y < xcoordinates.size(); y++)
    {
        for(int w = 0; w < xcoordinates.size(); w++)
        {
            if(w != y)
            {
                if(abs(xcoordinates[w] - xcoordinates[y]) < imstat->smallestDistx)
                {
                    imstat->smallestDistx = abs(xcoordinates[w] - xcoordinates[y]);
                }
            }
        }
    }
    imstat->xintervals = (dimensionVariables[0][1] - dimensionVariables[0][0])/imstat->smallestDistx;
    imstat->smallestDisty = 10000;
    for(int y = 0; y < ycoordinates.size(); y++)
    {
        for(int w = 0; w < ycoordinates.size(); w++)
        {
            if(w != y)
            {
                if(abs(ycoordinates[w] - ycoordinates[y]) < imstat->smallestDisty)
                {
                    imstat->smallestDisty = abs(ycoordinates[w] - ycoordinates[y]);
                }
            }
        }
    }
    imstat->yintervals = (dimensionVariables[1][1] - dimensionVariables[1][0])/imstat->smallestDisty;
    
    if(dim > 2)
    {
        imstat->smallestDistz = 10000;
        for(int y = 0; y < zcoordinates.size(); y++)
        {
            for(int w = 0; w < zcoordinates.size(); w++)
            {
                if(w != y)
                {
                    if(abs(zcoordinates[w] - zcoordinates[y]) < imstat->smallestDistz)
                    {
                        imstat->smallestDistz = abs(zcoordinates[w] - zcoordinates[y]);
                    }
                }
            }
        }
        imstat->zintervals = (dimensionVariables[2][1] - dimensionVariables[2][0])/imstat->smallestDistz;
    }

    vector<int> temp;
    std::vector <vector<int> > temp2;
    //cout << xcoordinates.size() << " " << ycoordinates.size() << " " << zcoordinates.size() << "\n";
    for(int x = 0; x < imstat->xintervals; x++)
    {
        temp.push_back(0);
    }
    for(int y = 0; y < imstat->yintervals; y++)
    {
        temp2.push_back(temp);
    }
    if(zcoordinates.size() > 0)
    {
        for(int z = 0; z < imstat->zintervals; z++)
        {
            returnable.push_back(temp2);
        }
    }
    else
    {
        returnable.push_back(temp2);
    }
    
    //vector <int> vectorHolder;
    for(int x = 0; x < numOfPoints; x++)
    {
    	//while(vectorHolder.size() > 0) vectorHolder.pop_back();
        int xindex = int((data[x][0] - dimensionVariables[0][0])/imstat->smallestDistx);
        //vectorHolder.push_back(xindex);
        int yindex = int((data[x][1] - dimensionVariables[1][0])/imstat->smallestDisty);
        //vectorHolder.push_back(yindex);
        int zindex = 0;
        if(dim > 2)
        {
            zindex = int((data[x][2] - dimensionVariables[2][0])/imstat->smallestDistz);
            //vectorHolder.push_back(zindex);
        }
        if(dim == 2)
        {
            //returnable[0][yindex][xindex] = returnable[0][yindex][xindex] + 1;
            returnable[0][yindex][xindex] = x + 1;
            //transistionList->push_back(vectorHolder);
        }
        if(dim == 3)
        {
            //returnable[zindex][yindex][xindex] = returnable[zindex][yindex][xindex] + 1;
            returnable[zindex][yindex][xindex] = x + 1;
            //transistionList->push_back(vectorHolder);
        }
    }
    return(returnable);
}

//*Maybe find a way to save this so we don't have to keep looking?*//
int findListIndex(std::vector <double> input, std::vector <std::vector<double> >  listOfPoints, int numOfPoints, int dim)
{
    int smallest = 0;
    int largest = numOfPoints;
    //std:: vector <int> recent;
    if(abs(input[0] - listOfPoints[0][0]) < 0.0001)
    {
        if(abs(input[1] - listOfPoints[0][1]) < 0.0001)
        {
            return(0);
        }
    }
    int compare = numOfPoints/2;
    while(compare < numOfPoints)
    {
        if(largest - smallest == 1)
        {
            return(-1);
        }
        //recent.push_back(compare);
        if(input[0] > listOfPoints[compare][0])
        {
            if(compare < largest)
            {
                largest = compare;
            }
            compare = (largest - smallest)/2 + smallest;
        }
        if((abs(input[0] - listOfPoints[compare][0]) < 0.0001) && (compare < numOfPoints))
        {
            if(dim == 2)
            {
                while((compare < numOfPoints) && abs(input[0] - listOfPoints[compare][0]) < 0.0001)
                {
                    if(largest - smallest == 1)
                    {
                        return(-1);
                    }
                    if((int)(input[1] > listOfPoints[compare][1]))
                    {
                        compare = compare - 1;
                    }
                    if(abs(input[1] - listOfPoints[compare][1]) < 0.0001)
                    {
                        if(abs(input[0] - listOfPoints[compare][0]) < 0.0001)
                        {
                            return(compare);
                        }
                    }
                    if(input[1] < listOfPoints[compare][1])
                    {
                        compare = compare + 1;
                    }
                    if(abs(input[1] - listOfPoints[compare][1]) < 0.0001)
                    {
                        if(abs(input[0] - listOfPoints[compare][0]) < 0.0001)
                        {
                            return(compare);
                        }
                    }
                }
            }
            if(dim == 3)
            {
                while((compare < numOfPoints) && abs(input[0] - listOfPoints[compare][0]) < 0.0001)
                {
                    if(largest - smallest == 1)
                    {
                        return(-1);
                    }
                    if(input[1] > listOfPoints[compare][1])
                    {
                        if(compare < largest)
                        {
                            largest = compare;
                        }
                        compare = (largest - smallest)/2 + smallest;
                    }
                    if((abs(input[0] - listOfPoints[compare][0]) < 0.0001) && (compare < numOfPoints))
                    {
                        while((compare < numOfPoints) && abs(input[0] - listOfPoints[compare][0]) < 0.0001 && abs(input[1] - listOfPoints[compare][1]) < 0.0001)
                        {
                            if(largest - smallest == 1)
                            {
                                return(-1);
                            }
                            if((int)(input[2] > listOfPoints[compare][2]))
                            {
                                if(compare < largest)
                                {
                                    largest = compare;
                                }
                                compare = compare - 1;
                            }
                            if(abs(input[2] - listOfPoints[compare][2]) < 0.0001)
                            {
                                if(abs(input[1] - listOfPoints[compare][1]) < 0.0001)
                                {
                                    if(abs(input[0] - listOfPoints[compare][0]) < 0.0001)
                                    {
                                        return(compare);
                                    }
                                }
                            }
                            if(input[2] < listOfPoints[compare][2])
                            {
                                if(compare > smallest)
                                {
                                    smallest = compare;
                                }
                                compare = compare + 1;
                            }
                            if(abs(input[2] - listOfPoints[compare][2]) < 0.0001)
                            {
                                if(abs(input[1] - listOfPoints[compare][1]) < 0.0001)
                                {
                                    if(abs(input[0] - listOfPoints[compare][0]) < 0.0001)
                                    {
                                        return(compare);
                                    }
                                }
                            }
                        }
                    }
                    if(input[1] < listOfPoints[compare][1])
                    {
                        if(compare > smallest)
                        {
                            smallest = compare;
                        }
                        compare = (largest - smallest)/2 + smallest;
                    }
                }
            }
        }
        if(input[0] < listOfPoints[compare][0])
        {
            if(compare > smallest)
            {
                smallest = compare;
            }
            compare = (largest - smallest)/2 + smallest;
        }
        if((abs(input[0] - listOfPoints[compare][0]) < 0.0001) && (compare < numOfPoints))
        {
            if(dim == 2)
            {
                while((compare < numOfPoints) && abs(input[0] - listOfPoints[compare][0]) < 0.0001)
                {
                    if(largest - smallest == 1)
                    {
                        return(-1);
                    }
                    if((int)(input[1] > listOfPoints[compare][1]))
                    {
                        if(compare < largest)
                        {
                            largest = compare;
                        }
                        compare = compare - 1;
                    }
                    if(abs(input[1] - listOfPoints[compare][1]) < 0.0001)
                    {
                        if(abs(input[0] - listOfPoints[compare][0]) < 0.0001)
                        {
                            return(compare);
                        }
                    }
                    if(input[1] < listOfPoints[compare][1])
                    {
                        if(compare > smallest)
                        {
                            smallest = compare;
                        }
                        compare = compare + 1;
                    }
                    if(abs(input[1] - listOfPoints[compare][1]) < 0.0001)
                    {
                        if(abs(input[0] - listOfPoints[compare][0]) < 0.0001)
                        {
                            return(compare);
                        }
                    }
                }
            }
        }
        if(dim == 3)
        {
            while((compare < numOfPoints) && abs(input[0] - listOfPoints[compare][0]) < 0.0001)
            {
                if(largest - smallest == 1)
                {
                    return(-1);
                }
                if(input[1] > listOfPoints[compare][1])
                {
                    if(compare < largest)
                    {
                        largest = compare;
                    }
                    compare = (largest - smallest)/2 + smallest;
                }
                if((abs(input[0] - listOfPoints[compare][0]) < 0.0001) && (compare < numOfPoints))
                {
                    while((compare < numOfPoints) && abs(input[0] - listOfPoints[compare][0]) < 0.0001 && abs(input[1] - listOfPoints[compare][1]) < 0.0001)
                    {
                        if(largest - smallest == 1)
                        {
                            return(-1);
                        }
                        if((int)(input[2] > listOfPoints[compare][2]))
                        {
                            if(compare < largest)
                            {
                                largest = compare;
                            }
                            compare = compare - 1;
                        }
                        if(abs(input[2] - listOfPoints[compare][2]) < 0.0001)
                        {
                            if(abs(input[1] - listOfPoints[compare][1]) < 0.0001)
                            {
                                if(abs(input[0] - listOfPoints[compare][0]) < 0.0001)
                                {
                                    return(compare);
                                }
                            }
                        }
                        if(input[2] < listOfPoints[compare][2])
                        {
                            if(compare > smallest)
                            {
                                smallest = compare;
                            }
                            compare = compare + 1;
                        }
                        if(abs(input[2] - listOfPoints[compare][2]) < 0.0001)
                        {
                            if(abs(input[1] - listOfPoints[compare][1]) < 0.0001)
                            {
                                if(abs(input[0] - listOfPoints[compare][0]) < 0.0001)
                                {
                                    return(compare);
                                }
                            }
                        }
                    }
                }
                if(input[1] < listOfPoints[compare][1])
                {
                    if(compare > smallest)
                    {
                        smallest = compare;
                    }
                    compare = (largest - smallest)/2 + smallest;
                }
            }
        }
    }
    if(largest - smallest == 1)
    {
        return(-1);
    }
    return(compare);
}

void removeDuplicates(std::vector <int> & stack)
{
	//std::vector <int> unique;
	//unique.push_back(stack[0]);
	if(stack.size() == 1)
	{
		return;
	}
	//for(int x = 0; x < stack.size(); x++)
	//{
		//cout << stack[x] << " ";
	//}
	//cout << "\n";
	for(int x = 0; x < stack.size(); x++)
	{
		//cout << stack[0] << "\n";
		//unique.push_back(stack[x]);
		for(int y = 0; y < stack.size(); y++)
		{
			if((x < stack.size()) && (y < stack.size()))
			{
				if(x != y)
				{
					if(stack[x] == stack[y])
					{
						//cout << "Thing " << stack[y] << "\n";
						stack.erase(stack.begin() + y);
					}
				}
				
			}
		}
	}
	//for(int x = 0; x < stack.size(); x++)
	//{
		//cout << stack[x] << " ";
	//}
	//cout << "\n";
	//cout << "Done this round\n";
}

void ASRGsingleCluster2D(int start, std::vector <int> * cluster, int currentCluster, int dim, std::vector <std::vector<double> > * listOfPoints, double ** dimensionVariables, struct imageStats imstat, std::vector <std::vector <std::vector <int> > > * dataGraph, int numOfPoints)
{
    std::vector <int> stack;
    std::vector <int> workingIndex;
    std::vector <int> up;
    std::vector <int> down;
    std::vector <int> left;
    std::vector <int> right;
    std::vector <int> upright;
    std::vector <int> downright;
    std::vector <int> upleft;
    std::vector <int> downleft;
    std::vector <int> back;
    std::vector <int> front;
    std::vector <int> backup;
    std::vector <int> backdown;
    std::vector <int> backleft;
    std::vector <int> backright;
    std::vector <int> backupright;
    std::vector <int> backdownright;
    std::vector <int> backupleft;
    std::vector <int> backdownleft;
    std::vector <int> frontup;
    std::vector <int> frontdown;
    std::vector <int> frontleft;
    std::vector <int> frontright;
    std::vector <int> frontupright;
    std::vector <int> frontdownright;
    std::vector <int> frontupleft;
    std::vector <int> frontdownleft;
    stack.push_back(start);
    int working;
    int zindex = 0;
    //from here, pop off the stack, look at adjacent, pop onto stack if adjacent is a 1
    working = (*stack.rbegin());
    //cout << working << "\n";
    //cout << "STARTING\n";
    int counter = 0;
    while(stack.size() > 0)
    {
    	//if(stack.size() > 1000)
    	//{
    		//if(counter < 10)
    		//{
    			//cout << "RemovingDuplicates " << counter << "\n";
    			//removeDuplicates(stack);
    		//}
    		//if(counter > 10 & counter%20 == 0)
    		//{
    			//cout << "RemovingDuplicates " << counter << "\n";
    			//removeDuplicates(stack);
    		//}
    		//counter++;
    	//}
    	//removeDuplicates(stack);
        //cout << stack.size() << "\n";
        //cout << (*cluster)[working] << "\n";
        if((*cluster)[working] > 0)
        {
            stack.pop_back();
        }
        if(working > -1 && (*cluster)[working] <= 0)
        {
            (*cluster)[working] = currentCluster;
            workingIndex = convertToIndex((*listOfPoints)[working], dimensionVariables, imstat, dim);
            //cout << workingIndex.size() << "\n";
//            if(dim == 3)
//            {
//                if(workingIndex[2] > 1)
//                {
//                    for(int x = workingIndex[1] - 1; x < workingIndex[1] + 2; x++)
//                    {
//                        for(int y = workingIndex[0] - 1; y < workingIndex[0] + 2; y++)
//                        {
//                            if(x > 0 && x < imstat.yintervals)
//                            {
//                                if(y > 0 && y < imstat.xintervals)
//                                {
//                                    cout << (*dataGraph)[workingIndex[2] -1][x][y];
//                                }
//                            }
//                        }
//                        cout << "\n";
//                    }
//                }
//                cout << "\n";
//                for(int x = workingIndex[1] - 1; x < workingIndex[1] + 2; x++)
//                {
//                    for(int y = workingIndex[0] - 1; y < workingIndex[0] + 2; y++)
//                    {
//                        if(x > 0 && x < imstat.yintervals)
//                        {
//                            if(y > 0 && y < imstat.xintervals)
//                            {
//                                cout << (*dataGraph)[workingIndex[2]][x][y];
//                            }
//                        }
//                    }
//                    cout << "\n";
//                }
//                cout << "\n";
//                if(workingIndex[2] < imstat.zintervals)
//                {
//                    for(int x = workingIndex[1] - 1; x < workingIndex[1] + 2; x++)
//                    {
//                        for(int y = workingIndex[0] - 1; y < workingIndex[0] + 2; y++)
//                        {
//                            if(x > 0 && x < imstat.yintervals)
//                            {
//                                if(y > 0 && y < imstat.xintervals)
//                                {
//                                    cout << (*dataGraph)[workingIndex[2] + 1][x][y];
//                                }
//                            }
//                        }
//                        cout << "\n";
//                    }
//                }
//                cout << "\n";
//            }
//            cout << "\n";
            if(workingIndex[0] > 0)
            {
                left.push_back(workingIndex[0] - 1);
                left.push_back(workingIndex[1]);
                if(dim > 2)
                {
                    left.push_back(workingIndex[2]);
                }
                if(dim > 3)
                {
                    left.push_back(workingIndex[3]);
                }
                if(workingIndex[1] > 0)
                {
                    upleft.push_back(workingIndex[0] - 1);
                    upleft.push_back(workingIndex[1] - 1);
                    if(dim > 2)
                    {
                        upleft.push_back(workingIndex[2]);
                    }
                    if(dim > 3)
                    {
                        upleft.push_back(workingIndex[3]);
                    }
                }
                //cout << upleft.size() << "\n";
                if((imstat.yintervals - workingIndex[1]) > 1)
                {
                    downleft.push_back(workingIndex[0] - 1);
                    downleft.push_back(workingIndex[1] + 1);
                    if(dim > 2)
                    {
                        downleft.push_back(workingIndex[2]);
                    }
                    if(dim > 3)
                    {
                        downleft.push_back(workingIndex[3]);
                    }
                }
            }
            if((imstat.xintervals - workingIndex[0]) > 1)
            {
                //cout << "There is a right\n";
                right.push_back(workingIndex[0] + 1);
                right.push_back(workingIndex[1]);
                if(dim > 2)
                {
                    right.push_back(workingIndex[2]);
                }
                if(dim > 3)
                {
                    right.push_back(workingIndex[3]);
                }
                if(workingIndex[1] > 0)
                {
                    //cout << "There is an up\n";
                    upright.push_back(workingIndex[0] + 1);
                    upright.push_back(workingIndex[1] - 1);
                    if(dim > 2)
                    {
                        upright.push_back(workingIndex[2]);
                    }
                    if(dim > 3)
                    {
                        upright.push_back(workingIndex[3]);
                    }
                }
                if((imstat.yintervals - workingIndex[1]) > 1)
                {
                    //cout << "There is a down\n";
                    downright.push_back(workingIndex[0] + 1);
                    downright.push_back(workingIndex[1]+1);
                    if(dim > 2)
                    {
                        downright.push_back(workingIndex[2]);
                    }
                    if(dim > 3)
                    {
                        downright.push_back(workingIndex[3]);
                    }
                }
            }
            if(dim > 2)
            {
                if(workingIndex[2] > 0)
                {
                    back.push_back(workingIndex[0]);
                    back.push_back(workingIndex[1]);
                    back.push_back(workingIndex[2] - 1);
                    if(workingIndex[0] > 0)
                    {
                        backleft.push_back(workingIndex[0] - 1);
                        backleft.push_back(workingIndex[1]);
                        backleft.push_back(workingIndex[2] - 1);
                        if(workingIndex[1] > 0)
                        {
                            backupleft.push_back(workingIndex[0] - 1);
                            backupleft.push_back(workingIndex[1] - 1);
                            backupleft.push_back(workingIndex[2] - 1);
                        }
                        if((imstat.yintervals - workingIndex[1]) > 1)
                        {
                            backdownleft.push_back(workingIndex[0] - 1);
                            backdownleft.push_back(workingIndex[1] + 1);
                            backdownleft.push_back(workingIndex[2] - 1);
                        }
                    }
                    if((imstat.xintervals - workingIndex[0]) > 1)
                    {
                        backright.push_back(workingIndex[0] + 1);
                        backright.push_back(workingIndex[1]);
                        backright.push_back(workingIndex[2] - 1);
                        if(workingIndex[1] > 0)
                        {
                            backupright.push_back(workingIndex[0] + 1);
                            backupright.push_back(workingIndex[1] - 1);
                            backupright.push_back(workingIndex[2] - 1);
                        }
                        if((imstat.yintervals - workingIndex[1]) > 1)
                        {
                            backdownright.push_back(workingIndex[0] + 1);
                            backdownright.push_back(workingIndex[1]+1);
                            backdownright.push_back(workingIndex[2] - 1);
                        }
                    }
                    if(workingIndex[1] > 0)
                    {
                        backup.push_back(workingIndex[0]);
                        backup.push_back(workingIndex[1] - 1);
                        backup.push_back(workingIndex[2] - 1);
                    }
                    if((imstat.yintervals - workingIndex[1]) > 1)
                    {
                        backdown.push_back(workingIndex[0]);
                        backdown.push_back(workingIndex[1]+1);
                        backdown.push_back(workingIndex[2] - 1);
                    }
                }
                if((imstat.zintervals - workingIndex[2]) > 1)
                {
                    front.push_back(workingIndex[0]);
                    front.push_back(workingIndex[1]);
                    front.push_back(workingIndex[2] + 1);
                    if(workingIndex[0] > 0)
                    {
                        frontleft.push_back(workingIndex[0] - 1);
                        frontleft.push_back(workingIndex[1]);
                        frontleft.push_back(workingIndex[2] + 1);
                        if(workingIndex[1] > 0)
                        {
                            frontupleft.push_back(workingIndex[0] - 1);
                            frontupleft.push_back(workingIndex[1] - 1);
                            frontupleft.push_back(workingIndex[2] + 1);
                        }
                        if((imstat.yintervals - workingIndex[1]) > 1)
                        {
                            frontdownleft.push_back(workingIndex[0] - 1);
                            frontdownleft.push_back(workingIndex[1] + 1);
                            frontdownleft.push_back(workingIndex[2] + 1);
                        }
                    }
                    if((imstat.xintervals - workingIndex[0]) > 1)
                    {
                        frontright.push_back(workingIndex[0] + 1);
                        frontright.push_back(workingIndex[1]);
                        frontright.push_back(workingIndex[2] + 1);
                        if(workingIndex[1] > 0)
                        {
                            frontupright.push_back(workingIndex[0] + 1);
                            frontupright.push_back(workingIndex[1] - 1);
                            frontupright.push_back(workingIndex[2] + 1);
                        }
                        if((imstat.yintervals - workingIndex[1]) > 1)
                        {
                            frontdownright.push_back(workingIndex[0] + 1);
                            frontdownright.push_back(workingIndex[1]+1);
                            frontdownright.push_back(workingIndex[2] + 1);
                        }
                    }
                    if(workingIndex[1] > 0)
                    {
                        frontup.push_back(workingIndex[0]);
                        frontup.push_back(workingIndex[1] - 1);
                        frontup.push_back(workingIndex[2] + 1);
                    }
                    if((imstat.yintervals - workingIndex[1]) > 1)
                    {
                        frontdown.push_back(workingIndex[0]);
                        frontdown.push_back(workingIndex[1]+1);
                        frontdown.push_back(workingIndex[2] + 1);
                    }
                }
            }
//            if(workingIndex[0] > 0)
//            {
//                left.push_back(workingIndex[0] - 1);
//                left.push_back(workingIndex[1]);
//                if(dim > 2)
//                {
//                    left.push_back(workingIndex[2]);
//                }
//                if(dim > 3)
//                {
//                    left.push_back(workingIndex[3]);
//                }
//                if(workingIndex[1] > 0)
//                {
//                    upleft.push_back(workingIndex[0] - 1);
//                    upleft.push_back(workingIndex[1] - 1);
//                    if(dim > 2)
//                    {
//                        upleft.push_back(workingIndex[2]);
//                    }
//                    if(dim > 3)
//                    {
//                        upleft.push_back(workingIndex[3]);
//                    }
//                }
//                if((imstat.yintervals - workingIndex[1]) > 1)
//                {
//                    downleft.push_back(workingIndex[0] - 1);
//                    downleft.push_back(workingIndex[1]+1);
//                    if(dim > 2)
//                    {
//                        downleft.push_back(workingIndex[2]);
//                    }
//                    if(dim > 3)
//                    {
//                        downleft.push_back(workingIndex[3]);
//                    }
//                }
//            }
//            if((imstat.xintervals - workingIndex[0]) > 1)
//            {
//                //cout << "There is a right\n";
//                right.push_back(workingIndex[0] + 1);
//                right.push_back(workingIndex[1]);
//                if(dim > 2)
//                {
//                    right.push_back(workingIndex[2]);
//                }
//                if(dim > 3)
//                {
//                    right.push_back(workingIndex[3]);
//                }
//                if(workingIndex[1] > 0)
//                {
//                    //cout << "There is an up\n";
//                    upright.push_back(workingIndex[0] + 1);
//                    upright.push_back(workingIndex[1] - 1);
//                    if(dim > 2)
//                    {
//                        upright.push_back(workingIndex[2]);
//                    }
//                    if(dim > 3)
//                    {
//                        upright.push_back(workingIndex[3]);
//                    }
//                }
//                if((imstat.yintervals - workingIndex[1]) > 1)
//                {
//                    //cout << "There is a down\n";
//                    downright.push_back(workingIndex[0] + 1);
//                    downright.push_back(workingIndex[1]+1);
//                    if(dim > 2)
//                    {
//                        downright.push_back(workingIndex[2]);
//                    }
//                    if(dim > 3)
//                    {
//                        downright.push_back(workingIndex[3]);
//                    }
//                }
           // }
            if(workingIndex[1] > 0)
            {
                //cout << "There is an up\n";
                up.push_back(workingIndex[0]);
                up.push_back(workingIndex[1] - 1);
                if(dim > 2)
                {
                    up.push_back(workingIndex[2]);
                }
                if(dim > 3)
                {
                    up.push_back(workingIndex[3]);
                }
            }
            if((imstat.yintervals - workingIndex[1]) > 1)
            {
                //cout << "There is a down\n";
                down.push_back(workingIndex[0]);
                down.push_back(workingIndex[1]+1);
                if(dim > 2)
                {
                    down.push_back(workingIndex[2]);
                }
                if(dim > 3)
                {
                    down.push_back(workingIndex[3]);
                }
            }
            stack.pop_back();
            //cout << "done up\n";
            if(dim > 2)
            {
                if(back.size() > 0)
                {
                    if((*dataGraph)[back[2]][back[1]][back[0]] > 0)
                    {
                        //start = findListIndex(convertToFrequancy(back, dimensionVariables, imstat, dim), (*listOfPoints), numOfPoints, dim);
                        start = (*dataGraph)[back[2]][back[1]][back[0]] - 1;
                        if(start > -1 && (*cluster)[start] == 0)
                        {
                            stack.push_back(start);
                            (*cluster)[start] = -1;
                        }
                    }
                    back.clear();
                }
                if(backleft.size() > 0)
                {
                    if((*dataGraph)[backleft[2]][backleft[1]][backleft[0]] > 0)
                    {
                        //start = findListIndex(convertToFrequancy(backleft, dimensionVariables, imstat, dim), (*listOfPoints), numOfPoints, dim);
                        start = (*dataGraph)[backleft[2]][backleft[1]][backleft[0]] - 1;
                        if(start > -1 && (*cluster)[start] == 0)
                        {
                           stack.push_back(start);
                           (*cluster)[start] = -1;
                        }
                    }
                    backleft.clear();
                }
                if(backright.size() > 0)
                {
                    if((*dataGraph)[backright[2]][backright[1]][backright[0]] > 0)
                    {
                        //start = findListIndex(convertToFrequancy(backright, dimensionVariables, imstat, dim), (*listOfPoints), numOfPoints, dim);
                        start = (*dataGraph)[backright[2]][backright[1]][backright[0]] - 1;
                        if(start > -1 && (*cluster)[start] == 0)
                        {
                            stack.push_back(start);
                            (*cluster)[start] = -1;
                        }
                    }
                    backright.clear();
                }
                if(backupleft.size() > 0)
                {
                    if((*dataGraph)[backupleft[2]][backupleft[1]][backupleft[0]] > 0)
                    {
                        //start = findListIndex(convertToFrequancy(backupleft, dimensionVariables, imstat, dim), (*listOfPoints), numOfPoints, dim);
                        start = (*dataGraph)[backupleft[2]][backupleft[1]][backupleft[0]] - 1;
                        if(start > -1 && (*cluster)[start] == 0)
                        {
                            stack.push_back(start);
                            (*cluster)[start] = -1;
                        }
                    }
                    backupleft.clear();
                }
                if(backupright.size() > 0)
                {
                    if((*dataGraph)[backupright[2]][backupright[1]][backupright[0]] > 0)
                    {
                        //start = findListIndex(convertToFrequancy(backupright, dimensionVariables, imstat, dim), (*listOfPoints), numOfPoints, dim);
                        start = (*dataGraph)[backupright[2]][backupright[1]][backupright[0]] - 1;
                        if(start > -1 && (*cluster)[start] == 0)
                        {
                            stack.push_back(start);
                            (*cluster)[start] = -1;
                        }
                    }
                    backupright.clear();
                }
                if(backdownleft.size() > 0)
                {
                    if((*dataGraph)[backdownleft[2]][backdownleft[1]][backdownleft[0]] > 0)
                    {
                        //start = findListIndex(convertToFrequancy(backdownleft, dimensionVariables, imstat, dim), (*listOfPoints), numOfPoints, dim);
                        start = (*dataGraph)[backdownleft[2]][backdownleft[1]][backdownleft[0]] - 1;
                        if(start > -1 && (*cluster)[start] == 0)
                        {
                            stack.push_back(start);
                            (*cluster)[start] = -1;
                        }
                    }
                    backdownleft.clear();
                }
                if(backdownright.size() > 0)
                {
                    if((*dataGraph)[backdownright[2]][backdownright[1]][backdownright[0]] > 0)
                    {
                        //start = findListIndex(convertToFrequancy(backdownright, dimensionVariables, imstat, dim), (*listOfPoints), numOfPoints, dim);
                        start = (*dataGraph)[backdownright[2]][backdownright[1]][backdownright[0]] - 1;
                        if(start > -1 && (*cluster)[start] == 0)
                        {
                            stack.push_back(start);
                            (*cluster)[start] = -1;
                        }
                    }
                    backdownright.clear();
                }
                if(backup.size() > 0)
                {
                    if((*dataGraph)[backup[2]][backup[1]][backup[0]] > 0)
                    {
                        //start = findListIndex(convertToFrequancy(backup, dimensionVariables, imstat, dim), (*listOfPoints), numOfPoints, dim);
                        start = (*dataGraph)[backup[2]][backup[1]][backup[0]] - 1;
                        if(start > -1 && (*cluster)[start] == 0)
                        {
                            stack.push_back(start);
                            (*cluster)[start] = -1;
                        }
                    }
                    backup.clear();
                }
                if(backdown.size() > 0)
                {
                    if((*dataGraph)[backdown[2]][backdown[1]][backdown[0]] > 0)
                    {
                        //start = findListIndex(convertToFrequancy(backdown, dimensionVariables, imstat, dim), (*listOfPoints), numOfPoints, dim);
                        start = (*dataGraph)[backdown[2]][backdown[1]][backdown[0]] - 1;
                        if(start > -1 && (*cluster)[start] == 0)
                        {
                            stack.push_back(start);
                            (*cluster)[start] = -1;
                        }
                    }
                    backdown.clear();
                }
                if(front.size() > 0)
                {
                    if((*dataGraph)[front[2]][front[1]][front[0]] > 0)
                    {
                        //start = findListIndex(convertToFrequancy(front, dimensionVariables, imstat, dim), (*listOfPoints), numOfPoints, dim);
                        start = (*dataGraph)[front[2]][front[1]][front[0]] - 1;
                        if(start > -1 && (*cluster)[start] == 0)
                        {
                            stack.push_back(start);
                            (*cluster)[start] = -1;
                        }
                    }
                    front.clear();
                }
                if(frontleft.size() > 0)
                {
                    if((*dataGraph)[frontleft[2]][frontleft[1]][frontleft[0]] > 0)
                    {
                        //start = findListIndex(convertToFrequancy(frontleft, dimensionVariables, imstat, dim), (*listOfPoints), numOfPoints, dim);
                        start = (*dataGraph)[frontleft[2]][frontleft[1]][frontleft[0]] - 1;
                        if(start > -1 && (*cluster)[start] == 0)
                        {
                            stack.push_back(start);
                            (*cluster)[start] = -1;
                        }
                    }
                    frontleft.clear();
                }
                if(frontright.size() > 0)
                {
                    if((*dataGraph)[frontright[2]][frontright[1]][frontright[0]] > 0)
                    {
                        //start = findListIndex(convertToFrequancy(frontright, dimensionVariables, imstat, dim), (*listOfPoints), numOfPoints, dim);
                        start = (*dataGraph)[frontright[2]][frontright[1]][frontright[0]] - 1;
                        if(start > -1 && (*cluster)[start] == 0)
                        {
                            stack.push_back(start);
                            (*cluster)[start] = -1;
                        }
                    }
                    frontright.clear();
                }
                if(frontupleft.size() > 0)
                {
                    if((*dataGraph)[frontupleft[2]][frontupleft[1]][frontupleft[0]] > 0)
                    {
                        //start = findListIndex(convertToFrequancy(frontupleft, dimensionVariables, imstat, dim), (*listOfPoints), numOfPoints, dim);
                        start = (*dataGraph)[frontupleft[2]][frontupleft[1]][frontupleft[0]] - 1;
                        if(start > -1 && (*cluster)[start] == 0)
                        {
                            stack.push_back(start);
                            (*cluster)[start] = -1;
                        }
                    }
                    frontupleft.clear();
                }
                if(frontupright.size() > 0)
                {
                    if((*dataGraph)[frontupright[2]][frontupright[1]][frontupright[0]] > 0)
                    {
                        //start = findListIndex(convertToFrequancy(frontupright, dimensionVariables, imstat, dim), (*listOfPoints), numOfPoints, dim);
                        start = (*dataGraph)[frontupright[2]][frontupright[1]][frontupright[0]] - 1;
                        if(start > -1 && (*cluster)[start] == 0)
                        {
                            stack.push_back(start);
                            (*cluster)[start] = -1;
                        }
                    }
                    frontupright.clear();
                }
                if(frontdownleft.size() > 0)
                {
                    if((*dataGraph)[frontdownleft[2]][frontdownleft[1]][frontdownleft[0]] > 0)
                    {
                        //start = findListIndex(convertToFrequancy(frontdownleft, dimensionVariables, imstat, dim), (*listOfPoints), numOfPoints, dim);
                        start = (*dataGraph)[frontdownleft[2]][frontdownleft[1]][frontdownleft[0]] - 1;
                        if(start > -1 && (*cluster)[start] == 0)
                        {
                            stack.push_back(start);
                            (*cluster)[start] = -1;
                        }
                    }
                    frontdownleft.clear();
                }
                if(frontdownright.size() > 0)
                {
                    if((*dataGraph)[frontdownright[2]][frontdownright[1]][frontdownright[0]] > 0)
                    {
                        //start = findListIndex(convertToFrequancy(frontdownright, dimensionVariables, imstat, dim), (*listOfPoints), numOfPoints, dim);
                        start = (*dataGraph)[frontdownright[2]][frontdownright[1]][frontdownright[0]] - 1;
                        if(start > -1 && (*cluster)[start] == 0)
                        {
                            stack.push_back(start);
                            (*cluster)[start] = -1;
                        }
                    }
                    frontdownright.clear();
                }
                if(frontup.size() > 0)
                {
                    if((*dataGraph)[frontup[2]][frontup[1]][frontup[0]] > 0)
                    {
                        //start = findListIndex(convertToFrequancy(frontup, dimensionVariables, imstat, dim), (*listOfPoints), numOfPoints, dim);
                        start = (*dataGraph)[frontup[2]][frontup[1]][frontup[0]] - 1;
                        if(start > -1 && (*cluster)[start] == 0)
                        {
                            stack.push_back(start);
                            (*cluster)[start] = -1;
                        }
                    }
                    frontup.clear();
                }
                if(frontdown.size() > 0)
                {
                    if((*dataGraph)[frontdown[2]][frontdown[1]][frontdown[0]] > 0)
                    {
                        //start = findListIndex(convertToFrequancy(frontdown, dimensionVariables, imstat, dim), (*listOfPoints), numOfPoints, dim);
                        start = (*dataGraph)[frontdown[2]][frontdown[1]][frontdown[0]] - 1;
                        if(start > -1 && (*cluster)[start] == 0)
                        {
                            stack.push_back(start);
                            (*cluster)[start] = -1;
                        }
                    }
                    frontdown.clear();
                }
            }
            if(upright.size() > 0)
            {
                if(dim == 2)
                {
                    zindex = 0;
                    upright.push_back(0.0);
                }
                if((*dataGraph)[upright[2]][upright[1]][upright[0]] > 0)
                {
                    //start = findListIndex(convertToFrequancy(upright, dimensionVariables, imstat, dim), (*listOfPoints), numOfPoints, dim);
                    start = (*dataGraph)[upright[2]][upright[1]][upright[0]] - 1;
                    if(start > -1 && (*cluster)[start] == 0)
                    {
                        stack.push_back(start);
                        (*cluster)[start] = -1;
                    }
                }
                upright.clear();
            }
            if(up.size() > 0)
            {
                if(dim == 2)
                {
                    zindex = 0;
                    up.push_back(0.0);
                }
                if((*dataGraph)[up[2]][up[1]][up[0]] > 0)
                {
                    //start = findListIndex(convertToFrequancy(up, dimensionVariables, imstat, dim), (*listOfPoints), numOfPoints, dim);
                    start = (*dataGraph)[up[2]][up[1]][up[0]] - 1;
                    if(start > -1 && (*cluster)[start] == 0)
                    {
                        stack.push_back(start);
                        (*cluster)[start] = -1;
                    }
                }
                up.clear();
            }
            //cout << "done upright\n";
            //cout << upleft.size() << "\n";
            if(upleft.size() > 0)
            {
                if(dim == 2)
                {
                    zindex = 0;
                    upleft.push_back(0.0);
                }
                if((*dataGraph)[upleft[2]][upleft[1]][upleft[0]] > 0)
                {
                    //start = findListIndex(convertToFrequancy(upleft, dimensionVariables, imstat, dim), (*listOfPoints), numOfPoints, dim);
                    start = (*dataGraph)[upleft[2]][upleft[1]][upleft[0]] - 1;
                    if(start > -1 && (*cluster)[start] == 0)
                    {
                        stack.push_back(start);
                        (*cluster)[start] = -1;
                    }
                }
                upleft.clear();
            }
            if(left.size() > 0)
            {
                if(dim == 2)
                {
                    zindex = 0;
                    left.push_back(0.0);
                }
                if((*dataGraph)[left[2]][left[1]][left[0]] > 0)
                {
                    //start = findListIndex(convertToFrequancy(left, dimensionVariables, imstat, dim), (*listOfPoints), numOfPoints, dim);
                    start = (*dataGraph)[left[2]][left[1]][left[0]] - 1;
                    if(start > -1 && (*cluster)[start] == 0)
                    {
                        stack.push_back(start);
                        (*cluster)[start] = -1;
                    }
                }
                left.clear();
            }
            if(downleft.size() > 0)
            {
                if(dim == 2)
                {
                    zindex = 0;
                    downleft.push_back(0.0);
                }
                if((*dataGraph)[downleft[2]][downleft[1]][downleft[0]] > 0)
                {
                    //start = findListIndex(convertToFrequancy(downleft, dimensionVariables, imstat, dim), (*listOfPoints), numOfPoints, dim);
                    start = (*dataGraph)[downleft[2]][downleft[1]][downleft[0]] - 1;
                    if(start > -1 && (*cluster)[start] == 0)
                    {
                        stack.push_back(start);
                        (*cluster)[start] = -1;
                    }
                }
                downleft.clear();
            }
            if(down.size() > 0)
            {
                if(dim == 2)
                {
                    zindex = 0;
                    down.push_back(0.0);
                }
                if((*dataGraph)[down[2]][down[1]][down[0]] > 0)
                {
                    //start = findListIndex(convertToFrequancy(down, dimensionVariables, imstat, dim), (*listOfPoints), numOfPoints, dim);
                    start = (*dataGraph)[down[2]][down[1]][down[0]] - 1;
                    if(start > -1 && (*cluster)[start] == 0)
                    {
                        stack.push_back(start);
                        (*cluster)[start] = -1;
                    }
                }
                down.clear();
            }
            if(downright.size() > 0)
            {
                if(dim == 2)
                {
                    zindex = 0;
                    downright.push_back(0.0);
                }
                if(dim > 2)
                {
                    zindex = downright[2];
                }
                if((*dataGraph)[downright[2]][downright[1]][downright[0]] > 0)
                {
                    //start = findListIndex(convertToFrequancy(downright, dimensionVariables, imstat, dim), (*listOfPoints), numOfPoints, dim);
                    start = (*dataGraph)[downright[2]][downright[1]][downright[0]] - 1;
                    if(start > -1 && (*cluster)[start] == 0)
                    {
                        stack.push_back(start);
                        (*cluster)[start] = -1;
                    }
                }
                downright.clear();
            }
            //cout << "done left\n";
            if(right.size() > 0)
            {
                if(dim == 2)
                {
                    zindex = 0;
                    right.push_back(0.0);
                }
                if((*dataGraph)[right[2]][right[1]][right[0]] > 0)
                {
                    //start = findListIndex(convertToFrequancy(right, dimensionVariables, imstat, dim), (*listOfPoints), numOfPoints, dim);
                    start = (*dataGraph)[right[2]][right[1]][right[0]] - 1;
                    if(start > -1 && (*cluster)[start] == 0)
                    {
                        stack.push_back(start);
                        (*cluster)[start] = -1;
                    }
                }
                right.clear();
            }
        }
        if(stack.size()>0)
        {
            working = (*stack.rbegin());
        }
        else
        {
            working = (*stack.end());
        }
    }
}

void ASRG(std::list <std::vector <double> > Seeds, std::vector <std::vector<double> >  listOfPoints, double ** dimensionVariables, struct imageStats imstat, std::vector <std::vector <std::vector <int> > > dataGraph, int numOfPoints, int dim)
{
    //keep track of which points you already looked at using an array that contains which cluster an index belongs to, and a stack to process the next round of adjacent pixels
    cout << "In ASRG\n";
    std::vector <int> cluster;
    std::vector <int> stack;
    std::list <std::vector<double> > :: iterator it;
    int currentCluster = 0;
    int start = -1;
    for(int x = 0; x < numOfPoints; x++)
    {
        cluster.push_back(0);
    }
    if(Seeds.size() > 0)
    {
        currentCluster++;
        it = Seeds.begin();
        start = findListIndex((*it), listOfPoints, numOfPoints, dim);
        ASRGsingleCluster2D(start, &cluster, currentCluster, dim, &listOfPoints, dimensionVariables, imstat, &dataGraph, numOfPoints);
        it++;
        while(it != Seeds.end())
        {
            start = findListIndex((*it), listOfPoints, numOfPoints, dim);
            if(cluster[start] < 1)
            {
                currentCluster++;
                ASRGsingleCluster2D(start, &cluster, currentCluster, dim, &listOfPoints, dimensionVariables, imstat, &dataGraph, numOfPoints);
            }
            it++;
        }
    }
    //cout <<"made it to here\n";
    string output_string = "";
    for(int y = 0; y < numOfPoints; y++)
    {
        if(cluster[y] == 0)
        {
        	output_string = "Index: " + to_string(y) + " of " + to_string(numOfPoints -1) + " " + to_string((double)(100.0*(y/(double)(numOfPoints-1)))) + " Percent Complete\r";
            //cout << "Index: " << y << " of "<< numOfPoints -1 << "\n";
            //cout <<  "\r" << (double)(100.0*(y/(double)(numOfPoints-1))) << " Percent Complete\n";
            cout << output_string;
            start = y;
            currentCluster++;
            ASRGsingleCluster2D(start, &cluster, currentCluster, dim, &listOfPoints, dimensionVariables, imstat, &dataGraph, numOfPoints);
            //cout << "Cluster " << currentCluster << " Complete\n";
            std::cout.flush();
        }
    }
    cout << "\n";
    for(int y = 0; y < numOfPoints; y++)
    {
        if(cluster[y] == 0)
        {
            start = y;
            currentCluster++;
            ASRGsingleCluster2D(start, &cluster, currentCluster, dim, &listOfPoints, dimensionVariables, imstat, &dataGraph, numOfPoints);
            cout << "Cluster " << currentCluster << " Complete\n";
        }
    }
    cout << "Final number of clusters: " << currentCluster << "\n";
    writeFile2(cluster, listOfPoints, currentCluster, numOfPoints, dim, 1);
} 

bool sortcol(const vector<int>& v1, const vector<int>& v2)
{
    return v1[0] < v2[0];
}

std::vector <int> sortList(std::vector <int> finalList, int numOfClusters, int numOfPoints)
{
	std::vector <std::vector <int> > lengths;
	std::vector <int> sortedFinalList;
	std::vector <int> temp;
	temp.push_back(0);
	temp.push_back(0);
	for(int x = 0; x < numOfClusters; x++)
	{
		temp[1] = (x+1);
		lengths.push_back(temp);
	}
	for(int x =0; x < numOfPoints; x++)
	{
		lengths[finalList[x]-1][0] = lengths[finalList[x]-1][0] + 1;
	}
	std::sort(lengths.begin(), lengths.end(),
	          sortcol);

	int current = 0;
	for(int x =0; x < numOfPoints; x++)
	{
		current = -1;
		for(int y =0; y < numOfClusters; y++)
		{
			if(finalList[x] == lengths[y][1])
			{
				current = y+1;
			}
		}
		if(current == -1)
		{
			cout << "ERROR!! cluster not found!\n";
		}
		sortedFinalList.push_back(current);
	}	
	return(sortedFinalList);
}

void writeFile2(std::vector <int> finalList, std::vector <std::vector<double> > listOfPoints, int numOfClusters, int numOfPoints, int dim, int flag)
{
	fstream aus;  //output file
	aus.precision(7);
    std::vector<double> holder;
    //std::vector<int>:: iterator it3;
    std::list <std::vector <int> >:: iterator it2;
    //int counter = 1;
	std::vector <int> finalList2 = sortList(finalList, numOfClusters, numOfPoints);
	int count = 0;
	string outfile = "";
	std::vector <std::vector <int> > lengths;
	std::vector <int> temp;
	temp.push_back(0);
	temp.push_back(0);
	for(int x = 0; x < numOfClusters; x++)
	{
		temp[1] = (x+1);
		lengths.push_back(temp);
	}
	for(int x =0; x < numOfPoints; x++)
	{
		//cout << x << "\n";
		lengths[finalList2[x]-1][0] = lengths[finalList2[x]-1][0] + 1;
	}
	
	/*for(int x =0; x < numOfClusters; x++)
	{
		cout << lengths[x][0] << " " << lengths[x][1] << "\n";
	}*/
	
	
    for(int x = 1; x <= numOfClusters; x++)
    {
		if(flag == 1)
		{
			outfile = "../Deconvolution/Datasets/outASRG" + to_string(x) + ".txt";
		}
        else
		{
			outfile = "../Deconvolution/Datasets/out" + to_string(x) + ".txt";
		}
        //fstream aus;  //output file
        //aus.open(outfile,ios::out);
        //aus.precision(7);
        string outstring;
        //double ind = 0.0;
        /*count = 0;
        for(int y = 0; y < numOfPoints; y++)
        {
            if(finalList2[y] == x)
            {
				count++;
			}
			//if(count > 0 & finalList2[y] != x)
			//{
				//break;
			//}
        }*/
        count = lengths[x-1][0];
        //cout << count << "\n";
        if(count > 10)
        {
        	aus.open(outfile,ios::out);
        	for(int y = 0; y < numOfPoints; y++)
        	{
            	if(finalList2[y] == x)
            	{
                	holder = listOfPoints[y];
                	if(dim == 2)
                	{
                    	outstring = to_string(holder.at(0)) + " " + to_string(holder.at(1)) + " " + to_string(holder.at(2)) + "\n";
                	}
                	if(dim == 3)
                	{
                    	outstring = to_string(holder.at(0)) + " " + to_string(holder.at(1)) + " " + to_string(holder.at(2)) + " " + to_string(holder.at(3)) + "\n";
                	}
                	if(dim == 4)
                	{
                    	outstring = to_string(holder.at(0)) + " " + to_string(holder.at(1)) + " " + to_string(holder.at(2)) + " " + to_string(holder.at(3)) + " " + to_string(holder.at(4)) + "\n";
                	}
                	aus << outstring;
            	}
        	}
        	aus.close();
        }
        /*for(int y = 0; y < numOfPoints; y++)
        {
            if(finalList2[y] == x)
            {
                holder = listOfPoints[y];
                if(dim == 2)
                {
                    outstring = to_string(holder.at(0)) + " " + to_string(holder.at(1)) + " " + to_string(holder.at(2)) + "\n";
                }
                if(dim == 3)
                {
                    outstring = to_string(holder.at(0)) + " " + to_string(holder.at(1)) + " " + to_string(holder.at(2)) + " " + to_string(holder.at(3)) + "\n";
                }
                if(dim == 4)
                {
                    outstring = to_string(holder.at(0)) + " " + to_string(holder.at(1)) + " " + to_string(holder.at(2)) + " " + to_string(holder.at(3)) + " " + to_string(holder.at(4)) + "\n";
                }
                aus << outstring;
            }
        }*/
        //aus.close();
    }
}

int main(int argc, char ** argv)
{
    clock_t time_req;
    time_req = clock();
    srand((int)time(NULL));
    std::vector <std::vector<double> >  listOfPoints;
    std::vector <std::vector<double> >  listOfPoints_temp;
    //ReadData(argc, argv, &listOfPoints, 1);
    vector < vector <int> > transistionList;
    
    ReadData_bin("Outlist.bin", &listOfPoints);
	
    pthread_mutex_init(&m, NULL);
    int dim = 0;
    dim = (int)listOfPoints.begin()->size() - 1;
	int numOfPoints = (int)listOfPoints.size();
	cout << "DIM: " << dim << " Length: " << numOfPoints << "\n";
    string runNum = argv[2];
    double ** dimensionVariables = NULL;
    double ** data; //= normalizeSet(numOfPoints, dim, &listOfPoints, &dimensionVariables);\

    //put a function here that does the seeded region Growing (ASRG) genetic algorithm. You can probably find a good one already made.
    //First set up the dataset as a 2D matrix use 0 for places that have no point, and a 1 for a location with a point
    //keep track of the indecies so you can revert to the points again. I made functions for this.
    //step 1 recreate data so its not normalized, and obtain the dimension variables. You need the dimension variables to convert to and from the matrix indecies
   data = (double**)malloc(numOfPoints*sizeof(double*));
   if(data != NULL)
   {
       for(int x = 0; x < numOfPoints; x++)
       {
           data[x] = (double*)malloc(dim*sizeof(double));
       }
   }
   std::vector <std::vector<double> > :: iterator it;
   std::vector<double>:: iterator it2;
   int count1 = 0;
   int count2 = 0;
   //find range of both dimensions
   for(it = listOfPoints.begin(); it != listOfPoints.end(); it++)
   {
       if(count1 < numOfPoints)
       {
           count2 = 0;
           for (it2 = it->begin(); it2 != it->end(); it2++)
           {
               if(count2 < dim)
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
	
    struct imageStats imstat; //this contains the variables that when conbined with dimension variables can convert frequency x and y coorindates to the matrix x and y coordinates and back again
    dimensionVariables = grabDimensions(dim, data, numOfPoints);
    std::vector <std::vector <std::vector <int> > > dataGraph;
    dataGraph = generateImageGraph(dimensionVariables, data, dim, numOfPoints, &imstat, &transistionList); //generate the image file suing 1s and 0s, 1s represent a data point, 0s represent a space
    //next collect the first point from each cluster produced by multi-k-means to use as seeds for the ASRG
    std::list <std::vector <double> > Seeds;
    //in this case we will skip the seeds. multiclustering takes too much memory
            //int counter = 1;
    //        for(int x = 0; x < numOfPoints; x++)
    //        {
    //            if(finalList[x] == counter)
    //            {
    //                Seeds.push_back(listOfPoints[x]);
    //                counter++;
    //            }
    //        }
            //good, so now we have the matrix, the seeds, the points, and conversion functions
            //next step is to write up the ASRG, the key is you only search the point directly adjacent to the current point, rather than linear search of comparisions of every point to every other point. you stop when there is a boarder of 0s around your cluster. You might eat up your other seeds in this process, so remember to keep track of which points are in which cluster. Maybe just an array.
    //int holder = findListIndex(listOfPoints[0], listOfPoints, numOfPoints, dim);
    //cout << holder << "\n";
    //holder = findListIndex(listOfPoints[28], listOfPoints, numOfPoints, dim);
    //cout << holder << "\n";
    //holder = findListIndex(listOfPoints[4000], listOfPoints, numOfPoints, dim);
    //cout << holder << "\n";
    //holder = findListIndex(listOfPoints[numOfPoints-1], listOfPoints, numOfPoints, dim);
    //cout << holder << "\n";
    ASRG(Seeds, listOfPoints, dimensionVariables, imstat, dataGraph, numOfPoints, dim);
    for(int x = 0; x < dim; x++)
    {
        free(dimensionVariables[x]);
    }
    free(dimensionVariables);
    for(int x = 0; x < numOfPoints; x++)
    {
        free(data[x]);
    }
    free(data);
    time_req = clock()- time_req;
    cout << "Program took " << (float)time_req/CLOCKS_PER_SEC << " seconds" << endl;
	return(0);
}

//Fix this later. It needs to be able to read in binary instad of text
void ReadData_bin(string filename, std::vector <std::vector<double> > * listOfPoints)
{
	std::cout.precision(17);
	fstream input;
	list <string> lines;
	double buf1;
	char thing[255];
    string buf2;
    std::vector<double> temp;
    input.open (filename, ios::in | ios::out | ios::binary);
    //input.read(thing, sizeof(char));
    //cout << thing << "\n";
    //input.read(reinterpret_cast<char*>(&buf1), sizeof(double));
    //cout << buf1 << "\n";
    int counter = 0;
    while (input.read(reinterpret_cast<char*>(&buf1), sizeof(double)))
    {
        //std::cout << buf1 << '\n';
        if(counter < 3)
        	temp.push_back(buf1);
        if(buf1 != 0 & counter > 2)
        	temp.push_back(buf1);
        if(buf1 == 0 & counter > 2)
        {
        	//temp.push_back(buf1);
        	//for(int x = 0; x < temp.size(); x++)
        	//cout << temp[x] << " ";
        	//cout << "\n";
        	(*listOfPoints).push_back(temp);
        	while(temp.size())
        		temp.pop_back();
        	counter = 0;
        }
        counter++;
    }

    input.close();
}
