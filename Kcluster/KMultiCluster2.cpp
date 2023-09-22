#include <iostream>
#include <fstream>
#include <cmath>
#include <list> 
#include <sstream>
#include <vector>
#include <algorithm>
#include <string.h>


using namespace std;
static int ** listOfRowsGraph;

//Reads in binary nmr spectral data output from thresholding algorithm
void ReadData_bin(string filename, std::vector <std::vector<double> > * listOfPoints);
//determines the max and min of each dimension in the spectrum
double ** grabDimensions(int dim, double ** data, int numOfPoints);

//from ASRG
//Main ASRG algorithm
void ASRG(std::list <std::vector <double> > Seeds, std::vector <std::vector<double> >  listOfPoints, double ** dimensionVariables, struct imageStats imstat, std::vector <std::vector <std::vector <int> > > dataGraph, int numOfPoints, int dim);
//Generates a single cluster using ASRG given a dataset and a starting data point.
void ASRGsingleCluster(int start, std::vector <int> * cluster, int currentCluster, int dim, std::vector <std::vector<double> > * listOfPoints, double ** dimensionVariables, struct imageStats imstat, std::vector <std::vector <std::vector <int> > > * dataGraph, int numOfPoints);
//Outputs the clusters as individual .txt files to be processed seperatly using multiprocessing.
void writeFile2(std::vector <int> finalList, std::vector <std::vector<double> > listOfPoints, int numOfClusters, int numOfPoints, int dim, int flag);
//impliments binary search to find the index of the datapoint passed in the "listOfPoints" list.
int findListIndex(std::vector <double> input, std::vector <std::vector<double> >  listOfPoints, int numOfPoints, int dim);
//Generates an array the same sizes as the spectra, where the locations with data contain their index in "listOfPoints" +1, and the locations without datapoints contain a 0. This graph is used in ASRG to determine the clusters.
std::vector <std::vector <std::vector <int> > > generateImageGraph(double ** dimensionVariables, double ** data, int dim, int numOfPoints, struct imageStats * imstat, vector < vector <int> > * transistionList); 
//Converts the corrdinates of a pixel in the graph to an index in the "listOfPoints" list.
std::vector <int> convertToIndex(vector<double> input, double ** dimensionVariables, struct imageStats imstat, int dim);

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
}; //this contains the variables that when conbined with dimension variables can convert frequency x and y coorindates to the matrix x and y coordinates and back again

double ** grabDimensions(int dim, double ** data, int numOfPoints)
{
    double ** dimensionVariables = NULL;
    dimensionVariables = (double**)malloc(dim*sizeof(double*));
    if(dimensionVariables == NULL)
    {
        cout << "ERROR\n";
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
    
    for(int x = 0; x < numOfPoints; x++)
    {
        int xindex = int((data[x][0] - dimensionVariables[0][0])/imstat->smallestDistx);
        int yindex = int((data[x][1] - dimensionVariables[1][0])/imstat->smallestDisty);
        int zindex = 0;
        if(dim > 2)
        {
            zindex = int((data[x][2] - dimensionVariables[2][0])/imstat->smallestDistz);
        }
        if(dim == 2)
        {
            returnable[0][yindex][xindex] = x + 1;
        }
        if(dim == 3)
        {
            returnable[zindex][yindex][xindex] = x + 1;
        }
    }
    return(returnable);
}

//*Maybe find a way to save this so we don't have to keep looking?*//
int findListIndex(std::vector <double> input, std::vector <std::vector<double> >  listOfPoints, int numOfPoints, int dim)
{
    int smallest = 0;
    int largest = numOfPoints;
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
	if(stack.size() == 1)
	{
		return;
	}
	for(int x = 0; x < stack.size(); x++)
	{
		for(int y = 0; y < stack.size(); y++)
		{
			if((x < stack.size()) && (y < stack.size()))
			{
				if(x != y)
				{
					if(stack[x] == stack[y])
					{
						stack.erase(stack.begin() + y);
					}
				}
				
			}
		}
	}
}

void ASRGsingleCluster(int start, std::vector <int> * cluster, int currentCluster, int dim, std::vector <std::vector<double> > * listOfPoints, double ** dimensionVariables, struct imageStats imstat, std::vector <std::vector <std::vector <int> > > * dataGraph, int numOfPoints)
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
    int counter = 0;
    while(stack.size() > 0)
    {
        if((*cluster)[working] > 0)
        {
            stack.pop_back();
        }
        if(working > -1 && (*cluster)[working] <= 0)
        {
            (*cluster)[working] = currentCluster;
            workingIndex = convertToIndex((*listOfPoints)[working], dimensionVariables, imstat, dim);
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
            if(workingIndex[1] > 0)
            {
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
            if(dim > 2)
            {
                if(back.size() > 0)
                {
                    if((*dataGraph)[back[2]][back[1]][back[0]] > 0)
                    {
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
                    start = (*dataGraph)[up[2]][up[1]][up[0]] - 1;
                    if(start > -1 && (*cluster)[start] == 0)
                    {
                        stack.push_back(start);
                        (*cluster)[start] = -1;
                    }
                }
                up.clear();
            }
            if(upleft.size() > 0)
            {
                if(dim == 2)
                {
                    zindex = 0;
                    upleft.push_back(0.0);
                }
                if((*dataGraph)[upleft[2]][upleft[1]][upleft[0]] > 0)
                {
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
                    start = (*dataGraph)[downright[2]][downright[1]][downright[0]] - 1;
                    if(start > -1 && (*cluster)[start] == 0)
                    {
                        stack.push_back(start);
                        (*cluster)[start] = -1;
                    }
                }
                downright.clear();
            }
            if(right.size() > 0)
            {
                if(dim == 2)
                {
                    zindex = 0;
                    right.push_back(0.0);
                }
                if((*dataGraph)[right[2]][right[1]][right[0]] > 0)
                {
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
        ASRGsingleCluster(start, &cluster, currentCluster, dim, &listOfPoints, dimensionVariables, imstat, &dataGraph, numOfPoints);
        it++;
        while(it != Seeds.end())
        {
            start = findListIndex((*it), listOfPoints, numOfPoints, dim);
            if(cluster[start] < 1)
            {
                currentCluster++;
                ASRGsingleCluster(start, &cluster, currentCluster, dim, &listOfPoints, dimensionVariables, imstat, &dataGraph, numOfPoints);
            }
            it++;
        }
    }
    string output_string = "";
    for(int y = 0; y < numOfPoints; y++)
    {
        if(cluster[y] == 0)
        {
        	output_string = "Index: " + to_string(y) + " of " + to_string(numOfPoints -1) + " " + to_string((double)(100.0*(y/(double)(numOfPoints-1)))) + " Percent Complete\r";
            cout << output_string;
            start = y;
            currentCluster++;
            ASRGsingleCluster(start, &cluster, currentCluster, dim, &listOfPoints, dimensionVariables, imstat, &dataGraph, numOfPoints);
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
            ASRGsingleCluster(start, &cluster, currentCluster, dim, &listOfPoints, dimensionVariables, imstat, &dataGraph, numOfPoints);
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
    std::list <std::vector <int> >:: iterator it2;
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
		lengths[finalList2[x]-1][0] = lengths[finalList2[x]-1][0] + 1;
	}

	
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
        string outstring;
        count = lengths[x-1][0];
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
    }
}

int main(int argc, char ** argv)
{
    clock_t time_req;
    time_req = clock();
    srand((int)time(NULL));
    std::vector <std::vector<double> >  listOfPoints;
    std::vector <std::vector<double> >  listOfPoints_temp;
    vector < vector <int> > transistionList;
    
    ReadData_bin("Outlist.bin", &listOfPoints);
	
    int dim = 0;
    dim = (int)listOfPoints.begin()->size() - 1;
	int numOfPoints = (int)listOfPoints.size();
	cout << "DIM: " << dim << " Length: " << numOfPoints << "\n";
    string runNum = argv[2];
    double ** dimensionVariables = NULL;
    double ** data; 

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
    std::list <std::vector <double> > Seeds;
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
    int counter = 0;
    while (input.read(reinterpret_cast<char*>(&buf1), sizeof(double)))
    {
        if(counter < 3)
        	temp.push_back(buf1);
        if(buf1 != 0 & counter > 2)
        	temp.push_back(buf1);
        if(buf1 == 0 & counter > 2)
        {
        	(*listOfPoints).push_back(temp);
        	while(temp.size())
        		temp.pop_back();
        	counter = 0;
        }
        counter++;
    }

    input.close();
}
