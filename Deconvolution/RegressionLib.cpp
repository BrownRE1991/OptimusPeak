#include "RegressionLib.h"

using namespace std;

void gene_inital(gene& g, int dim)
{
    g.height = 0.0;
    g.alpha = 0.0;
    if(g.center.size() < dim)
    {
        for(int x = 0; x < dim; x++)
        {
            g.center.push_back(0.0);
        }
    }
    if(g.width.size() < dim)
    {
       for(int x = 0; x < dim; x++)
        {
            g.width.push_back(0.0);
        }
    }
    if(g.center.size() == dim)
    {
        for(int x = 0; x < dim; x++)
        {
            g.center[x] = 0.0;
        }
    }
    if(g.width.size() == dim)
    {
        for(int x = 0; x < dim; x++)
        {
            g.width[x] = 0.0;
        }
    }
    g.dim = dim;
	g.fitness = 1000000000000000000;
}

double findSinglePointInModel(int dim, gene a, int numofPoints, vector < vector <double> > & field, int i)
{
    double val = 0.0;
    if(dim > 1)
    {
        val=(a).height*PF_Dual_quick(dim, field[i],(a).alpha,a.center,a.width);  //evaluate the hill
    }
	return(val);
}

double penalty(gene a, vector <double> edges)
{
	double results = 0.0;
	
	if(a.center[0] < edges[0])
	{
		results = results + edges[0] - a.center[0];
	}
	if(a.center[0] > edges[1])
	{
		results = results - edges[1] + a.center[0];
	}
	if(a.dim > 1)
	{
		if(a.center[1] < edges[2])
		{
			results = results + edges[2] - a.center[1];
		}
		if(a.center[1] > edges[3])
		{
			results = results - edges[3] + a.center[1];
		}
	}
	if(a.dim > 2)
	{
		if(a.center[2] < edges[4])
		{
			results = results + edges[4] - a.center[2];
		}
		if(a.center[2] > edges[5])
		{
			results = results - edges[5] + a.center[2];
		}
	}
	if(a.dim > 3)
	{
		if(a.center[3] < edges[6])
		{
			results = results + edges[6] - a.center[3];
		}
		if(a.center[3] > edges[7])
		{
			results = results - edges[7] + a.center[3];
		}
	}
	
	return(results);
}

double findFitness(double scale, int dim, gene a, int numofPoints, vector < vector <double> > & field, bool min, vector <double> edges)
{//evaluate fit
    int i;                  //loop index
    double val = 0.0;
    double delta = 0.0;
    double ttl = 0.0;   //value, different, total error
    ttl=0;  //zero the accumulator
    for(i=0;i<numofPoints;i++)
    {//loop over the data points
		val=(a).height*PF_Dual_quick(dim, field[i],(a).alpha,a.center,a.width);  //evaluate the hill
        delta=(val-field[i][dim])/scale;      //compute difference from the actual value
        ttl+=delta*delta;              //sum squared error
    }
    ttl/=((double)numofPoints);  //mean squared error
    ttl = sqrt(ttl);
    if(min == false)
    {
    	ttl = ttl - penalty(a, edges);
		return(1.0/ttl);
    }
    ttl = ttl + penalty(a, edges);
    return(ttl);   //return RMS error
}


void ReadData(string input1, vector <std::vector<double> > * listOfPoints)
{//read in the data

	fstream input;
	char buf1[256];
	string buf2;
    string line;
	
	string inputFile = input1;
	
    std::vector<double> out;
	input.open(inputFile,ios::in);
	input.getline(buf1, 255); //Get first line
    line = buf1;
    tokenize(line, ' ', &out);
    std::vector<double> buffer = std::vector<double> (out);
    out.erase(out.begin(), out.end());
    listOfPoints->push_back(buffer);
    
	while((strlen(buf1)>0))//while we have data
	{
        input.getline(buf1, 255);
        if(input.eof())
        {
            input.close();
            return;
        }
        line = buf1;
        tokenize(line, ' ', &out);
        std::vector<double> buffer = std::vector<double> (out);
        out.erase(out.begin(), out.end());
        listOfPoints->push_back(buffer);
		//lines.push_back(buf1);
	}
	input.close();
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

vector <double> findMax(vector <std::vector<double> > & listOfPoints, int numOfPoints, int dim)
{
	vector <double> currentMax;
	double MaxH = 0;
	vector <std::vector<double> > :: iterator it;
	for(int x = 0; x < dim+1; x++)
	{
		currentMax.push_back(0.0);
	}
	for(it = listOfPoints.begin(); it != listOfPoints.end(); it++)
	{
		if(dim == 2)
		{
			if(abs((*it).at(2)) > abs(MaxH))
			{
				(MaxH) = (*it).at(2);
				currentMax[0] = (*it).at(0);
				currentMax[1] = (*it).at(1);
				currentMax[2] = (*it).at(2);
			}
		}
		if(dim == 3)
		{
			if(abs((*it).at(3)) > abs(MaxH))
			{
				(MaxH) = (*it).at(3);
				currentMax[0] = (*it).at(0);
				currentMax[1] = (*it).at(1);
				currentMax[2] = (*it).at(2);
				currentMax[3] = (*it).at(3);
			}
		}
 		if(dim == 4)
 		{
 			if(abs((*it).at(4)) > abs(MaxH))
 			{
				(MaxH) = (*it).at(4);
				currentMax[0] = (*it).at(0);
				currentMax[1] = (*it).at(1);
				currentMax[2] = (*it).at(2);
				currentMax[3] = (*it).at(3);
				currentMax[4] = (*it).at(4);
 			}
 		}
	}
	return(currentMax);
}

double findMin(double ** listOfPoints, int numOfPoints, int dim)
{
	double min = 10000;
	for(int x = 0; x < numOfPoints; x++)
	{
		if(dim == 2)
		{
			if(listOfPoints[x][2] < min)
			{
				min = listOfPoints[x][2];
			}
		}
		if(dim == 3)
		{
			if(listOfPoints[x][3]< min)
			{
				min = listOfPoints[x][3];
			}
		}
 		if(dim == 4)
 		{
 			if(listOfPoints[x][4] < min)
 			{
				min = listOfPoints[x][4];
 			}
 		}
	}
	return(min);
}

double findMin(vector <std::vector<double> > & listOfPoints, int numOfPoints, int dim)
{
	double min = 10000;
	vector <std::vector<double> > :: iterator it;
	for(it = listOfPoints.begin(); it != listOfPoints.end(); it++)
	{
		if(dim == 2)
		{
			if((*it).at(2) < min)
			{
				min = (*it).at(2);
			}
		}
		if(dim == 3)
		{
			if((*it).at(3) < min)
			{
				min = (*it).at(3);
			}
		}
 		if(dim == 4)
 		{
 			if((*it).at(4) < min)
 			{
				min = (*it).at(4);
 			}
 		}
	}
	return(min);
}

vector <double> findApproxLinewidth(vector <std::vector<double> > & listOfPoints, int numOfPoints, int dim, vector <double> currentMax)
{
	int bestoOutOf = 5;
	int bestOutOfMod = 6;
	vector <double> linewidths;
	for(int x = 0; x < dim; x++)
	{
		linewidths.push_back(0.0);
	}
	vector <vector <double> * > set1;
	vector <vector <double> * > set2;
	vector <vector <double> * > set3;
	vector <vector <double> * > set4;

	vector <std::vector<double> > :: iterator it;
	it = listOfPoints.begin();

	for(it = listOfPoints.begin(); it != listOfPoints.end(); it++)
	{
		if(dim == 2)
		{
			if(abs((*it)[0] - currentMax[0]) < 0.00001)
			{
				set1.push_back(&(*it));
			}
			if((*it)[1] == currentMax[1])
			{
				set2.push_back(&(*it));
			}
		}
		if(dim == 3)
		{
			if(((*it)[0] == currentMax[0]) & ((*it)[2] == currentMax[2]))
			{
				set1.push_back(&(*it));
			}
			if(((*it)[1] == currentMax[1]) & ((*it)[2] == currentMax[2]))
			{
				set2.push_back(&(*it));
			}
			if(((*it)[0] == currentMax[0]) & ((*it)[1] == currentMax[1]))
			{
				set3.push_back(&(*it));
			}
		}
 		if(dim == 4)
 		{
			if(((*it)[0] == currentMax[0]) & ((*it)[2] == currentMax[2]) & ((*it)[3] == currentMax[3]))
			{
				set1.push_back(&(*it));
			}
			if(((*it)[1] == currentMax[1]) & ((*it)[2] == currentMax[2]) & ((*it)[3] == currentMax[3]))
			{
				set2.push_back(&(*it));
			}
			if(((*it)[0] == currentMax[0]) & ((*it)[1] == currentMax[1]) & ((*it)[3] == currentMax[3]))
			{
				set3.push_back(&(*it));
			}
			if(((*it)[0] == currentMax[0]) & ((*it)[1] == currentMax[1]) & ((*it)[2] == currentMax[2]))
			{
				set4.push_back(&(*it));
			}
 		}
	}

	vector <double> closest1;
	closest1.push_back(10000);
	closest1.push_back(0);
	vector <double> closest2;
	closest2.push_back(10000);
	closest2.push_back(0);
	vector <double> closest3;
	closest3.push_back(10000);
	closest3.push_back(0);
	vector <double> closest4;
	closest4.push_back(10000);
	closest4.push_back(0);
	vector <vector <double> > diflist;
	vector <double> diflist1;
	vector <double> diflist2;
	vector <double> diflist3;
	vector <double> diflist4;
	vector <double> closestTemp;

	//check every distance from from the point to half the height straight up and down.
	for(int y = 0; y < set1.size(); y++)
	{
		diflist1.push_back(abs((*set1[y])[dim] - (currentMax[dim]/2)));
		diflist.push_back((*set1[y]));
	}

	double temp;
	int temp2;

	//I look at the closest "bestoOutOf" points to the half height and choose the closest one to the center.
	if(diflist1.size() <= bestoOutOf*2)
	{
		bestoOutOf = int((diflist1.size()/bestOutOfMod)+2);
	}

	if(diflist1.size() > 1)
	{
		for(int x = 0; x < bestoOutOf; x++)
		{
			temp2 = int(std::min_element(diflist1.begin(), diflist1.end()) - diflist1.begin());
			temp = abs(diflist[temp2][1] - currentMax[1]);
			if(temp > 0)
			{
				closestTemp.push_back(temp);
			}
			diflist1[temp2] = 10000000000000;
		}
	}
	
	if(closestTemp.size() == 0)
	{
		closestTemp.push_back(0.04);
	}
	
	temp = (*std::min_element(closestTemp.begin(), closestTemp.end()));
	
	if(temp != 0)
	{
		linewidths[1] = 2*temp;
	}
	while(!diflist.empty())
	{
		diflist.pop_back();
	}
	while(!closestTemp.empty())
	{
		closestTemp.pop_back();
	}


	for(int y = 0; y < set2.size(); y++)
	{
		diflist2.push_back(abs((*set2[y])[dim] - (currentMax[dim]/2)));
		diflist.push_back((*set2[y]));
	}

	if(diflist2.size() <= bestoOutOf*2)
	{
		bestoOutOf = int((diflist2.size()/bestOutOfMod)+2);
	}
	
	if(diflist2.size() > 1)
	{
		for(int x = 0; x < bestoOutOf; x++)
		{
			temp2 = int(std::min_element(diflist2.begin(), diflist2.end()) - diflist2.begin());
			temp = abs(diflist[temp2][0] - currentMax[0]);
			if(temp > 0)
			{
				closestTemp.push_back(temp);
			}
			diflist2[temp2] = 10000000;
		}
	}
	if(closestTemp.size() == 0)
	{
		closestTemp.push_back(0.04);
	}
	bestoOutOf = 5;

	temp = (*std::min_element(closestTemp.begin(), closestTemp.end()));

	if(temp != 0)
	{
		linewidths[0] = 2*temp;
	}
	while(!diflist.empty())
	{
		diflist.pop_back();
	}
	while(!closestTemp.empty())
	{
		closestTemp.pop_back();
	}

	if(dim > 2)
	{
		for(int y = 0; y < set3.size(); y++)
		{
			diflist3.push_back(abs((*set3[y])[dim] - (currentMax[dim]/2)));
			diflist.push_back((*set3[y]));
		}
		if(diflist3.size() <= bestoOutOf*2)
		{
			bestoOutOf = int((diflist3.size()/bestOutOfMod)+2);
		}
		if(diflist3.size() > 1)
		{
			for(int x = 0; x < bestoOutOf; x++)
			{
				temp2 = int(std::min_element(diflist3.begin(), diflist3.end()) - diflist3.begin());
				temp = abs(diflist[temp2][2] - currentMax[2]);
				if(temp > 0)
				{
					closestTemp.push_back(temp);
				}
				diflist3[temp2] = 10000000;
			}
		}
		if(closestTemp.size() == 0)
		{
			closestTemp.push_back(0.04);
		}
		temp = (*std::min_element(closestTemp.begin(), closestTemp.end()));
		if(temp != 0)
		{
			linewidths[2] = 2*temp;
		}
		while(!diflist.empty())
		{
			diflist.pop_back();
		}
		while(!closestTemp.empty())
		{
			closestTemp.pop_back();
		}
	}
	if(dim > 3)
	{
		for(int y = 0; y < set4.size(); y++)
		{
			diflist4.push_back(abs((*set4[y])[dim] - (currentMax[dim]/2)));
			diflist.push_back((*set4[y]));
		}
		if(diflist4.size() <= bestoOutOf*2)
		{
			bestoOutOf = int((diflist4.size()/bestOutOfMod) + 2);
		}
		if(diflist4.size() > 1)
		{
			for(int x = 0; x < bestoOutOf; x++)
			{
				temp2 = int(std::min_element(diflist4.begin(), diflist4.end()) - diflist4.begin());
				temp = abs(diflist[temp2][3] - currentMax[3]);
				if(temp > 0)
				{
					closestTemp.push_back(temp);
				}
				diflist4[temp2] = 10000000;
			}
		}
		if(closestTemp.size() == 0)
		{
			closestTemp.push_back(0.04);
		}
				
		temp = (*std::min_element(closestTemp.begin(), closestTemp.end()));
		if(temp != 0)
		{
			linewidths[3] = 2*temp;
		}
		while(!diflist.empty())
		{
			diflist.pop_back();
		}
		while(!closestTemp.empty())
		{
			closestTemp.pop_back();
		}
	}
	linewidths.push_back((currentMax[dim]/2));
	return(linewidths);
}

double dist(vector <double> vec1, vector <double> vec2)
{
	double dist = 0;
	double sum = 0;
	
	for(int x = 0; x < vec1.size()-1; x++)
	{
		sum += pow((vec1[x] - vec2[x]),2);
	}
	dist = sum/(vec1.size()-1);
	return(dist);
}

vector <double> findMaxSecond(vector <std::vector<double> > & listOfPoints, int numOfPoints, int dim, vector <double> currentMax, vector <double> widths)
{
	double MaxH = 0;
	vector <double> max2;
	for(int x = 0; x < dim+1; x++)
	{
		max2.push_back(-10000000.0);
	}
	
	double w = 0;
	
	w = (*min_element(widths.begin(), widths.end()))/2; //I dont think I use this. Remove.
	double elipse = 0;
	double c = 1.0;
	double adjust = 1.0;
	
	vector <std::vector<double> > :: iterator it;
	for(it = listOfPoints.begin(); it != listOfPoints.end(); it++)
	{
		if(dim == 2)
		{
			elipse = (((*it)[0] - currentMax[0])/(adjust*widths[0]))*(((*it)[0] - currentMax[0])/(adjust*widths[0])) + (((*it)[1] - currentMax[1])/(adjust*widths[1]))*(((*it)[1] - currentMax[1])/(adjust*widths[1]));
		}
		if(dim == 3)
		{
			elipse = (((*it)[0] - currentMax[0])/(adjust*widths[0]))*(((*it)[0] - currentMax[0])/(adjust*widths[0])) + (((*it)[1] - currentMax[1])/(adjust*widths[1]))*(((*it)[1] - currentMax[1])/(adjust*widths[1])) + (((*it)[2] - currentMax[2])/(adjust*widths[2]))*(((*it)[2] - currentMax[2])/(adjust*widths[2]));
		}
		if(dim == 4)
		{
			elipse = (((*it)[0] - currentMax[0])/(adjust*widths[0]))*(((*it)[0] - currentMax[0])/(adjust*widths[0])) + (((*it)[1] - currentMax[1])/(adjust*widths[1]))*(((*it)[1] - currentMax[1])/(adjust*widths[1])) + (((*it)[2] - currentMax[2])/(adjust*widths[2]))*(((*it)[2] - currentMax[2])/(adjust*widths[2])) + (((*it)[3] - currentMax[3])/(adjust*widths[3]))*(((*it)[3] - currentMax[3])/widths[3]);
		}
		if(elipse > 1)
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
	it = listOfPoints.begin();
	if((max2[0] == -10000000.0) && (max2[1] == -10000000.0))
	{
		if(dim == 3)
		{
			if(max2[2] == -10000000.0)
			{
				max2[0] = (*it).at(0);
				max2[1] = (*it).at(1);
				max2[2] = (*it).at(2);
				max2[3] = (*it).at(3);
			}
		}
		else
		{
			max2[0] = (*it).at(0);
			max2[1] = (*it).at(1);
			max2[2] = (*it).at(2);
		}
	}
	return(max2);
}

void findStats4D(double * MaxH, double * MinH, double * SDx, double * SDy, double * SDz, double * SDw, vector <std::vector<double> > & listOfPoints, int numOfPoints, double * minC1, double * maxC1, double * minC2, double * maxC2, double * minC3, double * maxC3, double * minC4, double * maxC4, int dim)
{
	double sumx= 0;
	double sumy = 0;
	double sumz = 0;
	double sumw = 0;
	double meanx = 0;
	double meany = 0;
	double meanz = 0;
	double meanw = 0;
	list <std::vector<double> > :: iterator it;
	for(int x = 0; x < numOfPoints; x++)
	{
		sumx = sumx + listOfPoints[x][0];
		if(dim == 2)
		{
			if(listOfPoints[x][2] > *MaxH)
			{
				(*MaxH) = listOfPoints[x][2];
			}
			if(listOfPoints[x][2] < *MinH)
			{
				(*MinH) = listOfPoints[x][2];
			}
		}
		if(dim == 3)
		{
			if(listOfPoints[x][3] > *MaxH)
			{
				(*MaxH) = listOfPoints[x][3];
			}
			if(listOfPoints[x][3] < *MinH)
			{
				(*MinH) = listOfPoints[x][3];
			}
		}
 		if(dim == 4)
 		{
 			if(listOfPoints[x][4] > *MaxH)
 			{
				(*MaxH) = listOfPoints[x][4];
 			}
			if(listOfPoints[x][4] < *MinH)
			{
				(*MinH) = listOfPoints[x][4];
			}
 		}
		if(listOfPoints[x][0] > *maxC1)
		{
			(*maxC1) = listOfPoints[x][0];
		}
		if(listOfPoints[x][0] < *minC1)
		{
			(*minC1) = listOfPoints[x][0];
		}
		if(dim > 1)
		{
			if(listOfPoints[x][1] > *maxC2)
			{
				(*maxC2) = listOfPoints[x][1];
			}
			if(listOfPoints[x][1] < *minC2)
			{
				(*minC2) = listOfPoints[x][1];
			}
		}
		if(dim > 2)
		{
			if(listOfPoints[x][2] > *maxC3)
			{
				(*maxC3) = listOfPoints[x][2];
			}
			if(listOfPoints[x][2] < *minC3)
			{
				(*minC3) = listOfPoints[x][2];
			}
		}
		if(dim > 3)
		{
			if(listOfPoints[x][3] > *maxC4)
			{
				(*maxC4) = listOfPoints[x][3];
			}
			if(listOfPoints[x][3] < *minC4)
			{
				(*minC4) = listOfPoints[x][3];
			}
		}

		if(dim > 1)
		{
			sumy = sumy + listOfPoints[x][1];
		}
		if(dim > 2)
		{
			sumz = sumz + listOfPoints[x][2];
		}
 		if(dim > 3)
 		{
 			sumw = sumw + listOfPoints[x][3];
 		}
	}
	(meanx) = sumx/numOfPoints;
	if(dim > 1)
	{
		(meany) = sumy/numOfPoints;
	}
	if(dim > 2)
	{
		(meanz) = sumz/numOfPoints;
	}
	if(dim > 3)
	{
		(meanw) = sumw/numOfPoints;
	}
	double SDStartx = 0;
	double SDStarty = 0;
	double SDStartz = 0;
	double SDStartw = 0;
	for(int x = 0; x < numOfPoints; x++)
	{
		SDStartx = SDStartx + (listOfPoints[x][0]-(meanx))*(listOfPoints[x][0]-(meanx));
		if(dim > 1)
		{
			SDStarty = SDStarty + (listOfPoints[x][1]-(meany))*(listOfPoints[x][1]-(meany));
		}
		if(dim > 2)
		{
			SDStartz = SDStartz + (listOfPoints[x][2]-(meanz))*(listOfPoints[x][2]-(meanz));
		}
		if(dim > 3)
		{
			SDStartw = SDStartw + (listOfPoints[x][3]-(meanw))*(listOfPoints[x][3]-(meanw));
		}
	}
	(*SDx) = sqrt(SDStartx/numOfPoints);
	if(dim > 1)
	{
		(*SDy) = sqrt(SDStarty/numOfPoints);
	}
	if(dim > 2)
	{
		(*SDz) = sqrt(SDStartz/numOfPoints);
	}
	if(dim > 3)
	{
		(*SDz) = sqrt(SDStartz/numOfPoints);
	}
}


void Repair(gene& g)
{//run any needed repairs
    if(g.alpha>1.0) g.alpha=1.0;  //morph parameter cannot exceede 1
    if(g.alpha<0.0)g.alpha=0.0;  //morph parameter must be at least 0
    for(int x = 0; x < g.dim; x++)
    {
        if(g.width[x]<=0.0)g.width[x]=0.00001;
    }
}


void initializePopulation(vector <gene> & pop, double maxAlpha, double minAlpha, int dim, int popsize, vector <vector <double> > & field, int numofPoints, double mutsize, std::vector <double> max, std::vector <double> width, vector <double> edges)
{
	gene temp;
	gene_inital(temp, dim);
	for(int x = 0; x < popsize; x++)
	{
		if(pop.size() <= x)
		{
			gene_inital(temp, dim);
			pop.push_back(temp);
		}
		else
		{
			gene_inital(pop[x], dim);
		}
	}
	double ran;
	double min = findMin(field, numofPoints, dim);
	double mod = Gauss(mutsize);
	vector <double> startingLinewidth = width;
	for(int i = 0; i < popsize; i++)
	{
		mod = Gauss(mutsize);
		pop[i].height = max[max.size()-1] + mod;
		ran = ((double)rand()/RAND_MAX)*(maxAlpha - minAlpha);
		pop[i].alpha = minAlpha + ran;

		mod = Gauss(mutsize);
		pop[i].center[0] = max[0]  + mod;
		if(dim > 1)
		{
			mod = Gauss(mutsize);
			pop[i].center[1] = max[1]  + mod;
		}
		if(dim > 2)
		{
			mod = Gauss(mutsize);
			pop[i].center[2] = max[2]  + mod;
		}
		if(dim > 3)
		{
			mod = Gauss(mutsize);
			pop[i].center[3] = max[3]  + mod;
		}

		if(dim > 0)
		{
			mod = Gauss(mutsize);
			pop[i].width[0] = startingLinewidth[0]+ mod;
		}
		if(dim > 1)
		{
			mod = Gauss(mutsize);
			pop[i].width[1] = startingLinewidth[1] + mod;
		}
		if(dim > 2)
		{
			mod = Gauss(mutsize);
			pop[i].width[2] = startingLinewidth[2]  + mod;
		}
		if(dim > 3)
		{
			mod = Gauss(mutsize);
			pop[i].width[3] = startingLinewidth[3]  + mod;
		}
		pop[i].dim = dim;
		Repair(pop[i]);
		pop[i].fitness = findFitness(max[dim]*2 - min/2, dim, pop[i], numofPoints, field, true, edges);
	}
}

double vL(gene a)
{
	double volume = 0.0;
	double product = 1.0;
	for(int x = 0; x < a.dim; x++)
	{
		product = product*(a.width[x]);
	}
	double front = pow((M_PI/2), a.dim);
	volume = front*product;
	return(volume);
}

double vG(gene a)
{
	double volume = 0.0;
	double product = 1.0;
	for(int x = 0; x < a.dim; x++)
	{
		product = product*(a.width[x]);
	}
	double beta_G = M_PI/2.772588722239781;
	beta_G = pow(beta_G, a.dim/2);
	volume = beta_G*product;
	return(volume);
}

double peakVolume(gene a)
{
	double volume = 0.0;
	double product = 1.0;
	for(int x = 0; x < a.dim; x++)
	{
		product = product*(a.width[x]);
	}
	double other = a.alpha*pow((M_PI/2), a.dim);
	other = other + ((1-a.alpha)*pow((M_PI/2.772588722239781), a.dim/2));
	volume = a.height*product*other;
	return(volume);
}

void printGenestoFile(vector <gene> &a, string filename, double cutoff)
{
    fstream aus;  //output file
    aus.open(filename,ios::out);
    aus.precision(10);
	string line = "";
	for(int x = 0; x < a.size(); x++)
	{
		if(a[x].fitness < cutoff)
		{
			for(int i = 0; i < a[x].dim; i++)
			{
				line = line + to_string(a[x].center[i]) + ",";
			}
			for(int i = 0; i < a[x].dim; i++)
			{
				line = line + to_string(a[x].width[i]) + ",";
			}
			line = line + to_string(a[x].height) + ",";
			line = line + to_string(a[x].alpha) + ",";
			aus << line;
			aus << a[x].fitness;
			aus << ",";
			aus << peakVolume(a[x]);
			aus << "\n";
			line = "";
		}
	}
	aus.close();
}


