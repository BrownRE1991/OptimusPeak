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

void gene2hills_inital(gene2hills& g, int dim)
{
	g.fitness = 1000000000000000000;
	gene_inital(g.hill1, dim);
	gene_inital(g.hill2, dim);
	g.dim = dim;
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
	
	/*if(results > 0)
	{
		cout << results << "\n";
	}*/
	//cout << results << "\n";
	return(results);
}

//We want this to have  scalar penalty for going outside the range of the data. I should have the range somewhere from the initialization, so we just need to pass it here, and have the distance from the closest edge to the center be added to the fitness only if it it outside the range of the data. This should prevent run-aways.
double findFitness(double scale, int dim, gene a, int numofPoints, vector < vector <double> > & field, bool min, vector <double> edges)
{//evaluate fit
    int i;                  //loop index
    double val = 0.0;
    double delta = 0.0;
    double ttl = 0.0;   //value, different, total error
    ttl=0;  //zero the accumulator
	//cout << "got to here\n";
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
    //cout << "Fit: " << ttl << "\n";
    ttl = ttl + penalty(a, edges);
    return(ttl);   //return RMS error
}

double findFitness2hills(double scale, int dim, gene2hills a, int numofPoints, double ** field, bool min)
{//evaluate fit
    int i;                  //loop index
    double val = 0.0;
	double sum = 0.0;
    double ssRes = 0.0;
    double delta = 0.0;
    double ttl = 0.0;   //value, different, total error
	//double scale = maxH - minH;
    double * tempCen1 = (double *)malloc(dim*(sizeof(double)));
    double * tempWid1 = (double *)malloc(dim*(sizeof(double)));
    double * tempCen2 = (double *)malloc(dim*(sizeof(double)));
    double * tempWid2 = (double *)malloc(dim*(sizeof(double)));
    for(int x = 0; x < dim; x++)
    {
        tempCen1[x] = 0.0;
        tempWid1[x] = 0.0;
        tempCen2[x] = 0.0;
        tempWid2[x] = 0.0;
    }
    ttl=0;  //zero the accumulator
    for(i=0;i<numofPoints;i++)
    {//loop over the data points
        if(dim > 1)
        {
            tempCen1[0] = a.hill1.center[0];
            tempCen1[1] = a.hill1.center[1];
            tempWid1[0] = a.hill1.width[0];
            tempWid1[1] = a.hill1.width[1];
            tempCen2[0] = a.hill2.center[0];
            tempCen2[1] = a.hill2.center[1];
            tempWid2[0] = a.hill2.width[0];
            tempWid2[1] = a.hill2.width[1];
        }
        if(dim > 2)
        {
            tempCen1[2] = a.hill1.center[2];
            tempWid1[2] = a.hill1.width[2];
            tempCen2[2] = a.hill2.center[2];
            tempWid2[2] = a.hill2.width[2];
        }
        val=(a).hill1.height*PF_Dual(dim, field[i],(a).hill1.alpha,tempCen1,tempWid1);  //evaluate the hill
		val= val + (a).hill2.height*PF_Dual(dim, field[i],(a).hill2.alpha,tempCen2,tempWid2);  //evaluate the hill
		sum += field[i][dim];
        delta=(val-field[i][dim])/scale;      //compute difference from the actual value
        ssRes+=((val-field[i][dim])*(val-field[i][dim]));
        ttl+=delta*delta;              //sum squared error
    }
    double ssTot = 0.0;
    for(i=0;i<numofPoints;i++)
    {
        ssTot += ((field[i][dim] - sum)*(field[i][dim] - sum))/(double)numofPoints;
    }
    ttl/=((double)numofPoints);  //mean squared error
    ttl = sqrt(ttl);
    free(tempCen1);
    free(tempWid1);
    free(tempCen2);
    free(tempWid2);
	//r_squared = 1 - (ssRes/ssTot);
	//fit = ttl;
	if(min == false)
	{
		return(1.0/ttl);
	}
	//cout << "Fit: " << fit << "\n";
    return(ttl);   //return RMS error
}

void ReadData(string input1, vector <std::vector<double> > * listOfPoints)
{//read in the data

	fstream input;
	char buf1[256];
	string buf2;
	//list <string> lines;
    string line;
	
	string inputFile = input1;
	
	//cout << inputFile << "\n";
    std::vector<double> out;
	input.open(inputFile,ios::in);
	input.getline(buf1, 255); //Get first line
    line = buf1;
	//cout << line;
    tokenize(line, ' ', &out);
    std::vector<double> buffer = std::vector<double> (out);
    out.erase(out.begin(), out.end());
    listOfPoints->push_back(buffer);
    //char flag = 0;
    
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

// vector <double> findMax(vector < vector <double> >  & listOfPoints, int numOfPoints, int dim)
// {
// 	//cout << "Starting findMax\n";
// 	vector <double> currentMax;
// 	double MaxH = 0;
// 	list <std::vector<double> > :: iterator it;
// 	for(int x = 0; x < dim+1; x++)
// 	{
// 		currentMax.push_back(0.0);
// 	}
// 	for(int x = 0; x < numOfPoints; x++)
// 	{
// 		if(dim == 2)
// 		{
// 			if(abs(listOfPoints[x][2]) > abs(MaxH))
// 			{
// 				(MaxH) = listOfPoints[x][2];
// 				currentMax[0] = listOfPoints[x][0];
// 				currentMax[1] = listOfPoints[x][1];
// 				currentMax[2] = listOfPoints[x][2];
// 			}
// 		}
// 		if(dim == 3)
// 		{
// 			if(abs(listOfPoints[x][3]) > abs(MaxH))
// 			{
// 				(MaxH) = listOfPoints[x][3];
// 				currentMax[0] = listOfPoints[x][0];
// 				currentMax[1] = listOfPoints[x][1];
// 				currentMax[2] = listOfPoints[x][2];
// 				currentMax[3] = listOfPoints[x][3];
// 			}
// 		}
//  		if(dim == 4)
//  		{
//  			if(abs(listOfPoints[x][4]) > abs(MaxH))
//  			{
// 				(MaxH) = listOfPoints[x][4];
// 				currentMax[0] = listOfPoints[x][0];
// 				currentMax[1] = listOfPoints[x][1];
// 				currentMax[2] = listOfPoints[x][2];
// 				currentMax[3] = listOfPoints[x][3];
// 				currentMax[4] = listOfPoints[x][4];
//  			}
//  		}
// 	}
// 	//cout << "Done findMax\n";
// 	//cout << currentMax[2] << "\n";
// 	return(currentMax);
// }


vector <double> findMax(vector <std::vector<double> > & listOfPoints, int numOfPoints, int dim)
{
	//cout << "Starting findMax\n";
	vector <double> currentMax;
	double MaxH = 0;
	vector <std::vector<double> > :: iterator it;
	for(int x = 0; x < dim+1; x++)
	{
		currentMax.push_back(0.0);
	}
	for(it = listOfPoints.begin(); it != listOfPoints.end(); it++)
	{
		//cout << abs((*it).at(2)) << "\n";
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
	//cout << "Done findMax\n";
	//cout << currentMax[2] << "\n";
	return(currentMax);
}

double findMin(double ** listOfPoints, int numOfPoints, int dim)
{
	//cout << "Starting findMax\n";
	double min = 10000;
	//double MaxH = 0;
	//list <std::vector<double> > :: iterator it;
	// for(int x = 0; x < dim+1; x++)
	// {
	// 	currentMax.push_back(0.0);
	// }
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
	// for(it = listOfPoints.begin(); it != listOfPoints.end(); it++)
// 	{
// 		if(dim == 2)
// 		{
// 			if((*it).at(2) < min)
// 			{
// 				min = (*it).at(2);
// 			}
// 		}
// 		if(dim == 3)
// 		{
// 			if((*it).at(3) < min)
// 			{
// 				min = (*it).at(3);
// 			}
// 		}
//  		if(dim == 4)
//  		{
//  			if((*it).at(4) < min)
//  			{
// 				min = (*it).at(4);
//  			}
//  		}
// 	}
	//cout << "Done findMax\n";
	return(min);
}

double findMin(vector <std::vector<double> > & listOfPoints, int numOfPoints, int dim)
{
	//cout << "Starting findMax\n";
	double min = 10000;
	//double MaxH = 0;
	vector <std::vector<double> > :: iterator it;
	// for(int x = 0; x < dim+1; x++)
	// {
	// 	currentMax.push_back(0.0);
	// }
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
	//cout << "Done findMax\n";
	return(min);
}

vector <double> findApproxLinewidth(vector <std::vector<double> > & listOfPoints, int numOfPoints, int dim, vector <double> currentMax)
{
	// cout << "Starting Max: ";
// 	for(int x = 0; x < currentMax.size(); x++)
// 	{
// 		cout << currentMax[x] << " ";
// 	}
// 	cout << "\n";
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
			//cout << (*it)[0] << "\n";
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

	//cout << set1.size() << "\n";
	//cout << set2.size() << "\n";

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
	//double diff1 = 0;
	//double diff2 = 0;
	//double diff3 = 0;
	//double diff4 = 0;
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
	//cout << diflist1.size() << "\n";
	
	// for(int x = 0; x < diflist1.size(); x++)
// 	{
// 		cout << diflist1[x] << " ";
// 	}
// 	cout << "\n";

	double temp;
	//double temp1;
	int temp2;
	//int lastTemp;

	//I look at the closest "bestoOutOf" points to the half height and choose the closest one to the center.
	if(diflist1.size() <= bestoOutOf*2)
	{
		bestoOutOf = int((diflist1.size()/bestOutOfMod)+2);
		//cout << "bestoOutOf: " << bestoOutOf << "\n";
	}
	//cout << "got to here\n";
	//cout << "bestoOutOf: " << bestoOutOf << "\n";
	
	//Fixed! you need to turn the choice into unchoosable even if it is the center point.
	//Maybe set it to while size(clestestTemp < bestOutOf)
// 	while(closestTemp.size() < bestoOutOf)
// 	{
// 		cout << "got to here" << "\n";
// 		temp2 = std::min_element(diflist1.begin(), diflist1.end()) - diflist1.begin();
// 		//cout << temp2 << " " << diflist[temp2][1] << "\n";
// 		temp = abs(diflist[temp2][1] - currentMax[1]);
// 		if(temp > 0)
// 		{
// 			closestTemp.push_back(temp);
// 			//diflist1[temp2] = 10000000000000;
// 		}
// 		diflist1[temp2] = 10000000000000;
// 		// for(int y = 0; y < diflist1.size(); y++)
// // 		{
// // 			cout << diflist1[y] << " ";
// // 		}
// // 		cout << "\n";
// 	}
	if(diflist1.size() > 1)
	{
		for(int x = 0; x < bestoOutOf; x++)
		{
			//cout << "got to here" << x << "\n";
			temp2 = int(std::min_element(diflist1.begin(), diflist1.end()) - diflist1.begin());
			//cout << temp2 << " " << diflist[temp2][1] << "\n";
			temp = abs(diflist[temp2][1] - currentMax[1]);
			if(temp > 0)
			{
				closestTemp.push_back(temp);
				//diflist1[temp2] = 10000000000000;
			}
			diflist1[temp2] = 10000000000000;
			// for(int y = 0; y < diflist1.size(); y++)
			// {
			// 	cout << diflist1[y] << " ";
			// }
			// cout << "\n";
		}
	}
	
	if(closestTemp.size() == 0)
	{
		closestTemp.push_back(0.04);
	}
	
	//cout << "got to here2\n";
	temp = (*std::min_element(closestTemp.begin(), closestTemp.end()));
	
	/*for(int x = 0; x < closestTemp.size(); x++)
	{
		cout << closestTemp[x] << " ";
	}
	cout << "\n";*/
	//cout << temp << "\n";
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
	//cout << bestoOutOf << "\n";
	
	if(diflist2.size() > 1)
	{
		for(int x = 0; x < bestoOutOf; x++)
		{
			temp2 = int(std::min_element(diflist2.begin(), diflist2.end()) - diflist2.begin());
			//cout << temp2 << " " << diflist[temp2][0] << "\n";
			temp = abs(diflist[temp2][0] - currentMax[0]);
			if(temp > 0)
			{
				closestTemp.push_back(temp);
				//diflist2[temp2] = 10000000;
			}
			diflist2[temp2] = 10000000;
		}
	}
	if(closestTemp.size() == 0)
	{
		closestTemp.push_back(0.04);
	}
	
	/*for(int x = 0; x < closestTemp.size(); x++)
	{
		cout << closestTemp[x] << " ";
	}
	cout << "\n";*/
	bestoOutOf = 5;

	temp = (*std::min_element(closestTemp.begin(), closestTemp.end()));
	//cout << temp << "\n";

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
					//diflist2[temp2] = 10000000;
				}
				diflist3[temp2] = 10000000;
			}
		}
		if(closestTemp.size() == 0)
		{
			closestTemp.push_back(0.04);
		}
		temp = (*std::min_element(closestTemp.begin(), closestTemp.end()));
		//cout << temp << "\n";
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
					//diflist2[temp2] = 10000000;
				}
				diflist4[temp2] = 10000000;
			}
		}
		if(closestTemp.size() == 0)
		{
			closestTemp.push_back(0.04);
		}
				
		temp = (*std::min_element(closestTemp.begin(), closestTemp.end()));
		//cout << temp << "\n";
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
	
	//cout << vec1[0] << " " << vec1[1] << " " << vec1[2] << "\n";
	//cout << vec2[0]<< " " << vec2[1] << " " << vec2[2] << "\n";
	
	double sum = 0;
	//cout << vec1.size() << "\n";
	
	for(int x = 0; x < vec1.size()-1; x++)
	{
		sum += pow((vec1[x] - vec2[x]),2);
	}
	dist = sum/(vec1.size()-1);
	//cout << "dist:" << dist << "\n";
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
	//cout << w << "\n";
	double elipse = 0;
	double c = 1.0;
	double adjust = 1.0;
	
	vector <std::vector<double> > :: iterator it;
	for(it = listOfPoints.begin(); it != listOfPoints.end(); it++)
	{
		//cout << currentMax[0] << " " << currentMax[1] << " " << currentMax[2] << "\n";
		//cout << dist((*it), currentMax) << "\n";
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
		// if(dist((*it), currentMax) > w)
		if(elipse > 1)
		{
			//cout << "Outside smallest linewidth\n";
			
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
		//cout << "failed to find max2[0] and max2[1]\n";
		if(dim == 3)
		{
			if(max2[2] == -10000000.0)
			{
				//cout << "failed to find max2[2]\n";
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
	// cout << "second max\n";
// 	for(int x = 0; x < max2.size(); x++)
// 	{
// 		cout << max2[x] << " ";
// 	}
// 	cout << "\n";
	return(max2);
}

void findStats4D(double * MaxH, double * MinH, double * SDx, double * SDy, double * SDz, double * SDw, vector <std::vector<double> > & listOfPoints, int numOfPoints, double * minC1, double * maxC1, double * minC2, double * maxC2, double * minC3, double * maxC3, double * minC4, double * maxC4, int dim)
{
	//cout << "FS HERE1\n";
	//cout << listOfPoints.begin()->size() << "\n";
	//cout << listOfPoints.size() << "\n";
	//int dim = (int)(listOfPoints.begin()->size() - 1);
	double sumx= 0;
	double sumy = 0;
	double sumz = 0;
	double sumw = 0;
	double meanx = 0;
	double meany = 0;
	double meanz = 0;
	double meanw = 0;
	list <std::vector<double> > :: iterator it;
	//cout << dim << "\n";
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
		//cout << "Mins: " << minC1 << " " << minC2 << "\n";

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
	//cout << "FS HERE2\n";
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
	//cout << "FS HERE3\n";
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
	//cout << "FS HERE4\n";
}


void findStats4D(double * MaxH, double * MinH, double * SDx, double * SDy, double * SDz, double * SDw, vector <std::vector<double> > & listOfPoints, int numOfPoints, double * minC1, double * maxC1, double * minC2, double * maxC2, double * minC3, double * maxC3, double * minC4, double * maxC4)
{
	//cout << "FS HERE1\n";
	//cout << listOfPoints.begin()->size() << "\n";
	//cout << listOfPoints.size() << "\n";
	int dim = (int)(listOfPoints.begin()->size() - 1);
	double sumx= 0;
	double sumy = 0;
	double sumz = 0;
	double sumw = 0;
	double meanx = 0;
	double meany = 0;
	double meanz = 0;
	double meanw = 0;
	vector <std::vector<double> > :: iterator it;
	//cout << dim << "\n";
	for(it = listOfPoints.begin(); it != listOfPoints.end(); it++)
	{
		sumx = sumx + (*it).at(0);
		if(dim == 2)
		{
			if((*it).at(2) > *MaxH)
			{
				(*MaxH) = (*it).at(2);
			}
			if((*it).at(2) < *MinH)
			{
				(*MinH) = (*it).at(2);
			}
		}
		if(dim == 3)
		{
			if((*it).at(3) > *MaxH)
			{
				(*MaxH) = (*it).at(3);
			}
			if((*it).at(3) < *MinH)
			{
				(*MinH) = (*it).at(3);
			}
		}
 		if(dim == 4)
 		{
 			if((*it).at(4) > *MaxH)
 			{
				(*MaxH) = (*it).at(4);
 			}
			if((*it).at(4) < *MinH)
			{
				(*MinH) = (*it).at(4);
			}
 		}
		if((*it).at(0) > *maxC1)
		{
			(*maxC1) = (*it).at(0);
		}
		if((*it).at(0) < *minC1)
		{
			(*minC1) = (*it).at(0);
		}
		if(dim > 1)
		{
			if((*it).at(1) > *maxC2)
			{
				(*maxC2) = (*it).at(1);
			}
			if((*it).at(1) < *minC2)
			{
				(*minC2) = (*it).at(1);
			}
		}
		if(dim > 2)
		{
			if((*it).at(2) > *maxC3)
			{
				(*maxC3) = (*it).at(2);
			}
			if((*it).at(2) < *minC3)
			{
				(*minC3) = (*it).at(2);
			}
		}
		if(dim > 3)
		{
			if((*it).at(3) > *maxC4)
			{
				(*maxC4) = (*it).at(3);
			}
			if((*it).at(3) < *minC4)
			{
				(*minC4) = (*it).at(3);
			}
		}
		//cout << "Mins: " << minC1 << " " << minC2 << "\n";

		if(dim > 1)
		{
			sumy = sumy + (*it).at(1);
		}
		if(dim > 2)
		{
			sumz = sumz + (*it).at(2);
		}
 		if(dim > 3)
 		{
 			sumw = sumw + (*it).at(3);
 		}
	}
	//cout << "FS HERE2\n";
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
	for(it = listOfPoints.begin(); it != listOfPoints.end(); it++)
	{
		SDStartx = SDStartx + ((*it).at(0)-(meanx))*((*it).at(0)-(meanx));
		if(dim > 1)
		{
			SDStarty = SDStarty + ((*it).at(1)-(meany))*((*it).at(1)-(meany));
		}
		if(dim > 2)
		{
			SDStartz = SDStartz + ((*it).at(2)-(meanz))*((*it).at(2)-(meanz));
		}
		if(dim > 3)
		{
			SDStartw = SDStartw + ((*it).at(3)-(meanw))*((*it).at(3)-(meanw));
		}
	}
	//cout << "FS HERE3\n";
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
	//cout << "FS HERE4\n";
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

// void initializePopulation(vector <gene> & pop, double maxAlpha, double minAlpha, int dim, int popsize, double ** field, int numofPoints, vector <std::vector<double> > & listOfPoints, double mutsize)
// {
// 	gene temp;
// 	gene_inital(temp, dim);
// 	for(int x = 0; x < popsize; x++)
// 	{
// 		if(pop.size() <= x)
// 		{
// 			gene_inital(temp, dim);
// 			pop.push_back(temp);
// 		}
// 		else
// 		{
// 			gene_inital(pop[x], dim);
// 		}
// 	}
// 	double ran;
// 	vector <double> max = findMax(listOfPoints, numofPoints, dim);
// 	double min = findMin(listOfPoints, numofPoints, dim);
// 	double mod = Gauss(mutsize);
// 	vector <double> startingLinewidth = findApproxLinewidth(listOfPoints, numofPoints, dim, max);
// 	for(int i = 0; i < popsize; i++)
// 	{
// 		mod = Gauss(mutsize);
// 		//ran = ((double)rand()/RAND_MAX)*(boundHTop - boundHBot);
// 		pop[i].height = max[max.size()-1] + mod;
// 		ran = ((double)rand()/RAND_MAX)*(maxAlpha - minAlpha);
// 		pop[i].alpha = minAlpha + ran;
//
// 		mod = Gauss(mutsize);;
// 		//ran = ((double)rand()/RAND_MAX)*(maxC1- minC1);
// 		pop[i].center[0] = max[0]  + mod;
// 		if(dim > 1)
// 		{
// 			mod = Gauss(mutsize);
// 			//ran = ((double)rand()/RAND_MAX)*(maxC2 - minC2);
// 			pop[i].center[1] = max[1]  + mod;
// 		}
// 		if(dim > 2)
// 		{
// 			mod = Gauss(mutsize);
// 			//ran = ((double)rand()/RAND_MAX)*(maxC3 - minC3);
// 			pop[i].center[2] = max[2]  + mod;
// 		}
// 		if(dim > 3)
// 		{
// 			mod = Gauss(mutsize);
// 			//ran = ((double)rand()/RAND_MAX)*(maxC4 - minC4);
// 			pop[i].center[3] = max[3]  + mod;
// 		}
//
// 		if(dim > 0)
// 		{
// 			mod = Gauss(mutsize);
// 			//ran = ((double)rand()/RAND_MAX)*(maxW1 - minW1);
// 			pop[i].width[0] = startingLinewidth[0]+ mod;
// 		}
// 		if(dim > 1)
// 		{
// 			mod = Gauss(mutsize);
// 			//ran = ((double)rand()/RAND_MAX)*(maxW2 - minW2);
// 			pop[i].width[1] = startingLinewidth[1] + mod;
// 		}
// 		if(dim > 2)
// 		{
// 			mod = Gauss(mutsize);
// 			//ran = ((double)rand()/RAND_MAX)*(maxW3 - minW3);
// 			pop[i].width[2] = startingLinewidth[2]  + mod;
// 		}
// 		if(dim > 3)
// 		{
// 			mod = Gauss(mutsize);
// 			//ran = ((double)rand()/RAND_MAX)*(maxW4 - minW4);
// 			pop[i].width[3] = startingLinewidth[3]  + mod;
// 		}
// 		pop[i].dim = dim;
// 		Repair(pop[i]);
// 		pop[i].fitness = findFitness(max[dim]*2 - min/2, dim, pop[i], numofPoints, listOfPoints, true);
// 		//printGene(pop[i]);
// 		//Repair(pop[i]);
// 	}
// }

void initializePopulation(vector <gene> & pop, double maxAlpha, double minAlpha, int dim, int popsize, vector <vector <double> > & field, int numofPoints, double mutsize, std::vector <double> max, std::vector <double> width, vector <double> edges)
{
	//cout << "Made it to here\n";
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
	//cout << "Made it to here1\n";
	double ran;
	//vector <double> max = max;
	//cout << "Made it to here2\n";
	//cout << max[0] << " " << max[1] << "\n";
	double min = findMin(field, numofPoints, dim);
	//cout << "Made it to here3\n";
	double mod = Gauss(mutsize);
	vector <double> startingLinewidth = width;
	// cout << "Linewidths:\n";
// 	for(int x = 0; x < dim; x++)
// 	{
// 		cout << startingLinewidth[x] << " ";
// 	}
// 	cout << "\n";
	for(int i = 0; i < popsize; i++)
	{
		mod = Gauss(mutsize);
		//ran = ((double)rand()/RAND_MAX)*(boundHTop - boundHBot);
		//cout << max[max.size()-1] << "\n";
		//cout << max.size()-1 << "\n";
		pop[i].height = max[max.size()-1] + mod;
		//cout << pop[i].height << "\n";
		ran = ((double)rand()/RAND_MAX)*(maxAlpha - minAlpha);
		pop[i].alpha = minAlpha + ran;

		mod = Gauss(mutsize);
		//cout << max[0] << "\n";
		//ran = ((double)rand()/RAND_MAX)*(maxC1- minC1);
		pop[i].center[0] = max[0]  + mod;
		//cout << pop[i].center[0] << "\n";
		if(dim > 1)
		{
			mod = Gauss(mutsize);
			//ran = ((double)rand()/RAND_MAX)*(maxC2 - minC2);
			pop[i].center[1] = max[1]  + mod;
		}
		if(dim > 2)
		{
			mod = Gauss(mutsize);
			//ran = ((double)rand()/RAND_MAX)*(maxC3 - minC3);
			pop[i].center[2] = max[2]  + mod;
		}
		if(dim > 3)
		{
			mod = Gauss(mutsize);
			//ran = ((double)rand()/RAND_MAX)*(maxC4 - minC4);
			pop[i].center[3] = max[3]  + mod;
		}

		if(dim > 0)
		{
			mod = Gauss(mutsize);
			//ran = ((double)rand()/RAND_MAX)*(maxW1 - minW1);
			pop[i].width[0] = startingLinewidth[0]+ mod;
			//cout << pop[i].width[0] << "\n";
		}
		if(dim > 1)
		{
			mod = Gauss(mutsize);
			//ran = ((double)rand()/RAND_MAX)*(maxW2 - minW2);
			pop[i].width[1] = startingLinewidth[1] + mod;
		}
		if(dim > 2)
		{
			mod = Gauss(mutsize);
			//ran = ((double)rand()/RAND_MAX)*(maxW3 - minW3);
			pop[i].width[2] = startingLinewidth[2]  + mod;
		}
		if(dim > 3)
		{
			mod = Gauss(mutsize);
			//ran = ((double)rand()/RAND_MAX)*(maxW4 - minW4);
			pop[i].width[3] = startingLinewidth[3]  + mod;
		}
		pop[i].dim = dim;
		Repair(pop[i]);
		//cout << boundHTop << " " << max[dim]*2 << "\n";
		//cout << boundHBot << " " << min/2 << "\n";
		pop[i].fitness = findFitness(max[dim]*2 - min/2, dim, pop[i], numofPoints, field, true, edges);
		//printGene(pop[i]);
		//Repair(pop[i]);
	}
}

void initializePopulation2Hills(vector <gene2hills> & pop, double boundHTop, double boundHBot, double maxAlpha, double minAlpha, double maxC1, double minC1, double maxC2, double minC2, double maxC3, double minC3, double maxC4, double minC4,  double maxW1, double minW1, double maxW2, double minW2, double maxW3, double minW3, double maxW4, double minW4, int dim, int popsize, double ** field, int numofPoints, vector <std::vector<double> > & listOfPoints, double mutsize)
{
	//cout << "In initialzePopulation2Hills\n";
	gene2hills temp;
	gene2hills_inital(temp, dim);
	for(int x = 0; x < popsize; x++)
	{
		if(pop.size() <= x)
		{
			gene2hills_inital(temp, dim);
			pop.push_back(temp);
		}
		else
		{
			gene2hills_inital(pop[x], dim);
		}
	}
	double ran;
	vector <double> max = findMax(listOfPoints, numofPoints, dim);
	//cout << max[0] << " " << max[1] << " ";
// 	if(dim == 3)
// 	{
// 		cout << max[2] << " ";
// 		cout << max[3];
// 	}
//  	cout << "\n";
	
	double min = findMin(listOfPoints, numofPoints, dim);
	double mod = Gauss(mutsize);
	vector <double> startingLinewidth = findApproxLinewidth(listOfPoints, numofPoints, dim, max);
	// cout << "LW1:\n";
// 	for(int x = 0; x < dim; x++)
// 	{
// 		cout << startingLinewidth[x] << " ";
// 	}
// 	cout << "\n";
	
	vector <double> max2 = findMaxSecond(listOfPoints, numofPoints, dim, max, startingLinewidth);
	// cout << "Max2:\n";
// 	for(int x = 0; x < dim+1; x++)
// 	{
// 		cout << max2[x] << " ";
// 	}
// 	cout << "\n";
	vector <double> startingLinewidth2 = findApproxLinewidth(listOfPoints, numofPoints, dim, max2);
	// cout << "LW2:\n";
// 	for(int x = 0; x < dim; x++)
// 	{
// 		cout << startingLinewidth2[x] << " ";
// 	}
// 	cout << "\n";
	for(int i = 0; i < popsize; i++)
	{
		mod = Gauss(mutsize);
		//ran = ((double)rand()/RAND_MAX)*(boundHTop - boundHBot);
		pop[i].hill1.height = max[max.size()-1] + mod;
		ran = ((double)rand()/RAND_MAX)*(boundHTop - boundHBot);
		//cout << ran << " " << boundHTop << " " << boundHBot << "\n";
		pop[i].hill2.height = boundHBot + ran;
		ran = ((double)rand()/RAND_MAX)*(maxAlpha - minAlpha);
		pop[i].hill1.alpha = minAlpha + ran;
		ran = ((double)rand()/RAND_MAX)*(maxAlpha - minAlpha);
		pop[i].hill2.alpha = minAlpha + ran;
		
		mod = Gauss(mutsize);
		//ran = ((double)rand()/RAND_MAX)*(maxC1- minC1);
		pop[i].hill1.center[0] = max[0]  + mod;
		mod = Gauss(mutsize);
		//ran = ((double)rand()/RAND_MAX)*(maxC1- minC1);
		pop[i].hill2.center[0] = max2[0]  + mod;
		if(dim > 1)
		{
			mod = Gauss(mutsize);
			//ran = ((double)rand()/RAND_MAX)*(maxC2 - minC2);
			pop[i].hill1.center[1] = max[1]  + mod;
			mod = Gauss(mutsize);
			//ran = ((double)rand()/RAND_MAX)*(maxC2 - minC2);
			pop[i].hill2.center[1] = max2[1]  + mod;
		}
		if(dim > 2)
		{
			mod = Gauss(mutsize);
			//ran = ((double)rand()/RAND_MAX)*(maxC3 - minC3);
			pop[i].hill1.center[2] = max[2]  + mod;
			mod = Gauss(mutsize);
			//ran = ((double)rand()/RAND_MAX)*(maxC3 - minC3);
			pop[i].hill2.center[2] = max2[2]  + mod;
		}
		if(dim > 3)
		{
			mod = Gauss(mutsize);
			//ran = ((double)rand()/RAND_MAX)*(maxC4 - minC4);
			pop[i].hill1.center[3] = max[3]  + mod;
			mod = Gauss(mutsize);
			//ran = ((double)rand()/RAND_MAX)*(maxC4 - minC4);
			pop[i].hill2.center[3] = max2[3]  + mod;
		}

		if(dim > 0)
		{
			mod = Gauss(mutsize);
			//ran = ((double)rand()/RAND_MAX)*(maxW1 - minW1);
			pop[i].hill1.width[0] = startingLinewidth[0]+ mod;
			mod = Gauss(mutsize);
			//ran = ((double)rand()/RAND_MAX)*(maxW1 - minW1);
			pop[i].hill2.width[0] = startingLinewidth2[0]+ mod;
		}
		if(dim > 1)
		{
			mod = Gauss(mutsize);
			//ran = ((double)rand()/RAND_MAX)*(maxW2 - minW2);
			pop[i].hill1.width[1] = startingLinewidth[1]+ mod;
			mod = Gauss(mutsize);
			//ran = ((double)rand()/RAND_MAX)*(maxW2 - minW2);
			pop[i].hill2.width[1] = startingLinewidth2[1]+ mod;
		}
		if(dim > 2)
		{
			mod = Gauss(mutsize);
			//ran = ((double)rand()/RAND_MAX)*(maxW3 - minW3);
			pop[i].hill1.width[2] = startingLinewidth[2]+ mod;
			mod = Gauss(mutsize);
			//ran = ((double)rand()/RAND_MAX)*(maxW3 - minW3);
			pop[i].hill2.width[2] = startingLinewidth2[2]+ mod;
		}
		if(dim > 3)
		{
			mod = Gauss(mutsize);
			//ran = ((double)rand()/RAND_MAX)*(maxW4 - minW4);
			pop[i].hill1.width[3] = startingLinewidth[3]+ mod;
			mod = Gauss(mutsize);
			//ran = ((double)rand()/RAND_MAX)*(maxW4 - minW4);
			pop[i].hill2.width[3] = startingLinewidth2[3]+ mod;
		}
		pop[i].hill1.dim = dim;
		pop[i].hill2.dim = dim;
		Repair(pop[i].hill1);
		Repair(pop[i].hill2);
		pop[i].fitness = findFitness2hills(max[dim]*2 - min/2, dim, pop[i], numofPoints, field, true);
		//printGene(pop[i].hill1);
		//printGene(pop[i].hill2);
	}
}

void printGene(gene &a)
{
	cout << "Dim: " << a.dim << "\n";
	cout << "Height: " << a.height << " Alpha: " << a.alpha << "\n";
	for(int i = 0; i < a.dim; i++)
	{
		cout << a.center[i] << " ";
	}
	cout << "\n";
	for(int i = 0; i < a.dim; i++)
	{
		cout << a.width[i] << " ";
	}
	cout << "\n";
	cout << "Fitness: " << a.fitness << "\n";
	
}

double distanceFromTrue(gene a, gene b)
{
    double ret = 0.0;
    double dist = 0;
	dist = a.height - b.height;
	ret+= dist*dist;
	dist = a.alpha - b.alpha;
	ret+= dist*dist;
    for(int x = 0; x < a.dim; x++)
    {
        dist = a.center[x] - b.center[x];
        ret += dist * dist;
    }
	
    for(int x = 0; x < a.dim; x++)
    {
        dist = a.width[x] - b.width[x];
        ret += dist * dist;
    }
    return ret > 0.0 ? sqrt(ret) : 0.0;
}

void printtoFile(vector <double> a, string filename)
{
    fstream aus;  //output file
    aus.open(filename,ios::out);
    aus.precision(10);
       
    for(int i = 0; i < a.size(); i++)
    {
        aus << a[i] << "\n";
    }
    aus.close();
}

void printtoTwoFile(vector <double> a, vector <double> b, string filename)
{
    fstream aus;  //output file
    aus.open(filename,ios::out);
    aus.precision(10);
       
    for(int i = 0; i < a.size(); i++)
    {
        aus << a[i] << " " << b[i] << "\n";
    }
    aus.close();
}

void printGenetoFile(gene &a, string filename)
{
    fstream aus;  //output file
    aus.open(filename,ios::out);
    aus.precision(10);
	string line = "";
	for(int i = 0; i < a.dim; i++)
	{
		line = line + to_string(a.center[i]) + ",";
	}
	for(int i = 0; i < a.dim; i++)
	{
		line = line + to_string(a.width[i]) + ",";
	}
	line = line + to_string(a.height) + ",";
	line = line + to_string(a.alpha) + ",";
	line = line + to_string(a.fitness) + "\n";
	
	//cout << line << "\n";
	aus << line;
}

//See Allans hills.pdf in the Deconvolution folder for details
double vL(gene a)
{
	double volume = 0.0;
	double product = 1.0;
	for(int x = 0; x < a.dim; x++)
	{
		//cout << a.width[x] << " ";
		product = product*(a.width[x]);
	}
	//cout << "\n";
	//cout << product << "\n";
	double front = pow((M_PI/2), a.dim);
	volume = front*product;
	return(volume);
}

//See Allans hills.pdf in the Deconvolution folder for details
double vG(gene a)
{
	double volume = 0.0;
	double product = 1.0;
	for(int x = 0; x < a.dim; x++)
	{
		//cout << a.width[x] << " ";
		product = product*(a.width[x]);
	}
	//cout << "\n";
	//cout << product << "\n";
	double beta_G = M_PI/2.772588722239781;
	beta_G = pow(beta_G, a.dim/2);
	//cout << beta_G << "\n";
	//volume = sqrt(beta_G*product);
	volume = beta_G*product;
	return(volume);
}

//See Allans hills.pdf in the Deconvolution folder for details
//I'll do the sortcut
double peakVolume(gene a)
{
	double volume = 0.0;
	double product = 1.0;
	for(int x = 0; x < a.dim; x++)
	{
		//cout << a.width[x] << " ";
		product = product*(a.width[x]);
	}
	//cout << "\n";
	//cout << product << "\n";
	double other = a.alpha*pow((M_PI/2), a.dim);
	other = other + ((1-a.alpha)*pow((M_PI/2.772588722239781), a.dim/2));
	volume = a.height*product*other;
	//cout << "Volume: " << volume << "\n";
	//volume = a.height*((a.alpha*vL(a)) + ((1-a.alpha)*vG(a)));
	return(volume);
}

//put volume in here
void printGenestoFile(vector <gene> &a, string filename, double cutoff)
{
    fstream aus;  //output file
	//cout << filename << "\n";
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
			//line = line + to_string(a[x].fitness) + "\n";
			aus << line;
			aus << a[x].fitness;
			aus << ",";
			aus << peakVolume(a[x]);
			aus << "\n";
			line = "";
		}
		// for(int i = 0; i < a[x].dim; i++)
		// {
		// 	line = line + to_string(a[x].center[i]) + ",";
		// }
		// for(int i = 0; i < a[x].dim; i++)
		// {
		// 	line = line + to_string(a[x].width[i]) + ",";
		// }
		// line = line + to_string(a[x].height) + ",";
		// line = line + to_string(a[x].alpha) + ",";
		// //line = line + to_string(a[x].fitness) + "\n";
		// aus << line;
		// aus << a[x].fitness;
		// aus << ",";
		// aus << peakVolume(a[x]);
		// aus << "\n";
		// line = "";
	}
	aus.close();
}

void print2HillGenestoFile(vector <gene2hills> &a, string filename)
{
    fstream aus;  //output file
    aus.open(filename,ios::out);
    aus.precision(10);
	string line = "";
	for(int x = 0; x < a.size(); x++)
	{
		for(int i = 0; i < a[x].dim; i++)
		{
			line = line + to_string(a[x].hill1.center[i]) + ",";
		}
		for(int i = 0; i < a[x].dim; i++)
		{
			line = line + to_string(a[x].hill1.width[i]) + ",";
		}
		line = line + to_string(a[x].hill1.height) + ",";
		line = line + to_string(a[x].hill1.alpha) + ",";
		for(int i = 0; i < a[x].dim; i++)
		{
			line = line + to_string(a[x].hill2.center[i]) + ",";
		}
		for(int i = 0; i < a[x].dim; i++)
		{
			line = line + to_string(a[x].hill2.width[i]) + ",";
		}
		line = line + to_string(a[x].hill2.height) + ",";
		line = line + to_string(a[x].hill2.alpha) + ",";
		line = line + to_string(a[x].fitness) + "\n";
		aus << line;
		line = "";
	}
	aus.close();
}

