//#include <iostream>
//#include <fstream>
//#include <cstdio>
//#include <cstdlib>
//#include <cmath>
//#include <cstring>
//#include <vector>

//include the header for peakfields
#include "peaks.h"

/**********************************************************************/
//Lorentian and Guassian functions
double PF_Guass(int D, double * x,double *c,double *w){//Guassian in D dimensions at x
                                             //with center c and width w

double accu,dif;   //accumulator and difference calculation variable
int i;             //loop index

  accu=0.0; //zero the accumulator
  for(i=0;i<D;i++){//loop over dimensions
    dif=(x[i]-c[i])/w[i];  //normalized variable in dimension i
    accu+=(dif*dif);       //accumulate squared value
  }
  accu*=(-2.773);      //normalizing constant for half height
  return(exp(accu));   //return the height at x

}

double PF_Guass_quick(int D, vector <double> & x,std::vector <double> & c,std::vector <double> & w){//Guassian in D dimensions at x
                                             //with center c and width w

double accu,dif;   //accumulator and difference calculation variable
int i;             //loop index

  accu=0.0; //zero the accumulator
  for(i=0;i<D;i++){//loop over dimensions
    dif=(x[i]-c[i])/w[i];  //normalized variable in dimension i
    accu+=(dif*dif);       //accumulate squared value
  }
  accu*=(-2.773);      //normalizing constant for half height
  return(exp(accu));   //return the height at x

}

double PF_Lorentz(int D, double * x,double *c,double *w){//Lorentz in D dimensions at x
                                               //with center c and width w

double accu,dif;   //accumulator and difference calculation variable
int i;             //loop index

  accu=0.0; //zero the accumulator
  for(i=0;i<D;i++){//loop over dimensions
    dif=(x[i]-c[i])/w[i];  //normalized variable in dimension i
    accu+=(dif*dif);       //accumulate squared value
  }
  return(1/(4*accu+1));   //return the height at x

}

double PF_LorentzV2(int D, double * x,double *c,double *w){//Lorentz in D dimensions at x
                                               //with center c and width w

double accu,dif;   //accumulator and difference calculation variable
int i;             //loop index
double Ltotal = 1.0;

  accu=0.0; //zero the accumulator
  for(i=0;i<D;i++){//loop over dimensions
    dif=(x[i]-c[i])/w[i];  //normalized variable in dimension i
    accu+=(dif*dif);       //accumulate squared value
	Ltotal = Ltotal*(1/(4*accu+1));
  }
  return(Ltotal);   //return the height at x

}

double PF_LorentzV2_quick(int D, vector <double> & x, std::vector <double> & c,std::vector <double> & w){//Lorentz in D dimensions at x
                                               //with center c and width w

double accu,dif;   //accumulator and difference calculation variable
int i;             //loop index
double Ltotal = 1.0;

  accu=0.0; //zero the accumulator
  for(i=0;i<D;i++){//loop over dimensions
    dif=(x[i]-c[i])/w[i];  //normalized variable in dimension i
    accu+=(dif*dif);       //accumulate squared value
	Ltotal = Ltotal*(1/(4*accu+1));
  }
  return(Ltotal);   //return the height at x

}

//Combined hill with parameter alpha in D dim with point x, center c, width w
double PF_Dual(int D, double * x,double alpha,double *c,double *w){//morphed hill

double G,L; //values for the Gaussian and lorentian hills

  G=PF_Guass(D,x,c,w);    //Get the Guassian hill
  L=PF_LorentzV2(D,x,c,w);  //Get the Lorentzian hill
  return(alpha*L+(1-alpha)*G);  //compute and return value
  
}

//Combined hill with parameter alpha in D dim with point x, center c, width w
double PF_Dual_quick(int D, vector <double> & x,double alpha,std::vector <double> & c,std::vector <double> & w){//morphed hill
	
double G,L; //values for the Gaussian and lorentian hills

  G=PF_Guass_quick(D,x,c,w);    //Get the Guassian hill
  L=PF_LorentzV2_quick(D,x,c,w);  //Get the Lorentzian hill
  return(alpha*L+(1-alpha)*G);  //compute and return value
  
}

//constructors and destructors
peakfield::peakfield(){//default constructor

  dimv=0;  //signal that the peakfield is unallocated
  
}

peakfield::peakfield(int d,int n){//create a peakfield, dimension d, <=n peaks

  create(d,n);  //call the creation method

}

peakfield::peakfield(const peakfield &w){//copy constructor

   dimv=0;   //mark the peakfield as unallocated
   copy(w);  //copy it

}

peakfield::~peakfield(){//the destructor

  destroy();
  
}

//creation and deletion
void peakfield::create(int d,int n){//allocate the space and stuff

int i; //loop index

  if(dimv!=0)destroy();//Don't allocate over an existing structure
  dimv=d;  //assign the dimension
  npks=n;  //assign the maximum number of peaks
  cpks=0;  //start with no peaks
  actv=new int[npks];     //create the active peak boolean
  typv=new double[npks];  //create type vector
  hgt=new double[npks];   //create the height register
  ctr=new double*[npks];  //create the spine of the center array
  for(i=0;i<npks;i++)ctr[i]=new double[dimv];  //allocate center storage
  whh=new double*[npks];  //create the spine of the width array
  for(i=0;i<npks;i++)whh[i]=new double[dimv];  //allocate center storage

}

void peakfield::destroy(){//deallocate everything

int i; //loop index

  if(dimv!=0){//if it was allocated
    for(i=0;i<npks;i++)delete [] whh[i];   //delete width vectors
    delete [] whh;                         //delete width spine
    for(i=0;i<npks;i++)delete [] ctr[i];   //delete center vectors
    delete [] ctr;                         //delete center spine
    delete [] hgt;                         //delete heights
    delete [] typv;                        //delet peak type
    delete [] actv;                        //delete active variables

    //Everything's deleted so...
    dimv=0;                 //note no longer allocated
  }

}

void peakfield::copy(const peakfield &w){//make a copy

int i,j;

  if(dimv!=0){//if previously allocated
    if((dimv!=w.dimv)||(npks!=w.npks))destroy(); //destroy if wrong size
  }

  if(dimv==0){//if currently not allocated
    create(w.dimv,w.npks);   //allocate
  }

  cpks=w.cpks;  //copy number of peaks
  for(i=0;i<npks;i++){//loop over peaks
    actv[i]=w.actv[i];  //copy active status
    typv[i]=w.typv[i];  //copy peak type
    hgt[i]=w.hgt[i];    //copy
    for(j=0;j<dimv;j++){//loop over coordinates
      ctr[i][j]=w.ctr[i][j];  //copy center
      whh[i][j]=w.whh[i][j];  //copy width at half height
    }
  }
  
}

//useage

//Add returns 1 if it adds the peak and zero if it fails to add it
int peakfield::add(double *cntr,double *wdth,double pkh,double typ){//add a peak

int i; //loop index

 if(cpks<npks){//room?
    actv[cpks]=1;          //start the peak in its active state
    if(typ<=0)typ=0.0; else if(typ>=1.0)typ=1.0;  //force possible type
    typv[cpks]=typ;        //record type value
    hgt[cpks]=pkh;         //record peak type
    for(i=0;i<dimv;i++){//loop over coordinates
      ctr[cpks][i]=cntr[i];  //assign center
      whh[cpks][i]=wdth[i];  //assign widths
    }//done looping over coordinates
    cpks++;                //note new peak
    return(1);             //report peak was added
  } else return(0);        //report peak was not added

}

double dsqr(double *p,double *q,int dimv){//squared distance between centers

double delta,ttl;   //the difference and the total

  ttl=0.0;  //zero the accumulator
  for(int i=0;i<dimv;i++){//loop over coordinates
    delta=p[i]=q[i];  //find the difference
    ttl+=delta*delta; //accumulate the squared distance
  }
  return(ttl);  //return the total

}

int peakfield::closest(double *cntr){//find active peak closest to the center

int i,b;        //loop index, best pointer
double rd,nd;   //reference distance, new distance
 
  if(dimv!=0){//if the peakfield was implimented
    b=0;while((actv[b]==0)&&(b<cpks))b++;  //find first active peak
    if(b<cpks){//there WAS an active peak
      rd=dsqr(ctr[b],cntr,dimv);  //compute the initial reference distance
      for(i=b+1;i<cpks;i++)if(actv[i]==1){//loop over active peaks
	  nd=dsqr(ctr[i],cntr,dimv);  //get the distance to the trial peak
	if(nd<rd){//Better?
          b=i;    //assign new best peak index
	  rd=nd;  //update reference distance
	}
      }
    } else return(-1);  //no available active peak
  } else return(-1); //no available active peak
  return(0);
}

void peakfield::killpeak(int q){//inactivate peak q

  if((dimv==1)&&(q>=0)&&(q<cpks))actv[q]=0;  //turn off the peak

}

void peakfield::respeak(int q){//restore peak

  if((dimv==1)&&(q>=0)&&(q<cpks))actv[q]=1;  //turn on the peak

}

double peakfield::samfield(std::vector <double> pnt)
{//value of the peakfield at pnt
    int i;          //loop index
    double accu;    //value accumulator
    double GV,LV;   //Guassian and Lorentzian Values
    double * pnt2 = nullptr;
	pnt2 = (double *)malloc(pnt.size()*sizeof(double));
	// while(pnt2 == nullptr)
// 	{
// 		pnt2 = (double *)malloc(pnt.size()*sizeof(double));
// 	}
    for(int x = 0; x < pnt.size(); x++)
    {
        pnt2[x] = pnt[x];
    }
    accu=0.0;  //zero the accumulator
    if(dimv!=0){//if allocated
        for(i=0;i<cpks;i++)if(actv[i])
        {//loop over active peaks
            GV=PF_Guass(dimv,pnt2,ctr[i],whh[i]);    //get the Gaussian value
            LV=PF_Lorentz(dimv,pnt2,ctr[i],whh[i]);  //get the Lorentzian value
            accu+=hgt[i]*(typv[i]*LV+(1-typv[i])*GV);
        }
    }
    free(pnt2);
    return(accu);  //return the total of the active peaks
}

//I/O
void peakfield::write(ostream &aus){//write a peakfield description

int i,j; //loop index variables

 
  aus << dimv << endl;              //dimenasion
  if(dimv==0)return;                //not allocated!
  aus << npks << " -max" << endl;   //number of peaks, max
  aus << cpks << " -peaks" << endl; //number of peaks, extant
  for(i=0;i<cpks;i++){//loop over the peaks
    aus << actv[i] << " -on?" << endl;    //activity status
    aus << typv[i] << " -type" << endl;   //peak mixture type
    aus << hgt[i] << " -height" << endl;  //height
    for(j=0;j<dimv;j++){//loop over coordinates
      aus << ctr[i][j] << endl;   //center coordinate
      aus << whh[i][j] << endl;   //width at half height
    }
  }
}

void peakfield::read(istream &aus){//read a peakfield description

}
