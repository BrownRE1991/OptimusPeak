/*
 *   This code generates fields of peaks with mixed
 *   Lauentzian/Guassian hills at user specified locations.
 *
 *   A peak is specified by a center, a vector of widths at half
 *   height, a height, and a type.  The type is Gaussian if the
 *   parameters is zero, Lorentsian if the parameter is 1 and a blend
 *   according to the type otherwise.
 *
 */
 
#include <fstream>
#include <cmath>
#include <vector>

#ifndef	_PKFLD_H
#define	_PKFLD_H

using namespace std;

double PF_Dual(int D, double * x, double alpha,double *c,double *w); //morphed hill
double PF_Dual_quick(int D, vector <double> & x,double alpha,std::vector <double> & c,std::vector <double> & w);


class peakfield {

 public:

  //constructors and destructors
  peakfield();                    //default constructor
  peakfield(int d,int n);         //create a peakfield, dimension d, <=n peaks
  peakfield(const peakfield &w);  //copy constructor
  ~peakfield();                   //the destructor

  //creation and deletion
  void create(int d,int n);         //allocate the space and stuff
  void destroy();                   //deallocate everything
  void copy(const peakfield &w);    //deallocate everything

  //useage
  int add(double *cntr,double *wdth,double pkh,double typ); //add a peak
  int closest(double *cntr);           //find a peak closest to the center
  void killpeak(int q);                //inactivate peak q
  void respeak(int q);                 //restore peak                
  double samfield(std::vector <double>);        //value of the peakfield at pnt

  //I/O
  void write(ostream &aus);       //write a peakfield description
  void read(istream &aus);        //read a peakfield description

 private:

  int dimv;       //dimension of the data space   0=unallocate que
  int npks;       //maximum number of peaks
  int cpks;       //current number of peaks
  int *actv;      //active peak?  0=no 1=yes -1=never allocated
  double *typv;   //peak type
  double *hgt;    //peak heights
  double **ctr;   //peak centers
  double **whh;   //width at half height
  

};

#endif /* _PKFLD_H */
