
#include <fstream>
#include <iostream>
#include "time.h"
#include "math.h"
#include "qgstack.hh"
#ifndef QGSTRING_H
#define QGSTRING_H

class qgstring 
{
 friend class qgstack;
 static  const int Ngrid=3000;
 static  const double acc=1.0E-10;
 static const double alp = 1.;
 static const double gamma = 1.;
 public :
   qgstring(double *xio, double* dpio, int* N);//construction from external array
   qgstring(); // zero constructor
   qgstring(const qgstring & copy); // copy constructor
   qgstring(double s1,double s2,const qgstring & copy); // constructor for splitting
   qgstring(double p1, double p2, double bm); // constructor from p1+p2+impact_parameter
   qgstring(double p1, double p2, double x11, double x12, double x21, double x22); // p1+p2+coordinates of endpoints
   ~qgstring();
//   void Initialize(double * xio, double * dpio, int * N);
//   void GetX (double t, double * xout);
//   void GetdP(double t, double * dpout);
     void WriteX(char file[20]);
     void WriteX(); // Writes coordinates at zero time
     void WriteDir(); // Writes directrix
     double GetMass2(); // Invariant mass of the string squared in GeV^2
     double GetMass2(int N1, int N2); // Invariant mass of the string part from N1 to N2 points of the sigma grid, in GeV^2
     double GetMass2(double s1, double s2 ); // invariant mass as a function of sigma
     int GetNgrid(); // Returns number of points in the grid - 1 (The above Ngrid number)
     void TimeShift();//Evolvs directrix by time = 1 step
     void TimeShift(int Ntime); // Evolvs string by Ntime steps in time and shifts the time to zero. Sets directrix and momenta of the parts.  
     void TimeShift(double tfm); // Evolvs string by time tfm measured in fm
     void LR(); // switches ends (directrix <-> antidirectrix)
     void Boost3(double bb); // Boosts string directrix along x3 with velocity beta
     double GetBeta3(); // Returns string velocity in the lab frame, beta, assumed to be along x3
     double GetXi(int k, int i);
//     double GetdPi(int k, int i);
     double GetMomentum(int k, int N1, int N2);
     double GetMomentumLab(int k, int N1, int N2);    
     double GetSigmagrid (int i);
//     int GetDecayPoint(double lowlim);
//     int InstantGetDecayPoint(double lowlim);
     double GetSigmaDecay(double lowlim);
 private :
   double xdot[4][Ngrid+1]; double xprime[4][Ngrid+1]; double dpi[4][Ngrid]; double yd[4][2*Ngrid+1];
           double sigmagrid[Ngrid+1];
   /* x\dot, x\prime  -- time and sigma derivatives, yd -- directrix */
   double * xio; double* dpio; int* N;// double gamma; 
   double sigmat;
   /* gamma -- tension; sigmat -- total length*/
   double ltime; /*local time*/
   double beta; /*velocity along x3 in the lab frame*/
   double tmin; /*additional parameter to cut of absolute past after boost*/
   qgstring * NextPtr;
};

#endif
