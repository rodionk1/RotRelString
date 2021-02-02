#include "qgstring.hh"
#include "assert.h"
#include <iostream>

//DESTRUCTOR
qgstring::~qgstring()
{
};


//CONSTRUCTOR, RELOADED

//zero constructor
qgstring::qgstring()
{NextPtr=0;
   beta = 0;
      ltime = 0;
      tmin = 0.;
}

//constructor from p1+p2+impact parameter


qgstring::qgstring(double p1, double p2, double bm)
{
beta = 0;
   ltime = 0;
   tmin = 0.;
   //   gamma=1.; 
   
   sigmat=0.;
      NextPtr=0;   
// RAPIDITY OF STRING CM
  double beta = (p1+p2)/(fabs(p1)+fabs(p2));
  double ystring=0.5*log((1+beta)/(1-beta));
   std::cout << "string rapidity " << ystring << std::endl; 

// STRING ENERGY and SIGMA_TOT
  double Etot = fabs(p1)+fabs(p2);
  double sigmatot = Etot/gamma;
  double dsigma = sigmatot/Ngrid;
   std::cout << "sigma (total energy) " << sigmatot << std::endl;

   for (int i=0; i<=Ngrid; i++)
     {
	sigmagrid[i]=dsigma*i;
     };
   
   
// STRING ENERGY IN CMS (INVARIANT MASS) and SIGMA_MAX(CMS)   
  double Etotc = sqrt(2.*fabs(p1*p2)-2.*p1*p2);
  double sigmam = Etotc/gamma; 
 //  std::cout << "sigma in cms " << sigmam << std::endl;
   
// PARAMETER mpi, Lambert's W-function
  double x = (-0.5*bm/sigmam);
 // double mpi = -0.5*bm/(0.665*(1.+0.195*log(x+1.))*log(x+1.)+0.04);
   double al = -2.*(x-x*x+3*x*x*x/2-8*x*x*x*x/3+125*x*x*x*x*x/24)/bm;    
//   std::cout << "parameter al, preliminary "<< al << std::endl;
   
   double alt=-1.;
   
   while(fabs(fabs(al-alt)/alt)>1.e-9){ alt=al;
      al=2*log(alt*sigmam+exp(-0.5*alt*bm))/bm;
     };
//   std::cout <<"parameter al, after iterations " << al <<std::endl;
 
   
/************ MAKING DIRECTRIX *************/   

  /** temporary directrix - high precision 5x grid. **/
     
//double ydt[4][50*2*Ngrid+1];

   int prec = 500;
   double dsigmat = dsigma/prec;
  
   
   
    yd[1][Ngrid]=0.;
    yd[2][Ngrid]=0.5*bm;
    yd[3][Ngrid]=0.;
    //  double cc= 0.0308986;
    //   double al = log(sigmam/(2.*cc)+sqrt(pow(sigmam/(2*cc),2)+1))/(0.5*bm);
   
    //  yd[i][Ngrid]=xiold[i][0];
     for (int j=1; j<=Ngrid;j++)
               {  double bb = 0.5*(yd[2][Ngrid+j-1]+yd[2][Ngrid-j+1]);
for (int i=1; i<4; i++)
		    { 
		       
		  yd[i][Ngrid+j] = yd[i][Ngrid+j-1];
		  yd[i][Ngrid-j] = yd[i][Ngrid-j+1];
		    };
     //           xpp[2]=-1./cosh(2*log(sigmam/mpi)*bb/bm + ystring);
     //           xdd[3]=tanh(2*log(sigmam/mpi)*bb/bm + ystring);
                  for (int k=1;k<=prec; k++)
		  {		     
		    double xpp[4] = {0.,0.,0.,0.};
		    double xdd[4] = {0.,0.,0.,0.};
		       
          	  bb = 0.5*(yd[2][Ngrid+j]+yd[2][Ngrid-j]);     
		  xpp[2]=-1./cosh(al*bb+ystring);
                  bb+=0.5*dsigmat*xpp[2];
		  xpp[2]=-1./cosh(al*bb+ystring);
		  xdd[3]=tanh(al*bb+ystring);
     //           std::cout << bb <<" "<< xpp[2] << "  " << xdd[3] << std::endl;
                 for (int i=1;i<4;i++)
                   {
                  yd[i][Ngrid+j] +=  xpp[i]*dsigmat + xdd[i]*dsigmat;
                  yd[i][Ngrid-j] +=  xpp[i]*dsigmat - xdd[i]*dsigmat;
                    };
		    };
		  
            };
     
    
   
   yd[0][Ngrid]=0.;
   for (int j=1;j<=Ngrid;j++)
     {
   yd[0][Ngrid+j]=yd[0][Ngrid+j-1]+dsigma;
   yd[0][Ngrid-j]=yd[0][Ngrid-j+1]-dsigma;	
     };
   
  for (int i=0; i<Ngrid; i++)
   {
      dpi[0][i]=gamma*(sigmagrid[i+1]-sigmagrid[i]);
   };
   
   
/* Reworking table of momenta*/
   
   
  for (int k=0; k<4; k++)
     {	
	        for (int i = 0; i<Ngrid; i++)
	  {	     
	                  dpi[k][i]=0.5*gamma* ( yd[k][Ngrid+i+1] + yd[k][Ngrid-i] - yd[k][Ngrid-i-1] - yd[k][Ngrid+i] );
	  };
     };
   
   
};



// constructor from p1+p2+coordinates of endpoints


qgstring::qgstring(double p1, double p2, double x11, double x12, double x21, double x22)
{
//      gamma=1.; 
      
   beta = 0;
      ltime = 0;
      tmin = 0.;
   
   
      sigmat=0.;
      NextPtr=0;   
// RAPIDITY OF STRING CM
  double bm=sqrt((x11-x21)*(x11-x21)+(x12-x22)*(x12-x22));
   std::cout << "impact parameter" << bm << std::endl;
   
  double beta = (p1+p2)/(fabs(p1)+fabs(p2));
  double ystring=0.5*log((1+beta)/(1-beta));
   std::cout << "string rapidity " << ystring << std::endl; 

// STRING ENERGY and SIGMA_TOT
  double Etot = fabs(p1)+fabs(p2);
  double sigmatot = Etot/gamma;
  double dsigma = sigmatot/Ngrid;
   std::cout << "sigma (total energy) " << sigmatot << std::endl;

   for (int i=0; i<=Ngrid; i++)
     {
	sigmagrid[i]=dsigma*i;
     };
   
   
// STRING ENERGY IN CMS (INVARIANT MASS) and SIGMA_MAX(CMS)   
  double Etotc = sqrt(2.*fabs(p1*p2)-2.*p1*p2);
  double sigmam = Etotc/gamma; 
  std::cout << "sigma in cms " << sigmam << std::endl;
   
// PARAMETER mpi, Lambert's W-function
  double x = (-0.5*bm/sigmam);
 // double mpi = -0.5*bm/(0.665*(1.+0.195*log(x+1.))*log(x+1.)+0.04);
   double al = -2.*(x-x*x+3*x*x*x/2-8*x*x*x*x/3+125*x*x*x*x*x/24)/bm;    
   std::cout << "parameter al, preliminary "<< al << std::endl;
   
   double alt=-1.;
   
   while(fabs(fabs(al-alt)/alt)>1.e-9){ alt=al;
      al=2*log(alt*sigmam+exp(-0.5*alt*bm))/bm;
//      std::cout << al << "after iteration "<<std::endl;
   };
   std::cout <<"parameter al, after iterations " << al <<std::endl;
 
   
/************ MAKING DIRECTRIX *************/   

  /** temporary directrix - high precision 5x grid. **/
     
//double ydt[4][50*2*Ngrid+1];

   std::cout << "making dir" <<std::endl;
   
   int prec = 500;
   double dsigmat = dsigma/prec;
  
   
   
    yd[1][Ngrid]=0.;
    yd[2][Ngrid]=0.5*bm;
    yd[3][Ngrid]=0.;
    //  double cc= 0.0308986;
    //   double al = log(sigmam/(2.*cc)+sqrt(pow(sigmam/(2*cc),2)+1))/(0.5*bm);
   
    //  yd[i][Ngrid]=xiold[i][0];
     for (int j=1; j<=Ngrid;j++)
               {  double bb = 0.5*(yd[2][Ngrid+j-1]+yd[2][Ngrid-j+1]);
for (int i=1; i<4; i++)
		    { 
		       
		  yd[i][Ngrid+j] = yd[i][Ngrid+j-1];
		  yd[i][Ngrid-j] = yd[i][Ngrid-j+1];
		    };
     //           xpp[2]=-1./cosh(2*log(sigmam/mpi)*bb/bm + ystring);
     //           xdd[3]=tanh(2*log(sigmam/mpi)*bb/bm + ystring);
                  for (int k=1;k<=prec; k++)
		  {		     
		    double xpp[4] = {0.,0.,0.,0.};
		    double xdd[4] = {0.,0.,0.,0.};
		       
          	  bb = 0.5*(yd[2][Ngrid+j]+yd[2][Ngrid-j]);     
		  xpp[2]=-1./cosh(al*bb+(p1/fabs(p1))*ystring);
                  bb+=0.5*dsigmat*xpp[2];
		  xpp[2]=-1./cosh(al*bb+(p1/fabs(p1))*ystring);
		  xdd[3]= (p1/fabs(p1))*tanh(al*bb+(p1/fabs(p1))*ystring);
     //           std::cout << bb <<" "<< xpp[2] << "  " << xdd[3] << std::endl;
                 for (int i=1;i<4;i++)
                   {
                  yd[i][Ngrid+j] +=  xpp[i]*dsigmat + xdd[i]*dsigmat;
                  yd[i][Ngrid-j] +=  xpp[i]*dsigmat - xdd[i]*dsigmat;
                    };
		    };
		  
            };
     
    
   
   yd[0][Ngrid]=0.;
   for (int j=1;j<=Ngrid;j++)
     {
   yd[0][Ngrid+j]=yd[0][Ngrid+j-1]+dsigma;
   yd[0][Ngrid-j]=yd[0][Ngrid-j+1]-dsigma;	
     };

   std::cout << "Dir done, now rotation" << std::endl;
   
/* rotating and moving to the right position */
   /* Rotation matrix */
   
 double omega[3][3];  
   
   double cphi = (x12-x22)/bm;
   double sphi = (x11-x21)/bm;
   
   omega[1][1] = cphi;
   omega[2][2] = cphi;
   omega[1][2] = sphi;
   omega[2][1] = -sphi;
   
   std::cout << omega[1][1] << " " <<omega [1][2] << std::endl;
   std::cout << omega[2][1] << " " <<omega [2][2] << std::endl;
   
   for (int i=0; i<2*Ngrid+1; i++)
     {
	double ydold[3];
	ydold[1] = yd[1][i];
	ydold[2] = yd[2][i];	
for (int k=1; k<=2; k++)
	  {
	     yd[k][i]=omega[k][1]*ydold[1]+omega[k][2]*ydold[2];
	  };	
     };
/* rotation done, now the translation */   
   
   double P[3];
   P[1] = x11 - yd[1][Ngrid];
   P[2] = x12 - yd[2][Ngrid];
   
   for (int i=0; i<2*Ngrid+1; i++)
     {
	for(int k=1; k<=2; k++)
	  {
	   yd[k][i]+=P[k];  
	  };	
     };
   
   
   
  for (int i=0; i<Ngrid; i++)
   {
      dpi[0][i]=gamma*(sigmagrid[i+1]-sigmagrid[i]);
   };
   
   
/* Reworking table of momenta*/
   
   
  for (int k=0; k<4; k++)
     {	
	        for (int i = 0; i<Ngrid; i++)
	  {	     
	                  dpi[k][i]=0.5*gamma* ( yd[k][Ngrid+i+1] + yd[k][Ngrid-i] - yd[k][Ngrid-i-1] - yd[k][Ngrid+i] );
	  };
     };
   
   
};







qgstring::qgstring(double s1,double s2,const qgstring & copy) // constructor for splitting
{
   
   beta = copy.beta;
      ltime = copy.ltime;
      tmin = copy.tmin;
   
   
   NextPtr=0;
//   gamma=1.;
//copying directrix parts and making grids of sigma   
   double yp[4], ym[4],y0[4],dydp[4],dydm[4];
   double dyp[4][Ngrid+1], dym[4][Ngrid+1];
   for (int k=1; k<4; k++)
     { int i1 = (int) floor(fabs(s1-1.e-10)*Ngrid/copy.sigmagrid[Ngrid]);
	yp[k]=(copy.yd[k][Ngrid+i1]*(copy.sigmagrid[i1+1]-s1)+copy.yd[k][Ngrid+i1+1]*(s1-copy.sigmagrid[i1]))/(copy.sigmagrid[i1+1]-copy.sigmagrid[i1]);
	ym[k]=(copy.yd[k][Ngrid-i1]*(copy.sigmagrid[i1+1]-s1)+copy.yd[k][Ngrid-i1-1]*(s1-copy.sigmagrid[i1]))/(copy.sigmagrid[i1+1]-copy.sigmagrid[i1]);
        y0[k]=0.5*(yp[k]+ym[k]);
//	std::cout << s1<<" " << yp[k] <<" " <<ym[k] << std::endl;
      yd[k][Ngrid]=y0[k];
     };
   
   
   //   yd[0][Ngrid]=copy.yd[0][Ngrid];
   if (fabs(beta)>0.001) yd[0][Ngrid]=copy.yd[0][Ngrid]; 
   if (fabs(beta)<=0.001)  yd[0][Ngrid]=0;
   for (int i=1; i<=Ngrid; i++) {
      sigmagrid[i]= ((s2-s1)*i)/(Ngrid*1.0);
   //and directrix...
 int i1 = (int) floor(fabs(s1+sigmagrid[i]-1.e-15)*Ngrid/copy.sigmagrid[Ngrid]);      
// int i2 = (int) floor(fabs(s1+sigmagrid[i-1]-1.e-10)*Ngrid/copy.sigmagrid[Ngrid]);     
 for (int k=1; k<4; k++)
{ double s = s1+sigmagrid[i];
   
yd[k][Ngrid+i]=(copy.yd[k][Ngrid+i1]*(copy.sigmagrid[i1+1]-s)+copy.yd[k][Ngrid+i1+1]*(s-copy.sigmagrid[i1]))/(copy.sigmagrid[i1+1]-copy.sigmagrid[i1]);
yd[k][Ngrid-i]=(copy.yd[k][Ngrid-i1]*(copy.sigmagrid[i1+1]-s)+copy.yd[k][Ngrid-i1-1]*(s-copy.sigmagrid[i1]))/(copy.sigmagrid[i1+1]-copy.sigmagrid[i1]);    
yd[k][Ngrid+i]-=yp[k]; yd[k][Ngrid+i]+=y0[k];
yd[k][Ngrid-i]-=ym[k]; yd[k][Ngrid-i]+= y0[k];

   
/*   dydp[k]=(copy.yd[k][Ngrid+i1]*(copy.sigmagrid[i1+1]-s)+copy.yd[k][Ngrid+i1+1]*(s-copy.sigmagrid[i1]))/(copy.sigmagrid[i1+1]-copy.sigmagrid[i1]);
   dydp[k]-=(copy.yd[k][Ngrid+i2]*(copy.sigmagrid[i2+1]-s)+copy.yd[k][Ngrid+i2+1]*(s-copy.sigmagrid[i2]))/(copy.sigmagrid[i2+1]-copy.sigmagrid[i2]);
   dydm[k]=(copy.yd[k][Ngrid-i1]*(copy.sigmagrid[i1+1]-s)+copy.yd[k][Ngrid-i1-1]*(s-copy.sigmagrid[i1]))/(copy.sigmagrid[i1+1]-copy.sigmagrid[i1]);
   dydm[k]-=(copy.yd[k][Ngrid-i2]*(copy.sigmagrid[i2+1]-s)+copy.yd[k][Ngrid-i2-1]*(s-copy.sigmagrid[i2]))/(copy.sigmagrid[i2+1]-copy.sigmagrid[i2]);
*/
 /*  dyp[k][i]=yd[k][Ngrid+i]-yd[k][Ngrid+i-1];
   dym[k][i]=yd[k][Ngrid-i]-yd[k][Ngrid-i+1];
*/
 };
//      std::cout << "* "<<yd[3][Ngrid+i] << " " <<yd[3][Ngrid-i]<< std::endl;
       
      yd[0][Ngrid+i]=yd[0][Ngrid]+sigmagrid[i];
      yd[0][Ngrid-i]=yd[0][Ngrid]-sigmagrid[i];
   };
 /*  
   for (int i=1; i<=Ngrid; i++)
     {
              double dypn = sqrt(dyp[1][i]*dyp[1][i]+dyp[2][i]*dyp[2][i]+dyp[3][i]*dyp[3][i]);
	      double dymn = sqrt(dym[1][i]*dym[1][i]+dym[2][i]*dym[2][i]+dym[3][i]*dym[3][i]);
	      double ds=sigmagrid[i]-sigmagrid[i-1];
	for (int k=1; k<4; k++)
	  { dyp[k][i]=dyp[k][i]*ds/dypn;
	    dym[k][i]=dym[k][i]*ds/dymn;
	  
	yd[k][Ngrid+i]=yd[k][Ngrid+i-1]+dyp[k][i];
	yd[k][Ngrid-i]=yd[k][Ngrid-i+1]+dym[k][i];
	  };
     };
   
   */
   
//   std::cout << "***" <<std::endl;
   
}



//copy constructor
qgstring::qgstring(const qgstring & copy)
{
   
   beta = copy.beta;
      ltime = copy.ltime;
      tmin = copy.tmin;
   
   
   NextPtr=copy.NextPtr;
//gamma=1.;
   for(int i=0; i<Ngrid; i++)
       {
	  for(int j=0; j<=3; j++)
	    {
	     xdot[j][i]=copy.xdot[j][i];
	     xprime[j][i]=copy.xprime[j][i];
	      dpi[j][i]=copy.dpi[j][i];
	    }	  
	  sigmagrid[i]=copy.sigmagrid[i];
       };
   for(int i=0; i<=3; i++)
     {
	xdot[i][Ngrid]=copy.xdot[i][Ngrid];
        xprime[i][Ngrid]=copy.xprime[i][Ngrid];	
     };
   sigmagrid[Ngrid]=copy.sigmagrid[Ngrid];
       
   for(int i=0; i<=2*Ngrid; i++)
     {
	for(int j=0; j<4; j++)
	  {	     
	yd[j][i]=copy.yd[j][i];
	  };
	
     };
   for (int i=0; i<Ngrid; i++)
     { 
	dpi[0][i]=gamma*(sigmagrid[i+1]-sigmagrid[i]);
     };   
   
}


  
//constructor from external arrays  
qgstring::qgstring(double* xxi, double* dppi, int * NN)
  : xio(xxi), dpio(dppi), N(NN) 
{  
   beta = 0;
      ltime = 0;
      tmin = 0.;
   
   
   NextPtr=0;
//   gamma=1.; 
   sigmat=0.;
/* Computing xdot and xprime in old grid*/
/* Grids for different variables, notations
 x, sigma
  0    1   2   ...... N-1   N
 dp,x', xdot
    0    1   3           N-1     
 */   
   
   double xdold[4][*N];
   double xpold[4][*N];
   double xiold[4][*N+1];
   double dpiold[4][*N];
   double sigmaold[*N+1];
   for (int j=0; j<4; j++)
     {
	for (int i=0; i<=*N; i++)
	  {
           xiold[j][i]=*xio;
//	     std::cout << xiold[j][i] <<" ";
	   xio++;
	  }; 
//	std::cout << '\n';
     };

   for (int j=0; j<4; j++)
     {
	for (int i=0; i<*N; i++)
	  {
           dpiold[j][i]=*dpio;
	   dpio++;
	  };
     };
   sigmaold[0]=0.0;
    for (int j=0; j<*N; j++)
     {
	double dx=sqrt((xiold[1][j]-xiold[1][j+1])*(xiold[1][j]-xiold[1][j+1])+(xiold[2][j]-xiold[2][j+1])*(xiold[2][j]-xiold[2][j+1])+(xiold[3][j]-xiold[3][j+1])*(xiold[3][j]-xiold[3][j+1]));
	double dp=sqrt((dpiold[1][j])*(dpiold[1][j])+(dpiold[2][j])*(dpiold[2][j])+(dpiold[3][j])*(dpiold[3][j]));
	std::cout << dx << " " << dp << std::endl;
	double norm=0.0;
	for (int i=1; i<4; i++)
	  {
           xdold[i][j]=dpiold[i][j]/(sqrt(gamma*gamma*dx*dx+dp*dp));
	   xpold[i][j]=(xiold[i][j+1]-xiold[i][j])*gamma/(sqrt(gamma*gamma*dx*dx+dp*dp));
	   norm+=pow(xdold[i][j],2)+pow(xpold[i][j],2);
//	     std::cout << xdold[i][j]*1.e30 << " " << xpold[i][j]*1.e30 << " ";
	  };
//	std::cout << "norm = " << norm << std::endl;
// ensuring that xdot^2 + xprime^2 =1 
        for (int i=1;i<4;i++)
	  { 
	     xdold[i][j]/=sqrt(norm);
	     xpold[i][j]/=sqrt(norm);
	  };
	
	 std::cout << sqrt(gamma*gamma*dx*dx+dp*dp)/gamma <<'\n';
	sigmat+=sqrt(gamma*gamma*dx*dx+dp*dp)/gamma;
	sigmaold[j+1]=sigmat;
     };

/* Bringing to new grid in sigma*/

/*left point-> sigma_a <=> sigmagrid[i] ; xa   */
/*right point-> sigma_b <=> sigmagrid[i+1] ; xb */
   int num1=0, num2=0;
   sigmagrid[0]=0.;
   double dsigma=sigmat/Ngrid;
   
   //std::cout << dsigma << '\n';
   
  double xa[4], xb[4], dpab[4];
   xb[1]=xiold[1][0];
   xb[2]=xiold[2][0];
   xb[3]=xiold[3][0];
   xa[1]=xb[1]; xa[2]=xb[2]; xa[3]=xb[3];
   for(int i=0;i<Ngrid-1;i++)
     {	
     for(int j=1;j<4;j++)
     {
//	sigmagrid[i+1]=sigmagrid[i]+dsigma;
     sigmagrid[i+1]=sigmat*(i+1)/Ngrid;
     num2=num1;
     while(!((sigmaold[num2]<sigmagrid[i+1])&(sigmaold[num2+1]>=sigmagrid[i+1]))) num2++;
//	xa[j]=xiold[j][num1]+(sigmagrid[i]-sigmaold[num1])*xpold[j][num1];
	xb[j]=xiold[j][num2]+(-sigmaold[num2]+sigmagrid[i+1])*xpold[j][num2];
	dpab[j]=xdold[j][num1]*(sigmagrid[i+1]-sigmagrid[i]);
	dpab[j]+=(xdold[j][num1]*(sigmaold[num1+1]-sigmagrid[i+1])-xdold[j][num2]*(sigmaold[num2+1]-sigmagrid[i+1]));
	for(int k=num1+1; k<num2+1; k++)
	  {
	     dpab[j]+=xdold[j][k]*(sigmaold[k+1]-sigmaold[k]);
	  }; /* k */
	dpab[j]*=gamma;
/*	if (dpab[j]!=0.)
	  {
	     std::cout << sigmaold[num1+1]<<" "<<sigmagrid[i]<<" "<<sigmaold[num2+1]<<" "<<sigmagrid[i+1]<<'\n';
	  }; */
	
	dpi[j][i]=dpab[j]; /* table of dpi for new sigma grid */
     }; /* j */

	/*****  Corrections to satisfy xdot \perp xprime and xdot^2+xprime^2=1*****/
	/** conserving momentum and energy */
	double dxab[4]={0,xb[1]-xa[1],xb[2]-xa[2],xb[3]-xa[3]};
	
	
	double dx=sqrt((xa[1]-xb[1])*(xa[1]-xb[1])+(xa[2]-xb[2])*(xa[2]-xb[2])+(xa[3]-xb[3])*(xa[3]-xb[3]));
	double dp=sqrt((dpab[1])*(dpab[1])+(dpab[2])*(dpab[2])+(dpab[3])*(dpab[3]));

	double dpdx=dpab[1]*(xb[1]-xa[1])+dpab[2]*(xb[2]-xa[2])+dpab[3]*(xb[3]-xa[3]);
       // for (int j=1; j<4; j++) dxab[j]-=dpab[j]*dpdx/(dp*dp); // now momentum, energy and mass is conserved
//	for (int j=1; j<4; j++) dpab[j]-=(xb[j]-xa[j])*dpdx/(dx*dx);

/* NO CORRECTIONS: LOSS OF ACCURACY	
	if (gamma*dx>=1.e-10) {for (int j=1; j<4; j++) dpab[j]-=(xb[j]-xa[j])*dpdx/(dx*dx);}else{for (int j=1; j<4; j++) dxab[j]=0;}; //-=dpab[j]*dpdx/(dp*dp);};
	 dp=sqrt((dpab[1])*(dpab[1])+(dpab[2])*(dpab[2])+(dpab[3])*(dpab[3]));
	 dx=sqrt(pow(dxab[1],2)+pow(dxab[2],2)+pow(dxab[3],2));
*/		
        /****** now dpab \perp dxab ****/
	
   //	std::cout << "dsigma " << (sqrt(gamma*gamma*dx*dx+dp*dp)) << " "<<dsigma<<'\n';
        
/*	for (int j=1; j<4; j++)
	  {
           xdot[j][i]=dpab[j]/(sqrt(gamma*gamma*dx*dx+dp*dp));
	   xprime[j][i]=dxab[j]*gamma/(sqrt(gamma*gamma*dx*dx+dp*dp));  
     	  }; 
*/
      
	for (int j=1; j<4; j++)
	  {
           xdot[j][i]=dpab[j]/(gamma*dsigma);
	   xprime[j][i]=dxab[j]/dsigma;  
     	  }; 
//	std::cout << i << " " << xdot[3][i] << std::endl;
	
	/*****  now xdot^2+xprime^2=1  *****/
  //      std::cout<< xdot[1][i]*xprime[1][i]+xdot[2][i]*xprime[2][i]+xdot[3][i]*xprime[3][i] <<" ";
 //	std::cout << xdot[1][i]*xdot[1][i]+xprime[1][i]*xprime[1][i]+ xdot[2][i]*xdot[2][i]+xprime[2][i]*xprime[2][i]+ xdot[3][i]*xdot[3][i]+xprime[3][i]*xprime[3][i]<<std::endl;
/* ^^^ checking these conditions ^^^ */
	
	num1=num2; 
	for (int k=1; k<4; k++)
	  { xa[k]=xb[k];
	  };	
     }; /* i */
     /* Last point of the new grid. Sigmagrid and xb equal those of the opposite string end.*/
   sigmagrid[Ngrid]=sigmaold[*N]; num2=*N-1;
   for (int k=1; k<4; k++)
     {xb[k]=xiold[k][*N];
     };   
   for (int j=1; j<4; j++) /*j - vector components*/
     {	
     dpab[j]=xdold[j][num1]*(sigmagrid[Ngrid]-sigmagrid[Ngrid-1]);
     dpab[j]+=xdold[j][num1]*(sigmaold[num1+1]-sigmagrid[Ngrid])-xdold[j][num2]*(sigmaold[num2+1]-sigmagrid[Ngrid]);
	for(int k=num1; k<num2; k++)
	  {
	     dpab[j]+=xdold[j][k]*(sigmaold[k+1]-sigmaold[k]);
	  }; /* k */
	dpab[j]=gamma*dpab[j];
	dpi[j][Ngrid-1]=dpab[j]; /* table of dpi for new sigma grid, last point*/
     }; /*j- vector components*/


        double dxab[4]={0,xb[1]-xa[1],xb[2]-xa[2],xb[3]-xa[3]};
   
   
        double dx=sqrt((xa[1]-xb[1])*(xa[1]-xb[1])+(xa[2]-xb[2])*(xa[2]-xb[2])+(xa[3]-xb[3])*(xa[3]-xb[3]));
	double dp=sqrt(dpab[1]*dpab[1]+dpab[2]*dpab[2]+dpab[3]*dpab[3]);
   
   /* Correction to make xdot \perp xprime */
           double dpdx=dpab[1]*(xb[1]-xa[1])+dpab[2]*(xb[2]-xa[2])+dpab[3]*(xb[3]-xa[3]);
  
   
  /*  NO CORRECTIONS: LOSS of ACCURACY
    
           if (gamma*dx>1.e-10){for (int j=1; j<4; j++) dpab[j]-=(xb[j]-xa[j])*dpdx/(dx*dx);} else {for (int j=1; j<4; j++) dxab[j]=0.0 ;};//-=dpab[j]*dpdx/(dp*dp);
            dp=sqrt((dpab[1])*(dpab[1])+(dpab[2])*(dpab[2])+(dpab[3])*(dpab[3]));
            dx=sqrt(pow(dxab[1],2)+pow(dxab[2],2)+pow(dxab[3],2));
   
   
	for (int j=1; j<4; j++)
	  {
           xdot[j][Ngrid-1]=dpab[j]/(sqrt(gamma*gamma*dx*dx+dp*dp));
	   xprime[j][Ngrid-1]=(xb[j]-xa[j])*gamma/(sqrt(gamma*gamma*dx*dx+dp*dp));  
	  };
*/

   
   
	for (int j=1; j<4; j++)
	  {
           xdot[j][Ngrid-1]=dpab[j]/(gamma*dsigma);
	   xprime[j][Ngrid-1]=(xb[j]-xa[j])/dsigma;  
	  };
   
/* Tables of sigma, xprime, xdot, dpi done for new grid*/
   
//   std::cout << Ngrid -1 <<" " << xdot[3][Ngrid-1] << std::endl; 
/************ MAKING DIRECTRIX *************/   

   for (int i=1; i<4; i++)
     {
   yd[i][Ngrid]=xiold[i][0];
    for (int j=1; j<=Ngrid;j++)
	  {
	     yd[i][Ngrid+j] = yd[i][Ngrid+j-1] + xprime[i][j-1]*dsigma + xdot[i][j-1]*dsigma;
	     yd[i][Ngrid-j] = yd[i][Ngrid-j+1] + xprime[i][j-1]*dsigma - xdot[i][j-1]*dsigma;
	  };	
     };
   yd[0][Ngrid]=xiold[0][0];
   for (int j=1;j<=Ngrid;j++)
     {
   yd[0][Ngrid+j]=yd[0][Ngrid+j-1]+dsigma;
   yd[0][Ngrid-j]=yd[0][Ngrid-j+1]-dsigma;	
     };
   
  for (int i=0; i<Ngrid; i++)
   {
      dpi[0][i]=gamma*(sigmagrid[i+1]-sigmagrid[i]);
   };
   
   
/* Reworking table of momenta*/
   
   
  for (int k=0; k<4; k++)
     {	
	        for (int i = 0; i<Ngrid; i++)
	  {	     
	                  dpi[k][i]=0.5*gamma* ( yd[k][Ngrid+i+1] + yd[k][Ngrid-i] - yd[k][Ngrid-i-1] - yd[k][Ngrid+i] );
	  };
     };
   
   
   
/*   for (int i=0; i<Ngrid; i++)
     {
	std::cout << xprime[1][i]<< " "<<xprime[2][i] <<" " <<xdot[1][i]<< " " << xdot[2][i]<<'\n';
     };
 */  
   
 /*  std::ofstream outfile;
   outfile.open("dirxy",std::ios::app);
   for (int i=0; i<Ngrid+1; i++)
     {
	outfile << xdot[1][i]<<" "<< xprime[1][i]<<" "<<yd[1][Ngrid+i] <<" " <<  xdot[2][i]<<" " << xprime[2][i]  <<" "<< yd [2][Ngrid+i] <<'\n';
     };   
   outfile.close();*/
};


void qgstring::WriteX (char ofile [20])
{
   std::ofstream outfile;
   outfile.open(ofile,std::ios::app);
   if (! outfile)
     { 
	outfile.open(ofile,std::ios::out);
     };
   
//   outfile << "# String coordinates" << '\n';
//   outfile << "# time " << '\n';
//   outfile << "# x1" << '\n';
//   outfile << "# x2" << '\n';
//   outfile << "# x3" << '\n';
   for (int i=0; i<2*Ngrid+1; i++)
     {
	/* i - 'time', k - coordinate, j - proportional to 'length' */
	/* time */
	//	outfile  << 0.5*(yd[0][(Ngrid+i)%(2*Ngrid)] + ((Ngrid+i)/(2*Ngrid))*(yd[0][2*Ngrid]-yd[0][0]) + yd[0][(Ngrid+i)%(2*Ngrid)] + ((Ngrid+i)/(2*Ngrid))*(yd[0][2*Ngrid]-yd[0][0])) << '\n';
//	for (int k = 1; k<4; k++)	  {
	     for (int j=0;j<=Ngrid; j++)
	       {
		  for (int k=1; k<4;k++)
		    {		       
   outfile << " " << 0.5*(yd[k][(Ngrid+i-j)%(2*Ngrid)] + ((Ngrid+i-j)/(2*Ngrid))*(yd[k][(2*Ngrid)]-yd[k][0]) +yd[k][(Ngrid+i+j)%(2*Ngrid)]+((Ngrid+i+j)/(2*Ngrid))*(yd[k][(2*Ngrid)]-yd[k][0]));
	            }; outfile <<'\n';
	       }; outfile << '\n';	     
	  };	
   //  };      
}


void qgstring::WriteX ()
{
            /* k - coordinate, j - proportional to 'length' */
                for (int j=0;j<=Ngrid; j++)
                 {
                      for (int k=1; k<3;k++)
                        {
		 std::cout << 0.5*(yd[k][Ngrid+j]+yd[k][Ngrid-j]) <<" " ;
                        }; std::cout << '\n';
                 }; std::cout << '\n';
}





void qgstring::WriteDir()
{
            /* k - coordinate, j - proportional to 'length' */
                for (int j=0;j<=2*Ngrid; j++)
                 {
                      for (int k=0; k<=3;k++)
                        {
			  std::cout.precision(10);
		 std::cout << yd[k][j] <<" " ;
                        }; std::cout << '\n';
                 }; std::cout << '\n';
}



double qgstring::GetMass2()
{
  double Mass2, pp[3]={0,0,0};
 
 for (int k =1; k<4; k++)   pp[k-1]=0.5*gamma*(yd[k][2*Ngrid]-yd[k][0]); 
  /* for(int i=0;i<Ngrid;i++)
     {	
	for (int j=1;j<4;j++)
	  {
	     pp[j-1]+=dpi[j][i];
	  };
     };*/
   
  Mass2=gamma*gamma*sigmagrid[Ngrid]*sigmagrid[Ngrid]-pp[0]*pp[0]-pp[1]*pp[1]-pp[2]*pp[2];
  if (fabs(Mass2)<1.e-11) Mass2=0; // computing M to 10^-5 GeV accuracy 
  return  Mass2;
}

double qgstring::GetMass2(int N1, int N2)
{
  double Mass2, pp[3]={0,0,0};
/* Check if 0 < N1 < N2 < Ngrid */   
  
 for (int k =1; k<4; k++)   pp[k-1]=0.5*gamma*(yd[k][Ngrid+N2]-yd[k][Ngrid-N2])-0.5*gamma*(yd[k][Ngrid+N1]-yd[k][Ngrid-N1]);   
   
 /*  
   for(int i=N1;i<N2;i++)
     {
     for (int j=1;j<4;j++)
         {
	    pp[j-1]+=dpi[j][i];
	 };
     };
 */  
    Mass2=gamma*gamma*(sigmagrid[N2]-sigmagrid[N1])*(sigmagrid[N2]-sigmagrid[N1])-pp[0]*pp[0]-pp[1]*pp[1]-pp[2]*pp[2];	
  if ((Mass2)<1.e-11) Mass2=0.0; 
   return Mass2;
}



double qgstring::GetMass2 (double s1, double s2)
{
double Mass2, pp[3]={0,0,0};
   int N1 = (int) floor(s1*Ngrid/sigmagrid[Ngrid]);
   int N2 = (int) ceil(s2*Ngrid/sigmagrid[Ngrid]);
   for (int k =1; k<4; k++) pp[k-1]=0.5*gamma*(yd[k][Ngrid+N2]-yd[k][Ngrid-N2])-0.5*gamma*(yd[k][Ngrid+N1]-yd[k][Ngrid-N1]);
   for (int k =1; k<4; k++) pp[k-1]-=GetMomentum(k,N1,N1+1)*(s1-sigmagrid[N1])/(sigmagrid[N1+1]-sigmagrid[N1]);
   for (int k =1; k<4; k++) pp[k-1]-=GetMomentum(k,N2-1,N2)*(sigmagrid[N2]-s2)/(sigmagrid[N2]-sigmagrid[N2-1]);   
     Mass2=gamma*gamma*(s2-s1)*(s2-s1)-pp[0]*pp[0]-pp[1]*pp[1]-pp[2]*pp[2];
     if ((Mass2)<1.e-11) Mass2=0.0;
      return Mass2;
   
};



void qgstring::TimeShift()
{
   double ydt[4];
   for (int k=0;k<4;k++)
     {
	ydt[k]=yd[k][2*Ngrid]-yd[k][0];
      for (int i=0; i<2*Ngrid; i++)
	  {
	     yd[k][i]=yd[k][i+1];
	  };	
	yd[k][2*Ngrid]=yd[k][0]+ydt[k];
     };
ltime+=sigmagrid[1]-sigmagrid[0];   
}


void qgstring::TimeShift(int Ntime)
{
/*directrix*/   
   double ydt[4][2*Ngrid+1];
   for (int k =0; k<4; k++)
     {
	for (int i = 0; i<=2*Ngrid; i++)
	  {
	     ydt[k][i]=yd[k][((i+Ntime)%(2*Ngrid)+2*Ngrid)%(2*Ngrid)]+floor((i+Ntime)*1.0/(2.0*Ngrid))*(yd[k][2*Ngrid]-yd[k][0]);
	  };
     };   
    for (int k =0; k<4; k++)
     {	
	 for (int i = 0; i<=2*Ngrid; i++)
	  {
	        yd[k][i]=ydt[k][i];
	  };
     };
/*momentum*/   
/*
   for (int k=0; k<4; k++)
     {
	for (int i = 0; i<Ngrid; i++)
	  {
	     dpi[k][i]=0.5*gamma* ( yd[k][Ngrid+i+1] + yd[k][Ngrid-i] - yd[k][Ngrid-i-1] - yd[k][Ngrid+i] );
	  };
     };  */
 ltime+=Ntime*(sigmagrid[1]-sigmagrid[0]);  
   
}

void qgstring::TimeShift(double tfm)
{ 
   double p[4]; // Shift after period of directrix
   for (int k=0;k<4;k++) p[k]= yd[k][2*Ngrid]-yd[k][0];
   double ds  = yd[0][1]-yd[0][0];  
   int tup = (int) (ceil(tfm/ds));
   double tover = tup*ds - tfm;
   //   std::cout << "tfm " << tfm << " tup*ds " << tup*ds <<" tover " << tover << " ds "<< ds <<std::endl;
  TimeShift(tup); //overdoing shift 
   double ydt[4][2*Ngrid+1];
   for (int k = 0; k<4; k++)
     { for (int j = 1; j<=2*Ngrid; j++) 
       {    ydt[k][j] = (tover*yd[k][j-1]+(ds-tover)*yd[k][j])/ds;};
	ydt[k][0]=ydt[k][2*Ngrid]-p[k];
     } 
   for (int k = 0; k<4; k++)
     { for (int j = 0; j<=2*Ngrid; j++)
       {
	      yd[k][j] = ydt[k][j];
       };	
     }
   
}



int qgstring::GetNgrid()
{
   return Ngrid;
}







double qgstring::GetXi(int k, int i)
{
   double xi=0.5*(yd[k][Ngrid+i]+yd[k][Ngrid-i]);
return xi;
}





void qgstring::LR()
{
  double ydt[4][2*Ngrid+1];
  for (int k =1; k<4; k++)
  {
    for (int i = 0; i<=Ngrid; i++)
      {
	ydt[k][i] = yd[k][Ngrid+i]+0.5*(yd[k][0]-yd[k][2*Ngrid]);
        ydt[k][2*Ngrid-i]=yd[k][Ngrid-i]+0.5*(yd[k][2*Ngrid]-yd[k][0]);
      }; 
  };
  for (int k =1; k<4; k++)
    {
      for (int i = 0; i<=2*Ngrid; i++)
	{
	  yd[k][i] = ydt[k][i];
	};
    };


};



void qgstring::Boost3(double bb)
{
   /*saving boost velocity */
   beta = (beta+bb)/(1+beta*bb);
   /* obtaining  max x3, предполагая упорядочение по x3 - assuming ordering in x3 */
   double xmax;
     if (bb > 0) 
       {  if (GetXi(3, 0) > GetXi(3,Ngrid)) {xmax = GetXi(3,0);} else xmax = GetXi(3, Ngrid);
       } else
     { if (GetXi(3, 0) < GetXi(3,Ngrid)) {xmax = GetXi(3,0);} else xmax = GetXi(3, Ngrid);
     };
   
      
   /* obtaining parameter tmin*/
     tmin = (yd[0][Ngrid]- bb * xmax)/sqrt(1-bb*bb);
   
     //  std::cout <<  "x1 "<< GetXi(3,0) <<" x2 " << GetXi(3,Ngrid)<< " xmax " << xmax <<" tmin "<< tmin <<std::endl;
   
//   std::cout << "time "<< yd[0][0] << "  " << yd[0][2*Ngrid] << std::endl;
//   std::cout << "coord " << yd[3][0] << " " << yd[3][2*Ngrid] <<std::endl;

   

            /* k - coordinate, j - proportional to 'length' */
   /*
   std::cout << "dir before boost" << std::endl;
                for (int j=0;j<=2*Ngrid; j++)
                 {
                      for (int k=0; k<=3;k++)
                        {
		 std::cout << yd[k][j] <<" " ;
                        }; std::cout << '\n';
                 };
   std::cout << "********" <<std::endl;

   */
   
  /*Boosting direcrix*/ 
   double ydt[4][2*Ngrid+1];
   for (int i = 0; i<=2*Ngrid; i++)
     { ydt[0][i] = (yd[0][i] - bb*yd[3][i])/sqrt(1-bb*bb);
       ydt[3][i] = (yd[3][i] - bb*yd[0][i])/sqrt(1-bb*bb);
	ydt[1][i] = yd[1][i];
	ydt[2][i] = yd[2][i];
     };
//   std::cout << "time " << ydt[0][0] << " " << ydt[0][2*Ngrid] << std::endl;
//   std::cout << "coord " << ydt[3][0] << " " << ydt[3][2*Ngrid] <<std::endl;

   
   


            /* k - coordinate, j - proportional to 'length' */
   /*
   std::cout << "dir after boost" << std::endl;
                for (int j=0;j<=2*Ngrid; j++)
                 {
                      for (int k=0; k<=3;k++)
                        {
		 std::cout << ydt[k][j] <<" " ;
                        }; std::cout << '\n';
                 };
   std::cout << "********" <<std::endl;

   */
   
  sigmat = 0.5*(ydt[0][2*Ngrid]-ydt[0][0]);
//   std::cout << "total sigma " << sigmat << std::endl;
   double dsigma  = sigmat / (1.*Ngrid);
//   std::cout << "dsigma " << dsigma << std::endl;
/*reinterpolation*/
yd[0][0] = ydt[0][0];
yd[1][0] = ydt[1][0];
yd[2][0] = ydt[2][0];   
yd[3][0] = ydt[3][0];   
   for (int i=1; i< 2*Ngrid; i++)
		{
		  yd[0][i]=ydt[0][0]+i*dsigma;
                  int k=0;
		  while (yd[0][i]>ydt[0][k]) k++;
		   
		   k--;
              /*      k       k+1            */
	      /*          i                  */  
		   yd[3][i] = (ydt[3][k]*(ydt[0][k+1]-yd[0][i]) + ydt[3][k+1]*(yd[0][i]-ydt[0][k]))/(ydt[0][k+1]-ydt[0][k]);
   		   yd[2][i] = (ydt[2][k]*(ydt[0][k+1]-yd[0][i]) + ydt[2][k+1]*(yd[0][i]-ydt[0][k]))/(ydt[0][k+1]-ydt[0][k]);
		   yd[1][i] = (ydt[1][k]*(ydt[0][k+1]-yd[0][i]) + ydt[1][k+1]*(yd[0][i]-ydt[0][k]))/(ydt[0][k+1]-ydt[0][k]);
		}; // reinterpolation done
yd[0][2*Ngrid]= ydt[0][2*Ngrid];
yd[1][2*Ngrid]= ydt[1][2*Ngrid];
yd[2][2*Ngrid]= ydt[2][2*Ngrid];   
yd[3][2*Ngrid]= ydt[3][2*Ngrid];   
   sigmagrid[0]=0.;
   for (int i=0; i< Ngrid+1 ; i++) sigmagrid[i] = i*dsigma; 
//   std::cout << "time "<< yd[0][0] << "  " << yd[0][2*Ngrid] << std::endl;
//   std::cout << "coord " << yd[3][0] << " " << yd[3][2*Ngrid] <<std::endl;
   
   
/* now moving string to the past */		
   
//   std::cout << " yd[0][Ngrid] = " << yd[0][Ngrid] << " ydt[0][Ngrid] = "<< ydt[0][Ngrid] <<std::endl;
//   std::cout << "shift by " << (int) round((tmin - yd[0][Ngrid])/dsigma) <<  " " << round((tmin - yd[0][Ngrid])/dsigma) <<std::endl;
//	TimeShift((int) round((tmin - yd[0][Ngrid])/dsigma) );	
        TimeShift(tmin - yd[0][Ngrid]);
	//   std::cout << "after shift yd[0][Ngrid] = " << yd[0][Ngrid] << std::endl;
   /* setting local time to zero and sigma grid with y[0][...]  to 0...\sigma_\max */   
	ltime =0.;
	//	sigmagrid[0] = 0.;
/*	yd[0][Ngrid]=0.;	
	for (int i=1; i<=Ngrid; i++)
		{
		   yd[0][Ngrid-i]=-i*dsigma; yd[0][Ngrid+i]=i*dsigma;
		   sigmagrid[i] = i*dsigma;
		};*/
	TimeShift(0.5*dsigma);
	//   std::cout << "boost done with beta = " << bb << std::endl;
   //   std::cout << "fragment velocity in the lab frame betalab = " << beta << std::endl;
//WriteDir();
}


/*
double qgstring::GetdPi(int k, int i)
{
   return dpi[k][i];
}
*/

double qgstring::GetMomentum(int k, int N1, int N2)
{
   double pp=0.;
   pp=0.5*gamma*(yd[k][Ngrid+N2]-yd[k][Ngrid-N2])-0.5*gamma*(yd[k][Ngrid+N1]-yd[k][Ngrid-N1]);
  /* for (int i=N1; i<N2;i++)
     { 
	pp+=dpi[k][i];
     };*/
   return pp;
}


double qgstring::GetMomentumLab(int k, int N1, int N2)
{
   double pp[4],ppp[4];
   for (int j=0; j<4; j++) pp[j]=0.5*gamma*(yd[j][Ngrid+N2]-yd[j][Ngrid-N2])-0.5*gamma*(yd[j][Ngrid+N1]-yd[j][Ngrid-N1]);
  
   /* inverse boost along 3 axis */ 
   
   ppp[0] = (pp[0]+beta*pp[3])/sqrt(1-beta*beta);
   ppp[3] = (pp[3]+beta*pp[0])/sqrt(1-beta*beta);
   ppp[2] = pp[2];
   ppp[1] = pp[1];
   
   return ppp[k];
}


double qgstring::GetBeta3()
{
return beta;
}




double qgstring::GetSigmagrid(int i)
{
return sigmagrid[i];
}






double qgstring::GetSigmaDecay(double lowlim)
{   
       double sigmabreak = -1.;
       double dP[2*Ngrid+1], pno=1.,sigmaleft[2*Ngrid],sigmaright[2*Ngrid];
       int num=-1;
//   std::cout << yd[3][0] << " "<< yd[3][Ngrid]<<" "<< yd[3][2*Ngrid]<<std::endl;
 //  for (int it = 0; it<2*Ngrid; it++)
   int it = 0;
   int decayflag = 0;
   int sleft=0,sright=Ngrid;
   while((it < Ngrid)and(decayflag==0))  
 //  for (int it = 0; it<Ngrid; it++)
     { 
	dP[it]=0;

	if (GetMass2(0,sleft)<lowlim){do {sleft++;} while (GetMass2(0,sleft)<lowlim);} else {do {sleft--;} while (GetMass2(0,sleft)>lowlim); sleft++;};
        if (GetMass2(sright,Ngrid)<lowlim) {do {sright--;} while(GetMass2(sright,Ngrid)<lowlim);} else {do {sright++;} while (GetMass2(sright,Ngrid)>lowlim);sright--;}
	
	
	
	      /** point for 'right' is between sright and sright+1 **/
              double M2= sqrt(GetMass2(sright+1,Ngrid)), M1= sqrt(GetMass2(sright,Ngrid));
             // std::cout << t << "  sleft="<<sleft<<" Mass="<< M2  <<" Mass-1=" <<  M1  <<'\n';
             // std::cout << qgstring1->GetMass2(0,sright+1) << " " << M2break<<" "<<qgstring1->GetMass2(0,sright) << std::endl;
              double sigma1=GetSigmagrid(sright),sigma2=GetSigmagrid(sright+1);
            //  double sigmabreak = sigma1+(sqrt(mp2)-M1)*(sigma2-sigma1)/(M2-M1);
             
              double dels=sigma2-sigma1;
              double dM2=GetMass2(sright,sright+1);
              double pdp=0.;
              if (dM2>0){
	      double EdelE=GetMomentum(0,sright,sright+1)*GetMomentum(0,sright+1,Ngrid);
              for (int i=1;i<4;i++) pdp+=GetMomentum(i,sright,sright+1)*GetMomentum(i,sright+1,Ngrid);
              double dsigma=dels*(-(EdelE-pdp)+sqrt(pow(EdelE-pdp,2)+(lowlim-M2*M2)*dM2))/dM2;
           //   std::cout << "Sbold = "<< sigmabreak << " Sbnew = " <<sigma2-dsigma<<std::endl;
	      sigmaright[it]=sigma2-dsigma;
	      }else{
		 double EdelE=GetMomentum(0,sright,sright+1)*GetMomentum(0,sright+1,Ngrid);
		  for (int i=1;i<4;i++) pdp+=GetMomentum(i,sright,sright+1)*GetMomentum(i,sright+1,Ngrid);
		  double dsigma = (lowlim- M2*M2)*dels/(2*EdelE-2*pdp); //sigmaright[it]=sigma1+(sqrt(lowlim)-M1)*(sigma2-sigma1)/(M2-M1);
	      sigmaright[it]= sigma2-dsigma; };
	
	
		 
	   /* the same for sleft */ 
		 
		 
             M2= sqrt(GetMass2(0,sleft));
	     M1= sqrt(GetMass2(0,sleft-1));
		             /* std::cout << t << "  sleft="<<sleft<<" Mass="<< M2  <<" Mass-1=" <<  M1  <<'\n'; */
//		 std::cout << qgstring1->GetMass2(sleft,Ng) << " " << M2break<<" "<<qgstring1->GetMass2(sleft-1,Ng) << std::endl;
		              sigma1=GetSigmagrid(sleft-1); sigma2=GetSigmagrid(sleft);
		           //   double sigmabreak = sigma1+ (sqrt(mp2)-M1)*(sigma2-sigma1)/(M2-M1);
		              //double xbreak[4],pbreak[4];
		  
	                         dels=sigma2-sigma1;
		              dM2=GetMass2(sleft-1,sleft);
		              if(dM2>0)
		   {               double pdp=0.;
		                   double EdelE=GetMomentum(0,sleft-1,sleft)*GetMomentum(0,0,sleft-1);
		                   for (int i=1;i<4;i++) pdp+=GetMomentum(i,sleft-1,sleft)*GetMomentum(i,0,sleft-1);
		                   double dsigma=dels*(-(EdelE-pdp)+sqrt(pow(EdelE-pdp,2)+(lowlim-M1*M1)*dM2))/dM2;
		      // std::cout << "Sbold = "<< sigmabreak << " Sbnew = " <<sigma1+dsigma<<std::endl;
		                   sigmaleft[it]=sigma1+dsigma;} else { 
				       double pdp=0.;
				       double EdelE=GetMomentum(0,sleft-1,sleft)*GetMomentum(0,0,sleft-1);
				       for (int i=1;i<4;i++) pdp+=GetMomentum(i,sleft-1,sleft)*GetMomentum(i,0,sleft-1);
				       double dsigma = (lowlim-M1*M1)*dels/(2*EdelE-2*pdp); //sigmaleft[it]= sigma1+ (sqrt(lowlim)-M1)*(sigma2-sigma1)/(M2-M1);
				   sigmaleft[it]=sigma1+dsigma;};
		 
		 


	 /* Now corrections for local time after transformation, local time = yd[0][Ngrid] */
	 double xc, xl = GetXi(3,0), xr = GetXi(3,Ngrid);
	 if (fabs(beta)>0) xc = - yd[0][Ngrid]/beta; 

	 if(((beta < 0)and(xl>=xc)and(xr>=xc))or((beta>0)and(xl<=xc)and(xr<=xc))){sigmaleft[it]=sigmaright[it];};//no breaking	 
	 if(((beta < 0)and(xl>xc)and(xr<xc))or((beta>0)and(xl<xc)and(xr>xc))){
	   sleft = 0;
         do {sleft++;} while (beta*GetXi(3,sleft)<beta*xc);
	   double st=(sigmagrid[sleft-1]*(GetXi(3,sleft)-xc)+sigmagrid[sleft]*(xc-GetXi(3,sleft-1)))/(GetXi(3,sleft)-GetXi(3,sleft-1));
           if (st>sigmaleft[it]) sigmaleft[it]=st;
         };//sigmaleft modification if needed	 
	 if(((beta < 0)and(xl<xc)and(xr>xc))or((beta>0)and(xl>xc)and(xr<xc))){
           sright = Ngrid;
         do {sright--;} while (beta*GetXi(3,sright)<beta*xc);	
	   double st =(sigmagrid[sright]*(GetXi(3,sright+1)-xc)+sigmagrid[sright+1]*(xc-GetXi(3,sright)))/(GetXi(3,sright+1)-GetXi(3,sright)); 
	    if (st < sigmaright[it]) sigmaright[it] = st;
         };//sigmaright modification if needed 

	
         if (sigmaleft[it] >=sigmaright[it])
           { num=-1; dP[it] = 0;// splitting impossible
           }else
           { 
              for (int N=sleft; N<=sright+1; N++)
              {
                double xxl[4], xxr[4];
                   for (int i=1; i<4; i++){xxl[i] = (yd[i][Ngrid+N-1]+yd[i][Ngrid-N+1]-yd[i][Ngrid+N]-yd[i][Ngrid-N]);};
		 
                   dP[it] += 0.25*alp*(xxl[1]*xxl[1] + xxl[2]*xxl[2] + xxl[3]*xxl[3]);
	      };
	   double xxl[4];
	      int N=sleft;
	      for (int i=1; i<4; i++){  xxl[i] = (yd[i][Ngrid+N-1]+yd[i][Ngrid-N+1]-yd[i][Ngrid+N]-yd[i][Ngrid-N]); };
	      
      dP[it]-=(sigmaleft[it]-sigmagrid[sleft-1])*0.25*alp*(xxl[1]*xxl[1] + xxl[2]*xxl[2] + xxl[3]*xxl[3])/(sigmagrid[sleft]-sigmagrid[sleft-1]);

              N=sright+1;
	                    for (int i=1; i<4; i++)
		{xxl[i] = (yd[i][Ngrid+N-1]+yd[i][Ngrid-N+1]-yd[i][Ngrid+N]-yd[i][Ngrid-N]); };
	      
          dP[it]-=(-sigmaright[it]+sigmagrid[sright+1])*0.25*alp*(xxl[1]*xxl[1] + xxl[2]*xxl[2] + xxl[3]*xxl[3])/(sigmagrid[sright+1]-sigmagrid[sright]);
	      
	      
	   };
	TimeShift(); //std::cout << dP[it] << " ";
//	std::cout << 0.5*(yd[3][2*Ngrid]+yd[3][0])-yd[3][Ngrid] << std::endl;
//	std::cout << sigmaleft[it] << " " << sigmaright[it] << std::endl;	
       /* Checking whether break occurs now, simulating random number */
	double dpno = exp(-dP[it]);
	double z = (rand()*1.0)/(RAND_MAX*1.0);
	if (z > dpno) decayflag = it; // decay at this time
	it++;
     };    
      //  TimeShift(-2*Ngrid);
      //  TimeShift(-Ngrid);
     TimeShift(-it);
   
  //   double P[2*Ngrid]; //2*Ngrid steps is the full cycle
      double P[Ngrid]; //Ngrid steps is the full cycle
   P[0]=1.-exp(-dP[0]); pno = exp(-dP[0]);
 //    for ( int it = 1; it < 2*Ngrid; it++)
     for ( int itt = 1; itt < Ngrid; itt++)     
     { 
	P[itt] = P[itt-1]+pno*(1-exp(-dP[itt]));
	pno=pno*exp(-dP[itt]);
     };
   	 
     //   std::cout << " Probabilities of decay after step 0, 1, 2, 3: "<< P[0] << " "<< P[1] << " "<< P[2] << " " << P[3] << std::endl;
 //  std::cout << P[2*Ngrid-1] << std::endl;
  
 //  std::cout << yd[3][0] <<" "<<yd[3][Ngrid]<<" "<< yd[3][2*Ngrid]<<std::endl;
   
  //     double z=(rand()*1.0*P[2*Ngrid-1])/RAND_MAX;
     double z=(rand()*1.0*P[Ngrid-1])/RAND_MAX;   
   
 //  if (P[2*Ngrid-1]<1.e-10) 
       if ((P[Ngrid-1]<1.e-10)and(decayflag==0)) 
     {std::cout << "break impossible" <<std::endl;
	std::cout << "momentum " <<  sqrt(gamma*gamma*GetSigmagrid(Ngrid)*GetSigmagrid(Ngrid)-GetMass2()) << std::endl;
	std::cout << "sigma " << GetSigmagrid(Ngrid) << std::endl;
        sigmabreak = -1.;
     }else
     {
if (decayflag==0)
	  {	     
     it =0;
     while (z>P[it]) it++; 
 //  now 'it' gives the time of decay
	std::cout << "decay time it = "<< it << std::endl; 
	
   TimeShift(it);} else {TimeShift(decayflag); it = decayflag;}; // moving string to the right time
 
//   std::cout << it << " " << lowlim<<" " << sigmaleft[it] <<" "<< GetMass2(0.,sigmaleft[it])<<" "<<GetMass2(sigmaright[it],sigmagrid[Ngrid]) <<" "<< sigmaright[it]<<std::endl;
    

   
          int sleft=0,sright=Ngrid;
           do 
     {
	sleft++;
     }
    while (sigmagrid[sleft]<sigmaleft[it]);
            do 
     {
	sright--;
     }
    while (sigmagrid[sright]>sigmaright[it]);
   
   
    
//   std::cout << sigmaleft[it]<<" "<< sleft << " " <<sright <<" " << sigmaright[it] << std::endl;
    
   if (sleft>sright)
     {
     sigmabreak = sigmaleft[it] + ((sigmaright[it] - sigmaleft[it])*rand())/RAND_MAX;
     }
   else
     {	
   double dp[Ngrid+1]; double iP[Ngrid+1], xr[4];
   for (int ii=0; ii<Ngrid+1; ii++) {dp[ii]=0.; iP[ii]=0;};
     int iN=sleft-1;
   for (int ii=1; ii<4; ii++) xr[ii]=  (yd[ii][Ngrid+iN]+yd[ii][Ngrid-iN]-yd[ii][Ngrid+iN+1]-yd[ii][Ngrid-iN-1]);
      dp[sleft-1]=-(sigmaleft[it]-sigmagrid[sleft-1])*0.25*alp*(xr[1]*xr[1] + xr[2]*xr[2] + xr[3]*xr[3])/(sigmagrid[sleft]-sigmagrid[sleft-1]);
   
   
       iN=sright;
       for (int ii=1; ii<4; ii++)  xr[ii]=  (yd[ii][Ngrid+iN]+yd[ii][Ngrid-iN]-yd[ii][Ngrid+iN+1]-yd[ii][Ngrid-iN-1]);
       dp[sright]=-(-sigmaright[it]+sigmagrid[sright+1])*0.25*alp*(xr[1]*xr[1] + xr[2]*xr[2] + xr[3]*xr[3])/(sigmagrid[sright+1]-sigmagrid[sright]);
   
           for (int N=sleft-1; N<=sright; N++)
               {	          // double xr[4];
	              for (int i=1; i<4; i++)
	         {   // xl[i] = (yd[i][Ngrid+N-1]+yd[i][Ngrid-N+1]-yd[i][Ngrid+N]-yd[i][Ngrid-N]);
	            xr[i] =  (yd[i][Ngrid+N]+yd[i][Ngrid-N]-yd[i][Ngrid+N+1]-yd[i][Ngrid-N-1]);};
	            dp[N]+= 0.25*alp*(xr[1]*xr[1] + xr[2]*xr[2] + xr[3]*xr[3] );
	         iP[N]=iP[N-1]+dp[N];
//		  std::cout << N << " " << dp[N]<<" " << sigmagrid[N] << std::endl;
	       };
	
               z=(rand()*1.0)/RAND_MAX;
//	std::cout << z << std::endl;
               z=z*iP[sright];
	       num=sleft-1;
               while(iP[num]<z) num++;
              if (num == sleft-1) 
	  { 
	     sigmabreak = sigmaleft[it] + ((sigmagrid[sleft]-sigmaleft[it])*rand())/RAND_MAX;
	  }else
	     {
		if(num == sright)
		  {
		     sigmabreak = sigmagrid[sright] +  ((sigmaright[it]-sigmagrid[sright])*rand())/RAND_MAX;
		  }else
		     { //break is between num and num+1
		     sigmabreak = sigmagrid[num]+((sigmagrid[num+1]-sigmagrid[num])*rand())/RAND_MAX;
		     };				
	      };
	
	
              //std::cout << "break at " << num << "; probability was "<< iP[sright]<< std::endl;
   
     };
//      sigmabreak = 0.5*(sigmaleft[it]+sigmaright[it]); 
//   std::cout << "it="<<it << " sigmaleft "<< sigmaleft[it] << " decay " << sigmabreak << " right " <<sigmaright[it] << std::endl;
//   std::cout << " m2left "<< GetMass2(0.,sigmaleft[it]) << " decay " << sigmabreak << " right " << GetMass2(sigmaright[it],sigmagrid[Ngrid]) << std::endl;
//	for (int i=sleft-1; i<=sright; i++) std::cout << dP[i] << std::endl;
    
//          sigmabreak = sigmaleft[it] + ((sigmaright[it] - sigmaleft[it])*rand())/RAND_MAX;


     };
   

   
          return sigmabreak;
     }
      
