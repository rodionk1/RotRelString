#include <fstream>
#include <iostream>
//#include <string>
#include <iomanip>
#include "stdlib.h"
#include "time.h"
#include "math.h"
#include "qgstring.hh"
#include "qgstring.cc"
#include "qgstack.hh"

//#include "RandGauss.h"

//using namespace std;

double ff(double z);
double fdistr(double z);
double Fdistr(double z);
double RandGauss();

int main()
{   
   srand(time(NULL));
   
   int Nsim=5000;
    char of[8] = "frag100";
    char off[17] = "fragfile_100_b=0";
     double b=0.5;
     double ppp=100,mpp=-100.;
     double R0 = 0.5, alprime = 0.01; // Regge-inspired parameters, vertex radius and pomeron slope
   

     double xin[4][3] = {{0,0,0},{0,0,0},{0,0,0},{0,0,0}};
     //   double pin[4][3] = {0,0,0,0,0,0,1,0,1,0,0,0};                                                                                                        

     double pin[4][2]={{0,0},{0,0},{0,0},{0,0}};

     pin[3][0] = mpp; pin[3][1]=ppp;

     double * xxi = &xin[0][0]; double *  ppi = &pin[0][0];
     int Ni=2;

   
 for (int nsim=1; nsim<=Nsim; nsim++) //simulations cycle
     {  	std::cout << "******** Sim " << nsim << "*******" <<std::endl;
  qgstring * qgstring1 = new qgstring (xxi, ppi, &Ni);
	 qgstack *stack1 = new qgstack();
	 qgstack *stack2 = new qgstack(); // storage of primary fragments
	stack1->put(*qgstring1);
	delete qgstring1;
//	std::cout << "put done " << std::endl; 	
  do{	
//     std::cout << "taking from stack "<<std::endl;
      qgstring1 = stack1->get();     
   std::cout << "get done, mass2 =  " << qgstring1->GetMass2() <<std::endl;
     
      int Ng=qgstring1->GetNgrid();
      int N=Ng;
   
      double mp2=pow(0.144,2); 
if(qgstring1->GetMass2()<2.*mp2) 
       { //WRITING TO FILE     
//	  std::cout << "writing low mass string m = " << qgstring1->GetMass2() <<std::endl;	  
   	  std::ofstream outfile;
		  
	  
	  outfile.open(of,std::ios::app);
	                  if (! outfile)
			    {outfile.open(of,std::ios::out);};
	  outfile.precision(19);
	  for (int k=0; k<4; k++) outfile << qgstring1->GetMomentum(k,0,N) << " ";
	               outfile <<'\n';
	  
	               double pout[4];
	  
	  outfile.close();
	  delete qgstring1; // Freeing memory
	  
       }else
       {// making a breakpoint   
     
   if (qgstring1->GetMass2()<100.25*mp2) //low mass string
       {  
	  std::cout << "doing low mass string m = " << sqrt(qgstring1->GetMass2())/(0.144) <<" mpi" <<std::endl;
	  
	  // Breaking in exactly 2 pieces with pi-meson masses, writing them both into the file and freeing memory
	  
	  
	 /* Boost matrix */
	  double lam[4][4];
	  
	  
	  
	  double pstring[4],ptot=0.,p1,p2, M2=qgstring1->GetMass2(), E=qgstring1->GetMomentumLab(0,0,N);
	  for (int i=1; i<4; i++){ pstring[i]=qgstring1->GetMomentumLab(i,0,N); ptot+=pstring[i]*pstring[i];};
      	  ptot = sqrt(ptot);

	  
	  double beta[4];
	  beta[0]=1.;
	  beta[1]=pstring[1]/E;
	  beta[2]=pstring[2]/E;
	  beta[3]=pstring[3]/E;
	  double b2=beta[1]*beta[1]+beta[2]*beta[2]+beta[3]*beta[3];
	  double gam = 1./sqrt(1-b2);
	  
	  lam[0][0]=gam;
	  for (int i=1;i<4; i++) { 
	     lam[0][i] = beta[i]*gam; lam[i][0]=lam[0][i]; lam[i][i] = 1 + beta[i]*beta[i]*(gam-1)/b2;
               };
	  
	  lam[1][2]=beta[1]*beta[2]*(gam-1)/b2;
	  lam[2][1]=lam[1][2];
	  lam[2][3]=beta[2]*beta[3]*(gam-1)/b2;
	  lam[3][2]=lam[2][3];
	  lam[1][3]=beta[1]*beta[3]*(gam-1)/b2;
	  lam[3][1]=lam[1][3];
	  
	  double pdecay = sqrt(M2/4-mp2);
	  double p1c[4], p2c[4];
	  double z1 = (2.*rand())/RAND_MAX-1.;
	  double z2 = (8.*atan(1.)*rand())/RAND_MAX;
	  p1c[3] = pdecay*z1;
	  p2c[3] = -p1c[3];
	  p1c[2] = pdecay*cos(z2)*sqrt(1-z1*z1);
	  p2c[2] = -p1c[2];
	  p1c[1] = pdecay*sin(z2)*sqrt(1-z1*z1);
	  p2c[1] = -p1c[1];
	  
	  p1c[0] = sqrt (mp2+p1c[1]*p1c[1]+p1c[2]*p1c[2]+p1c[3]*p1c[3]);
	  p2c[0] = sqrt (mp2+p2c[1]*p2c[1]+p2c[2]*p2c[2]+p2c[3]*p2c[3]);
	  
//	  p1 = (ptot*M2+sqrt(E*E*M2*(M2-4.*mp2)))/(2.*M2);
//	  p2 = (ptot*M2-sqrt(E*E*M2*(M2-4.*mp2)))/(2.*M2);
	  
//	  std::cout << p1 << " " << p2 << " " << ptot<<std::endl;
	  
	  
/*  boost */
	  
	  double p1o[4], p2o[4];
	  
	  for (int k=0;k<4;k++)
	    { p1o[k]=0.; p2o[k]=0.;
	    for (int j=0;j<4; j++)
	      {
		 p1o[k]+=lam[k][j]*p1c[j];
		 p2o[k]+=lam[k][j]*p2c[j];
	      };	       	       
	    };
	  
	  std::ofstream outfile;
	  
	  outfile.open(of,std::ios::app);
	                  if (! outfile)
	    {outfile.open(of,std::ios::out);};
	              outfile.precision(19);
	  outfile << p1o[0] <<" ";
	  for (int k=1; k<4; k++) outfile << p1o[k] << " ";
	               outfile <<'\n';
	  
	 //            double pout[4];
	  
                 	  
	              outfile << p2o[0] << " " << p2o[1] << " " << p2o[2] << " " << p2o[3] << std::endl;
	        
                      outfile.close();	  
                delete qgstring1;
       } //ALL FOR LOW MASS
	  else
       { //INTERMEDIATE AND HIGH
			 
//		 std::cout << "doing string of mass m = " << sqrt(qgstring1->GetMass2()/mp2) <<" mpi" << std::endl;
		           int tt=0; double sDecay = -1.; int iDecay;
		        //  do{
		      
		      
		             sDecay = qgstring1->GetSigmaDecay(4.01*mp2);
		
			     			  
			      if ( sDecay < 0)
			        {  
				// we break exactly into 2 pi-mesons
				// 
	  
				   
				   
				   
	

				   
	 /* Boost matrix */
	  double lam[4][4];
	  
	  
	  
	  double pstring[4],ptot=0.,p1,p2, M2=qgstring1->GetMass2(), E=qgstring1->GetMomentumLab(0,0,N);
	  for (int i=1; i<4; i++){ pstring[i]=qgstring1->GetMomentumLab(i,0,N); ptot+=pstring[i]*pstring[i];};
      	  ptot = sqrt(ptot);

	  
	  double beta[4];
	  beta[0]=1.;
	  beta[1]=pstring[1]/E;
	  beta[2]=pstring[2]/E;
	  beta[3]=pstring[3]/E;
	  double b2=beta[1]*beta[1]+beta[2]*beta[2]+beta[3]*beta[3];
	  double gam = 1./sqrt(1-b2);
	  
				   
	if (b2>0){			   
	  lam[0][0]=gam;
	  for (int i=1;i<4; i++) { 
	     lam[0][i] = beta[i]*gam; lam[i][0]=lam[0][i]; lam[i][i] = 1 + beta[i]*beta[i]*(gam-1)/b2;
               };
	  
	  lam[1][2]=beta[1]*beta[2]*(gam-1)/b2;
	  lam[2][1]=lam[1][2];
	  lam[2][3]=beta[2]*beta[3]*(gam-1)/b2;
	  lam[3][2]=lam[2][3];
	  lam[1][3]=beta[1]*beta[3]*(gam-1)/b2;
	  lam[3][1]=lam[1][3];
	}else { 
					
	   for (int i = 0; i<4; i++) lam[i][i]=1.;
	   for (int i = 1; i<4; i++) {lam[0][i]=0.; lam[i][0]=0.;};
	             lam[1][2]=0;
	             lam[2][1]=0;
	             lam[2][3]=0;
	             lam[3][2]=0;
	             lam[1][3]=0;
	             lam[3][1]=0;	   
	};
				   
				   
				   
				   
	  double pdecay = sqrt(M2/4-mp2);
	  double p1c[4], p2c[4];
	  double z1 = (2.*rand())/RAND_MAX-1.;
	  double z2 = (8.*atan(1.)*rand())/RAND_MAX;
	  p1c[3] = pdecay*z1;
	  p2c[3] = -p1c[3];
	  p1c[2] = pdecay*cos(z2)*sqrt(1-z1*z1);
	  p2c[2] = -p1c[2];
	  p1c[1] = pdecay*sin(z2)*sqrt(1-z1*z1);
	  p2c[1] = -p1c[1];
	  
	  p1c[0] = sqrt (mp2+p1c[1]*p1c[1]+p1c[2]*p1c[2]+p1c[3]*p1c[3]);
	  p2c[0] = sqrt (mp2+p2c[1]*p2c[1]+p2c[2]*p2c[2]+p2c[3]*p2c[3]);
	  
//	  p1 = (ptot*M2+sqrt(E*E*M2*(M2-4.*mp2)))/(2.*M2);
//	  p2 = (ptot*M2-sqrt(E*E*M2*(M2-4.*mp2)))/(2.*M2);
	  
//	  std::cout << p1 << " " << p2 << " " << ptot<<std::endl;
	  
	  
/*  boost */
	  
	  double p1o[4], p2o[4];
	  
	  for (int k=0;k<4;k++)
	    { p1o[k]=0.; p2o[k]=0.;
	    for (int j=0;j<4; j++)
	      {
		 p1o[k]+=lam[k][j]*p1c[j];
		 p2o[k]+=lam[k][j]*p2c[j];
	      };	       	       
	    };
	  
	  std::ofstream outfile;
	  
	  outfile.open(of,std::ios::app);
	                  if (! outfile)
	    {outfile.open(of,std::ios::out);};
         outfile.precision(19);	              
	  outfile << p1o[0] <<" ";
	  for (int k=1; k<4; k++) outfile << p1o[k] << " ";
	               outfile <<'\n';
	  
                 	  
	              outfile << p2o[0] << " " << p2o[1] << " " << p2o[2] << " " << p2o[3] << std::endl;
	  
                   outfile.close();	  
				 				   
	     delete qgstring1;
								
				}
			      else
				{
								   
		               /** now sDecay lies between iDecay and iDecay+1 **/
			      iDecay=0;
		              while (qgstring1->GetSigmagrid(iDecay)<sDecay) iDecay++;
		              iDecay--;
		//	      std::cout << "decay between  " << iDecay << " and "<< iDecay+1 << " of "<< Ng << std::endl;

		  //         std::cout << "left part mass2 "<<qgstring1->GetMass2(0.,sDecay)<< " right part mass2 " << qgstring1->GetMass2(sDecay,qgstring1->GetSigmagrid(N)) << std::endl;
			   
	      
	         qgstring* qgsr = new qgstring(sDecay, qgstring1->GetSigmagrid(N), *qgstring1);
		 //                 qgstring1->LR();          
		 qgstring* qgsl = new qgstring(0.,sDecay,*qgstring1); //qgsl->LR();
		 //   std::cout << "after break ml = " << qgsl->GetMass2() <<" mr = " << qgsr->GetMass2() <<std::endl; 
				   
			double betal = qgsl->GetMomentum(3,0,Ng)/qgsl->GetMomentum(0,0,Ng);
			double betar = qgsr->GetMomentum(3,0,Ng)/qgsr->GetMomentum(0,0,Ng);	   
//	qgsl->WriteDir();
//        qgsr->WriteDir();				   
            if (qgsl->GetMass2()<100.26*mp2) {stack2->put(*qgsl); stack1->put(*qgsl);} else 
	      { if((qgsl->GetBeta3()==0)and(fabs(betal)>0.985)){qgsl->Boost3(betal); std::cout << "betal = " << betal << std::endl;} stack1->put(*qgsl); };
				                                         //        std::cout << "   put left" <<std::endl;
	    if (qgsr->GetMass2()<100.26*mp2) {stack2->put(*qgsr); stack1->put(*qgsr);} else
	      { if((qgsr->GetBeta3()==0)and(fabs(betar)>0.985)){qgsr->Boost3(betar); std::cout << "betar = " << betar << std::endl;} stack1->put(*qgsr);};
//qgsl->WriteDir(); qgsr->WriteDir();							   
                 delete qgsl; delete qgsr; delete qgstring1;
                                                    			      
	      };

       }; //ALL FOR INTERMEDIATE AND HIGH
	  
       }; // ALL FOR BREAKPOINT
     
       }while(not(stack1->isEmpty())); // ALL FOR EMPTY STACK

	

	std::ofstream ofile;
	ofile.open(of,std::ios::app);
	                if (! ofile)
	  {ofile.open(of,std::ios::out);
          };
	
	 ofile << "# ***** end of event *******\n";

	ofile.close();
	
	       
	std::ofstream fragfile;
	fragfile.open(off,std::ios::app);
	if (! fragfile)
	  {	     
	     fragfile.open(off,std::ios::out);
	  };
	fragfile.precision(19);
	      do
	  {	     
	     qgstring1 = stack2->get();	     
	     fragfile << qgstring1->GetMomentumLab(0,0,qgstring1->GetNgrid()) << " " << qgstring1->GetMomentumLab(1,0,qgstring1->GetNgrid()) << " ";
	     fragfile  << qgstring1->GetMomentumLab(2,0,qgstring1->GetNgrid()) << " " << qgstring1->GetMomentumLab(3,0,qgstring1->GetNgrid()) << std::endl;
   	     delete qgstring1;
	  }
	while(not(stack2->isEmpty())); // ALL FOR EMPTY STACK
	
	
	fragfile << "#******** End of event **********" << std::endl;
	        fragfile.close();
	
	
	
	
	
	
     
     };
   //end simulation cycle
   
};




/* Regge-form distribution of valence quark */

double fdistr(double z)
{   
      double ff;
      if (z<0.01) ff = 2*sqrt(0.01)/0.01 ; else ff=pow(sqrt(1-z),3)/sqrt(z);
      return ff;
};

double Fdistr(double z)
{   
      double FF;
      if(z<0.01) FF=2*sqrt(z);
      return FF;
};










double RandGauss()
{   
        double r;
        double v1,v2,fac;
        do
     {	
	            v1 = 2.0 * rand()/RAND_MAX - 1.0;
	            v2 = 2.0 * rand()/RAND_MAX - 1.0;
	            r = v1*v1 + v2*v2;
     }   
       while ( r > 1.0 );
   
        fac = sqrt(-2.0*log(r)/r);
        return v2*fac;   
};

