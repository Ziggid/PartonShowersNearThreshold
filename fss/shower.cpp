/* -----------------------------------------------------------------------------
 * This program generates final state showers. 
 * The program creates an instance of the quark-class as defined in partons.h
 * which generates a gluon and another quark, the latter creating yet another
 * quark-gluon pair etc., untill the virtuality of a new quark drops below the
 * cutoff value t0. The new quark is then put on mass shell, and the result is a
 * generated (simple) parton shower. The program keeps doing this until it 
 * generates a shower with at least 2 branchings.
 * At the moment the generated gluons are not allowed to branch, and are put on
 * mass shell.
 * After a shower with at leas 2 branchings has been generated its particle's 4-
 * momenta and other properties are printed, and a momentum conservation check
 * is done.
 *
 * Name:    shower.cpp
 * Author:  Bart Verouden
 * Version: 1.01 - 20 October 2009
 * Uses:    partons.cpp (and header files)
 *
 * Includes:
 *  -   int main()          : Creates the parton showers.
 *  -   double P(double z)  : The splitting function.
 *	-	double A(double z0)	: (a_s/2PI) times The integrand over z of the 
 *							  splitting function from z0 to 1-z0.
 *	- 	double z0(double t, double E)
 *							: The cutoff value of z as a function of the virtual 
 *							  mass squared and the energy
 * -----------------------------------------------------------------------------
 */

#include <stdio.h>
#include <math.h>
#include "shower.h"
#include "partons.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
using namespace std;

bool debug2 = false;

#define PI 3.1415926535897932385

/*
 * Global variables used as external variables in partons.cpp.
 */
double tmin =1.01, Q2 = 1000.,binsize = 0.01;
int nobmax=0, quark::nob=0;
int nbins =(int)(0.25*Q2/binsize)+1;

int main () {
	quark * godfather = NULL;    // The initial quark.
	int nsim =100000;

	cout << "Q2 = "<<Q2<<", nsim = "<<nsim<<endl;    

	DynBinDat *pt2= new DynBinDat((int)Q2);
	pt2->add(Q2/4);

	for(int n=0; n<nsim;n++)
	{	
		if((n*100)%nsim==0)
			cout <<(n*100)/nsim<<"%\n";
		delete godfather;
		godfather=new quark(Q2,NULL);
		godfather->trackpt2(pt2);
	}
	
	stringstream ss, ss2;
	ss<<"pt2_cul_as="<<a_s;
	pt2->plotCum(ss.str());
	ss2<<"pt2_as="<<a_s;
	pt2->plotSec(ss2.str());

	// code for generating mathematica 3d plot of the shower.
	// godfather->graphplot(1);

    /*
     * Call the print function of the initial quark, which in turn calls
     * the print function of the children, hence, effectively, printing the
     * shower.
     */
	godfather->print();


	return 0;
}

/* -----------------------------------------------------------------------------
 * The q->qg splitting function
 * -----------------------------------------------------------------------------
 */
double P(double z) {
    return (1+z*z)/(1-z);
}

/*
 * A= (a_s/2*PI)Integrate[P(z),{z,z0,1-z0}]
 */
double A(double z0){
	return (a_s/(2.*PI))*(-3./2.+3.*z0+2.*log(-1.+1./z0));
}

/*
 * z0(t) gives the cutoff value of z as a function of the virtual mass squared, 
 * and the energy.
 */
double z0(double t, double E) {
	return fmax(0.5*(1.0-sqrt(1.-t/(E*E))), 1.0e-6);
}
