/* -----------------------------------------------------------------------------
 * This file contains the constructors and functions for the parton class and
 * its subclasses.
 * Name:    partons.cpp
 * Author:  Bart Verouden
 * Version: 1.01- 20 October 2009
 * Uses:    rng.cpp (and header file)
 *
 * Includes:
 *  -   class parton
 *  -   class quark (subclass of parton)
 *  -   class gluon (subclass of parton)
 * -----------------------------------------------------------------------------
 */
#include "partons.h"
#include "rng.h"
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "vec4d.h"
using namespace std;

#define PI 3.1415926535897932385
#define DEFAULTSEED 123456789L

// set to true to enable debugging
bool debug=false;

ofstream graph("shower.nb");

/* -----------------------------------------------------------------------------
 * Constructors and destructors
 * -----------------------------------------------------------------------------
 */
 
/*
 * -----------------------------------------------------------------------------
 * The parton contructor requires as input the maximum virtuality available to
 * the parton, as well as a pointer to the parent parton, which must be set NULL
 * if the parton is an initial parton. The constructor initializes all variables
 * to zero (exept for t), sets the pointer to parent, and initializes the RNG if 
 * not done so already. 
 * -----------------------------------------------------------------------------
 */
parton::parton(double tmax, parton *parent)
{
	long testseed;
	
	/* Initialize variables.*/	 
	phi=theta=z=0.0;
	p = *new vec4d();   
	
	/* Set parent if the parton isn't an initial parton*/
    if(parent) this->parent = parent;
    else {
		this->parent = NULL;
        /* Initialize random number generator if not initialized yet */
        GetSeed(&testseed);
		if(testseed == DEFAULTSEED){
			if(debug) cout << "initializing RNG" << endl; 
			PutSeed(-1); // The seed -1 generates a seed depending on the time.
		}
    }
	
	/* initializing child pointers to NULL */
	qchild=NULL;
	gchild=NULL;
}

/* -----------------------------------------------------------------------------
 * The gluon constructor calls the parton constructor with the same values. No
 * futher branching from a gluon is allowed at this time, hence the values of
 * the pointers to the children are set to NULL
 * -----------------------------------------------------------------------------
 */
gluon::gluon(double tmax, parton *parent) :
	parton(tmax, parent)
{
	gchild = qchild = NULL; // Gluons are not allowed to branch
	t=0.;					// and are therefore chosen to be on shell
}


/* -----------------------------------------------------------------------------
 * The quark constructor is the core of the final state shower and calls the 
 * parton constructor with the same values. A virtual mass t<tmax and the
 * splitting variable z:z0<z<1-z0 are generated by veto techniques. After that
 * the kinematic variables are calculated, as well as the gluon sisters and if
 * t>tmin, a new quark-gluon pair is created.
 * -----------------------------------------------------------------------------
 */
quark::quark(double tmax, parton *parent) :
	parton(tmax, parent)
{
	extern double A(double z0), z0(double t, double E), P(double z);
	extern double tmin, Q2;
	double told, E, z00;
	int snob;
	
	/* Keep track of number of branches in the shower.*/
	if(parent){	
		snob = nob;
		nob = snob + 1; // Increase number of branches
	}
	else 
		nob = 0; // Reset number of branches
	
	if(debug) 
		cout << "Generating new quark\n";
	
	/* Start the Veto-Algorithm.*/
	double scale = (parent ? parent->getz()*parent->getE() : tmax);
	t = tmax;
	do{
		if(tmax<tmin){ //This happens if the particles energy is less then tmin.
			t = 0.;
			break;
		}
		told = t; // Set the (new) upper limit for the allowed t.
		do 	t = Ginv(log(Random()) + G(told, scale), scale);   // Pick a t.....
		while(t >= told);							 		   // .......<told.
		if(t < tmin) break;
		/* Calculate energy (corresponding to the selected t) */
		if(!parent)
			E = (Q2 + t) / (2. * sqrt(Q2));
		else
			{E = parent->getz() * parent->getE();}
		/* 
		 * Calculate splitting variable z in the interval (z0,1-z0) 
		 * corresponding to the selected t and E. z is randomly choosen, and
		 * accepted with a chance proportional to the splitting function P(z).
		 */
		z00 = z0(t,E);
		do 
			z= z00 + Random()*(1.-2.*z00);
		while (P(z) / P(1.0 - z00) < Random());
	/* Veto if f/g < R, or if the parent's z falls out of its new boundaries. */
	} while(f(t, tmax) / g(t, scale)<Random() || (parent ? !checkz() : false));

	if(debug) cout << "generated t = "<< t << endl;
	
	if(t < tmin)
		t = 0.;	// put the quark on shell if the virtuality is to low.
	
	/* Adjust parent's splitting variable if the daughter quark got a mass. */
	else if(parent){
		if(debug) cout << "z:" << parent->getz();
		parent->setz(parent->transformz(t));
		if(debug) 
			cout<< " -> "<< parent->getz()<< ", "<< "z0 = "<< z0(t, E)<< endl;
	}
	
	/* Set Kinematics */
	if(parent) {
		double zeff = parent->getz();// Get parent's final z.
		p.E = zeff * parent->getE(); // Set quark's final E.
		p.z = sqrt(p.E*p.E-t);		 // Set quarks 3 momentum in it's own frame.
				
		double Es = (1.0 - zeff) * parent->getE();	// Sister gluon's Energy.
		double ps = Es;								// Sister gluon's |p|.
		
        phi = 2.0*PI*Random();	// Pick phi random from flat distribution.
		
		getsister()->setphi(phi); // Set sister's phi the same as this phi.
		
		/* 
		 * Calculate the angle between the new quark and gluon from 4(or 3)-
		 * momentum conservation, in order to calculate the quark's theta 
		 * and the gluon's theta from transverse momentum conservation.
		 */
        double costhn = 2.0*p.E*Es-parent->gett()+t; 
		double costhd = 2.0*p.z*ps;
		double costh = costhn/costhd;
		if(debug){
			cout <<" costh = " << costhn << "/ " << costhd << " = ";
			cout <<costh<<endl;
		}
		double sinth = sqrt(1.0-costh*costh);
		double xi=(1.0/acos(costh))*atan(ps*sinth/(p.z+ps*costh));
		
		/* Set the quark and the sister gluon's theta */
		theta = fabs(xi*acos(costh)); 
		getsister()->settheta(-(1.0-xi)*acos(costh));
		getsister()->setE(Es); // Set the sister's energy.
		getsister()->setpz(ps);// Set the sister's 3 momentum in it's own frame.
	} else {
		p.E = (Q2+t)/(2.*sqrt(Q2));
		p.z =(Q2-t)/(2.*sqrt(Q2));
		theta = 0.0;
	}
	
	/* Transform the 3-momenta to the COM frame.*/
	transformcoordinates();	
	if(parent)
		getsister()->transformcoordinates();
	
	/* Create a new quark-gluon pair if t>tmin */
	if (t>tmin) {
		gchild = new gluon(fmin(t,(1.-z)*(1.-z)*E*E), this);		
		qchild = new quark(fmin(t,z*z*E*E), this);
	} else {
		qchild = NULL;
		gchild = NULL;	
	}
}

quark::~quark() {
	delete qchild;
	delete gchild;
}

gluon::~gluon() {}

/* -----------------------------------------------------------------------------
 * The implementations of the parton's virtual getsister() function that returns
 * it's sister parton, which is a quark for the gluon and vica versa.
 * -----------------------------------------------------------------------------
 */
parton *quark::getsister(void) {
	return parent->getgchild();
}

parton *gluon::getsister(void) {
	return parent->getqchild();
}

/* -----------------------------------------------------------------------------
 * The implementations of the parton's overloaded virtual trackpt2 function
 * which tracks the distribution of the transverse momenum of the emitted gluons
 * either in an external vector, or in a DynBinDat object that is passed to the
 * function by reference.
 * -----------------------------------------------------------------------------
 */
void quark::trackpt2(DynBinDat *pt) {
	if(qchild) {
		qchild->trackpt2(pt);
		gchild->trackpt2(pt);
	}
}

void gluon::trackpt2(DynBinDat *pt) {
	pt->add(pt2());
}

/* -----------------------------------------------------------------------------
 * The implementations of the partons virtual print() function. It prints
 * information about the quark/gluon and calls the print function of it's
 * daughters.
 * -----------------------------------------------------------------------------
 */
void quark::print(void){
	cout<<(parent?(t>0?"A virtual ":"A final "):"\nThe initial ");
	cout<<"quark was emitted with ";
	printf("E=%f, phi=%f, theta=%f, t=%f\n",p.E, phi, theta, t);
	printf("Momentum: ( %9f, %9f, %9f )\n", p.x,p.y,p.z);
    if (gchild) gchild->print();
    if (qchild) qchild->print();
}

void gluon::print(void){
	printf("A gluon was emitted with: E=%10f, phi=%f, theta=%f\n"
		,p.E, phi, theta);
	printf("Momentum: ( %9f, %9f, %9f )\n", p.x,p.y,p.z);
    if (gchild) gchild->print();
    if (qchild) qchild->print();
}

/* -----------------------------------------------------------------------------
 * The parton::transformcoordinates function transforms coordinates (px, py, pz)
 * in its own frame to it's parents frame, an then calls the
 * transformcoordinates function of the parent, using the new coordinates as
 * input. So effectively, transformcoordinates transforms the coordinates to the
 * initial (C.O.M.) frame.
 * -----------------------------------------------------------------------------
 */
void parton::transformcoordinates(double *px, double *py, double *pz)
{
    double xt, yt, zt;

    xt= cos(theta)*cos(phi)**px-cos(theta)*sin(phi)**py+sin(theta)*cos(phi)**pz;
    yt= sin(phi)**px + cos(phi)**py + sin(theta)*sin(phi)**pz;
    zt= -sin(theta)*cos(phi)**px + sin(theta)*sin(phi)**py + cos(theta)**pz;

    *px=xt;
    *py=yt;
    *pz=zt;

    if(parent) parent->transformcoordinates(px,py,pz);
}

/* -----------------------------------------------------------------------------
 * The double parton::transformz function returns the transformed z-value where
 * the masses of its daughters have been taken into account.
 * -----------------------------------------------------------------------------
 */
double parton::transformz(double td)
{
	return (z-0.5)*(t-td)/t+(t+td)/(2.*t);
}

/* -----------------------------------------------------------------------------
 * The bool quark::checkz function checks if for a generated t, the generated z
 * still lies between z0 and 1-z0 which are dependent on the energery, which
 * gets changed when the quark gets a mass.
 * -----------------------------------------------------------------------------
 */
bool quark::checkz(void)
{
	extern double z0(double t, double E);
	double z00 = z0(t,parent->transformz(t)*parent->getE());
	return(z >= z00 && z<=1.-z00);
}


/* -----------------------------------------------------------------------------
 * Returns the transverse momentum squared of the sub-shower initiated by the 
 * parton, where the z-axis aligns with the direction of the first parton.
 * -----------------------------------------------------------------------------
 */
double parton::pt2(void)
{
	if (getqchild())
		return (getqchild()->pt2()+getgchild()->pt2());
	else return (p.x*p.x+p.y*p.y);
}

/* -----------------------------------------------------------------------------
 * double parton::f(double t, double tmax) represents the integrand of the 
 * Sudakov form factor. The implementation of z0 used in the splitting variable
 * integration is different for the first quark then for later quarks. This is
 * because the energy of the first parton is not fixed, but depends on t.
 * -----------------------------------------------------------------------------
 */
double parton::f(double t, double tmax)
{
	extern double z0(double t, double E), A(double z0);
	if(parent)
		return (A(z0(t,parent->getz()*parent->getE())) / t);
	else
		return(A(t / (tmax+t)) / t);
}

/* -----------------------------------------------------------------------------
 * double g(double t, double tmax) serves as an overestimation of the function f
 * used in the veto algorithm. g is not dimensionless (like f) and is to rough 
 * an overestimation for Q2 of order 10000 or higher, hence a better fitting
 * (dimensionless) overestimation should be looked for.
 * -----------------------------------------------------------------------------
 */
double g(double t, double tmax) 
{
	return(0.07*sqrt(tmax)/(t*log(t)));
}

/* -----------------------------------------------------------------------------
 * double G(double t, double tmax) is the primitive of g(t,tmax).
 * -----------------------------------------------------------------------------
 */
double G(double t, double tmax)
{
	extern double tmin;
	return(0.07*sqrt(tmax)*log(log(t)/log(tmin)));
}

/* -----------------------------------------------------------------------------
 * double Ginv(double x, double tmax) is the inverse of G(t, tmax).
 * -----------------------------------------------------------------------------
 */
double Ginv(double x, double tmax)
{
	extern double tmin;
	return pow(tmin, exp(x/(0.07*sqrt(tmax))));
}

/* -----------------------------------------------------------------------------
 * Generate Mathematica code that creates a 3d representation of the shower
 * -----------------------------------------------------------------------------
 */
void parton::graphplot(int i)
{	
	extern int nobmax;
	if(!parent){
		i=1;
		graph << "q0={0,0,0}\n";
		graph << "q1={"<<p.x<<","<<p.y<<","<<p.z<<"}\n";
	} else {
		if(this==parent->qchild){
			graph<<"q"<<i<<"=";
			graph<<"q"<<i-1<<"+"<<"{"<<p.x<<","<<p.y<<","<<p.z<<"}\n";
		} else {
			graph<<"g"<<i<<"=";
			graph<<"q"<<i-1<<"+"<<"{"<<(nobmax-i+2)*p.x<<","<<(nobmax-i+2)*p.y;
			graph<<","<<(nobmax-i+2)*p.z<<"}\n";
		}

	}

	if(qchild){
		gchild->graphplot(i+1);
		qchild->graphplot(i+1);
	} else if(!parent || this==parent->qchild){
		graph<<"tree={a0->a1";
		for(int j=1;j<i;j++)
			graph<<",a"<<j<<"->b"<<j+1<<",a"<<j<<"->a"<<j+1;
		graph<<"}\n";
		graph<<"GraphPlot3D[tree, VertexCoordinateRules ->{a0->q0, a1->q1";
		for(int j=2;j<=i;j++)
			graph<<",a"<<j<<"->q"<<j<<",b"<<j<<"->g"<<j;
		graph<<"}]\n";
		graph.close();
	}
}