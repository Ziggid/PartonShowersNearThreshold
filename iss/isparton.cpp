/* -----------------------------------------------------------------------------
 * This file contains the constructors and functions for the isparton class and
 * its subclasses.
 * Name:    ispartons.cpp
 * Author:  Bart Verouden
 * Version: 0.2 - 4 December 2009
 * Uses: 	isparton.h (header file)
 *			mstwpdf.h (mstw parton distribution functions)
 *			rng.h (randon number generator)
 *
 * Includes:
 *  -   class isparton
 *  -   class isquark (subclass of isparton)
 *  -   class isgluon (subclass of isparton)
 *	-	double A(double x, double zmax) 	
 *			: A represents the z-integral of the overestimation from x to zmax. 
 *	-	double G(double t, double x, double zmax)
 *			: G is the primitive of the overestimation of the t-integrand.
 *	-	double Ginv(double y, double x, double zmax)
 *			: Ginv is the inverse of G.
 * -----------------------------------------------------------------------------
 */
#include "isparton.h"
#include "mstwpdf.h"
#include "rng.h"
#include <cmath>
#include <vector>

double sfty=1.	;// A safety factor to ensure that the overestimation of the
				 // exponent in the Sudakov form factor is indeed an 
				 // overestimation. Should be adjusted according to s_hat and
				 // stot.

/*
 * Declaration of funtions used in the veto-algorithm
 */
double A(double x, double zmax);
//double P(double z);
//double g(double t, double x, double zmax);
double G(double t, double x, double zmax);
double Ginv(double y, double x, double z);

/* -----------------------------------------------------------------------------
 * Constructors
 * -----------------------------------------------------------------------------
 */

/* -----------------------------------------------------------------------------
 * The isquark constructor requires as input the maximum virtuality, the x value
 * and a pointer to it's daughter quark (NULL if absent), and generates a quark 
 * with t and z, 0=<-m^2=t<tmax and x<z<zmax, by means of a veto-algorithm. If a
 * daughter quark was present, a (timelike) daughter gluon must be created as 
 * well with virtuality t4=m^2<t4_max. The gluon virtuality is set to zero here 
 * (so no timelike sideshowering is generated), but the maximum available 
 * virtuality is stored to check for kinematical consistency.
 *
 * If the quarks generated (and accepted) t is non-zero, the (backward) 
 * showering continues, and a mother quark is created, following the same 
 * procedure. Hence the isquark constructor is the kernel of the backwards
 * showering.
 * -----------------------------------------------------------------------------
 */
isquark::isquark(double tmax, double xx, isquark *d){
	extern c_mstwpdf *pdf;
	
	extern double maxfactor;
	
	z = 1.234;
	x = xx;
	t = tmax;
	daughter = d;
	s = 0;
	mother = NULL;
	gdaughter = NULL;
	isresolved = true;
	
	
	/* Keep track of number of branches in the shower.*/
	int snob;
	if(daughter) {	
		snob = nob;
		nob = snob + 1;
	}
	else { 
		nob = 0; // (re)set the number of branches to zero.
		isparton::kinbreak = false; // reset kinematical breakdown boolean.

		/* 
		 * Initialize the random number generator if this is the first shower 
		 * generated.
		 */
		long testseed;
        GetSeed(&testseed);
        if(testseed == DEFAULT)
			PutSeed(-1); // The seed -1 generates a seed depending on the time.
	}
	
	/* Veto part */
	double z_max = zmax();
	double test;	
	do {
		t = selectt(t); 	// Pick a t.
		if(t == 0.) break; 	// No further showering allowed when t=0.
		double R = Random();		
		/* 
	 	 *We solve 
	 	 *Integrate[(as/2*PI)(1+zmax^2)/(sqrt(z')*(1-z')),{z',x,z}]==R*A(x,zmax)
	 	 */
		// Pick a z.
		z = pow2(tanh((1 - R)*atanh(sqrt(x)) + R * atanh(sqrt(z_max)))); 
	/* 
	 * Veto if the actual z-integrand is smaller than a random fraction R of 
	 * the overestimation.
	 */
	 	test=(1 / sfty ) * ((1 + pow2(z)) / (1 + pow2(z_max))) * sqrt(z) *
		(pdf->parton(2, x / z, t) / pdf->parton(2, x, t));
		maxfactor = fmax(maxfactor, test);
	} while (test < Random());


	if(t > 0) 
		mother = new isquark(t, x / z, this); // Create a new mother quark.
}

/*
 * The destructor of the initial state quark.
 */
isquark::~isquark()
{
	delete mother;
	delete gdaughter;
	mother=NULL;
	gdaughter=NULL;
}



/* -----------------------------------------------------------------------------
 * The isgluon constructor needs as input the maximum virtuality available for a
 * timelike shower initiated for the gluon (which we won't do, but we store the
 * virtuality available), and a pointer to the mother quark.
 * The 4-momentum is set by the mother quark, so nothing needs to be done here 
 * except for storing the two input values.
 * -----------------------------------------------------------------------------
 */
isgluon::isgluon(double tm, isquark *mom)
{
	tmax = tm;
	mother = mom;
}

/* -----------------------------------------------------------------------------
 * The isquark::construct(isquark *q2) function fills in the kinematics of the
 * isquark instance it is a member from. The isquark q2, passed by reference, is
 * the quark from the other quark line that is still resolved with the highest
 * virtuality. The system is thus in the COM frame of q2 and the daughter quark
 * of the quark currenctly being constructed.
 * -----------------------------------------------------------------------------
 */
int isquark::construct(isquark *q2)
{
	extern double stot;
	
	double t2 = q2->gett();
	double x2 = q2->getx();
	s = x * x2 * stot;
	
	/* Fill in the kinematics */
	double t4_max;
	if(daughter == NULL) {
		// This is the situation where the quark partook in the hard scattering
		double pz = 0.5 * sqrt((pow2(s + t + t2) - 4 * t * t2) / s);
		double E  = 0.5 * (s + t2 - t) / sqrt(s);
		if(!(q2->isresolved)) 
			pz = -pz;
		p=Vec4(0.,0.,pz,E);
	} else {
		// Variables used in kinematical reconstruction
		double t1 = daughter->t;
		double s1 = daughter->z*s + t1 +t2;
		double s3 = s + t + t2;
		double r1 = sqrt(s1 * s1 - 4 * t2 * t1);
		double r3 = sqrt(s3 * s3 - 4 * t2 * t);
		
		// calculate maximum virtuality gluon.
		if(t2 == 0.) {
			t4_max = ((t1 / daughter->z) - t) * daughter->z*s *
			(1. / (daughter->z*s + t1) - 1. / (s + t));
		} else {
			t4_max = (s1 * s3 - r1 * r3) / (2 * t2) - t1 - t;
		}
		if(t4_max < 0.) {	// Kinematical breakdown. Note, and abort shower.
			isparton::kinbreak = true;
			return 1;	// Return value for kinematical failure.	
		}

		double E  = 0.5 * (1. / sqrt(daughter->z*s)) * (s - t1 + t2);
		double pz = 0.5 * (1. / daughter->p.pz()) * (s3 - 2 * q2->getp().e()*E);
		double pt = sqrt(t4_max * (0.5 * (s1*s3 + r1*r3) -  t2 * (t1 + t)) 
					/ (r1*r1));
		double phi = 2.*M_PI*Random(); // Generate a random azimuthal angle phi.
		
		p = Vec4(pt * sin(phi), pt * cos(phi), pz, E);	// Set 4-momentum
		
		gdaughter = new isgluon(t4_max, this);	// Create new gluon
		gdaughter->setp(p - daughter->getp());	// Set it's 4-momentum
	}
	isresolved = false;
	return 0;	// Return value for succes.
}

/* -----------------------------------------------------------------------------
 * isquark::selectt(tmax) returns a virtuality t<tmax, with probality distribu-
 * tion g(t)*Integrate(g(t'), {t', t, zmax})
 * -----------------------------------------------------------------------------
 */
double isquark::selectt(double tmax)
{	
	extern double tmin;
	double t = tmax, told, z_max = zmax();
	
	if(z_max < x)
		return 0;

	told = t; // Set the (new) upper limit for the allowed t
	do {	
		t = Ginv(log(Random())+G(told,x, z_max), x, z_max);	// Pick a t.....
	}while(t >= told);							 		// ........<told.
	if(t < tmin) {
		return 0.;
	}
	return t;	
}


/* -----------------------------------------------------------------------------
 * The isquark::current() function returns a pointer to the quark in the shower
 * closest to the hard scattering, that is still resolved.
 * -----------------------------------------------------------------------------
 */
isquark * isquark::current(void)
{
	if(isresolved)
		return daughter;
	else if(t > 0)
		return mother->current();
	else return this;
}

/* -----------------------------------------------------------------------------
 * isquark::zmax() returns the maximum value of z such that if z<zmax, the
 * energy of the gluon is larger than the resolution, res, in the C.O.M. frame 
 * of the two quarks that partook in the hard scattering.
 * -----------------------------------------------------------------------------
 */
double isquark::zmax(void)
{
	extern double stot, res;
	double gamma = history.gamma();
	double xe = 2.* res/ (sqrt(stot) * gamma);
	return(x / (x + xe));
}

/* -----------------------------------------------------------------------------
 * The isquark::ptdist(*dat) adds the absolute value squared of the transverse
 * momenta of the gluons in the shower to the DynBinDat object dat.
 * -----------------------------------------------------------------------------
 */
void isquark::ptdist(DynBinDat *dat)
{
	if(this->daughter != NULL) {
		dat->add(gdaughter->getp().pT2());
		daughter->ptdist(dat);
		++gluoncount;	
	}
	return;
}

/* -----------------------------------------------------------------------------
 * The isquark::Edist(*dat, &maxE) adds the absolute value squared of the 
 * energy of the gluons in the shower to the DynBinDat object z, and updates the
 * maximum emitted energy maxE.
 * -----------------------------------------------------------------------------
 */
void isquark::Edist(DynBinDat *dat, double &maxE)
{
	if(this->daughter != NULL) {
		maxE = fmax(maxE,gdaughter->getp().e());
		dat->add(gdaughter->getp().e());
		daughter->Edist(dat, maxE);
	}
	return;
}



/* -----------------------------------------------------------------------------
 * isquark::pt2max() returns the pt^2 of the gluon with the highest transverse
 * momentum in this quark line.
 * -----------------------------------------------------------------------------
 */
double isquark::pt2max(void)
{
	if(this->daughter != NULL) {
		return(fmax(gdaughter->getp().pT2(), daughter->pt2max()));
	} else
	return 0.;
}


/* -----------------------------------------------------------------------------
 * RotnBoost(q, g, M) applies the matrix M to q and g, and to all q's offspring.
 * -----------------------------------------------------------------------------
 */
void RotnBoost(isquark *q, isgluon *g, RotBstMatrix M)
{
	q->p.rotbst(M);
	if(g != NULL)
		g->p.rotbst(M);
	if(q->daughter != NULL){
		RotnBoost(q->daughter, q->gdaughter, M);
	}	
}

/* -----------------------------------------------------------------------------
 * Quark and gluon's output stream operators.
 * -----------------------------------------------------------------------------
 */
ostream& operator << (ostream& os, const isquark& q) {
	os << "q: p = " << q.p << ", x=" << q.x << ", t=" << q.t << endl; 
	return os;
}

ostream& operator << (ostream& os, const isgluon& g) {
	os << "g: p = " << g.p << ", pt^2 = " << g.p.pT2() << endl; 
	return os;
}

/* -----------------------------------------------------------------------------
 * Functions needed for the veto-algortihm.
 * -----------------------------------------------------------------------------
 */
 
/* -----------------------------------------------------------------------------
 * We overestimate the z-integrand appearing in the Sudakov form-factor:
 * [(1+z^2)/(1-z)]*(1/z)*[pdf(x',t)/pdf(x,t)] by
 * [(1+z_max^2)/(1-z)]*[1/sqrt(z)]. "A" represents the z-integral of this 
 * overestimation from x to zmax.
 * -----------------------------------------------------------------------------
 */
double A(double x, double zmax)
{	
	extern double as;
	return (sfty * (as / M_PI) * (1 + pow2(zmax)) * (atanh(sqrt(zmax)) - 
		atanh(sqrt(x))));
}

/* -----------------------------------------------------------------------------
 * P is the splitting function from the DGLAP equations.
 * -----------------------------------------------------------------------------
 * 
	double P(double z) {
	return ((1 + z * z) / (1 - z));
}*/

/* ----------------------------------------------------------------------------- 
 * g is the overestimation of the t-integrand in the Sudakov form factor.
 * -----------------------------------------------------------------------------
 *
double g(double t, double x, double zmax)
{
	return((1 / t) * A(x, zmax));
}*/

/* -----------------------------------------------------------------------------
 * G is the primitive of the overestimation of the t-integrand, g. 
 * -----------------------------------------------------------------------------
 */
double G(double t, double x, double zmax)
{
	return(A(x, zmax) * log(t));
}

/* -----------------------------------------------------------------------------
 * Ginv is the inverse of G.
 * -----------------------------------------------------------------------------
 */
double Ginv(double y, double x, double zmax)
{
	return(exp(y / A(x, zmax)));
}
