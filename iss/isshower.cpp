/* -----------------------------------------------------------------------------
 * This file generates initial state (spacelike) showers by the use of the
 * isparton class and its subclasses. The varialbles defined in the section
 * Configurable variables can be changed by creating a configuration file and
 * passing its path as an argument when executing the program. These values need
 * to be assigned. Also an global pointer to a c_mstwpdf instance named
 * pdf must be declared and initialized here, as well as all static members of
 * isparton and it's subclasses.
 * This specific file uses a file bindat.h, which contains the classes BinDat
 * and DynBinDat that I created, and enables the user to quickly and easily
 * distibute data over bins and export and/or plot the distributions.
 * Name:    isshower.cpp
 * Author:  Bart Verouden
 * Version: 1.1 - 28 May 2010
 * Uses: 	isparton.h 		(contains the parton classes)
 *			mstwpdf.h 		(mstw parton distribution functions)
 *			bindat.h 		(bin-tool)
 *			ConfigFile.h	(handles configuration files)
 *
 * Includes:
 * 	-	int main(int argc, char* argv[])
 *	-	class isshower
 *	-	void selectz(void)
 *	-	void outconf(void)
 * -----------------------------------------------------------------------------
 */
#include <cmath>
#include <sstream>
#include <ctime>
#include <fstream>
#include "isparton.h"
#include "isshower.h"
#include "mstwpdf.h"
#include "bindat.h"
#include "vec4.h"
#include "rng.h"
#include "ConfigFile.h"


using namespace std;

/* Configurable variables (Default values entered) */
double as=0.2;		// The (fixed) strong coupling constant
double tmin = 1.;	// The cutoff virtuality
double stot=1e8; 	// Beam center of mass energy squared in GeV^2
double s_hat=1e4;	// Center of mass energy of the interacting partons.
double res=2.;		// The minimum energy for emitted gluons to have to be
					// considered resolvable.
bool do_iss = true;	// simulate the initial state shower.
bool do_thr = false;// simulate the threshold shower.
int nos = 100000;	// number of showers to simulate.
bool trackpt = true;// track the transverse momenta of the  gluons.
bool trackE = false;// track the energy of the gluons in the lab frame.
string output = "out";// set the output folder for plots and data.

/* Global variables */
double maxfactor = 0;// Track the correctness of the overestimation for the
					// Veto algorithm
c_mstwpdf *pdf; 	// MSTW2008 parton distribution functions
double maxE;		// The maximum gluon-energy emitted by the shower.



/* Declaration of static class members */
int 	isquark::nob=0;
double 	isquark::x2;
Vec4	isparton::p_1;
Vec4 	isparton::p_2;
bool 	isparton::kinbreak;
bool	isshower::kinbreak=false;
RotBstMatrix isparton::history;
int isparton::gluoncount = 0;


/* -----------------------------------------------------------------------------
 * Main function - As an argument the user can specify the path to the 
 * configuration he wants to use.
 * -----------------------------------------------------------------------------
 */
int main(int argc, char* argv[]) 
{
	extern double sfty;
	bool redone = false;
		
	if(argc==1){
		cout<<"No configuration file specified. Using default values:\n";
		outconf();
	}
	else {
		try{
			ConfigFile input(argv[1]);
			cout << "Using \'" << argv[1] <<"\' as configuration file.\n";
			input.readInto(as, "as");
			input.readInto(tmin, "tmin");
			input.readInto(stot, "stot");
			input.readInto(s_hat, "s_hat");
			input.readInto(res, "res");
			input.readInto(do_iss, "do_iss");
			input.readInto(do_thr, "do_thr");
			input.readInto(nos, "nos");
			input.readInto(trackpt, "trackpt");
			input.readInto(trackE, "trackE");
			input.readInto(output, "output");
		}
		catch(ConfigFile::file_not_found) {
			cout<< "Error: Config file not found at \""<< argv[1] <<"\"\n";
			cout<< "Using default configuration:\n";
			outconf();
		}
	}
	
	/* Initialize the parton distribution functions */
	char filename[] = "mstw2008nlo.00.dat";
	pdf = new c_mstwpdf(filename);	
	
	isshower *iss = NULL;
	isshower *news= NULL;	
	
	DynBinDat *pt = NULL;
	DynBinDat *ptn= NULL;
	DynBinDat *y_nor = NULL;
	DynBinDat *y_thr = NULL;	
	DynBinDat *z_nor = NULL;
	DynBinDat *z_thr = NULL;	
	
	double x1, x2, z;
	x1 = x2 = sqrt(s_hat / stot);	
	
	do{
		if(do_iss) {
			if(trackpt)
				pt = new DynBinDat(1000);
			if(trackE) {
				y_nor = new DynBinDat(100);
				z_nor = new DynBinDat(100);
			}
		}
		if(do_thr) {
			if(trackpt)
				ptn= new DynBinDat(1000);
			if(trackE) {
				y_thr = new DynBinDat(100);
				z_thr = new DynBinDat(100);
			}
		}		
	
		maxfactor=0;
		
		cout << endl;
		for(int i = 0; i < nos; i++) {
 			if((i * 100) % nos == 0) {
 				cout << "\r" << (i * 100) / nos << "%"<<flush;
			}
			if(do_iss) {
				/* Do a 'normal' initial state shower. */
				do{
					delete iss;
					iss = new isshower(s_hat, x1, x2);
				}
				while(isshower::kinbreak);
				if(trackpt)
					iss->ptdist(pt);
				if(trackE)
					iss->Edist(z_nor,y_nor);
				if(i == nos-1){
					cout << "\nInitial state shower\n\n";
 					cout <<endl<< *iss << endl<<endl;
				}
				delete iss;
				iss = NULL;
			}
			
			if(do_thr) {
				/* Do a threshold adapted initial state shower. */
				do{
					delete news;
					news = NULL;
					z=selectz();
					news = new isshower(s_hat*pow2(1-z), x1, x2);
				}
				while(isshower::kinbreak);
				if(trackpt)
					news->ptdist(ptn);
				if(trackE)
					news->Edist(z_thr, y_thr);
				if(i == nos-1){
					cout << "\nThreshold shower:\n\n";
 					cout <<endl<< *news<<endl;
				}
				delete news;
				news = NULL;				
			}
		}
		
		if(maxfactor>1) {
			if(redone) {
				cout<<"Encountering too large momentum fractions. Settings out";
				cout<<" of scope of program.\n";
				cout<<"Exiting program. Results may be unreliable because the ";
				cout<<"overestimation in veto-\nalgorithm is not valid at all ";
				cout<<"times.\n\n";
				goto prints;
			}		
			delete pt;
			pt = NULL;	
			delete ptn;
			ptn = NULL;
			cout << "Safety factor to low, increasing sfty from " << sfty;
			sfty *= (4./3. * maxfactor - 1./3.);
			cout << ", to " << sfty << endl;
			redone = true;

		}
	} while(maxfactor>1);

	/* output statements */
	// Switching working dir to output folder
	stringstream ss;
	ss << "mkdir "<<output;
	if(chdir(output.c_str())!=0) {
		system(ss.str().c_str());	
		chdir(output.c_str());
	}
	
	// Generating output
	if(do_iss && do_thr) {
		/* 
		 * Set the sizes of the bins for the data of the normal and the initial
		 * state shower equal, so they can be compared straighforwardly.
		 */
		if(trackpt)
			eqRange(pt, ptn);	
		if(trackE) {
			eqRange(y_nor, y_thr);
			eqRange(z_nor, z_thr);
		}
	}
	if(do_iss) {
		if(trackpt)
			pt->plotSec("pt_nor");
		if(trackE) {
			y_nor->plotSec("y_nor");	
			z_nor->plotSec("z_nor");
		}
	}	
	if(do_thr) {
		if(trackpt)		
			ptn->plotSec("pt_thr");
		if(trackE) {
			y_thr->plotSec("y_thr");
			z_thr->plotSec("z_thr");
		}
	}

	delete pt;
	pt = NULL;
	delete ptn;
	ptn = NULL;
	
	delete pdf;
	pdf = NULL;
}

/* -----------------------------------------------------------------------------
 * Selectz() generates a 'z' according to the appropriate threshold resummed 
 * distribution.
 * -----------------------------------------------------------------------------
 */
double selectz()
{
	long testseed;
	GetSeed(&testseed);
	if(testseed == DEFAULT)
		PutSeed(-1); // The seed -1 generates a seed depending on the time.
	
	double z;
	do
		z = Random();
	while(exp(-.5*pow2(log(1-z))) < Random());
	return z; 
}

/* -----------------------------------------------------------------------------
 * Print the current configuration to standerd output.
 * -----------------------------------------------------------------------------
 */
void outconf(void)
{
	cout << endl;
	cout << "a_s = " << as<<endl;
	cout << "tmin = "<<tmin<<endl;
	cout << "s_hadron = " <<stot<<endl;
	cout << "Q^2 = " << s_hat << endl;
	cout << "generating "<< nos <<" events\n";
	cout << "resolution = "<< res << endl;
	if(trackpt)
		cout << "tracking transverse momenta\n";
	if(trackE)
		cout << "tracking energies\n";
	if(do_iss && do_thr)
		cout << "simulating normal and threshold shower\n";
	else {
		if(do_iss)
			cout<< "simulating initial state shower\n";
		if(do_thr)
			cout<< "simulating threshol shower\n";
	}
	cout << endl;
}

/* -----------------------------------------------------------------------------
 * The isshower constructor backwardly generates an initial state shower of
 * quark - antiquark pair coming from a proton - antiproton collision. Both
 * quark lines are evolved seperately from scale tmax down to tmin. 
 * The kinematics, however, are filled in later, where the order of kinematical
 * evaluation is equal to the decreasing order in virtuality.
 * -----------------------------------------------------------------------------
 */
isshower::isshower(double tmax, double x1, double x2)
{
	isshower::kinbreak = false; // Track kinematical breakdown.
	
	K=0;
		
	/* 
	 * Generate sequences of virtualities at which branchings should take place
	 * for both quark lines.
	 */
	q1=new isquark(tmax, x1, NULL);
	q2=new isquark(tmax, x2, NULL);
	
	// Construct the kinematics for the two quarks entering the hard scattering
	q1->construct(q2);
	q2->construct(q1);
	
	// Reset the history of rotations and boosts.
	isparton::history.reset();
	
	double t1,t2;
	t1 = q1->current()->gett();
	t2 = q2->current()->gett();
	
	/*
	 * Pick the (next) branching that is closest to the hard scattering and
	 * construct the kinematics.
	 */
	while(t1 > 0 || t2 > 0) {
		if(t1 > t2){
			if(q1->current()->mother->construct(q2->current())) 
			{
				isshower::kinbreak = true;// Shower was kinematically incorrect,
				return;					  // abort.
			}
		}
		else
			if(q2->current()->mother->construct(q1->current())) 
			{
				isshower::kinbreak = true;// Shower was kinematically incorrect.
				return;					  // abort.
			}
			
		/*
		 * After constructing the kinematics of a branching, the system is
		 * boosted and rotated to the COM-frame of the quarks in both quark 
		 * lines that are to become unresolved. This is done using the 
		 * RotBstMatrix class,.as defined in vec4.h, from PYTHIA 8.1, and 
		 * applying that matrix to all 4-vectors.
		 */
		RotBstMatrix M; 	// declare matrix.
		// initialize matrix
		M.toCMframe(q1->current()->getp(), q2->current()->getp());	
		// apply matrix to all other 4-momenta.
		RotnBoost(q1->current(), NULL, M);
		RotnBoost(q2->current(), NULL, M);
		// Update the history of boosts and rotations.
		isparton::history.rotbst(M);
		
		t1 = q1->current()->gett();
		t2 = q2->current()->gett();
	}
	
	/*
	 * After the shower is fully unresolved, boost back to the COM frame of the
	 * two colliding protons. This will mean that the orientation of the two 
	 * colliding partons, and thus that of the outgoing photon, will not be 
	 * trivial, and wil not be oriented along the z axis.
	 */	
	// Create the four vector of the protons	
	Vec4 pp1 = q1->current()->getp()/q1->current()->getx();
	Vec4 pp2 = q2->current()->getp()/q2->current()->getx();
	RotBstMatrix M;				     // declare boost matrix
 	M.toCMframe(pp1, pp2);	         // initialize boost matrix
	RotnBoost(q1->current(), NULL, M); // apply matrix to all other 4-momenta.
	RotnBoost(q2->current(), NULL, M); 
	isparton::history.rotbst(M);  // Update the history of boosts and rotations.
}

/* -----------------------------------------------------------------------------
 * isshower destructor.
 * -----------------------------------------------------------------------------
 */
isshower::~isshower()
{
	delete q1;
	delete q2;
	q1=q2=NULL;
}

/* -----------------------------------------------------------------------------
 * Returns the four momentum of the outgoing photon.
 * -----------------------------------------------------------------------------
 */
Vec4 isshower::photon(void)
{
	return q1->getp()+q2->getp();
}

/* -----------------------------------------------------------------------------
 * Get the transverse momentum squared of all gluons in the shower and add
 * them to the supplied DynbinDat object to generate a distribution.
 * -----------------------------------------------------------------------------
 */
void isshower::ptdist(DynBinDat *dat)
{
	q1->current()->ptdist(dat);
	q2->current()->ptdist(dat);
}

void isshower::Edist(DynBinDat *z, DynBinDat *y)
{
	double maxE=0;
	q1->current()->Edist(z, maxE);
	q2->current()->Edist(z, maxE);
	
	if(maxE>0)
	 	y->add(maxE);
}

/* -----------------------------------------------------------------------------
 * Print a textual representation of the initial state shower to the ostream os.
 * -----------------------------------------------------------------------------
 */
ostream& operator << (ostream& os, const isshower& iss){
	isquark *temp = iss.q1->current();
	os << "Quark from proton\n";
	os<<"___________________________________________________________________\n";
	do{
		os << *temp << endl;
		if(temp->daughter != NULL)
			os << *temp->gdaughter << endl;
	} while((temp = temp->daughter) != NULL);
	os<<"xXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXx\n";
	os<<".........................Hard Scattering...........................\n";
	os<<"xXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXx\n";
	os<<"Anti-quark from anti-proton\n";
	os<<"___________________________________________________________________\n";
	temp = iss.q2->current();
	do{
		os << *temp << endl;
		if(temp->daughter != NULL)
			os << *temp->gdaughter << endl;
	} while((temp = temp->daughter) != NULL);	
	return os;
}
