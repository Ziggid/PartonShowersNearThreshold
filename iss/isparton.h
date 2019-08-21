/* -----------------------------------------------------------------------------
 * The class isparton is the base class of the subclasses isquark and isgluon. 
 * These classes are made to generate an initial state parton shower.
 * Name:    ispartons.h (header file of ispartons.cpp)
 * Author:  Bart Verouden
 * Version: 1.0 - 21 May 2010
 * Uses: 	vec4.h (four-vector handling from Pythia)
 *			bindat.h (binning tool)
 *
 * Includes:
 *  -   class isparton
 *  -   class isquark (subclass of isparton)
 *  -   class isgluon (subclass of isparton)
 * -----------------------------------------------------------------------------
 */
#ifndef _ISPARTON_
#define _ISPARTON_

// Define pi if not yet done.
#ifndef M_PI
#define M_PI 3.1415926535897932385
#endif //M_PI

//#include "isshower.h"
#include "vec4.h"
#include "bindat.h"
#include <iostream>
#include <vector>
using namespace std;


class isquark; // Reference to the isquark class needed for isparton.
class isgluon; // Reference to the isgluon class needed for isparton.

class isparton
{
	protected:
	Vec4 p;					// The parton's four momentum (px, py, pz, E)
	isquark *mother;		// Pointer to the parton's mother quark.
	
	static Vec4 p_1;		// The 4-momentum of the final quark.
	static Vec4 p_2; 		// The momentum of the 'other' quark in the hard 
							// scattering.
	
	public:
	static bool kinbreak;	// Boolean representing the kinematical beakdown.
	static RotBstMatrix history;	// RotBstMatrix tracking the history of
									// the rotations and boosts made.
	static int gluoncount;	// Keeps track of the total number of of gluons 
							// emitted.
	
	/* Getters and setters */
	Vec4 getp(void){return p;}
	void setp(Vec4 pp){p=pp;}
	Vec4 getp_2(void){return p_2;}	
		
	/* Friends */
	friend void RotnBoost(isquark *q, isgluon *g, RotBstMatrix M);
	friend class isshower;
};

class isgluon : public isparton
{
	private :
	double tmax;	// The maximum virtual mass available for timelike shower.
					// (timelike showering off gluons not implemented)
		
	public :
	/* Constructors */	
	isgluon(double tmax, isquark *mom);
	
	/* Friends */
	friend ostream& operator<<(ostream&, const isgluon&);	
};

class isquark : public isparton
{
	private :
	double x;	// The fractional COM energy (x*x_2*s=s^hat)
	double z;	// The splitting variable
	double t;	// t=-M^2>=0
	double s;	// the COM energy of the subsystem (s^hat)
	
	
	/* private functions */
 	double selectt(double tmax);
	void restart(void);
	double zmax(void);
	int construct(isquark *q2);	
	
	
	public :
	bool isresolved;
	static int nob;		// The number of branches in the shower.
	static double x2;	// the x variable of the 'other' hard parton in the hard
						// scattering
	isquark *daughter;	// The quarks daugther quark
	isgluon *gdaughter; // The quarks daughter gluon
	
	static int noq;
	
	/* Constructors */
	isquark(double tmax, double xx, isquark *d);
	~isquark();
	
	/* Getters */
	double gett(void){return t;}
	double getz(void){return z;}
	double getx(void){return x;}
	double gets(void){return s;}
	
	/* Setters */
	void sett(double tt){t=tt;}
	void setz(double zz){z=zz;}
	void setx(double xx){x=xx;}
	void sets(double ss){s=ss;}
	
	/* Public functions */
	void ptdist(DynBinDat *dat);
	void Edist(DynBinDat *z, double &maxE);
	isquark * current(void);	
	double pt2max(void);
	
	/* Friends */
	//friend ostream& operator<<(ostream&, const isshower&);
	friend ostream& operator<<(ostream&, const isquark&);
	friend class isshower;	
};

#endif // _ISPARTON_
