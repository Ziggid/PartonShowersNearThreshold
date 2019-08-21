/* -----------------------------------------------------------------------------
 * The class isshower is the base of the initial state showering algorithm for
 * p-pbar collisions.
 * Name:    isshower.h (header file of isshower.cpp)
 * Author:  Bart Verouden
 * Version: 1.1 - 28 May 2010
 * Uses: 	vec4.h (four-vector handling from Pythia)
 *			bindat.h (binning tool)
 *			isparton.h (contains initial state partons)
 *
 * Includes:
 *  -   class isparton
 *  -   class isquark (subclass of isparton)
 *  -   class isgluon (subclass of isparton)
 * -----------------------------------------------------------------------------
 */
#ifndef _ISSHOWER_
#define _ISSHOWER_

#include "isparton.h"
#include "bindat.h"
#include "vec4.h"
#include <iostream>
using namespace std;

class isshower
{
	protected:	
	Vec4 K;		// Sum of all emitted gluon four vectors.
	
	public:
	/* Constructors and destructors */
	isshower(double tmax, double x1, double x2);
	~isshower();
	
	/* General public functions */
	void ptdist(DynBinDat *dat);
	void Edist(DynBinDat *z, DynBinDat *y);
	Vec4 photon(void);

	isquark *q1, *q2; // The two colliding parton lines.	
	
	/* Static variables */
	static bool kinbreak;	
	
	/* Friends */
	friend ostream& operator<<(ostream&, const isshower&);
};
#endif // _ISSHOWER_

double selectz(void);
void outconf(void);
