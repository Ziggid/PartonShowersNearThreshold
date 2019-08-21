/* -----------------------------------------------------------------------------
 * The class parton is the abstract base class of the subclasses quark and
 * gluon. It contains the pure virtual function print(), making it an abstract
 * class. These classes are made to generate a final state parton shower.
 * Name:    partons.h (header file of partons.cpp)
 * Author:  Bart Verouden
 * Version: 1.01 - 20 October 2009
 *
 * Includes:
 *  -   class parton
 *  -   class quark (subclass of parton)
 *  -   class gluon (subclass of parton)
 * -----------------------------------------------------------------------------
 */
#ifndef _PARTON_
#define _PARTON_
#define PI 3.1415926535897932385
#include "vec4d.h"
#include "bindat.h"

double g(double t, double tmax);
double G(double t, double tmax);
double Ginv(double x, double tmax);

class parton
{
    protected:
        double phi;             // The azimuthal angle w.r.t. the parent frame
        double theta;           // The polar angle w.r.t. the parent frame
        double t;               // The virtuality of the parton
		double z;



        /*
         * Pointers to child and parent partons. These are NULL if nonexisting.
         */
        class parton *qchild;
        class parton *gchild;
        class parton *parent;

    public:
		vec4d p;				// the 4-vector in the COM-frame	
		/*
		 * Constructors and destructors
		 */
		parton(double tmax, parton *parent);
		virtual ~parton(){};
		
		/*
		 * General public functions
		 */
		void transformcoordinates(double *px, double *py, double *pz);
		void transformcoordinates(void){ transformcoordinates(&p.x, &p.y,&p.z);}
		double transformz(double td);
        virtual void print(void) =0;
		double f(double t, double tmax);
		double pt2(void);
		void graphplot(int i);
		virtual void trackpt2(DynBinDat *pt){};
		
		/*
		 * Getters
		 */
		double gett(void){ return this->t;}
		double getE(void){ return p.E;}
		double getz(void){ return z;}
		parton *getgchild(){ return gchild;}
		parton *getqchild(){ return qchild;}
		virtual parton *getsister()=0;
		/*
		 * Setters
		 */
		void setphi(double phi){ this->phi=phi;}
		void setz(double z){ this->z=z;}
		void settheta(double theta){ this->theta=theta;}
		void setE(double E){p.E=E;}
		void setpz(double pz){p.z=pz;}
};

class quark : public parton
{
    public :
		static int nob; // nob is the number of branches in the shower.
	
		/*
		 * Constructors and destructors.
		 */	
		quark(double tmax, parton *parent);
		~quark();
		
		/*
		 * Implementations of virtual functions inhereted from parton class.
		 */
		parton *getsister();
    	void print(void);
		void trackpt2(DynBinDat *pt);
		
		/*
		 * General public functions.
		 */
		bool checkz(void);
};

class gluon : public parton
{
    public :
		/*
		 * Constructors and destructors.
		 */
		gluon(double tmax, parton *parent);
		~gluon();
		
		/*
		 * Implementations of virtual functions inhereted from parton class.
		 */
		parton *getsister();
    	void print(void);
		void trackpt2(DynBinDat *pt);
};

#endif // _PARTON_

