/* -----------------------------------------------------------------------------
 * The BinDat and DynBinDat class provides a tool to automatically and easily 
 * bin and manipulate generated data. The BinDat constructor requires as input
 * a minimum and maximum value, as well as the number of required bins, whereas
 * the DynBinDat constructer only requires the number of bins as an input. Since
 * DynBinDat rebins the data when an input value falls out of the current 
 * minimum-maximum range, the dynamic version is somewhat slower, which might
 * become a problem for large sets of data.
 *
 * BinDat allows you to 
 *	- add data without actively having to bin it. (add)
 *	- get the vector of bins as an output (get)
 *	- generate a text file with frequencies (out)
 *	- generate and run a gnuplot file of the data (plot)
 * All output methods (get, out and plot) can be called returning either the 
 * plain frequencies, or a normalized distribution by adding the "Sec" resp. 
 * "Norm" suffix.
 *
 * Name:    BinDat.h
 * Author:  Bart Verouden
 * Version: 1.3 - 28 April 2010
 *
 * Includes:
 *  -   class BinDat
 *  -   class DynBinDat (subclass of BinDat)
 * -----------------------------------------------------------------------------
 */
#ifndef _BINDAT_
#define _BINDAT_
#include <vector>
#include <string>
#include <fstream>
using namespace std;

class BinDat
{
protected:
	vector<int> data;
	double 		binsize;
	int 		nentries;
	double 		xmin, xmax;
	int 		nbins;	
	int 		missed;
	
	// Default constructer - Should not manually be called, hence protected.
	BinDat(void){nentries = 0; xmin=xmax=0;}	

private:
	ofstream output;

public:
	// Constructor
	BinDat(int nbins, double xmi, double xma);

	
	/* public functions */
	void add(double x);
	
	void outSec(string filename);
	void outNorm(string filename);
	void outCum(string filename);
	void outWeighted(string filename, double (*weight)(double));
	void outSmooth(string filename, int n);
	
	vector <int> getSec(void);
	vector <double> getNorm(void);
	vector <int> getCum(void);
	vector <double> getWeighted(double (*weight)(double));
	vector <double> getSmooth(int n);
	
	void plotSec(string name);
	void plotNorm(string name);	
	void plotCum(string name);
	void plotWeighted(string name, double (*weight)(double));
	void plotSmooth(string name, int d);
};

class DynBinDat : public BinDat
{
private:
	vector<double> ddata;	
	void refill(void);
	
public:
	// Constructor
	DynBinDat(int nb){nbins=nb; xmin=xmax=0;}
	
	/* public functions */	
	void add(double x); // Overloads the BinDat function.
	vector<double> getEntries(void){return ddata;} 	// Return the original input
													// as a vector.
	
	void rescale(int nbins, double xmi, double xma);
	
	/* Friends */
	friend void eqRange(DynBinDat *, DynBinDat *);
};

#endif //_BINDAT_
