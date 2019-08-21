/* -----------------------------------------------------------------------------
 * This file contains the constructors and functions for the BinDat class and
 * its subclass DynBinDat. For more info see the header file.
 * Name:    BinDat.cpp
 * Author:  Bart Verouden
 * Version: 1.3 - 28 April 2010
 * Uses:    header file
 *
 * Includes:
 *  -   class BinDat
 *  -   class DynBinDat (subclass of BinDat)
 * -----------------------------------------------------------------------------
 */
#include "bindat.h"
#include <iostream>
#include <sstream>
#include <cmath>

/* -----------------------------------------------------------------------------
 * Constructor - The BinDat constructor requires a range of of x-values, and
 * the number of bins
 * -----------------------------------------------------------------------------
 */
BinDat::BinDat(int nbins, double xmi, double xma)
{
	xmin=xmi; xmax=xma;
	binsize = (xmax-xmin)/nbins;
	data.assign(nbins,0);
	nentries=0;
	missed = 0;
}

/* -----------------------------------------------------------------------------
 * Add 1 to the bin where x falls in.
 * -----------------------------------------------------------------------------
 */
void BinDat::add(double x) 
{
	if(x>xmax || x<xmin) {
		++missed;
		return;
	}
	if(x==xmax) 
		data.at(nbins-1)++;
	else
		data.at((int)((x-xmin)/binsize))++;
	nentries++;
}

/* -----------------------------------------------------------------------------
 * Write the data to the file <filename>
 * -----------------------------------------------------------------------------
 */
void BinDat::outSec(string filename) 
{
	output.open(filename.c_str());
	if(!output.is_open()) {
		cout << "Could not open file: " << filename;
		cout <<", check for write protection.\n";
	} else
		for(int i=0;i<(int)data.size();i++) {
			output<<(0.5+i)*binsize+xmin<<" "<<data.at(i)<<endl;
		}
	output.close();	
}

/* -----------------------------------------------------------------------------
 * Normalize the data such that (number of entries) x binsize = 1, and write the
 * normalized data to the file <filename>
 * -----------------------------------------------------------------------------
 */
void BinDat::outNorm(string filename)
{
	output.open(filename.c_str());
	if(!output.is_open()){
		cout << "Could not open file: " << filename;
		cout <<", check for write protection.\n";
	} else {
		for(int i=0; i<(int)data.size(); i++) {
			output<<(0.5+i)*binsize+xmin<<" ";
			output<<(double)(data.at(i)/((double)nentries*binsize))<<endl;
		}
	}
	output.close();	
}

/* -----------------------------------------------------------------------------
 * Write the cumulative data to the file <filename>
 * -----------------------------------------------------------------------------
 */
void BinDat::outCum(string filename)
{
	output.open(filename.c_str());
	if(!output.is_open()){
		cout << "Could not open file: " << filename;
		cout <<", check for write protection.\n";
	} else {
		int temp =0;
		for(int i=0;i<(int)data.size();i++) {
			output<<(0.5+i)*binsize+xmin<<" "<<(temp+=data.at(i))<<endl;
		}
	}
	output.close();		
}

/* -----------------------------------------------------------------------------
 * Write the weighed data to the file <filename>
 * -----------------------------------------------------------------------------
 */
void BinDat::outWeighted(string filename, double (*weight)(double))
{
	output.open(filename.c_str());
	if(!output.is_open()){
		cout << "Could not open file: " << filename;
		cout <<", check for write protection.\n";
	} else{
		for(int i=0;i<(int)data.size();i++) {
			double x = (0.5+i)*binsize+xmin;
			output<< x <<" "<<data.at(i) * weight(x)<<endl;
		}
	}
	output.close();		
}


/* -----------------------------------------------------------------------------
 * Write the data to the file <filename> where each bin is avaraged with its n
 * neighbors on either side.
 * -----------------------------------------------------------------------------
 */
void BinDat::outSmooth(string filename, int n)
{
	output.open(filename.c_str());
	if(!output.is_open()){
		cout << "Could not open file: " << filename;
		cout <<", check for write protection.\n";
	} else{
		for(int i=0; i<(int)data.size();++i) {
			double out=0;
			int d=0;
			for(int j=-n; j<=n; ++j) {
				if(i+j>=0 && i+j<(int)data.size()) {
					out +=data.at(i+j);
					++d;
				}
			}
			double x = (0.5+i)*binsize+xmin;
			output<< x <<" "<<out/d<<endl;
		}
	}
	output.close();
}

/* -----------------------------------------------------------------------------
 * Return the data as a vector.
 * -----------------------------------------------------------------------------
 */
vector <int> BinDat::getSec(void){return data;}

/* -----------------------------------------------------------------------------
 * Normalize the data such that (number of entries) x binsize = 1, and return 
 * the normalized data as a vector.
 * -----------------------------------------------------------------------------
 */
vector <double> BinDat::getNorm(void)
{
	vector <double> result;
	
	for(int i=0; i<(int)data.size(); ++i	)
		result.push_back((double)(data.at(i)/((double)nentries*binsize)));
	return result;
}

/* -----------------------------------------------------------------------------
 * Return the cumulative data as a vector.
 * -----------------------------------------------------------------------------
 */
vector <int> BinDat::getCum(void)
{
	vector <int> result;
	int temp=0;
	for(int i=0; i<(int)data.size(); ++i	)
		result.push_back((temp+=data.at(i)));
	return result;
}

/* -----------------------------------------------------------------------------
 * Return the weighted data as a vector.
 * -----------------------------------------------------------------------------
 */
vector <double> BinDat::getWeighted(double (*weight)(double))
{
	vector <double> result;
	for(int i=0; i<(int)data.size(); ++i	) {
		double x = (0.5+i)*binsize+xmin;
		result.push_back(data.at(i) * weight(x));
	}
	return result;
}

/* -----------------------------------------------------------------------------
 * Return the data as a vector, where each bin is avaraged with its n
 * neighbors on either side.
 * -----------------------------------------------------------------------------
 */
vector <double> BinDat::getSmooth(int n)
{
	vector <double> result;

	for(int i=0; i<(int)data.size();++i) {
		double out=0;
		int d=0;
		for(int j=-n; j<=n; ++j) {
			if(i+j>=0 && i+j<(int)data.size()) {
				out +=data.at(i+j);
				++d;
			}
		}
		result.push_back(out/d);
	}
	return result;	
}

/* -----------------------------------------------------------------------------
 * Write the data to the file <name>.dat, generate a gnuplot file <name>.plot
 * and run gnuplot on that file creating a plot <name>.eps
 * -----------------------------------------------------------------------------
 */
void BinDat::plotSec(string name)
{
	stringstream datfile, plotfile, command;
	
	datfile<<name<<".dat";
	plotfile<<name<<".plot";
	command<<"gnuplot "<<plotfile.str();
	
	outSec(datfile.str().c_str());
	ofstream plot(plotfile.str().c_str());
	if(plot.is_open()){
		plot<<"set style data histeps\n";
		plot<<"set term postscript eps enhanced color\n";
		plot<<"set output \""<<name<<".eps\"\n";
		plot<<"plot \""<<datfile.str()<<"\" w lines\n";
		plot.close();
	
		system(command.str().c_str());	
	} else {
		cout << "Could not open file: " << plotfile.str().c_str();
		cout <<", check for write protection.\n";
	}
}

/* -----------------------------------------------------------------------------
 * Normalize the data such that (number of entries) x binsize = 1, write the 
 * normalized data to the file <name>.dat, generate a gnuplot file <name>.plot
 * and run gnuplot on that file creating a plot <name>.eps
 * -----------------------------------------------------------------------------
 */
void BinDat::plotNorm(string name)
{
	stringstream datfile, plotfile, command;
	
	datfile<<name<<".dat";
	plotfile<<name<<".plot";
	command<<"gnuplot "<<plotfile.str();
	
	outNorm(datfile.str().c_str());
	ofstream plot(plotfile.str().c_str());
	if(plot.is_open()){
		plot<<"set style data histeps\n";
		plot<<"set term postscript eps enhanced color\n";
		plot<<"set output \""<<name<<".eps\"\n";
		plot<<"plot \""<<datfile.str()<<"\" w lines\n";
		plot.close();
		system(command.str().c_str());	
	} else {
		cout << "Could not open file: " << plotfile.str().c_str();
		cout <<", check for write protection.\n";
	}
}

/* -----------------------------------------------------------------------------
 * Write the culumative data to the file <name>.dat, generate a gnuplot file 
 * <name>.plot and run gnuplot on that file creating a plot <name>.eps
 * -----------------------------------------------------------------------------
 */
void BinDat::plotCum(string name)
{
	stringstream datfile, plotfile, command;
	
	datfile<<name<<".dat";
	plotfile<<name<<".plot";
	command<<"gnuplot "<<plotfile.str();
	
	outCum(datfile.str().c_str());
	ofstream plot(plotfile.str().c_str());
	if(plot.is_open()){
		plot<<"set style data histeps\n";
		plot<<"set term postscript eps enhanced color\n";
		plot<<"set output \""<<name<<".eps\"\n";
		plot<<"plot \""<<datfile.str()<<"\" w lines\n";
		plot.close();
	
		system(command.str().c_str());	
	} else {
		cout << "Could not open file: " << plotfile.str().c_str();
		cout <<", check for write protection.\n";
	}
}

/* -----------------------------------------------------------------------------
 * Write the weighted data to the file <name>.dat, generate a gnuplot file 
 * <name>.plot and run gnuplot on that file creating a plot <name>.eps
 * -----------------------------------------------------------------------------
 */
void BinDat::plotWeighted(string name, double (*weight)(double))
{
	stringstream datfile, plotfile, command;
	
	datfile<<name<<".dat";
	plotfile<<name<<".plot";
	command<<"gnuplot "<<plotfile.str();
	
	outWeighted(datfile.str().c_str(),weight);
	ofstream plot(plotfile.str().c_str());
	if(plot.is_open()){
		plot<<"set style data histeps\n";
		plot<<"set term postscript eps enhanced color\n";
		plot<<"set output \""<<name<<".eps\"\n";
		plot<<"plot \""<<datfile.str()<<"\" w lines\n";
		plot.close();
	
		system(command.str().c_str());	
	} else {
		cout << "Could not open file: " << plotfile.str().c_str();
		cout <<", check for write protection.\n";
	}
}

/* -----------------------------------------------------------------------------
 * Write the data to the file <name>.dat, where each bin is avaraged with its n
 * neighbors on either side, generate a gnuplot file <name>.plot and run gnuplot
 * on that file creating a plot <name>.eps.
 * -----------------------------------------------------------------------------
 */
void BinDat::plotSmooth(string name, int d)
{
	stringstream datfile, plotfile, command;
	
	datfile<<name<<".dat";
	plotfile<<name<<".plot";
	command<<"gnuplot "<<plotfile.str();
	
	outSmooth(datfile.str().c_str(),d);
	ofstream plot(plotfile.str().c_str());
	if(plot.is_open()){
		plot<<"set style data histeps\n";
		plot<<"set term postscript eps enhanced color\n";
		plot<<"set output \""<<name<<".eps\"\n";
		plot<<"plot \""<<datfile.str()<<"\" w lines\n";
		plot.close();
	
		system(command.str().c_str());	
	} else {
		cout << "Could not open file: " << plotfile.str().c_str();
		cout <<", check for write protection.\n";
	}
}

/* -----------------------------------------------------------------------------
 * Functions for the dynamic subclass DynBinDat
 * -----------------------------------------------------------------------------
 */
 
 /* -----------------------------------------------------------------------------
 * Add one to the bin where x belongs. If x is out of the current range
 * (xmin, xmax), adjust range, and update current bin-vector data, befor adding 
 * x.
 * -----------------------------------------------------------------------------
 */
void DynBinDat::add(double x)
{	
	if(ddata.size()==0) {
		xmin=xmax=x;
		ddata.push_back(x);
	}else if(x==xmin && xmin==xmax)
		ddata.push_back(x);
	else {
		if(x<xmin){
			xmin=x;
			binsize=(xmax-xmin)/nbins;
			refill();
		} else if (x>xmax) {
			xmax=x;
			binsize=(xmax-xmin)/nbins;
			refill();
		}
		this->BinDat::add(x);
		ddata.push_back(x);
	}
	nentries = ddata.size();
}

/* -----------------------------------------------------------------------------
 * Bin all elements from the input history ddata into data. This is a private
 * function which is meant to be called after a rescaling of the bin-range.
 * -----------------------------------------------------------------------------
 */
void DynBinDat::refill(void)
{
	vector<double>::iterator it;
	data.assign(nbins,0);
	
	for(it=ddata.begin(); it<ddata.end(); it++) {
		if(*it==xmax)
			data.at(nbins-1)++;
		else
			data.at((int)((*it-xmin)/binsize))++;
	}
}

/* -----------------------------------------------------------------------------
 * Rescale the DynBinDat object to the new values of xmin, xmax and nbins. The 
 * new range must enclose the old range to prevent existing data to fall out of 
 * the new range.
 * -----------------------------------------------------------------------------
 */
void DynBinDat::rescale(int nb, double xmi, double xma)
{
	if(xmi>xmin || xma<xmax) {
		cout <<  "Incorrect rescaling: Data out of new bounds.\n";
		return;
	}
	if(nb<=0) {
		cout << "Incorrect rescaling: Invalid new number of bins\n";
		return;
	}
	nbins = nb;
	xmin = xmi;
	xmax = xma;
	
	binsize=(xmax - xmin) / nbins;
	
	refill();
}

/* -----------------------------------------------------------------------------
 * eqRange(DynBinDat *a, DynBinDat *b) rescales the two DynBinDat objects such
 * that they have the same range and number of bins. The maximum of the two 
 * number of bins is taken for the rescaled versions, and no data may fall out
 * of the new ranges, stricly increasing the ranges.
 * -----------------------------------------------------------------------------
 */
void eqRange(DynBinDat *a, DynBinDat *b)
{
	double xmin = fmin(a->xmin, b->xmin);
	double xmax = fmax(a->xmax, b->xmax);
	int nbins  = max(a->nbins, b->nbins);
	a->rescale(nbins, xmin, xmax);
	b->rescale(nbins, xmin, xmax);
}
