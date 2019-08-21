/* -----------------------------------------------------------------------------
 * This file contains the constructors and functions for the vec4d class.
 * Name:    vec4d.cpp
 * Author:  Bart Verouden
 * Version: 1.0 - 14 October 2009
 * Uses: vec4d.h
 *
 * Includes:
 *  -   class vec4d
 * -----------------------------------------------------------------------------
 */
#include "vec4d.h"

/* -----------------------------------------------------------------------------
 * Constructors and destructors
 * -----------------------------------------------------------------------------
 */

/*
 * The default constructor initializes all values at zero.
 */ 
vec4d::vec4d()
{
	E=0.;
	x=0.;
	y=0.;
	z=0.;
}

/*
 * The vec4d(E, pz) constructor creates the 4-vector (E,0,0,p).
 */
vec4d::vec4d(double E, double pz)
{
	this->E=E;
	x=0.;
	y=0.;
	z=pz;
}

/*
 * The vec4d(E, px, py, pz) constructor creates the 4-vector (E,px,py,pz).
 */
vec4d::vec4d(double E, double px, double py, double pz)
{
	this->E=E;
	x=px;
	y=py;
	z=pz;
}

/*
 * The destructor.
 */
vec4d::~vec4d()
{
}

/* -----------------------------------------------------------------------------
 * Operators
 * -----------------------------------------------------------------------------
 */	

/*
 * The '=' operator sets all elements of the lhs 4-vector equal to the elements
 * of the rhs 4-vector.
 */
vec4d & vec4d::operator =(const vec4d &p)
{
	if(this != &p) {
		E=p.E;
		x=p.x;
		y=p.y;
		z=p.z;
	}
	return *this;
}

/*
 * The '+' operator.
 */
const vec4d vec4d::operator+(const vec4d &p)
{
	return vec4d(E+ p.E,x+p.x, y+p.y, z+p.z);;
}

/*
 * The '-' operator.
 */
const vec4d vec4d::operator-(const vec4d &p)
{
	return vec4d(E-p.E,x-p.x, y-p.y, z-p.z);

}

/*
 * The '* double r' operator returns a 4-vector of which all elements are
 * multiplied by r.
 */
const vec4d vec4d::operator*(const double r)
{
	return vec4d(r*E,r*x,r*y,r*z);

}

/*
 * vec4d p * vec4d q returns the Minkosvkian dot product of 4-vectors p and q.
 */
const double vec4d::operator*(vec4d p)
{
	return E*p.E -x*p.x -y*p.y-z*p.z;
}

/*
 * double r * vec4d p = p * r.
 */
vec4d operator*(double r, vec4d& p)
{
	return p*r;
}

