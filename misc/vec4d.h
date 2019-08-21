/* -----------------------------------------------------------------------------
 * Name:    vec4d.h (header file of vec4d.cpp)
 * Author:  Bart Verouden
 * Version: 1.0 - 14 October 2009
 *
 * Includes:
 *  -   class vec4d
 * -----------------------------------------------------------------------------
 */

#ifndef VEC4D_H
#define VEC4D_H

class vec4d{
public:
	double E;
	double x;
	double y;
	double z;
	
	/*
	 * Constructors and destructors
	 */
    vec4d();
	vec4d(double E, double p);
	vec4d(double E, double px, double py, double pz);
		
    ~vec4d();
	
	/*
	 * Operators
	 */	 	
	vec4d & operator=(const vec4d &p);
	const vec4d operator+( const vec4d& p);
	const vec4d operator-( const vec4d& p);
	const vec4d operator*(const double r);
	const double operator*(const vec4d p);
	
	/*
	 * General public functions
	 */
	double abs2(void){ return (*this)*(*this);}

};

vec4d operator*(double r, vec4d& p);

#endif
