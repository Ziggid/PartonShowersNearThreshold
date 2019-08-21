/**************************************************************************
* This file is copied from Pythia 8.1 with a few minor adjustments
*/
#ifndef _VEC4_
#define _VEC4_

// Stdlib header files for mathematics.
#include <cmath>
#include <cstdlib>

// Stdlib header files for strings and containers.
#include <string>
#include <vector>
#include <map>
#include <deque>

// Stdlib header file for input and output.
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

// Define pi if not yet done.
#ifndef M_PI
#define M_PI 3.1415926535897932385
#endif
using namespace std;

// Powers of small integers - for balance speed/code clarity.
inline double pow2(const double& x) {return x*x;}
inline double pow3(const double& x) {return x*x*x;}
inline double pow4(const double& x) {return x*x*x*x;}
inline double pow5(const double& x) {return x*x*x*x*x;}
inline double pow6(const double& x) {return x*x*x*x*x*x;}

//**************************************************************************

// Forward reference to RotBstMatrix class.
class RotBstMatrix;

//**************************************************************************

// Vec4 class.
// This class implements four-vectors, in energy-momentum space.
// (But can equally well be used to hold space-time four-vectors.)

class Vec4 {

public:

  // Constructors.
  Vec4(double xIn = 0., double yIn = 0., double zIn = 0., double tIn = 0.)
    : xx(xIn), yy(yIn), zz(zIn), tt(tIn) { }
  Vec4(const Vec4& v) : xx(v.xx), yy(v.yy), zz(v.zz), tt(v.tt) { }
  Vec4& operator=(const Vec4& v) { if (this != &v) { xx = v.xx; yy = v.yy; 
    zz = v.zz; tt = v.tt; } return *this; }
  Vec4& operator=(double value) { xx = value; yy = value; zz = value; 
    tt = value; return *this; }
      
  // Member functions for input.
  void reset() {xx = 0.; yy = 0.; zz = 0.; tt = 0.;}
  void p(double xIn, double yIn, double zIn, double tIn) 
    {xx = xIn; yy = yIn; zz = zIn; tt = tIn;}
  void p(Vec4 pIn) {xx = pIn.xx; yy = pIn.yy; zz = pIn.zz; tt = pIn.tt;} 
  void px(double xIn) {xx = xIn;}
  void py(double yIn) {yy = yIn;}
  void pz(double zIn) {zz = zIn;}
  void e(double tIn) {tt = tIn;}

  // Member functions for output.
  double px() const {return xx;}
  double py() const {return yy;}
  double pz() const {return zz;}
  double e() const {return tt;}
  double mCalc() const {double temp = tt*tt - xx*xx - yy*yy - zz*zz;
    return (temp >= 0.) ? sqrt(temp) : -sqrt(-temp);}
  double m2Calc() const {return tt*tt - xx*xx - yy*yy - zz*zz;}
  double pT() const {return sqrt(xx*xx + yy*yy);}
  double pT2() const {return xx*xx + yy*yy;}
  double pAbs() const {return sqrt(xx*xx + yy*yy + zz*zz);}
  double pAbs2() const {return xx*xx + yy*yy + zz*zz;}
  double eT() const {double temp = xx*xx + yy*yy;
    return tt * sqrt( temp / (temp + zz*zz) );}
  double eT2() const {double temp = xx*xx + yy*yy;
    return tt*tt * temp / (temp + zz*zz);}
  double theta() const {return atan2(sqrt(xx*xx + yy*yy), zz);}
  double phi() const {return atan2(yy,xx);}
  double thetaXZ() const {return atan2(xx,zz);}
  double pPos() const {return tt + zz;}
  double pNeg() const {return tt - zz;}

  // Member functions that perform operations.
  void rescale3(double fac) {xx *= fac; yy *= fac; zz *= fac;}
  void rescale4(double fac) {xx *= fac; yy *= fac; zz *= fac; tt *= fac;}
  void flip3() {xx = -xx; yy = -yy; zz = -zz;}
  void flip4() {xx = -xx; yy = -yy; zz = -zz; tt = -tt;}
  void rot(double thetaIn, double phiIn); 
  void rotaxis(double phiIn, double nx, double ny, double nz); 
  void rotaxis(double phiIn, const Vec4& n);
  void bst(double betaX, double betaY, double betaZ); 
  void bst(double betaX, double betaY, double betaZ, double gamma); 
  void bst(const Vec4& pIn); 
  void bst(const Vec4& pIn, double mIn); 
  void bstback(const Vec4& pIn); 
  void bstback(const Vec4& pIn, double mIn); 
  void rotbst(const RotBstMatrix& M); 

  // Operator overloading with member functions
  Vec4 operator-() {Vec4 tmp; tmp.xx = -xx; tmp.yy = -yy; tmp.zz = -zz; 
    tmp.tt = -tt; return tmp;}
  Vec4& operator+=(const Vec4& v) {xx += v.xx; yy += v.yy; zz += v.zz; 
    tt += v.tt; return *this;}
  Vec4& operator-=(const Vec4& v) {xx -= v.xx; yy -= v.yy; zz -= v.zz; 
    tt -= v.tt; return *this;}
  Vec4& operator*=(double f) {xx *= f; yy *= f; zz *= f; 
    tt *= f; return *this;}
  Vec4& operator/=(double f) {xx /= f; yy /= f; zz /= f; 
    tt /= f; return *this;}

  // Operator overloading with friends
  friend Vec4 operator+(const Vec4& v1, const Vec4& v2);
  friend Vec4 operator-(const Vec4& v1, const Vec4& v2);
  friend Vec4 operator*(double f, const Vec4& v1);
  friend Vec4 operator*(const Vec4& v1, double f);
  friend Vec4 operator/(const Vec4& v1, double f);
  friend double operator*(const Vec4& v1, const Vec4& v2);

  // Invariant mass of a pair and its square.
  friend double m(const Vec4& v1, const Vec4& v2);
  friend double m2(const Vec4& v1, const Vec4& v2);

  // Scalar and cross product of 3-vector parts.
  friend double dot3(const Vec4& v1, const Vec4& v2);
  friend Vec4 cross3(const Vec4& v1, const Vec4& v2);

  // theta is polar angle between v1 and v2.
  friend double theta(const Vec4& v1, const Vec4& v2);
  friend double costheta(const Vec4& v1, const Vec4& v2);

  // phi is azimuthal angle between v1 and v2 around z axis.
  friend double phi(const Vec4& v1, const Vec4& v2);  
  friend double cosphi(const Vec4& v1, const Vec4& v2);

  // phi is azimuthal angle between v1 and v2 around n axis.
  friend double phi(const Vec4& v1, const Vec4& v2, const Vec4& n);
  friend double cosphi(const Vec4& v1, const Vec4& v2, const Vec4& n);

  // Print a four-vector
  friend ostream& operator<<(ostream&, const Vec4& v) ;

private:

  // Constants: could only be changed in the code itself.
  static const double TINY;

  // The four-vector data members.
  double xx, yy, zz, tt;

};

// Implementation of operator overloading with friends.

inline Vec4 operator+(const Vec4& v1, const Vec4& v2) 
  {Vec4 v = v1 ; return v += v2;}

inline Vec4 operator-(const Vec4& v1, const Vec4& v2) 
  {Vec4 v = v1 ; return v -= v2;}

inline Vec4 operator*(double f, const Vec4& v1) 
  {Vec4 v = v1; return v *= f;}

inline Vec4 operator*(const Vec4& v1, double f) 
  {Vec4 v = v1; return v *= f;}

inline Vec4 operator/(const Vec4& v1, double f) 
  {Vec4 v = v1; return v /= f;}

inline double operator*(const Vec4& v1, const Vec4& v2)
  {return v1.tt*v2.tt - v1.xx*v2.xx - v1.yy*v2.yy - v1.zz*v2.zz;}  

//**************************************************************************

// RotBstMatrix class.
// This class implements 4 * 4 matrices that encode an arbitrary combination
// of rotations and boosts, that can be applied to Vec4 four-vectors.

class RotBstMatrix {

public:

  // Constructors.
  RotBstMatrix() {for (int i = 0; i < 4; ++i) { for (int j = 0; j < 4; ++j) 
    { M[i][j] = (i==j) ? 1. : 0.; } } } 
  RotBstMatrix(const RotBstMatrix& Min) {
    for (int i = 0; i < 4; ++i) { for (int j = 0; j < 4; ++j) {
    M[i][j] = Min.M[i][j]; } } }
  RotBstMatrix& operator=(const RotBstMatrix& Min) {if (this != &Min) {
    for (int i = 0; i < 4; ++i) { for (int j = 0; j < 4; ++j) {
    M[i][j] = Min.M[i][j]; } } } return *this; }

  // Member functions.
  void rot(double = 0., double = 0.);
  void rot(const Vec4& p);
  void bst(double = 0., double = 0., double = 0.);
  void bst(const Vec4&);
  void bstback(const Vec4&);
  void bst(const Vec4&, const Vec4&);
  void toCMframe(const Vec4&, const Vec4&);
  void fromCMframe(const Vec4&, const Vec4&);
  void rotbst(const RotBstMatrix&);
  void invert();
  void reset();
  
  double gamma(){ return M[0][0];} // Added by Bart Verouden.

  // Crude estimate deviation from unit matrix.
  double deviation() const;

  // Print a transformation matrix.
  friend ostream& operator<<(ostream&, const RotBstMatrix&) ;

  // Private members to be accessible from Vec4. 
  friend class Vec4;

private:

  // Constants: could only be changed in the code itself.
  static const double TINY;

  // The rotation-and-boost matrix data members.
  double M[4][4];

};
#endif //_VEC4_
