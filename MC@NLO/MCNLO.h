/* -----------------------------------------------------------------------------
 * Header file for MCNLO.c and logterm.c
 * -----------------------------------------------------------------------------
 */

// Arrays for data storage.
double NLO[50][2], MC[4][50][2], EXCL[4][50][2], INCL[4][50][2];
double UNWZ1[4][50][2], UNWZ2[4][50][2], UNWY1[4][50][2], UNWY2[4][50][2];
double YMC[4][50][2], YNLO[4][50][2], YCT[4][50][2];
double COUNTMC[4][50][2], COUNTNLO[4][50][2];

double Q(double x);
double R(double x);
double xm(double x);
double G(double x, int nmbr);
double F1(double x);
double F2(double x);
void doMCy(double xm, double *y);
void doMCz(double xm, double *z);
void f(double z, double R, double xm, double *f, double *df);
double rtsafe(void (*funcd)(double, double, double, double *, double *),
    double x1, double x2, double xacc, double R);
 
