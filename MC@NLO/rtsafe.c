#include <math.h>
#include <stdio.h>
#define MAXIT 100



double rtsafe(void (*funcd)(double, double, double, double *, double *),
    double x1, double x2, double xacc, double R)
{
	int j;
	double temp,xh,xl,rts;
	double df,dx,dxold,f,fh,fl;

	(*funcd)(x1,R,x2,&fl,&df);
   	(*funcd)(x2,R,x2,&fh,&df);
	if (fl > 0.0 && fh > 0.0){
        // root lies below x0, returning bogus result smaller then x0
        return x1/2;
    }
    if(fl < 0.0 && fh < 0.0){
        printf("something wrong\n");
        return 0.0;
    }
	if (fl == 0.0) return x1;
	if (fh == 0.0) return x2;
	if (fl < 0.0) {
		xl=x1;
		xh=x2;
	} else {
		xh=x1;
		xl=x2;
	}
	rts=0.5*(x1+x2);
	dxold=fabs(x2-x1);
	dx=dxold;
	(*funcd)(rts,R,x2,&f,&df);
	for (j=1;j<=MAXIT;j++) {
		if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0)
			|| (fabs(2.0*f) > fabs(dxold*df))) {
			dxold=dx;
			dx=0.5*(xh-xl);
			rts=xl+dx;
			if (xl == rts) return rts;
		} else {
			dxold=dx;
			dx=f/df;
			temp=rts;
			rts -= dx;
			if (temp == rts) return rts;
		}
		if (fabs(dx) < xacc) return rts;
		(*funcd)(rts,R,x2,&f,&df);
		if (f < 0.0)
			xl=rts;
		else
			xh=rts;
	}
	printf("Maximum number of iterations exceeded in rtsafe");
	return 0.0;
}
#undef MAXIT
