#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "veto_toy.h"
#include "rng.h"
#define PI 3.1415926535897932385

/*
 * This progam gives a toy model for the Veto Algorithm as described in Pythia
 * 6.4 Physics and Manual.
 */
double tmin=0, told, tmax=1, t, tot=0, ttest, scale=4, z, zold,n;
int nhits=100000, nbins=100;


int main(void)
{
    int i,j,k, dat[101];
    double d, R, r, test, a[20];
    FILE *data1= fopen("data1.dat", "wt"), *data2=fopen("data2.dat", "wt"),
        *data3=fopen("data3.dat", "wt");

    //initializing the data
    for(i=0;i<nbins;i++){
        dat[i]=0;
    }

    //for(i=0;i<7;i++) a[i]=pow(2,i-1);
    for(i=0;i<20;i++) a[i]=.25+.5*(i);
    j=k=0;
    /* initialize random generator */
    PutSeed(-1);

    /*
     * First we are going to do a simple decay algorithm for a function g(t)
     * which has a known primitive that has a well defined inverse. To test that
     * this algorithm works we are going to do nhits=100000 calculations and see
     * what fraction 'r' of the resulting t-values are smaller then ttest. If
     * our code is correct we should see that r is to a good approximation equal
     * to the theoretical value  1-exp(-ttest*ttest/2).
     *
     * As an additional test we are going to see what the average is of all
     * generated numbers (we expect this to be sqrt(0.5*pi)) and see if the
     * deviation is of the order of the standard deviation sqrt((2-0.5*pi)/N).
     *
     * A more visual test is to devide our region of t's into bins and count
     * how many of the generated t's belong in each bin. If we normalize and
     * plot this we can compare the graph with the expected distribution.
     */
    ttest = Random();
    printf("t_test = %f\n\n", ttest);
    printf("\nSimple decay algorithm for g(t) = t\n");
    test = 0;
    for(i=0;i<nhits;i++){
        R=Random();
        t=Ginv(G(tmin)-log(R));
        if(t<scale)                       //refusing t larger then our set scale
            dat[(int)(nbins/scale*t)]++;  //save data
        if(t<ttest)
            k++;
        test = test+t/nhits;
    }
    /*
     * Writing data to file. To normalize we must first devide by nhits*dt =
     * nhits*scale/nbins
     */
    for(i=0;i<=nbins;i++){
        fprintf(data1,"%f %f\n",(double)i*scale/nbins,
            (double)dat[i]*nbins/(nhits*scale));
        dat[i]=0;
    }
    r=(float)k/nhits;
    printf("Simulation:  left=%f, right=%f \n",r,1-r);
    printf("Theoretical: left=%f, right=%f \n", 1-exp(-ttest*ttest/2),
        exp(-ttest*ttest/2));
    printf("\nExpected average: %f, actual average: %f\n",sqrt(0.5*PI),test);
    printf("Standard deviation: %f,  deviation from average: %f\n",
        sqrt((2-0.5*PI)/nhits), fabs(sqrt(0.5*PI)-test));
    k=0;
    test=0;

    /*
     * Now suppose g(t) was not nice enough, i.e. it has no primitive that has
     * a well-defined inverse. Now we can try to find a function h(t) that has a
     * primitive with a well-defined inverse, such that h(t) >= g(t) for t >= 0.
     * If we find such a function we can us the veto algorithm. First we are
     * again going to calculate 100.000 t-values and see if the ratio r matches
     * the theory again. Then we will again look at the standard deviation of
     * the data. At last we will make another plot. We are going to pretend that
     * g(t) = t is not nice enough, but h(t) = t^2 + a is, and since h(t)>=g(t)
     * if a>1/4 we can use the veto algorithm.
     */
    printf("\n\nVeto algorithm for t^2+a > t for t>0\n");
    for(i=0;i<nhits;++i){
        t=j=0;
        do {
            told = t;
            j++;
            do {
                t = Hinv(H(told,a[0])-log(Random()),a[0]);
            } while (t <= told);
        } while (g(t)/h(t,a[0])<= Random());
        if(t<scale)                       //refusing t larger then our set scale
            dat[(int)(nbins/scale*t)]++;  //save data
        if(t<ttest)
            k++;
        test = test+t/nhits;
    }
    /*
     * Writing data to file. To normalize we must first devide by nhits*dt =
     * nhits*scale/nbins
     */
    for(i=0;i<=nbins;i++){
        fprintf(data2,"%f %f\n",(double)i*scale/nbins,
            (double)dat[i]*nbins/(nhits*scale));
        dat[i]=0;
    }
    printf("Simulation:  left=%f, right=%f \n",(float)k/nhits,1-(float)k/nhits);
    printf("Theoretical: left=%f, right=%f \n", 1-exp(-ttest*ttest/2),
        exp(-ttest*ttest/2));
    printf("\nExpected average: %f, actual average: %f\n",sqrt(0.5*PI),test);
    printf("Standard deviation: %f,  deviation from average: %f\n",
        sqrt((2-0.5*PI)/nhits), fabs(sqrt(0.5*PI)-test));

    /*
     * Now we are going to look at how the efficiency changes as h(t) moves
     * further away from g(t). As a measure of the efficiency we are going to
     * use the number of calls made to Random()
     *
     * This test is relatively time consuming, and hence commented out. The test
     * results in a lineair relation between the number of calls and a, for
     * large a. More interesting results migh be found for real singular
     * functions when h is very close to g.
     *
     printf("\nEfficiency test\n");
     for(j=0;j<20;j++){
        printf("*");
        n=0;
        for(i=0;i<nhits;++i){
            t=0;
            do {
                n++
                told = t;
                do {
                    t = Hinv(H(told,a[j])-log(Random()),a[j]);
                    n++;
                } while (t <= told);
            } while (g(t)/h(t,a[j])<= Random());
        }
        fprintf(data3,"%f %e\n",a[j], n);
     }


    /*
     * Now we do a simple particle branching simulation where the distribution
     * of z is x^2 on the interval (0,1)
     */
    printf("\n\nThe next simulation gives a list of (t,z) branchings up until");
    printf("l\nit reaches the scale = %.1f\n\n", scale);
    t=told=0;
    z=1;
    while(t <= scale)
        {
         printf("parton is at (%f, %f)\n", t, z);
         do {
            told = t;
            do {
                t = Hinv(H(told,1)-log(Random()),1);
            } while (t <= told);
         } while (g(t)/h(t,1)<= Random());
         z=pow(3*(Random()*(1./3.)),1./3.)*z; // Some arbitrary distribution to
    }                                         // choose a new z(<z_old))
    return 0;
}

/*
 * g(t) original decay probability
 */
double g(double t) {
    return t;
}

/*
 * the primitive of g(t)
 */
double G(double t) {
    return t*t/2;
}

/*
 * the inverse of the primitive of g(t)
 */
double Ginv(double x) {
    return pow(2*x, 1./2.);
}

/*
 * h(t) = t^2+a >= g(t) for a >= 1/4  (t>0)
 */
double h(double t, double a) {
    return(t*t+a);
}

/*
 * the primitive of h(t)
 */
double H(double t, double a) {
    return(t*t*t/3 + a*t);
}

/*
 * the inverse of the primitive of h(t)
 */
double Hinv(double x, double a) {
    return(-pow(2*a*a*a/(3*x + sqrt(4*a*a*a+9*x*x)),1./3.)
        + pow((3*x+sqrt(4*a*a*a+9*x*x))/2, 1./3.));
}
