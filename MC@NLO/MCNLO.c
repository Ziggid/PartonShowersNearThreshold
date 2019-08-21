#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "rng.h"
#include "MCNLO.h"

int nbin = 50, nsim = 100000, j;
double a =0.3, B =2, V =1, x0 =0.02, xdead =0.6, sigma = 4.675, binsize;

int main(void){
    int i,k,p, count1, count2, procnr=1;
    double x, y, z[20], temp;

    FILE *monte[4], *excl[4], *incl[4], *unwy[4], *unwz[4], *ymc[4], *ynlo[4],
        *yct[4], *countmc[4], *countnlo[4], *nlo = fopen("nlo.dat", "wt");
    
	/*
     * Initialisation
     */
    
	binsize = (1-x0)/nbin;

    /* initializing files */
    monte[0] = fopen("monte1.dat", "wt");
    monte[1] = fopen("monte2.dat", "wt");
    monte[2] = fopen("monte3.dat", "wt");
    monte[3] = fopen("monte4.dat", "wt");

    excl[0] = fopen("excl1.dat", "wt");
    excl[1] = fopen("excl2.dat", "wt");
    excl[2] = fopen("excl3.dat", "wt");
    excl[3] = fopen("excl4.dat", "wt");

    incl[0] = fopen("incl1.dat", "wt");
    incl[1] = fopen("incl2.dat", "wt");
    incl[2] = fopen("incl3.dat", "wt");
    incl[3] = fopen("incl4.dat", "wt");

    unwy[0] = fopen("unwy1.dat", "wt");
    unwy[1] = fopen("unwy2.dat", "wt");
    unwy[2] = fopen("unwy3.dat", "wt");
    unwy[3] = fopen("unwy4.dat", "wt");

    unwz[0] = fopen("unwz1.dat", "wt");
    unwz[1] = fopen("unwz2.dat", "wt");
    unwz[2] = fopen("unwz3.dat", "wt");
    unwz[3] = fopen("unwz4.dat", "wt");

    ymc[0] = fopen("ymc1.dat", "wt");
    ymc[1] = fopen("ymc2.dat", "wt");
    ymc[2] = fopen("ymc3.dat", "wt");
    ymc[3] = fopen("ymc4.dat", "wt");

    ynlo[0] = fopen("ynlo1.dat", "wt");
    ynlo[1] = fopen("ynlo2.dat", "wt");
    ynlo[2] = fopen("ynlo3.dat", "wt");
    ynlo[3] = fopen("ynlo4.dat", "wt");

    yct[0] = fopen("yct1.dat", "wt");
    yct[1] = fopen("yct2.dat", "wt");
    yct[2] = fopen("yct3.dat", "wt");
    yct[3] = fopen("yct4.dat", "wt");

    countmc[0] = fopen("countmc1.dat", "wt");
    countmc[1] = fopen("countmc2.dat", "wt");
    countmc[2] = fopen("countmc3.dat", "wt");
    countmc[3] = fopen("countmc4.dat", "wt");

    countnlo[0] = fopen("countnlo1.dat", "wt");
    countnlo[1] = fopen("countnlo2.dat", "wt");
    countnlo[2] = fopen("countnlo3.dat", "wt");
    countnlo[3] = fopen("countnlo4.dat", "wt");


    /* initialize random generator */
    PutSeed(-1);

    /* initializing arrays */
    for(j=0; j<4; j++)
        for(i=0; i<nbin; i++) {
            NLO[i][1] = MC[j][i][1] = EXCL[j][i][1] = INCL[j][i][1] =
                UNWY1[j][i][1] = UNWY2[j][i][1] = UNWZ1[j][i][1] =
                UNWZ2[j][i][1] = YMC[j][i][1]= YNLO[j][i][1]= YCT[j][i][1]=
                COUNTMC[j][i][1] = COUNTNLO[j][i][1] = 0;
            NLO[i][0] = MC[j][i][0] = EXCL[j][i][0] = INCL[j][i][0] =
                UNWY1[j][i][0] = UNWY2[j][i][0] = UNWZ1[j][i][0] =
                 UNWZ2[j][i][0] = YMC[j][i][0]= YNLO[j][i][0]= YCT[j][i][0]=
                 COUNTMC[j][i][0] = COUNTNLO[j][i][0] =  x0+(i+0.5)*binsize;
    }

    /*
     * Monte Carlo 1 (Q(x) = 1)
     */
    printf("simulation %d/11\n", procnr++);
    j=0; // j=0 => Q(x) = 1;
    for(i=0; i<nsim; i++) {
        doMCy(1, &y);
        if(y>x0)
            MC[0][(int)((y-x0)/binsize)][1]++;
    }

    /*
     * NLO hit and miss. Not necessary. We might as well just plot it. 
     *
    printf("simulation %d/11\n", procnr++);
    for(i=0; i<nsim; i++) {
        do {
            while((y=Random()) < x0){}
            } while(R(y)*(1-log(y))/(y*(1-log(x0))* R(x0)/x0) <= Random ());
        NLO[(int)((y-x0)/binsize)][1]++;
    }*/


    for(i=0; i<nbin; i++){
        y=NLO[i][0];
        NLO[i][1]= a*R(y)/y;
    }

    /* Renormalising MC and writing to files  */
    for(k=0; k<nbin; k++) {
        MC[0][k][1]*=(sigma*nbin/nsim);
        fprintf(monte[0], "%f %f\n", MC[0][k][0], MC[0][k][1]);
        fprintf(nlo, "%f %f\n", NLO[k][0], NLO[k][1]);
    }

    /*
     * MC for all Q(x) with dead-zone at 0.6
     */
    for(j=1; j<4; j++) {
        printf("simulation %d/11\n", procnr++);
        for(i=0; i<nsim; i++) {
            doMCy(1,&y);
            if(y > x0)
                MC[j][(int)((y-x0)/binsize)][1]++;
        }
        // renormalising, setting bin 0 to unity.
        for(k=nbin-1; k>=0; k--)
            MC[j][k][1] /= MC[j][0][1];
        // writing data to file
        for(k=0; k<nbin; k++)
            fprintf(monte[j],"%f %f\n", MC[j][k][0], MC[j][k][1]);
    }

    /*
     * MC@NLO for the exclusive variable y = max{x1, x2,...,xn) where x1...xn
     * are the emitted photons including the NLO photon.
     */
    for(j=1; j<4; j++) {
        printf("simulation %d/11\n", procnr++);
        for(i=0; i<nsim; i++) {
            x = Random(); // represents the NLO-photon
            doMCy(xm(x),&temp);
            y = max(x,temp);
            if(x > x0 && y > x0){   // add weight
                EXCL[j][(int)((y-x0)/binsize)][1] += F1(x);
                if(x==y){
                    YNLO[j][(int)((y-x0)/binsize)][1] += F1(x);
                    COUNTNLO[j][(int)((y-x0)/binsize)][1]++;
                }
                else {
                    YMC[j][(int)((y-x0)/binsize)][1] += F1(x);
                    COUNTMC[j][(int)((y-x0)/binsize)][1]++;
                }
            }
            doMCy(1,&y);    // counter term
            if(x > x0 && y > x0){    // add counter-weight
                EXCL[j][(int)((y-x0)/binsize)][1] += F2(x);
                YCT[j][(int)((y-x0)/binsize)][1] += F2(x);            
            }
        }
        /* writing data to files */
        for(k=0; k<nbin; k++) {
            EXCL[j][k][1] *= (double)nbin/nsim;
            YMC[j][k][1] *= (double)nbin/nsim;
            YNLO[j][k][1] *= (double)nbin/nsim;
            YCT[j][k][1] *= (double)nbin/nsim;

            fprintf(excl[j], "%f %f\n", EXCL[j][k][0], EXCL[j][k][1]);
            fprintf(ymc[j], "%f %f\n", YMC[j][k][0], YMC[j][k][1]);
            fprintf(ynlo[j], "%f %f\n", YNLO[j][k][0], YNLO[j][k][1]);
            fprintf(yct[j], "%f %f\n", YCT[j][k][0], YCT[j][k][1]);

            fprintf(countmc[j], "%f %f\n", COUNTMC[j][k][0], COUNTMC[j][k][1]);
            fprintf(countnlo[j], "%f %f\n", COUNTNLO[j][k][0], COUNTNLO[j][k][1]);
        }
    }


    /*
     * MC@NLO for the inclusive variable z={x1, x2,...,xn} where x1...xn are the
     * emitted photons including the NLO photon. Hence z is a spectrum of
     * photons that should al contribute to the histogram.
     */
    for(i=0; i<20; i++) z[i]=0; // we allow a maximum of 20 branchings
    for(j=1; j<4; j++) {
        printf("simulation %d/11\n", procnr++);
        for(i=0; i<nsim; i++) {
            x = Random(); // represents the NLO-photon
            doMCz(xm(x),z);
            k=0;
            while(z[k]>x0){   // add the MC photons weighted to the data.
                INCL[j][(int)((z[k]-x0)/binsize)][1] += F1(x);
                k++;
            }
            if(x>x0){         // add the NLO photon weighted to the data.
                INCL[j][(int)((x-x0)/binsize)][1] += F1(x);
            }
            doMCz(1,z);
            k=0;
            while(z[k]>x0){  // add the MC counter-term weighted to the data.
               INCL[j][(int)((z[k]-x0)/binsize)][1] += F2(x);
                k++;
            }
        }
        /* writing data to files */
        for(k=0; k<nbin; k++) {
            INCL[j][k][1] *= (double)nbin/nsim;
            fprintf(incl[j], "%f %f\n", INCL[j][k][0], INCL[j][k][1]);
        }
    }

    /*
     * The unweighted version of the exclusive variable y. The program selects a
     * random x with probability proportional to F1(x)+F2(x). Then it chooses to
     * add a normal event or a counter event with probabilities F1/(F1+F2) and
     * F2/(F1+F2).
     */
    for(j=1; j<4; j++){
        printf("simulation %d/11\n", procnr++);
        count1=count2=0;
        for(i=0; i<nsim; i++) {
            do{ x = Random(); }
            while((F1(x)+F2(x))/(F1(1)+ F2(1))< Random() );
            if((F1(x)/(F1(x)+F2(x))) >= Random()) {
                count1++;
                doMCy(xm(x), &temp);
                y=max(x,temp);
                if(x > x0 && y > x0)
                    UNWY1[j][(int)((y-x0)/binsize)][1]++;
            } else {
                count2++;
                doMCy(1, &y);
                if(y > x0)
                    UNWY1[j][(int)((y-x0)/binsize)][1]++;
            }
        }
        for(k=0; k<nbin; k++){ //normalise to the weigted data and write to file
            UNWY1[j][k][1] *=EXCL[j][nbin-1][1]/UNWY1[j][nbin-1][1];
            fprintf(unwy[j], "%f %f\n", UNWY1[j][k][0], UNWY1[j][k][1]);
        }
    }

    /*
     * The unweighted version of the inclusive variable z. Same method as in the
     * y-case.
     */
     for(i=0; i<20; i++) z[i]=0; // we allow a maximum of 20 branchings
     for(j=1; j<4; j++){
        printf("simulation %d/11\n", procnr++);
        count1=count2=0;
        for(i=0; i<nsim; i++) {
            do {x = Random();}
            while( (F1(x)+F2(x))/(F1(1)+ F2(1))< Random() );
            if((F1(x)/(F1(x)+F2(x))) >= Random()) {
                count1++;
                doMCz(xm(x), z);
                k=0;
                while(z[k]>x0){   // add the MC photons unweighted to the data.
                    UNWZ1[j][(int)((z[k]-x0)/binsize)][1]++;
                    k++;
                }
                if(x>x0){         // add the NLO photon unweighted to the data.
                    UNWZ1[j][(int)((x-x0)/binsize)][1]++;
                }
            } else {
                count2++;
                doMCz(1,z);
                k=0;
                while(z[k]>x0){  // add the MC counter-term unweighted to the data.
                    UNWZ1[j][(int)((z[k]-x0)/binsize)][1]++;
                    k++;
                }
            }

        }
        for(k=0; k<nbin; k++){ //normalise to the weigted data and write to file
            UNWZ1[j][k][1] *=INCL[j][nbin-1][1]/UNWZ1[j][nbin-1][1];
            fprintf(unwz[j], "%f %f\n", UNWZ1[j][k][0], UNWZ1[j][k][1]);
        }
     }

    return 0;
}

/*
 * Q characterizes the sudakov form-factor in the MC generator.
 */
double Q(double x) {
    switch (j) // j should be between 0 and 3 to distinguish the different
               // choices for Q(x). j is defined globally.
    {
        case 0:
        { return 1; }
        case 1:
        { if(xdead - x > 0) return 1; }
        case 2:
        { if(xdead - x > 0) return(G(x/xdead,1)); }
        case 3:
        { if(xdead - x > 0) return(G(x/xdead,2)); }
    }
    return 0;
}


double R(double x) {
    return (B + x * (1 + .5*x + 20*x*x));
}

double xm(double x) {
    return (1-x);
}

double G(double x, int nmbr) {
    if(nmbr == 1)
        return(pow(1-x, 2)/(x*x + pow(1-x,2)));
    return(64*pow(1-x, 2)/(pow(x, 4) + 64*pow(1-x, 2)));
}

double F1(double x)
{
    return (a*(R(x)-B*Q(x))/x);
}

double F2(double x)
{
    return (B+a*(V+B*(Q(x)- 1)/x));
}

/*
 * MC simulation for the variable y. Since the MC photon's are strictly
 * decreasing in energy, only the first photon is a candidate to be the hardest
 * photon. Hence we stop after one branching.
 */

void doMCy(double xm, double *y) {
    while(Q(*y = xm*pow(Random(),1/a)) < Random())
        xm = *y;
    return;

}

/*
 * MC simulation for the variable z. Here we need the full distribution, so we
 * let the system to keep on brancing until it generates a photon with energy
 * less then x0.
 */
void doMCz(double xm, double *z) {
    int i=0;
    while(xm > x0 && i<19) { // we accept a maximum of 19 branchings.
        while(Q(z[i] = xm*pow(Random(),1/a)) < Random())
            xm = z[i];
        xm=z[i];
        i++;
    }
    // resetting al unused positions in the array, so that no old data is used.
    while(i<19 && z[i] != 0) {
    z[i]=0;
    i++;
    }
    return;
}

