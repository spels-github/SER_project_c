/*
 ============================================================================
 Name        : integrand.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : c library for python SER project
 ============================================================================
 */

#include <math.h>
#include <stdio.h>

static double bessi0( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) and n=0.  */
/*------------------------------------------------------------*/
{
   double ax,ans;
   double y;

   if ((ax=fabs(x)) < 3.75) {
      y=x/3.75,y=y*y;
      ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
         +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
   } else {
      y=3.75/ax;
      ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
         +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
         +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
         +y*0.392377e-2))))))));
   }
   return ans;
}


double i_c(double arg, void *user_data) {
	double x = ((double)arg);
	double R = ((double *)user_data)[0];
    double r0 = ((double *)user_data)[1];
    double Da = ((double *)user_data)[2];
    double t = ((double *)user_data)[3];
    return (x / sqrt(pow(R,2) - pow(x,2)))*exp(-((pow(r0,2) + pow(x,2)) / (4*Da*t))) * bessi0(x*r0 / (2*Da*t));

}

double i_b(void *user_data) {
	double L = ((double *)user_data)[0];
    double z0 = ((double *)user_data)[1];
    double Da = ((double *)user_data)[2];
    double t = ((double *)user_data)[3];
    return exp(-(pow(z0,2)/(4*Da*t))-(t*Da*pow((3.1415/L),2)));

}
