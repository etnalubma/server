#include <math.h>
#include <stdlib.h>
#include <stdarg.h>

#include "randomgen.h"
#include "densities.h"
#define PI 3.14159265

int factorial(int i){return i == 0 ? 1 : i*factorial(i-1);}

double gamma(double x, double alpha, double beta, double gamma_fn){
    if(x<=0)
        return 0;
    else
        return (pow(beta, -alpha)*pow(x, alpha-1)*exp(-x/beta))/gamma_fn;
}

double normal(double x, double u, double s2){
    return (1./sqrt(2*PI*s2))*exp(-(pow(x-u, 2))/2*s2);
}

double lognormal(double x, double u, double s2){
    if(x<=0)
        return 0;
    else
        return (1./x*sqrt(2*PI*s2))*exp(-(pow(log(x)-u, 2))/2*s2);
}

double cumulative_normal(double x, double u, double s2){
    /*
    Abramowitz & Stegun (1964) give the approximation for Φ(x) with the absolute error |ε(x)| < 7.5·10−8
    http://www.math.sfu.ca/~cbm/aands/page_932.htm
    */
    int i;
    double z, t, bp;
    double b[6] = {0.2316419, 0.319381530, -0.356563782, 1.781477937, -1.821255978, 1.330274429};
    z = (x - u)/sqrt(s2);
    
    if(z<0)
        return 1-cumulative_normal(-z, 0, 1);    
    
    t = 1./(1. + b[0]*z);
    for(i=1;i<6;i++)
        bp += b[i]*pow(t, i);
        
    return 1 - normal(z, 0, 1)*bp;
}

double cumulative_gamma(double x, double alpha, double beta){
    double sum=0;
    int i;
    
    if(x<0)
        return 0;

    for(i=0; i<alpha-1; i++){
        sum += pow(x/beta, i)/factorial(i);
    }
    
    return 1 - exp(-x/beta)*sum;
    
}
