#include <math.h>
#include <stdlib.h>
#include <stdarg.h>

#include "randomgen.h"
#include "densities.h"
#define PI 3.14159265

double gamma(double x, double alpha, double beta, double gamma_fn){
    if(x<=0)
        return 0;
    else
        return (pow(beta, -alpha)*pow(x, alpha-1)*exp(-x/beta))/gamma_fn;
}

double normal(double x, double u, double s2){
    return (1./sqrt(2*PI*s2))*exp(-(pow(x-u, 2))/(2*s2));
}

double lognormal(double x, double u, double s2){
    if(x<=0)
        return 0;
    else
        return (1./(x*sqrt(2*PI*s2)))*exp(-(pow(log(x)-u, 2))/(2*s2));
}

