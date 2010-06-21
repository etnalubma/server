#ifndef DENSITIES_H
#define DENSITIES_H

#include "randomgen.h"

double gamma(double x, double alpha, double beta, double gamma_fn);
double normal(double x, double u, double s2);
double lognormal(double x, double u, double s2);

#endif
