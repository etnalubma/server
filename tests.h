#ifndef TESTS_H
#define TESTS_H
#include <stdarg.h>
#include "randomgen.h"

double chi_square(double a, double x);

double calculate_t(int k, int m, int *n, double *p);
double simulate_t(randgen rg, int s, int m, int k, double *p, double t_zero);
double calculate_t_sim(randgen rg, int m, int k, double *p);
double calculate_d(int n, double *y, double (*f)(double));
double simulate_d(randgen rg, int s, int n, double d_zero);

#endif
