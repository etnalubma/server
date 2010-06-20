#ifndef ESTIMATE_H
#define ESTIMATE_H

#include <stdarg.h>
#include "randomgen.h"

#define Z95 1.96
#define Z99 2.576

typedef struct sestimate *estimate;

estimate estimateInit(randgen rg, float (*generator)(randgen), int n, float d, int (*cont)(float, float, int));

estimate estimateDestroy(estimate e);

void sample_mean_sigma(estimate e, float *mean, float *sigma, int *samples);

void sample_prob(estimate e, float *prob, int *samples);

int confidence_interval(float sigma, float z, float delta);

double estimate_mean(int n, double *sample);

double estimate_s2(int n, double *sample);

double estimate_median(int n, double *sample);

double estimate_skewness(int n, double *sample);

double estimate_ln_mean(int n, double *sample);

double estimate_ln_s2(int n, double *sample);

double estimate_gamma_t(int n, double *sample);

#endif
