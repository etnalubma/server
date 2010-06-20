#include <math.h>
#include <stdlib.h>

#include "estimate.h"

struct sestimate{
    float (*generator)(randgen);
    int n;
    float d;
    int (*cont)(float, float, int);
    randgen rg;
};

estimate estimateInit(randgen rg, float (*generator)(randgen), int n, float d, int (*cont)(float, float, int)){
    estimate e;
    
    e = (estimate)calloc(1, sizeof(struct sestimate));
    e->generator = generator;
    e->n = n;
    e->d = d;
    e->cont = cont;
    e->rg = rg;
    
    return e;    
}

estimate estimateDestroy(estimate e){
    e->generator = NULL;
    free(e);
    e = NULL;
    return e;
}

void sample_mean_sigma(estimate e, float *mean, float *sigma, int *samples){
    float x, s2, m, a;
    int i;
    x = e->generator(e->rg);
    m = x;
    s2 = 0;
    
    for(i=2; i<e->n; i++){
        x = e->generator(e->rg);
        a = m;
        m = m + (x-m)/(float)i;
        s2 = (1. - 1./(i-1))*s2 + i*pow(m-a, 2);            
    }
    while(e->cont(e->d, sqrt(s2), i)){
        i++;
        x = e->generator(e->rg);
        a = m;
        m = m + (x-m)/(float)i;
        s2 = (1. - 1./(i-1))*s2 + i*pow(m-a, 2);            
    }
    *mean = m;
    *sigma = sqrt(s2);
    *samples = i;
}

void sample_prob(estimate e, float *prob, int *samples){
    float x, p;
    int i;
    x = e->generator(e->rg);
    p = x;
    
    for(i=1; i<e->n; i++){
        x = e->generator(e->rg);
        p = p + (x-p)/(float)i;
    }
    while(e->cont(e->d, p, i)){
        i++;
        x = e->generator(e->rg);
        x = e->generator(e->rg);
        p = p + (x-p)/(float)i;
    }
    *prob = p;
    *samples = i;    
}

int confidence_interval(float sigma, float z, float delta) {   
    return delta < 2*z*sigma;
}

double estimate_mean(int n, double *sample){
    int i=0;
    double sum = 0;
    
    for(i=0; i<n; i++){
        sum += sample[i];
    }
    return sum/(double)n;
}

double estimate_s2(int n, double *sample){
    int i = 0;
    double sum = 0;
    double mean = estimate_mean(n, sample);
    
    for(i=0; i<n; i++){
        sum += pow(sample[i] - mean, 2);
    }
    
    return sum/(double)(n-1);
}

double estimate_median(int n, double *sample){    
    return n%2 == 0 ? (sample[n/2 - 1] + sample[n/2])/2. : sample[(n+1)/2 - 1];
}

double estimate_skewness(int n, double *sample){
    int i = 0;
    double sum = 0;
    double mean = estimate_mean(n, sample);
    double s2 = estimate_s2(n, sample);
    
    for(i=0; i<n; i++){
        sum += pow(sample[i] - mean, 3) / (double)n;
    }
    
    return sum/pow(s2, 3./2.);
}


double estimate_ln_mean(int n, double *sample){
    int i=0;
    double sum = 0;
    
    for(i=0; i<n; i++){
        sum += log(sample[i]);
    }
    return sum/(double)n;
}

double estimate_ln_s2(int n, double *sample){
    int i = 0;
    double sum = 0;
    double mean = estimate_ln_mean(n, sample);
    
    for(i=0; i<n; i++){
        sum += pow(log(sample[i]) - mean, 2);
    }
    
    return sum/(double)n;

}

double estimate_gamma_t(int n, double *sample){
    double ln_mean = estimate_ln_mean(n, sample);
    double mean = estimate_mean(n, sample);
    
    double t = 1./(log(mean) - ln_mean);
    return t;
}
