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
