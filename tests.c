#include <math.h>
#include <stdlib.h>

#include "special.h"
#include "tests.h"
#include "sort.h"

#define MAX(A,B) A > B ? A : B

int inv_gen(randgen rg, double *prob);
double f_uni(double y);

double chi_square(double a, double x){
    return gammq(a/2, x/2.);
}

double calculate_t(int k, int m, int *n, double *p){
    int i;
    double t, tmp=0;
    
    for(i=0;i<k;i++){
        tmp = pow(n[i] - m*p[i], 2)/(m*p[i]);
        t += tmp;
/*        printf("%i \t %f \t %f \n", n[i], m*p[i], pow(n[i] - m*p[i], 2)/(m*p[i]));*/

    }
    
    return t;
}

double calculate_t_sim(randgen rg, int m, int k, double *p){
    int *n;
    int i;
    double t;
    
    n = calloc(k, sizeof(int));    
    for(i=0; i<m; i++){
        n[inv_gen(rg, p)]++;
    }

    t = calculate_t(k, m, n, p);
    free(n);
    return t;
}

double simulate_t(randgen rg, int s, int m, int k, double *p, double t_zero){
    int j, e=0;
    double t;

    for(j=0; j<s; j++){
        t = calculate_t_sim(rg, m, k, p);
        if(t > t_zero)
            e++;
            
    }
    return e/(double)s;
}

double calculate_d(int n, double *y, double (*f)(double)){
    int i=0;
    double d1=0, d2=0, tmp;

    bublesort(n, y);
    
    for(i=0; i<n; i++){       
        tmp = (i+1)/(double)n - (*f)(y[i]);
        d1 = MAX(d1, tmp);
    }

    for(i=0; i<n; i++){
        tmp = (*f)(y[i]) - i/(double)n;
        d2 = MAX(d2, tmp);        
    }    

    return MAX(d1, d2);
}

double simulate_d(randgen rg, int s, int n, double d_zero){
    double d;
    int i, j, e=0;
    double *u;
    
    for(j=0; j<s; j++){
        u = calloc(n, sizeof(double));

        for(i=0; i<n; i++){
            u[i] = random_gen(rg);
        }

        d = calculate_d(n, u, &f_uni);

        if(d > d_zero)
            e++;
            
        free(u);
    }
    return e/(double)s;
}

double f_uni(double y){
    return y;
}


int inv_gen(randgen rg, double *prob){
    int i = 0;
    double rand, p;
    
    rand = random_gen(rg);
    p = prob[i];
    while(rand >= p){
        i = i+1;
        p = p + prob[i];
    }
    return i;
}
